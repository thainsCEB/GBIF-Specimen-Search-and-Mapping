#!/usr/bin/env Rscript
# ============================================================
# GBIF phylogenetic sampling  ->  one specimen per taxon
#
# A specialised sibling of gbif_population_sampling.R for building a
# phylogenetic sampling sheet: it returns ONE best specimen per
# species / subspecies, preferring tissue over skin, and searching
# collections in a priority order (e.g. FMNH first, then ANSP, then
# anywhere). Eggs are excluded. Synonyms are covered automatically
# (GBIF occurrences are indexed under the accepted taxon key).
#
# Selection, per taxon (best -> worst), default order:
#   1. material            tissue > liquid nitrogen > skin > (other)
#   2. collection order    the --collections list, then anywhere
#   3. coordinate quality  GPS > province centroid > country centroid
#   4. has a collection date
#   5. smaller coordinate uncertainty
# Use --collection_first to make the collection cascade outrank material
# (FMNH-anything before ANSP-anything, with tissue>skin within a museum).
#
# Packages: optparse, rgbif, dplyr, readr, stringr
#   --geocode also needs: sf, rnaturalearth
#
# Examples
# --------
#   # One tissue-first sample per taxon, FMNH then ANSP then anywhere:
#   Rscript gbif_phylo_sampling.R \
#     --taxon "Struthio camelus,Struthio molybdophanes" \
#     --collections "FMNH,ANSP" \
#     --output ostrich_phylo.csv
#
#   # A single genus, one per subspecies, native range only, with coords:
#   Rscript gbif_phylo_sampling.R --taxon Pyrrhura \
#     --rank_level subspecies --has_coordinate true \
#     --collections "FMNH,ANSP,AMNH" --output pyrrhura_phylo.csv
# ============================================================
suppressPackageStartupMessages({
  library(optparse); library(rgbif); library(dplyr); library(readr); library(stringr)
})

OSTRICH_NATIVE <- c("ZA","KE","BW","TZ","SD","EG","NA","CD","SL","DZ","SO","ET",
                    "ZW","UG","ML","LY","GH","TD","SN","NG","NE","MA",
                    "YE","SA","SY","IL","JO","OM","IQ")
PRESETS <- list("ostrich-native" = list(countries = OSTRICH_NATIVE))

option_list <- list(
  make_option(c("-t", "--taxon"), type = "character", default = NULL,
              help = "Scientific name or comma-separated list (synonyms covered) [required unless --taxon_key]."),
  make_option("--taxon_key", type = "character", default = NULL,
              help = "Comma-separated GBIF taxonKeys (skips name resolution)."),
  make_option(c("-k", "--keep_taxa"), type = "character", default = NULL,
              help = paste("Comma-separated taxa to KEEP; all others are dropped. Synonyms",
                           "are resolved (via the GBIF backbone) so records under synonymous",
                           "names are kept too. e.g. 'Crocodylus niloticus,Crocodylus suchus'.")),
  make_option(c("-K", "--keep_taxa_file"), type = "character", default = NULL,
              help = "File (one taxon per line, or a CSV with a taxon/species column) to KEEP; combined with --keep_taxa."),
  make_option("--restrict_taxon", action = "store_true", default = TRUE,
              help = "Keep only records whose genus matches --taxon/--keep_taxa (drops off-target --include hits) [default: on]."),
  make_option("--no_restrict_taxon", action = "store_false", dest = "restrict_taxon",
              help = "Allow --include specimens of any taxon (full filter bypass)."),
  make_option(c("-r", "--rank_level"), type = "character", default = "taxon",
              help = "One sample per: 'taxon' (species+subspecies), 'species', or 'subspecies' [default: taxon]."),
  make_option("--per_taxon", type = "integer", default = 1L,
              help = "How many specimens to keep per taxon [default: 1]."),
  make_option(c("-c", "--collections"), type = "character", default = NULL,
              help = "Institution search order, e.g. 'FMNH,ANSP'. Always searched in this order (collection order is the primary ranking). By default only these collections are used."),
  make_option(c("-f", "--collection_first"), action = "store_true", default = FALSE,
              help = "Search the listed --collections in order first, THEN fall back to other collections not in the list (turns off the default restriction)."),
  make_option("--restrict_collections", action = "store_true", default = TRUE,
              help = "Use only the listed --collections [default: on]. --collection_first or --no_restrict_collections also considers others (after the list)."),
  make_option("--no_restrict_collections", action = "store_false", dest = "restrict_collections",
              help = "Also consider collections not in --collections (ranked after the listed ones)."),
  make_option(c("-p", "--preparations"), type = "character", default = "tissue,liquid nitrogen,skin",
              help = "Material preference/filter, best first [default: 'tissue,liquid nitrogen,skin']."),
  make_option("--keep_unknown_prep", action = "store_true", default = TRUE,
              help = "Allow specimens with no recorded preparation as a last resort [default: on]."),
  make_option("--drop_unknown_prep", action = "store_false", dest = "keep_unknown_prep",
              help = "Never use specimens with no recorded preparation."),
  make_option(c("-x", "--exclude"), type = "character", default = "",
              help = "Per-institution catalogue-number exclusions 'CODE:regex;...' (off by default; e.g. 'NHMUK:E' for NHMUK eggs)."),
  make_option("--no_exclude_eggs", action = "store_false", dest = "exclude_eggs", default = TRUE,
              help = "Do not drop egg collections/preparations (eggs are dropped by default)."),
  make_option("--countries", type = "character", default = NULL,
              help = "Comma-separated ISO-2 country codes to restrict to."),
  make_option(c("-e", "--preset"), type = "character", default = NULL,
              help = "Named country preset (e.g. 'ostrich-native')."),
  make_option("--basis", type = "character", default = "PRESERVED_SPECIMEN",
              help = "Comma-separated basisOfRecord ['any' for all] [default: PRESERVED_SPECIMEN]."),
  make_option("--has_coordinate", type = "character", default = "any",
              help = "Restrict to records with/without coordinates: any|true|false [default: any]."),
  make_option(c("-g", "--geocode"), type = "character", default = "none",
              help = "Geocode coordinate-less records: none|country|state [default: none]."),
  make_option("--no_jitter", action = "store_true", default = FALSE),
  make_option(c("-C", "--country_centroid"), action = "store_true", default = FALSE,
              help = paste("Write a SEPARATE CSV in which specimens with no GPS are filled",
                           "with their country centroid (main output keeps GPS only).")),
  make_option("--centroid_output", type = "character", default = NULL,
              help = "Path for the country-centroid CSV [default: <output>_centroid.csv]."),
  make_option("--limit_total", type = "integer", default = 0L),
  make_option(c("-I", "--include"), type = "character", default = NULL,
              help = paste("CSV of must-include specimens (columns: institution, species",
                           "[+ optional subspecies], catalogue). Added after selection",
                           "(looked up on GBIF to fill columns, or from the given fields).")),
  make_option(c("-y", "--include_only"), action = "store_true", default = FALSE,
              help = "Output ONLY the --include specimens (no gap-fill even if --keep_taxa is given)."),
  make_option("--include_add_all", action = "store_true", default = FALSE,
              help = "With --include, also add the full automatic one-per-taxon genus selection (old behavior)."),
  make_option("--summary", type = "character", default = NULL,
              help = "Summary-counts TSV path [default: <output>_summary.tsv]."),
  make_option(c("-o", "--output"), type = "character", default = "gbif_phylo_samples.csv"),
  make_option("--quiet", action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "GBIF -> one specimen per taxon for phylogenetics."))
say <- function(...) if (!opt$quiet) message(...)
splitcsv <- function(x) if (is.null(x)) NULL else trimws(strsplit(x, ",")[[1]])
`%||%` <- function(a, b) if (is.null(a) || length(a)==0 || (length(a)==1 && is.na(a))) b else a

# Retry wrapper for transient GBIF errors (e.g. "Service Unavailable" 503).
gbif_try <- function(fn, tries = 6, quiet = FALSE) {
  for (a in seq_len(tries)) {
    r <- tryCatch(fn(), error = function(e) e)
    if (!inherits(r, "error")) return(r)
    if (a == tries) stop("GBIF request failed after ", tries, " tries: ", conditionMessage(r))
    wait <- min(60, 2^a)
    if (!quiet) message("  GBIF unavailable (", conditionMessage(r), "); retry ", a, "/", tries - 1, " in ", wait, "s ...")
    Sys.sleep(wait)
  }
}

# ── Resolve taxon keys (use accepted key so synonyms are covered) ──
if (!is.null(opt$taxon_key)) {
  taxon_keys <- as.integer(splitcsv(opt$taxon_key))
} else if (!is.null(opt$taxon)) {
  say("Resolving taxa (synonyms are searched via the accepted key) ...")
  taxon_keys <- integer(0)
  for (nm in splitcsv(opt$taxon)) {
    m <- gbif_try(function() name_backbone(name = nm), quiet = opt$quiet)
    key <- if (isTRUE(m$status == "SYNONYM") && !is.null(m$acceptedUsageKey)) m$acceptedUsageKey else m$usageKey
    if (is.null(key) || isTRUE(m$matchType == "NONE")) { warning("No GBIF match for '", nm, "'"); next }
    taxon_keys <- c(taxon_keys, key)
    say("  '", nm, "' -> taxonKey ", key, " (", m$rank, " ", m$scientificName,
        if (isTRUE(m$status=="SYNONYM")) " [synonym -> accepted]" else "", ")")
  }
  if (!length(taxon_keys)) stop("No taxa resolved.")
} else stop("Provide --taxon or --taxon_key.")
taxon_keys <- unique(taxon_keys)

countries <- NULL
if (!is.null(opt$preset)) { if (is.null(PRESETS[[opt$preset]])) stop("Unknown --preset"); countries <- PRESETS[[opt$preset]]$countries }
if (!is.null(opt$countries)) countries <- toupper(splitcsv(opt$countries))
basis     <- if (tolower(opt$basis) == "any") NULL else splitcsv(opt$basis)
has_coord <- switch(opt$has_coordinate, any = NULL, "true" = TRUE, "false" = FALSE,
                    stop("--has_coordinate must be any|true|false"))

# ── Fetch ─────────────────────────────────────────────────────
gbif_search <- function(params, limit_total, pause = 0.3) {
  got <- 0L; offset <- 0L; page <- 300L; acc <- list()
  repeat {
    q <- c(params, list(limit = page, start = offset))
    res <- gbif_try(function() do.call(occ_search, q)); d <- res$data
    if (is.null(d) || !nrow(d)) break
    acc[[length(acc)+1]] <- d; got <- got + nrow(d); offset <- offset + page
    if (isTRUE(res$meta$endOfRecords) || offset >= 100000L) break
    if (limit_total > 0 && got >= limit_total) break
    Sys.sleep(pause)
  }
  if (length(acc)) bind_rows(acc) else NULL
}
say("Querying GBIF ...")
grid <- expand.grid(tk = taxon_keys, cc = if (is.null(countries)) NA else countries, stringsAsFactors = FALSE)
raw <- list()
for (i in seq_len(nrow(grid))) {
  p <- list(taxonKey = grid$tk[i], basisOfRecord = basis, hasCoordinate = has_coord,
            fields = "all")
  if (!is.na(grid$cc[i])) p$country <- grid$cc[i]
  d <- gbif_search(p, opt$limit_total)
  if (!is.null(d)) raw[[length(raw)+1]] <- d
}
if (!length(raw)) stop("No records returned.")
occ <- bind_rows(raw)
say("  raw records: ", nrow(occ))

# ── Normalise ─────────────────────────────────────────────────
gc <- function(df, name) if (name %in% names(df)) df[[name]] else NA
subsp_of <- function(sci, ie) {
  ie <- as.character(ie); if (length(ie)==1) ie <- rep(ie, length(sci))
  out <- ifelse(!is.na(ie) & nzchar(ie), ie, NA_character_); need <- is.na(out)
  if (any(need)) out[need] <- vapply(sci[need], function(s) {
    toks <- str_extract_all(s %||% "", "[A-Za-z]+")[[1]]
    ep <- toks[-1][toks[-1] == tolower(toks[-1])]; if (length(ep) >= 2) ep[2] else NA_character_
  }, character(1))
  out
}
gid <- gc(occ, "gbifID"); if (all(is.na(gid))) gid <- gc(occ, "key")
df <- tibble(
  catalog_number = as.character(gc(occ,"catalogNumber")),
  Collection     = as.character(gc(occ,"institutionCode")),
  collectionCode = as.character(gc(occ,"collectionCode")),
  species        = as.character(gc(occ,"species")),
  scientificName = as.character(gc(occ,"scientificName")),
  Country        = as.character(gc(occ,"country")),
  countryCode    = as.character(gc(occ,"countryCode")),
  stateProvince  = as.character(gc(occ,"stateProvince")),
  locality       = as.character(gc(occ,"locality")),
  preparations   = as.character(gc(occ,"preparations")),
  Lat            = suppressWarnings(as.numeric(gc(occ,"decimalLatitude"))),
  Lon            = suppressWarnings(as.numeric(gc(occ,"decimalLongitude"))),
  coordinateUncertaintyInMeters = suppressWarnings(as.numeric(gc(occ,"coordinateUncertaintyInMeters"))),
  date           = as.character(gc(occ,"eventDate")),
  year           = suppressWarnings(as.integer(gc(occ,"year"))),
  sex            = as.character(gc(occ,"sex")),
  basisOfRecord  = as.character(gc(occ,"basisOfRecord")),
  gbifID         = as.character(gid),
  occurrenceID   = as.character(gc(occ,"occurrenceID"))
)
df$subspecies <- subsp_of(df$scientificName, gc(occ,"infraspecificEpithet"))
df$taxon <- ifelse(!is.na(df$subspecies) & nzchar(df$subspecies), paste(df$species, df$subspecies),
                   ifelse(!is.na(df$species) & nzchar(df$species), df$species, "(unidentified)"))
df$coord_precision <- ifelse(!is.na(df$Lat), "exact", NA_character_)
df$ID <- ifelse(!is.na(df$Collection) & nzchar(df$Collection) & !is.na(df$catalog_number),
                paste0(df$Collection, "_", df$catalog_number), df$gbifID)
df <- df %>% mutate(.k = ifelse(!is.na(catalog_number) & nzchar(catalog_number),
                                paste(Collection, catalog_number), paste0("_g_", gbifID))) %>%
             distinct(.k, .keep_all = TRUE) %>% select(-.k)
say("  unique specimens: ", nrow(df))

# ── Keep-list: retain only listed taxa (synonyms resolved) ───
keep_taxa <- splitcsv(opt$keep_taxa)
if (!is.null(opt$keep_taxa_file)) {
  tf <- suppressWarnings(readLines(opt$keep_taxa_file, warn = FALSE))
  if (length(tf) && grepl(",", tf[1]) && any(grepl("taxon|species|name", tolower(tf[1])))) {
    d0 <- suppressWarnings(readr::read_csv(opt$keep_taxa_file, show_col_types = FALSE))
    col <- intersect(c("taxon","species","scientificname","name"), tolower(names(d0)))
    if (length(col)) tf <- as.character(d0[[names(d0)[tolower(names(d0)) == col[1]][1]]])
  }
  keep_taxa <- unique(c(keep_taxa, trimws(tf[nzchar(trimws(tf))])))
}
canon <- function(x) {
  t <- str_extract_all(tolower(x %||% ""), "[a-z]+")[[1]]
  if (length(t) >= 3) paste(t[1], t[2], t[3]) else if (length(t) >= 2) paste(t[1], t[2]) else if (length(t)) t[1] else ""
}
if (!is.null(keep_taxa) && length(keep_taxa)) {
  say("Resolving keep-list (", length(keep_taxa), " taxa) and their synonyms ...")
  keepset <- character(0)
  for (nm in keep_taxa) {
    m <- gbif_try(function() name_backbone(name = nm), quiet = opt$quiet)
    acc_key <- if (isTRUE(m$status == "SYNONYM") && !is.null(m$acceptedUsageKey)) m$acceptedUsageKey else m$usageKey
    acc_can <- m$canonicalName %||% m$scientificName
    keepset <- c(keepset, canon(nm), canon(acc_can %||% nm))
    if (!is.null(acc_key) && !is.na(acc_key)) {
      syn <- tryCatch(gbif_try(function() name_usage(key = acc_key, data = "synonyms", limit = 1000)$data, quiet = opt$quiet),
                      error = function(e) NULL)
      if (!is.null(syn) && nrow(syn)) {
        cn <- if ("canonicalName" %in% names(syn)) syn$canonicalName else syn$scientificName
        keepset <- c(keepset, vapply(cn, canon, character(1)))
      }
    }
  }
  keepset <- unique(keepset[nzchar(keepset)])
  rec_bin <- tolower(trimws(df$species))                 # accepted species binomial
  rec_sci <- vapply(df$scientificName, canon, character(1))
  rec_tax <- vapply(df$taxon, canon, character(1))
  hit <- rec_bin %in% keepset | rec_sci %in% keepset | rec_tax %in% keepset
  say("  keep-list: ", sum(hit), " of ", nrow(df), " records kept (",
      length(keepset), " accepted+synonym names).")
  df <- df[hit, ]
  if (!nrow(df)) stop("Keep-list removed every record — check the names in --keep_taxa.")
}

# ── Exclude eggs & per-institution catalogue rules ───────────
if (isTRUE(opt$exclude_eggs)) {
  egg <- grepl("egg", ifelse(is.na(df$collectionCode), "", df$collectionCode), ignore.case = TRUE) |
         grepl("egg", ifelse(is.na(df$preparations),  "", df$preparations),  ignore.case = TRUE)
  egg[is.na(egg)] <- FALSE
  if (any(egg)) say("  excluded ", sum(egg), " egg records (collection/preparation)")
  df <- df[!egg, ]
}
if (!is.null(opt$exclude) && nzchar(opt$exclude)) {
  for (rule in trimws(strsplit(opt$exclude, ";")[[1]])) {
    if (!nzchar(rule)) next
    pp <- strsplit(rule, ":", fixed = TRUE)[[1]]; if (length(pp) < 2) next
    code <- trimws(pp[1]); pat <- trimws(paste(pp[-1], collapse=":"))
    hit <- !is.na(df$Collection) & toupper(df$Collection) == toupper(code) &
           !is.na(df$catalog_number) & grepl(pat, df$catalog_number)
    hit[is.na(hit)] <- FALSE
    if (any(hit)) say("  exclude ", code, " /", pat, "/: dropped ", sum(hit))
    df <- df[!hit, ]
  }
}

# ── Material preference / filter ─────────────────────────────
# 'any' among the tokens disables filtering (keep every preparation) while the
# remaining tokens still set the preference order. This matters for groups like
# reptiles where specimens are 'whole animal'/'alcohol', not 'skin'.
prep_tokens <- tolower(splitcsv(opt$preparations))
filter_prep <- !("any" %in% prep_tokens)
kw <- setdiff(prep_tokens, "any")
if (!length(kw)) kw <- c("tissue","liquid nitrogen","skin")   # default preference order
prep_low <- tolower(ifelse(is.na(df$preparations), "", df$preparations))
df$.prep_pref <- Inf
for (j in seq_along(kw)) { hit <- grepl(kw[j], prep_low, fixed = TRUE) & is.infinite(df$.prep_pref); df$.prep_pref[hit] <- j }
has_info <- nzchar(prep_low)
if (filter_prep) {
  keep <- is.finite(df$.prep_pref) | (opt$keep_unknown_prep & !has_info)
  df <- df[keep, ]
}
# unknown / non-preferred material ranks after the listed materials
df$.prep_pref[is.infinite(df$.prep_pref)] <- length(kw) + 1L
say("  after material step: ", nrow(df), " specimens",
    if (!filter_prep) " (no material filter; 'any')" else "")

# ── Collection search order ──────────────────────────────────
# Collections are always searched in the listed order (primary ranking).
# By default only the listed collections are used (--restrict_collections);
# --collection_first (or --no_restrict_collections) also considers others,
# ranked after the listed ones.
coll <- splitcsv(opt$collections)
restrict_eff <- isTRUE(opt$restrict_collections) && !isTRUE(opt$collection_first)
if (!is.null(coll) && restrict_eff) {
  df <- df[toupper(df$Collection) %in% toupper(coll), ]
  say("  restricted to collections ", paste(coll, collapse = ","), ": ", nrow(df), " specimens")
} else if (!is.null(coll)) {
  say("  collection order ", paste(coll, collapse = " > "), " then others")
}
df$.coll_rank <- if (is.null(coll)) 1L else {
  m <- match(toupper(df$Collection), toupper(coll)); ifelse(is.na(m), length(coll) + 1L, m)
}

# ── Optional geocoding ───────────────────────────────────────
if (opt$geocode != "none" && requireNamespace("sf", quietly = TRUE) &&
    requireNamespace("rnaturalearth", quietly = TRUE)) {
  suppressPackageStartupMessages({ library(sf); library(rnaturalearth) })
  world <- ne_countries(scale = "medium", returnclass = "sf")
  wpt <- suppressWarnings(st_point_on_surface(world)); wxy <- st_coordinates(wpt)
  iso <- world$iso_a2_eh; iso[is.na(iso)|iso=="-99"] <- world$iso_a2[is.na(iso)|iso=="-99"]
  cpt <- setNames(lapply(seq_along(iso), function(i) c(wxy[i,2], wxy[i,1])), iso)
  spt <- list()
  if (opt$geocode == "state") {
    st <- tryCatch(ne_states(returnclass="sf"), error=function(e) NULL)
    if (!is.null(st)) { sp <- suppressWarnings(st_point_on_surface(st)); sxy <- st_coordinates(sp)
      for (i in seq_len(nrow(st))) spt[[paste(st$iso_a2[i], tolower(trimws(st$name[i])))]] <- c(sxy[i,2], sxy[i,1]) }
  }
  set.seed(42)
  for (i in which(is.na(df$Lat))) {
    cc <- df$countryCode[i]; s <- df$stateProvince[i]; lat <- NA; lon <- NA; pr <- "none"
    if (opt$geocode=="state" && !is.na(cc) && !is.na(s) && nzchar(s)) { h <- spt[[paste(cc, tolower(trimws(s)))]]; if (!is.null(h)) { lat<-h[1];lon<-h[2];pr<-"state" } }
    if (is.na(lat) && !is.na(cc) && !is.null(cpt[[cc]])) { h <- cpt[[cc]]; lat<-h[1];lon<-h[2];pr<-"country" }
    if (!is.na(lat)) { if (!opt$no_jitter) { a <- if (pr=="state") 0.15 else 0.8; lat<-lat+runif(1,-a,a); lon<-lon+runif(1,-a,a) }
      df$Lat[i]<-round(lat,4); df$Lon[i]<-round(lon,4); df$coord_precision[i]<-pr }
  }
} else if (opt$geocode != "none") warning("--geocode needs sf + rnaturalearth; skipped.")

df$coord_source <- ifelse(is.na(df$coord_precision), "none",
                   ifelse(df$coord_precision=="exact","GPS (recorded)",
                   ifelse(df$coord_precision=="state","province centroid",
                   ifelse(df$coord_precision=="country","country centroid", df$coord_precision))))
df$.prec_rank <- ifelse(is.na(df$coord_precision), 3L,
                        c(exact=0L, state=1L, country=2L)[df$coord_precision])
df$.date_rank <- ifelse(!is.na(df$date) & nzchar(df$date), 0L, 1L)
df$.unc <- ifelse(is.na(df$coordinateUncertaintyInMeters), Inf, df$coordinateUncertaintyInMeters)

# ── Select one (per_taxon) best specimen per taxon ───────────
rank_col <- switch(opt$rank_level, taxon = "taxon", species = "species", subspecies = "subspecies",
                   stop("--rank_level must be taxon|species|subspecies"))
# Collection order is always the primary key, then material, then coord/date.
ord <- c(".coll_rank", ".prep_pref", ".prec_rank", ".date_rank", ".unc")
sel <- df %>% group_by(across(all_of(rank_col))) %>%
  arrange(across(all_of(ord)), .by_group = TRUE) %>%
  slice_head(n = opt$per_taxon) %>% ungroup()
say("Selected ", nrow(sel), " specimen(s) across ",
    length(unique(sel[[rank_col]])), " ", opt$rank_level, "(s).")

# ── Must-include specimens from a list (bypass all filters) ──
apply_include_list <- function(final, path, quiet = FALSE, genus = NA) {
  say2 <- function(...) if (!quiet) message(...)
  final$forced_include <- FALSE
  if (is.null(path)) return(final)
  inc <- suppressWarnings(readr::read_csv(path, show_col_types = FALSE))
  fx <- function(cands) { h <- intersect(tolower(cands), tolower(names(inc)))
    if (length(h)) names(inc)[tolower(names(inc)) == h[1]] else NA_character_ }
  ci_i <- fx(c("institution","institutioncode","collection","inst","institution_code"))
  ci_c <- fx(c("catalogue","catalog","catalog_number","catalognumber","catalogue_number","cat","catno","catalog_no"))
  ci_t  <- fx(c("taxon","scientificname","scientific_name"))            # full name column (optional)
  ci_s  <- fx(c("species","sp","name"))                                 # species (epithet or binomial)
  ci_ss <- fx(c("subspecies","subsp","ssp","subsp.","race","infraspecificepithet"))  # optional
  if (is.na(ci_i) || is.na(ci_c)) stop("--include needs an institution column and a catalogue column")
  inst <- trimws(as.character(inc[[ci_i]])); cat_ <- trimws(as.character(inc[[ci_c]]))
  # Build a taxon name per row: prefer a full taxon column; else Species (+optional
  # subspecies), prefixing the run's genus when Species is a bare epithet.
  if (!is.na(ci_t)) {
    tax <- trimws(as.character(inc[[ci_t]]))
  } else if (!is.na(ci_s)) {
    spv <- trimws(as.character(inc[[ci_s]]))
    ssv <- if (!is.na(ci_ss)) trimws(as.character(inc[[ci_ss]])) else rep(NA_character_, nrow(inc))
    tax <- vapply(seq_along(spv), function(k) {
      s <- spv[k]; if (is.na(s) || !nzchar(s)) return(NA_character_)
      if (!grepl("\\s", s) && !is.na(genus) && nzchar(genus)) s <- paste(genus, s)  # bare epithet -> add genus
      ss <- ssv[k]
      if (!is.na(ss) && nzchar(ss)) s <- paste(s, tail(strsplit(trimws(ss), "\\s+")[[1]], 1))
      s
    }, character(1))
  } else {
    tax <- rep(NA_character_, nrow(inc))
  }
  keyf <- function(a, b) paste(toupper(trimws(ifelse(is.na(a),"",a))), toupper(trimws(ifelse(is.na(b),"",b))))
  inc_key <- keyf(inst, cat_); have <- keyf(final$Collection, final$catalog_number)
  final$forced_include <- have %in% inc_key
  need <- which(!(inc_key %in% have) & nzchar(ifelse(is.na(cat_), "", cat_)))
  say2("Include list: ", nrow(inc), " listed; ", length(need), " to add (",
       sum(inc_key %in% have), " already present).")
  gen1 <- function(x) { s <- tolower(trimws(x %||% "")); if (!nzchar(s)) return(""); strsplit(s, "\\s+")[[1]][1] }
  rows <- list()
  for (i in need) {
    want_gen <- if (!is.na(tax[i]) && nzchar(tax[i])) gen1(tax[i]) else NA_character_
    o <- tryCatch(gbif_try(function() occ_search(institutionCode = inst[i], catalogNumber = cat_[i], limit = 20, fields = "all")$data, quiet = quiet),
                  error = function(e) NULL)
    if (is.null(o) || !nrow(o))
      o <- tryCatch(gbif_try(function() occ_search(catalogNumber = cat_[i], limit = 50, fields = "all")$data, quiet = quiet),
                    error = function(e) NULL)
    o1 <- NULL
    if (!is.null(o) && nrow(o)) {
      icv <- toupper(trimws(ifelse(is.na(o$institutionCode), "", o$institutionCode))) == toupper(inst[i])
      ccv <- toupper(trimws(as.character(o$catalogNumber))) == toupper(cat_[i])
      cand <- which(icv & ccv); if (!length(cand)) cand <- which(ccv)
      # keep only candidates whose genus matches the row's taxon (guards against
      # catalogue numbers reused across a museum's other collections)
      if (length(cand) && !is.na(want_gen)) {
        gok <- vapply(cand, function(k)
          gen1(o$species[k]) == want_gen || gen1(o$scientificName[k]) == want_gen, logical(1))
        cand <- cand[gok]
      }
      if (length(cand)) o1 <- o[cand[1], , drop = FALSE]
    }
    if (!is.null(o1)) {
      gidv <- gc(o1, "gbifID"); if (all(is.na(gidv))) gidv <- gc(o1, "key")
      r <- tibble(
        catalog_number = as.character(gc(o1,"catalogNumber")), Collection = as.character(gc(o1,"institutionCode")),
        collectionCode = as.character(gc(o1,"collectionCode")), species = as.character(gc(o1,"species")),
        subspecies = subsp_of(as.character(gc(o1,"scientificName")), gc(o1,"infraspecificEpithet")),
        scientificName = as.character(gc(o1,"scientificName")), Country = as.character(gc(o1,"country")),
        countryCode = as.character(gc(o1,"countryCode")), stateProvince = as.character(gc(o1,"stateProvince")),
        locality = as.character(gc(o1,"locality")), preparations = as.character(gc(o1,"preparations")),
        Lat = suppressWarnings(as.numeric(gc(o1,"decimalLatitude"))), Lon = suppressWarnings(as.numeric(gc(o1,"decimalLongitude"))),
        coordinateUncertaintyInMeters = suppressWarnings(as.numeric(gc(o1,"coordinateUncertaintyInMeters"))),
        date = as.character(gc(o1,"eventDate")), year = suppressWarnings(as.integer(gc(o1,"year"))),
        sex = as.character(gc(o1,"sex")), basisOfRecord = as.character(gc(o1,"basisOfRecord")),
        gbifID = as.character(gidv), occurrenceID = as.character(gc(o1,"occurrenceID")))
      say2("  + ", inst[i], " ", cat_[i], " (from GBIF)")
    } else {
      sp <- NA_character_; ss <- NA_character_
      if (!is.na(tax[i]) && nzchar(tax[i])) { ss <- subsp_of(tax[i], NA)
        tk <- strsplit(tax[i], "\\s+")[[1]]; sp <- if (length(tk) >= 2) paste(tk[1], tk[2]) else tax[i] }
      r <- tibble(catalog_number = cat_[i], Collection = inst[i], collectionCode = NA_character_, species = sp,
        subspecies = ss, scientificName = ifelse(is.na(tax[i]), NA_character_, tax[i]), Country = NA_character_,
        countryCode = NA_character_, stateProvince = NA_character_, locality = NA_character_, preparations = NA_character_,
        Lat = NA_real_, Lon = NA_real_, coordinateUncertaintyInMeters = NA_real_, date = NA_character_,
        year = NA_integer_, sex = NA_character_, basisOfRecord = NA_character_, gbifID = NA_character_, occurrenceID = NA_character_)
      say2("  + ", inst[i], " ", cat_[i], " (no taxon-matching GBIF record; from list)")
    }
    rows[[length(rows) + 1]] <- r
  }
  if (length(rows)) {
    a <- bind_rows(rows)
    a$taxon <- ifelse(!is.na(a$subspecies) & nzchar(a$subspecies), paste(a$species, a$subspecies),
                ifelse(!is.na(a$species) & nzchar(a$species), a$species, "(unidentified)"))
    a$coord_precision <- ifelse(!is.na(a$Lat), "exact", NA_character_)
    a$coord_source <- ifelse(is.na(a$coord_precision), "none", "GPS (recorded)")
    a$ID <- ifelse(!is.na(a$Collection) & nzchar(a$Collection) & !is.na(a$catalog_number),
                   paste0(a$Collection, "_", a$catalog_number), a$gbifID)
    a$forced_include <- TRUE
    final <- bind_rows(final, a)
  }
  final
}
# genus context for --include (used to expand bare species epithets)
inc_genus <- {
  g <- vapply(splitcsv(opt$taxon) %||% character(0),
              function(s) strsplit(trimws(s), "\\s+")[[1]][1], character(1))
  if (length(g) && length(unique(tolower(g))) == 1) g[1] else NA_character_
}
sel <- apply_include_list(sel, opt$include, quiet = opt$quiet, genus = inc_genus)

# ── Restrict everything (incl. --include) to the searched taxon ──
if (isTRUE(opt$restrict_taxon)) {
  gsrc <- c(splitcsv(opt$taxon), keep_taxa)
  allowed_gen <- unique(tolower(vapply(gsrc, function(s) {
    t <- strsplit(trimws(s %||% ""), "\\s+")[[1]]; if (length(t)) t[1] else "" }, character(1))))
  allowed_gen <- allowed_gen[nzchar(allowed_gen)]
  if (length(allowed_gen)) {
    g1 <- tolower(vapply(strsplit(trimws(ifelse(is.na(sel$species), "", sel$species)), "\\s+"),
                         function(t) if (length(t)) t[1] else "", character(1)))
    g2 <- tolower(vapply(strsplit(trimws(ifelse(is.na(sel$scientificName), "", sel$scientificName)), "\\s+"),
                         function(t) if (length(t)) t[1] else "", character(1)))
    keeprow <- g1 %in% allowed_gen | g2 %in% allowed_gen
    if (any(!keeprow))
      say("restrict_taxon: dropped ", sum(!keeprow), " off-target record(s): ",
          paste(head(sel$ID[!keeprow], 10), collapse = ", "),
          if (sum(!keeprow) > 10) ", ..." else "")
    sel <- sel[keeprow, ]
  }
}

# ── Assemble the final panel when an --include list is given ─
# Default: the include specimens, PLUS one auto-picked specimen for each
# species in --keep_taxa not already covered by an include, filled from the
# priority --collections in order, then by preparation. --include_only keeps
# only the includes; --include_add_all restores the full genus selection.
if (!is.null(opt$include) && !isTRUE(opt$include_add_all)) {
  binname <- function(x) { t <- str_extract_all(tolower(x %||% ""), "[a-z]+")[[1]]
    if (length(t) >= 2) paste(t[1], t[2]) else if (length(t)) t[1] else "" }
  inc_rows <- sel[sel$forced_include %in% TRUE, ]
  inc_sp   <- unique(vapply(inc_rows$species, binname, character(1)))
  gap_rows <- sel[0, ]
  if (isTRUE(opt$include_only) || is.null(keep_taxa) || !length(keep_taxa)) {
    say("include: keeping only the ", nrow(inc_rows), " include specimen(s)",
        if (is.null(keep_taxa) || !length(keep_taxa)) " (no --keep_taxa gap-fill)." else " (--include_only).")
  } else {
    target_sp <- unique(vapply(keep_taxa, binname, character(1))); target_sp <- target_sp[nzchar(target_sp)]
    gap_sp <- setdiff(target_sp, inc_sp)
    if (length(gap_sp)) {
      pool <- df
      pool$.spbin <- vapply(pool$species, binname, character(1))
      pool <- pool[pool$.spbin %in% gap_sp, ]
      if (nrow(pool)) {
        gap_ord <- c(".coll_rank", ".prep_pref", ".prec_rank", ".date_rank", ".unc")  # priority, then preparation
        gap_rows <- pool %>% group_by(.spbin) %>%
          arrange(across(all_of(gap_ord)), .by_group = TRUE) %>%
          slice_head(n = 1) %>% ungroup() %>% select(-.spbin)
        gap_rows$forced_include <- FALSE
      }
    }
    filled <- unique(vapply(gap_rows$species, binname, character(1)))
    say("include: ", nrow(inc_rows), " from list + ", nrow(gap_rows), " gap-filled species (",
        length(setdiff(gap_sp, filled)), " of ", length(target_sp),
        " target species had no specimen available).")
  }
  sel <- bind_rows(inc_rows, gap_rows)
}

# ── Write CSV ─────────────────────────────────────────────────
cols <- c("ID","catalog_number","Collection","collectionCode","species","subspecies","taxon",
          "scientificName","Country","countryCode","stateProvince","locality","preparations",
          "Lat","Lon","coordinateUncertaintyInMeters","coord_precision","coord_source",
          "date","year","sex","basisOfRecord","gbifID","occurrenceID","forced_include")
write_csv(sel[, cols], opt$output)
message("Wrote ", opt$output, ": ", nrow(sel), " specimens (", sum(!is.na(sel$Lat)), " with coordinates).")

# ── Summary ───────────────────────────────────────────────────
mk <- function(v, categ) { t <- sort(table(v), decreasing = TRUE)
  data.frame(category = categ, value = names(t), n = as.integer(t), stringsAsFactors = FALSE) }
prep_match <- function(p) { pl <- tolower(ifelse(is.na(p), "", p))
  for (j in seq_along(kw)) if (grepl(kw[j], pl, fixed = TRUE)) return(kw[j]); "unknown/none" }
prep_types <- vapply(sel$preparations, prep_match, character(1))
S <- rbind(mk(ifelse(is.na(sel$Collection)|!nzchar(sel$Collection),"(none)",sel$Collection), "institution"),
           mk(sel$taxon, "taxon"),
           mk(prep_types, "material"),
           mk(sel$coord_source, "coord_source"))
for (ct in c("institution","taxon","material","coord_source")) {
  b <- S[S$category==ct, c("value","n")]
  cat(sprintf("\n== per %s (n=%d) ==\n", ct, sum(b$n))); print(b, row.names = FALSE)
}
summary_path <- if (!is.null(opt$summary)) opt$summary else paste0(sub("\\.[^.]*$","",opt$output), "_summary.tsv")
write_tsv(S, summary_path); message("\nWrote ", summary_path, ": summary counts (", nrow(S), " rows).")

# ── Optional: separate country-centroid CSV ──────────────────
# Same rows as --output, but specimens with no GPS get their country
# centroid so the whole set can be mapped. The main output is untouched.
if (isTRUE(opt$country_centroid)) {
  fill_country_centroids <- function(d) {
    if (!requireNamespace("sf", quietly = TRUE) || !requireNamespace("rnaturalearth", quietly = TRUE)) {
      warning("--country_centroid needs sf + rnaturalearth; skipped."); return(NULL) }
    suppressPackageStartupMessages({ library(sf); library(rnaturalearth) })
    w  <- ne_countries(scale = "medium", returnclass = "sf")
    wp <- suppressWarnings(st_point_on_surface(w)); xy <- st_coordinates(wp)
    iso <- w$iso_a2_eh; iso[is.na(iso) | iso == "-99"] <- w$iso_a2[is.na(iso) | iso == "-99"]
    by_iso <- setNames(seq_len(nrow(w)), toupper(iso))
    nmA <- tolower(trimws(w$admin)); nmL <- tolower(trimws(w$name_long)); nmN <- tolower(trimws(w$name))
    miss <- which(is.na(d$Lat) | is.na(d$Lon)); filled <- 0L
    for (i in miss) {
      idx <- NA_integer_; cc <- d$countryCode[i]
      if (!is.na(cc) && nzchar(cc)) { k <- toupper(trimws(cc)); if (k %in% names(by_iso)) idx <- by_iso[[k]] }
      if ((is.null(idx) || is.na(idx)) && !is.na(d$Country[i])) {
        cn <- tolower(trimws(d$Country[i])); hit <- which(nmA == cn | nmL == cn | nmN == cn)
        if (length(hit)) idx <- hit[1]
      }
      if (!is.null(idx) && !is.na(idx)) {
        d$Lat[i] <- round(xy[idx, 2], 4); d$Lon[i] <- round(xy[idx, 1], 4)
        d$coord_precision[i] <- "country"; d$coord_source[i] <- "country centroid"; filled <- filled + 1L
      }
    }
    attr(d, "filled") <- filled; d
  }
  dc <- fill_country_centroids(sel)
  if (!is.null(dc)) {
    cpath <- if (!is.null(opt$centroid_output)) opt$centroid_output
             else paste0(sub("\\.[^.]*$", "", opt$output), "_centroid.csv")
    write_csv(dc[, cols], cpath)
    message("Wrote ", cpath, ": ", nrow(dc), " specimens with country centroids filled (",
            attr(dc, "filled"), " added; ", sum(!is.na(dc$Lat)), " now mappable).")
  }
}

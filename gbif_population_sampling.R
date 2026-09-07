#!/usr/bin/env Rscript
# ============================================================
# GBIF specimen search  ->  sample sheet for specimen_map.R
#
# Queries the GBIF occurrence API (via rgbif) for specimen records
# of one or more taxa, optionally restricted to a set of countries
# and/or holding institutions, and writes a tidy CSV whose column
# names are auto-detected by specimen_map.R (ID, Lat, Lon, species,
# subspecies, Collection ...).
#
# Records with no coordinates can optionally be geocoded to the
# province or country centroid (Natural Earth) so they still plot;
# every row is tagged in `coord_precision` (exact|state|country).
#
# Packages: optparse, rgbif, dplyr, readr, stringr
#   --geocode / --range_shapefile additionally need: sf, rnaturalearth
#   (rnaturalearthhires improves province coverage; falls back to
#    country centroids if absent)
#
# Examples
# --------
#   # Ostrich preserved specimens across the native range, coords only:
#   Rscript gbif_population_sampling.R \
#     --taxon "Struthio camelus,Struthio molybdophanes" \
#     --preset ostrich-native \
#     --output ostrich_gbif.csv
#
#   # Also place coordinate-less records at province/country centroids
#   # and emit a dissolved range shapefile for the map extent:
#   Rscript gbif_population_sampling.R \
#     --taxon Struthio --preset ostrich-native \
#     --geocode state \
#     --priority "FMNH,AMNH,NHMUK,MVZ,MSU,MNHN" \
#     --range_shapefile range/ostrich_native_range.shp \
#     --output ostrich_gbif.csv
#
#   # Any taxon / countries by ISO-2 code:
#   Rscript gbif_population_sampling.R --taxon "Panthera leo" \
#     --countries KE,TZ,ZA --output lion.csv
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(rgbif)
  library(dplyr)
  library(readr)
  library(stringr)
})

# Native (extant) range of Struthio, as ISO 3166-1 alpha-2 codes.
OSTRICH_NATIVE <- c("ZA","KE","BW","TZ","SD","EG","NA","CD","SL","DZ","SO","ET",
                    "ZW","UG","ML","LY","GH","TD","SN","NG","NE","MA",   # Africa
                    "YE","SA","SY","IL","JO","OM","IQ")                  # Middle East
PRESETS <- list("ostrich-native" = list(countries = OSTRICH_NATIVE))

# ── Flags ─────────────────────────────────────────────────────
option_list <- list(
  make_option(c("-t", "--taxon"), type = "character", default = NULL,
              help = "Comma-separated scientific names to search (resolved via GBIF backbone)."),
  make_option("--taxon_key", type = "character", default = NULL,
              help = "Comma-separated GBIF taxonKeys (skips name resolution)."),
  make_option(c("-c", "--countries"), type = "character", default = NULL,
              help = "Comma-separated ISO-2 country codes to restrict to."),
  make_option(c("-e", "--preset"), type = "character", default = NULL,
              help = "Named country preset (e.g. 'ostrich-native')."),
  make_option(c("-b", "--basis"), type = "character", default = "PRESERVED_SPECIMEN",
              help = "Comma-separated basisOfRecord ['any' for all] [default: PRESERVED_SPECIMEN]."),
  make_option(c("-i", "--institutions"), type = "character", default = NULL,
              help = paste("Limit output to these institutionCodes (comma-separated,",
                           "case-insensitive; multiple allowed), searched/ranked in the",
                           "listed order. Regions come from --countries/--preset.")),
  make_option(c("-f", "--collection_first"), action = "store_true", default = FALSE,
              help = "Prefer --institutions in order but also keep other collections (ranked after the listed ones in --minimize) instead of restricting to the list."),
  make_option(c("-p", "--preparations"), type = "character", default = "tissue,liquid nitrogen,skin",
              help = paste("Keep only records whose preparations match these keywords",
                           "(comma-separated, case-insensitive substring, so all variations",
                           "like 'study skin', 'Mounted skin', 'frozen tissue' match), in",
                           "preference order for --minimize",
                           "[default: 'tissue,liquid nitrogen,skin']. 'any' keeps all types.")),
  make_option("--keep_unknown_prep", action = "store_true", default = TRUE,
              help = "Keep records with no preparation recorded (the 'none' bucket) [default: on]."),
  make_option("--drop_unknown_prep", action = "store_false", dest = "keep_unknown_prep",
              help = "Drop records that have no preparation recorded."),
  make_option(c("-x", "--exclude"), type = "character", default = "",
              help = paste("Institution-specific exclusion rules as 'CODE:regex' on the",
                           "catalogue number, separated by ';'. Drops matching records.",
                           "Off by default. e.g. 'NHMUK:E;AMNH:^EGG' (NHMUK eggs carry an E).")),
  make_option("--priority", type = "character", default = NULL,
              help = "Comma-separated institutionCodes to flag in an is_priority column."),
  make_option("--has_coordinate", type = "character", default = "any",
              help = "Restrict to records with/without coordinates: any|true|false [default: any]."),
  make_option("--no_geospatial_issue", action = "store_true", default = FALSE,
              help = "Drop records GBIF flags with a geospatial issue."),
  make_option("--limit_total", type = "integer", default = 0L,
              help = "Cap total records fetched per taxon x country (0 = no cap)."),
  make_option(c("-g", "--geocode"), type = "character", default = "none",
              help = "Place coordinate-less records at centroid: none|country|state [default: none]."),
  make_option("--no_jitter", action = "store_true", default = FALSE,
              help = "Do not jitter centroid-placed points."),
  make_option(c("-C", "--country_centroid"), action = "store_true", default = FALSE,
              help = paste("Write a SEPARATE CSV in which specimens with no GPS are filled",
                           "with their country centroid (main output keeps GPS only).")),
  make_option("--centroid_output", type = "character", default = NULL,
              help = "Path for the country-centroid CSV [default: <output>_centroid.csv]."),
  make_option(c("-r", "--range_shapefile"), type = "character", default = NULL,
              help = "Write a dissolved shapefile of the countries in the result (map extent)."),
  # --- minimisation: smallest panel that still covers the diversity ---
  make_option(c("-m", "--minimize"), action = "store_true", default = FALSE,
              help = "Reduce to a minimal representative panel (one best specimen per stratum)."),
  make_option("--strata", type = "character", default = "taxon,countryCode,stateProvince,locality",
              help = paste("Comma-separated columns defining a coverage cell.",
                           "Default keeps one specimen per unique locality:",
                           "taxon,countryCode,stateProvince,locality. Coarsen to",
                           "'taxon,countryCode' for a smaller panel.")),
  make_option("--per_stratum", type = "integer", default = 1L,
              help = "How many specimens to keep per stratum [default: 1]."),
  make_option("--collection_priority", type = "character", default = "",
              help = "Ranked institutionCodes preferred when choosing a stratum's representative (searched in this order). Defaults to the --institutions order."),
  make_option("--require_coords", action = "store_true", default = FALSE,
              help = "When minimizing, only consider records that have (or were geocoded to) coordinates."),
  make_option("--output_full", type = "character", default = NULL,
              help = "When minimizing, also write the complete (un-minimized) set here."),
  make_option(c("-I", "--include"), type = "character", default = NULL,
              help = paste("CSV of must-include specimens (columns: institution, taxon,",
                           "catalogue). After the search, any not already in the final",
                           "sheet are looked up on GBIF to fill all columns and appended",
                           "(or added from the given fields if GBIF has no record).",
                           "These bypass all filters.")),
  make_option(c("-u", "--summary"), type = "character", default = NULL,
              help = "Summary-counts TSV path [default: <output>_summary.tsv]."),
  make_option(c("-o", "--output"), type = "character", default = "gbif_samples.csv",
              help = "Output CSV path [default: gbif_samples.csv]."),
  make_option("--quiet", action = "store_true", default = FALSE,
              help = "Suppress progress messages.")
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Search GBIF for specimens -> specimen_map.R sample sheet."))
say <- function(...) if (!opt$quiet) message(...)

splitcsv <- function(x) if (is.null(x)) NULL else trimws(strsplit(x, ",")[[1]])

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

# ── Resolve taxon keys ────────────────────────────────────────
if (!is.null(opt$taxon_key)) {
  taxon_keys <- as.integer(splitcsv(opt$taxon_key))
} else if (!is.null(opt$taxon)) {
  say("Resolving taxa ...")
  taxon_keys <- integer(0)
  for (nm in splitcsv(opt$taxon)) {
    m <- gbif_try(function() name_backbone(name = nm), quiet = opt$quiet)
    if (is.null(m$usageKey) || isTRUE(m$matchType == "NONE")) {
      warning("No GBIF match for '", nm, "' — skipped"); next
    }
    taxon_keys <- c(taxon_keys, m$usageKey)
    say("  '", nm, "' -> taxonKey ", m$usageKey, " (", m$rank, " ", m$scientificName, ")")
  }
  if (!length(taxon_keys)) stop("No taxa resolved; nothing to search.")
} else {
  stop("Provide --taxon or --taxon_key.")
}

# ── Assemble country list ─────────────────────────────────────
countries <- NULL
if (!is.null(opt$preset)) {
  if (is.null(PRESETS[[opt$preset]])) stop("Unknown --preset '", opt$preset, "'")
  countries <- PRESETS[[opt$preset]]$countries
}
if (!is.null(opt$countries)) countries <- toupper(splitcsv(opt$countries))

basis     <- if (tolower(opt$basis) == "any") NULL else splitcsv(opt$basis)
insts     <- splitcsv(opt$institutions)
priority  <- splitcsv(opt$priority)
has_coord <- switch(opt$has_coordinate, any = NULL, "true" = TRUE, "false" = FALSE,
                    stop("--has_coordinate must be any|true|false"))

# ── Fetch (loop taxon x country, page via start offset) ───────
say("Querying GBIF ...")
grid <- expand.grid(tk = taxon_keys,
                    cc = if (is.null(countries)) NA else countries,
                    stringsAsFactors = FALSE)
raw <- list()
for (i in seq_len(nrow(grid))) {
  tk <- grid$tk[i]; cc <- grid$cc[i]
  start <- 0L; page <- 300L
  repeat {
    res <- gbif_try(function() occ_search(
      taxonKey          = tk,
      country           = if (is.na(cc)) NULL else cc,
      basisOfRecord     = basis,
      hasCoordinate     = has_coord,
      hasGeospatialIssue = if (opt$no_geospatial_issue) FALSE else NULL,
      limit = page, start = start, fields = "all"
    ), quiet = opt$quiet)
    d <- res$data
    if (is.null(d) || !nrow(d)) break
    raw[[length(raw) + 1]] <- d
    got_here <- sum(vapply(raw, nrow, 0L))            # rough; fine for small taxa
    start <- start + page
    if (isTRUE(res$meta$endOfRecords) || start >= 100000L) break
    if (opt$limit_total > 0 && nrow(d) < page) break
  }
}
if (!length(raw)) stop("No records returned for the given query.")
occ <- bind_rows(raw)
say("  raw records: ", nrow(occ))

# ── Normalise columns (tolerate missing fields) ──────────────
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
gc <- function(df, name) if (name %in% names(df)) df[[name]] else NA
subsp_of <- function(sci, ie) {
  ie <- as.character(ie)
  if (length(ie) == 1) ie <- rep(ie, length(sci))   # recycle a missing field to full length
  out <- ifelse(!is.na(ie) & nzchar(ie), ie, NA_character_)
  need <- is.na(out)
  if (any(need)) {
    parsed <- vapply(sci[need], function(s) {
      toks <- str_extract_all(s %||% "", "[A-Za-z]+")[[1]]
      ep <- toks[-1][toks[-1] == tolower(toks[-1])]   # lowercase epithets after genus
      if (length(ep) >= 2) ep[2] else NA_character_
    }, character(1))
    out[need] <- parsed
  }
  out
}

gid <- gc(occ, "gbifID"); if (all(is.na(gid))) gid <- gc(occ, "key")

df <- tibble(
  catalog_number = as.character(gc(occ, "catalogNumber")),
  Collection     = as.character(gc(occ, "institutionCode")),
  collectionCode = as.character(gc(occ, "collectionCode")),
  species        = as.character(gc(occ, "species")),
  scientificName = as.character(gc(occ, "scientificName")),
  Country        = as.character(gc(occ, "country")),
  countryCode    = as.character(gc(occ, "countryCode")),
  stateProvince  = as.character(gc(occ, "stateProvince")),
  locality       = as.character(gc(occ, "locality")),
  preparations   = as.character(gc(occ, "preparations")),
  Lat            = suppressWarnings(as.numeric(gc(occ, "decimalLatitude"))),
  Lon            = suppressWarnings(as.numeric(gc(occ, "decimalLongitude"))),
  coordinateUncertaintyInMeters = suppressWarnings(as.numeric(gc(occ, "coordinateUncertaintyInMeters"))),
  date           = as.character(gc(occ, "eventDate")),
  year           = suppressWarnings(as.integer(gc(occ, "year"))),
  sex            = as.character(gc(occ, "sex")),
  basisOfRecord  = as.character(gc(occ, "basisOfRecord")),
  gbifID         = as.character(gid),
  occurrenceID   = as.character(gc(occ, "occurrenceID"))
)
df$subspecies <- subsp_of(df$scientificName, gc(occ, "infraspecificEpithet"))
# combined taxon label: species + subspecies (so a single colour scale can
# distinguish both species and subspecies, with no NA/grey category)
df$taxon <- ifelse(!is.na(df$subspecies) & nzchar(df$subspecies),
                   paste(df$species, df$subspecies),
                   ifelse(!is.na(df$species) & nzchar(df$species), df$species, "(unidentified)"))
df$coord_precision <- ifelse(!is.na(df$Lat), "exact", NA_character_)
df$ID <- ifelse(!is.na(df$Collection) & nzchar(df$Collection) & !is.na(df$catalog_number),
                paste0(df$Collection, "_", df$catalog_number), df$gbifID)
df$is_priority <- if (!is.null(priority)) df$Collection %in% priority else NA

# dedupe by specimen identity (institution + catalog no.), else gbifID
df <- df %>%
  mutate(.key = ifelse(!is.na(catalog_number) & nzchar(catalog_number),
                       paste(Collection, catalog_number), paste0("_gbif_", gbifID))) %>%
  distinct(.key, .keep_all = TRUE) %>% select(-.key)
say("  unique specimens: ", nrow(df))

# ── Institution-specific exclusion rules (e.g. NHMUK eggs) ───
# Each rule is CODE:regex applied to the catalogue number; records from that
# institution whose catalogue number matches are dropped. Handles quirks like
# NHMUK egg registrations carrying an 'E' in the number.
if (!is.null(opt$exclude) && nzchar(opt$exclude)) {
  for (rule in trimws(strsplit(opt$exclude, ";")[[1]])) {
    if (!nzchar(rule)) next
    parts <- strsplit(rule, ":", fixed = TRUE)[[1]]
    if (length(parts) < 2) { warning("Ignoring malformed --exclude rule: '", rule, "'"); next }
    code <- trimws(parts[1]); pat <- trimws(paste(parts[-1], collapse = ":"))
    hit <- !is.na(df$Collection) & toupper(df$Collection) == toupper(code) &
           !is.na(df$catalog_number) & grepl(pat, df$catalog_number)
    hit[is.na(hit)] <- FALSE
    if (any(hit)) say("  exclude ", code, " /", pat, "/: dropped ", sum(hit), " records")
    df <- df[!hit, ]
  }
}

# ── Limit to requested museum collections (client-side) ──────
# By default the sheet is restricted to --institutions (searched/ranked in the
# listed order in --minimize). With --collection_first the listed institutions
# are preferred first but others are also kept (ranked after them in --minimize).
if (!is.null(insts) && !isTRUE(opt$collection_first)) {
  keep <- toupper(df$Collection) %in% toupper(insts)
  say("  collection filter (", paste(insts, collapse = ", "), "): ",
      sum(keep), " of ", nrow(df), " kept")
  df <- df[keep, ]
} else if (!is.null(insts)) {
  say("  collection order ", paste(insts, collapse = " > "), " then others (--collection_first)")
}

# ── Limit by preparation (default: tissue, then skin) ────────
# Preference order = the --preparations list; used both to filter and to
# rank in --minimize. Matching is case-insensitive substring, so every
# variation ("study skin", "skin; skeleton", "frozen tissue", …) matches.
# 'any' among the tokens keeps every preparation (no filter) while the other
# tokens still set the preference order — needed for reptiles etc. where
# specimens are 'whole animal'/'alcohol', not 'skin'.
prep_tokens <- tolower(splitcsv(opt$preparations))
filter_prep <- !("any" %in% prep_tokens)
kw <- setdiff(prep_tokens, "any"); if (!length(kw)) kw <- c("tissue","liquid nitrogen","skin")
prep_low <- tolower(ifelse(is.na(df$preparations), "", df$preparations))
df$.prep_pref <- Inf
for (j in seq_along(kw)) {
  hit <- grepl(kw[j], prep_low, fixed = TRUE) & is.infinite(df$.prep_pref)
  df$.prep_pref[hit] <- j                          # first (most preferred) match wins
}
has_info <- nzchar(prep_low)
if (filter_prep) {
  keep <- is.finite(df$.prep_pref) | (opt$keep_unknown_prep & !has_info)
  say("  preparation filter (", paste(kw, collapse = " > "), "): ", sum(keep),
      " of ", nrow(df), " kept",
      if (opt$keep_unknown_prep) paste0(" (incl. ", sum(!has_info & keep), " unknown-prep)") else "")
  df <- df[keep, ]
}

# ── Optional geocoding of coordinate-less rows ───────────────
if (opt$geocode != "none") {
  ok <- requireNamespace("sf", quietly = TRUE) &&
        requireNamespace("rnaturalearth", quietly = TRUE)
  if (!ok) {
    warning("--geocode needs sf + rnaturalearth; leaving coordinate-less rows blank.")
  } else {
    suppressPackageStartupMessages({ library(sf); library(rnaturalearth) })
    world <- ne_countries(scale = "medium", returnclass = "sf")
    wpt   <- suppressWarnings(st_point_on_surface(world))
    wxy   <- st_coordinates(wpt)
    iso   <- world$iso_a2_eh; iso[is.na(iso) | iso == "-99"] <- world$iso_a2[is.na(iso) | iso == "-99"]
    cpt   <- setNames(lapply(seq_along(iso), function(i) c(wxy[i, 2], wxy[i, 1])), iso)

    spt <- list()
    if (opt$geocode == "state") {
      states <- tryCatch(ne_states(returnclass = "sf"), error = function(e) NULL)
      if (!is.null(states)) {
        spt_pt <- suppressWarnings(st_point_on_surface(states))
        sxy <- st_coordinates(spt_pt)
        for (i in seq_len(nrow(states))) {
          key <- paste(states$iso_a2[i], tolower(trimws(states$name[i])))
          spt[[key]] <- c(sxy[i, 2], sxy[i, 1])
        }
      } else warning("Province layer unavailable (install rnaturalearthhires); using country centroids.")
    }

    set.seed(42)
    need <- which(is.na(df$Lat))
    filled <- 0L
    for (i in need) {
      cc <- df$countryCode[i]; stt <- df$stateProvince[i]; lat <- NA; lon <- NA; prec <- "none"
      if (opt$geocode == "state" && !is.na(cc) && !is.na(stt) && nzchar(stt)) {
        hit <- spt[[paste(cc, tolower(trimws(stt)))]]
        if (!is.null(hit)) { lat <- hit[1]; lon <- hit[2]; prec <- "state" }
      }
      if (is.na(lat) && !is.na(cc) && !is.null(cpt[[cc]])) {
        hit <- cpt[[cc]]; lat <- hit[1]; lon <- hit[2]; prec <- "country"
      }
      if (!is.na(lat)) {
        if (!opt$no_jitter && prec %in% c("state", "country")) {
          amt <- if (prec == "state") 0.15 else 0.8
          lat <- lat + runif(1, -amt, amt); lon <- lon + runif(1, -amt, amt)
        }
        df$Lat[i] <- round(lat, 4); df$Lon[i] <- round(lon, 4); df$coord_precision[i] <- prec
        filled <- filled + 1L
      }
    }
    say("  geocoded ", filled, " records; precision mix: ",
        paste(names(table(df$coord_precision)), table(df$coord_precision),
              sep = "=", collapse = ", "))
  }
}

# Human-readable coordinate source: real GPS vs a geocoded centroid.
df$coord_source <- ifelse(is.na(df$coord_precision), "none",
                   ifelse(df$coord_precision == "exact",   "GPS (recorded)",
                   ifelse(df$coord_precision == "state",   "province centroid",
                   ifelse(df$coord_precision == "country", "country centroid",
                          df$coord_precision))))

df_full <- df   # keep the complete set for the range shapefile / --output_full
df_full$forced_include <- FALSE

# ── Must-include specimens from a list (bypass all filters) ──
# Reads a CSV with institution/taxon/catalogue columns; for each listed
# specimen not already in the final sheet, looks it up on GBIF to fill every
# column, or (if GBIF has no record) adds it from the provided fields.
apply_include_list <- function(final, path, priority = NULL, quiet = FALSE) {
  say2 <- function(...) if (!quiet) message(...)
  final$forced_include <- FALSE
  if (is.null(path)) return(final)
  inc <- suppressWarnings(readr::read_csv(path, show_col_types = FALSE))
  fx <- function(cands) { h <- intersect(tolower(cands), tolower(names(inc)))
    if (length(h)) names(inc)[tolower(names(inc)) == h[1]] else NA_character_ }
  ci_i <- fx(c("institution","institutioncode","collection","inst","institution_code"))
  ci_c <- fx(c("catalogue","catalog","catalog_number","catalognumber","catalogue_number","cat","catno","catalog_no"))
  ci_t <- fx(c("taxon","scientificname","scientific_name","species","name"))
  if (is.na(ci_i) || is.na(ci_c)) stop("--include needs an institution column and a catalogue column")
  inst <- trimws(as.character(inc[[ci_i]])); cat_ <- trimws(as.character(inc[[ci_c]]))
  tax  <- if (!is.na(ci_t)) trimws(as.character(inc[[ci_t]])) else rep(NA_character_, nrow(inc))
  keyf <- function(a, b) paste(toupper(trimws(ifelse(is.na(a),"",a))), toupper(trimws(ifelse(is.na(b),"",b))))
  inc_key <- keyf(inst, cat_); have <- keyf(final$Collection, final$catalog_number)
  final$forced_include <- have %in% inc_key
  need <- which(!(inc_key %in% have) & nzchar(ifelse(is.na(cat_), "", cat_)))
  say2("Include list: ", nrow(inc), " listed; ", length(need), " to add (",
       sum(inc_key %in% have), " already present).")
  rows <- list()
  for (i in need) {
    o <- tryCatch(gbif_try(function() occ_search(institutionCode = inst[i], catalogNumber = cat_[i], limit = 20, fields = "all")$data, quiet = quiet),
                  error = function(e) NULL)
    if (is.null(o) || !nrow(o))
      o <- tryCatch(gbif_try(function() occ_search(catalogNumber = cat_[i], limit = 50, fields = "all")$data, quiet = quiet),
                    error = function(e) NULL)
    if (!is.null(o) && nrow(o)) {
      icv <- toupper(trimws(ifelse(is.na(o$institutionCode), "", o$institutionCode))) == toupper(inst[i])
      ccv <- toupper(trimws(as.character(o$catalogNumber))) == toupper(cat_[i])
      pick <- which(icv & ccv); if (!length(pick)) pick <- which(ccv); if (!length(pick)) pick <- 1L
      o1 <- o[pick[1], , drop = FALSE]
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
      say2("  + ", inst[i], " ", cat_[i], " (not found on GBIF; from list)")
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
    if (!is.null(priority)) a$is_priority <- a$Collection %in% priority
    a$forced_include <- TRUE
    final <- bind_rows(final, a)
  }
  final
}

# ── Optional minimisation: smallest representative panel ──────
# One (or --per_stratum) best specimen per coverage cell. A cell is the
# combination of --strata columns (default subspecies x country). Within a
# cell, specimens are ranked by: (1) collection priority, (2) coordinate
# quality (exact > state > country > none), (3) has a date, (4) smaller
# coordinate uncertainty. This keeps the fewest samples that still span
# every subspecies-in-country present in the data.
if (opt$minimize) {
  strata <- splitcsv(opt$strata)
  bad <- setdiff(strata, names(df))
  if (length(bad)) stop("--strata column(s) not found: ", paste(bad, collapse = ", "))
  prio <- splitcsv(opt$collection_priority)
  if (is.null(prio) || !length(prio)) prio <- insts   # default to the --institutions order

  cand <- df
  if (opt$require_coords) cand <- cand[!is.na(cand$Lat), ]

  prec_rank <- c(exact = 0, state = 1, country = 2)
  cand <- cand %>%
    mutate(
      .coll_rank = { m <- match(Collection, prio); ifelse(is.na(m), length(prio) + 1L, m) },
      .prec_rank = ifelse(is.na(coord_precision), 3L, prec_rank[coord_precision]),
      .date_rank = ifelse(!is.na(date) & nzchar(date), 0L, 1L),
      .unc       = ifelse(is.na(coordinateUncertaintyInMeters), Inf, coordinateUncertaintyInMeters)
    ) %>%
    group_by(across(all_of(strata))) %>%
    arrange(.prep_pref, .coll_rank, .prec_rank, .date_rank, .unc, .by_group = TRUE) %>%
    slice_head(n = opt$per_stratum) %>%
    ungroup() %>%
    select(-.coll_rank, -.prec_rank, -.date_rank, -.unc)

  n_strata <- nrow(unique(cand[, strata, drop = FALSE]))
  say("  minimized: ", nrow(df_full), " -> ", nrow(cand),
      " specimens across ", n_strata, " strata (", paste(strata, collapse = " x "), ")")
  df <- cand
}

# apply the must-include list to the final sheet (after minimisation)
df <- apply_include_list(df, opt$include, priority = priority, quiet = opt$quiet)

# ── Write CSV ─────────────────────────────────────────────────
cols <- c("ID","catalog_number","Collection","collectionCode","species","subspecies","taxon",
          "scientificName","Country","countryCode","stateProvince","locality","preparations",
          "Lat","Lon","coordinateUncertaintyInMeters","coord_precision","coord_source",
          "date","year","sex","basisOfRecord","gbifID","occurrenceID","is_priority","forced_include")
write_csv(df[, cols], opt$output)
message("Wrote ", opt$output, ": ", nrow(df), " specimens (",
        sum(!is.na(df$Lat)), " with coordinates).")
if (opt$minimize && !is.null(opt$output_full)) {
  write_csv(df_full[, cols], opt$output_full)
  message("Wrote ", opt$output_full, ": ", nrow(df_full), " specimens (complete set).")
}

# ── Summary counts: per institution / taxon / preparation ────
# Describes the written sample sheet (the minimized panel when --minimize).
# Printed to stdout and written as a tidy TSV (category, value, n).
# A specimen with several preparations (e.g. "skin; skeleton") is counted
# once per preparation type, so preparation totals are prep-instances.
summarise_counts <- function(d) {
  taxon <- ifelse(!is.na(d$subspecies) & nzchar(d$subspecies),
                  paste(d$species, d$subspecies),
                  ifelse(!is.na(d$species) & nzchar(d$species), d$species, "(unidentified)"))
  inst  <- ifelse(is.na(d$Collection) | !nzchar(d$Collection), "(none)", d$Collection)
  prep_list <- strsplit(ifelse(is.na(d$preparations), "", d$preparations), "\\s*[;,|/]\\s*")
  prep <- unlist(lapply(prep_list, function(x) { x <- x[nzchar(x)]; if (!length(x)) "(none)" else x }))
  mk <- function(v, categ) {
    t <- sort(table(v), decreasing = TRUE)
    data.frame(category = categ, value = names(t), n = as.integer(t), stringsAsFactors = FALSE)
  }
  src <- ifelse(is.na(d$coord_source), "none", d$coord_source)
  rbind(mk(inst, "institution"), mk(taxon, "taxon"),
        mk(prep, "preparation"), mk(src, "coord_source"))
}
S <- summarise_counts(df)
for (ct in c("institution", "taxon", "preparation", "coord_source")) {
  blk <- S[S$category == ct, c("value", "n")]
  cat(sprintf("\n== counts per %s (n=%d%s) ==\n", ct, sum(blk$n),
              if (ct == "preparation") " prep-instances" else ""))
  print(blk, row.names = FALSE)
}
summary_path <- if (!is.null(opt$summary)) opt$summary else paste0(sub("\\.[^.]*$", "", opt$output), "_summary.tsv")
write_tsv(S, summary_path)
message("\nWrote ", summary_path, ": summary counts (", nrow(S), " rows).")

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
  dc <- fill_country_centroids(df)
  if (!is.null(dc)) {
    cpath <- if (!is.null(opt$centroid_output)) opt$centroid_output
             else paste0(sub("\\.[^.]*$", "", opt$output), "_centroid.csv")
    write_csv(dc[, cols], cpath)
    message("Wrote ", cpath, ": ", nrow(dc), " specimens with country centroids filled (",
            attr(dc, "filled"), " added; ", sum(!is.na(dc$Lat)), " now mappable).")
  }
}

# ── Optional dissolved range shapefile ───────────────────────
if (!is.null(opt$range_shapefile)) {
  ok <- requireNamespace("sf", quietly = TRUE) &&
        requireNamespace("rnaturalearth", quietly = TRUE)
  if (!ok) {
    warning("--range_shapefile needs sf + rnaturalearth; skipped.")
  } else {
    suppressPackageStartupMessages({ library(sf); library(rnaturalearth) })
    world <- ne_countries(scale = "medium", returnclass = "sf")
    iso   <- world$iso_a2_eh; iso[is.na(iso) | iso == "-99"] <- world$iso_a2[is.na(iso) | iso == "-99"]
    present <- unique(na.omit(df_full$countryCode))   # extent from the full pull, not the minimized set
    sel <- world[iso %in% present, ]
    if (nrow(sel)) {
      dir.create(dirname(opt$range_shapefile), showWarnings = FALSE, recursive = TRUE)
      diss <- st_sf(name = "range", geometry = st_union(st_geometry(sel)))
      st_write(diss, opt$range_shapefile, delete_dsn = TRUE, quiet = TRUE)
      message("Wrote ", opt$range_shapefile, ": dissolved range of ", nrow(sel), " countries.")
    } else warning("range_shapefile: no countries matched; skipped.")
  }
}

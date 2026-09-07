#!/usr/bin/env Rscript
# ============================================================
# Merge & filter IUCN Red List range maps
#
# Takes one or more IUCN range downloads (.zip, or already-unzipped
# .shp/.gpkg), unzips as needed, merges them into one layer, and keeps
# only the ranges you want by IUCN `presence` and `origin` codes.
# By default it keeps EXTANT + HISTORICAL native range and drops
# INTRODUCED (and vagrant / assisted-colonisation) polygons.
#
# The result is written as a shapefile or GeoPackage suitable for
# specimen_map.R --shapefile (map extent and optional range overlay).
#
# IUCN attribute codes
#   presence: 1 extant · 2 probably extant · 3 possibly extant
#             4 possibly extinct · 5 extinct · 6 presence uncertain
#   origin:   1 native · 2 reintroduced · 3 introduced · 4 vagrant
#             5 origin uncertain · 6 assisted colonisation
#   seasonal: 1 resident · 2 breeding · 3 non-breeding · 4 passage
#             5 seasonal occurrence uncertain
#
# Packages: optparse, sf, dplyr
#
# Examples
# --------
#   # One or many zips (or a folder of zips) -> filtered GeoPackage:
#   Rscript merge_iucn_ranges.R \
#     --input "Struthio_camelus.zip,Struthio_molybdophanes.zip" \
#     --output ostrich_iucn_range.gpkg
#
#   # A directory of IUCN zips, dissolved to one polygon for a map extent:
#   Rscript merge_iucn_ranges.R --input iucn_downloads/ \
#     --dissolve all --output ostrich_range_dissolved.shp
#
#   # Keep only strictly native, extant range (no historical/extinct):
#   Rscript merge_iucn_ranges.R --input ranges/ \
#     --presence extant --origin native --output extant_native.gpkg
# ============================================================
suppressPackageStartupMessages({
  library(optparse)
  library(sf)
  library(dplyr)
})

# keyword -> IUCN code(s); grouped keywords expand to several codes
PRESENCE_KW <- list(extant = c(1,2,3), historical = c(4,5),
                    "probably_extant"=2, "possibly_extant"=3,
                    "possibly_extinct"=4, extinct = 5, uncertain = 6)
ORIGIN_KW   <- list(native = 1, reintroduced = 2, introduced = 3,
                    vagrant = 4, uncertain = 5, assisted = 6)
SEASONAL_KW <- list(resident = 1, breeding = 2, non_breeding = 3,
                    passage = 4, uncertain = 5)

option_list <- list(
  make_option(c("-i", "--input"), type = "character",
              help = "Comma-separated IUCN .zip / .shp / .gpkg files and/or folders [required]."),
  make_option(c("-o", "--output"), type = "character", default = "iucn_ranges_merged.gpkg",
              help = "Output .gpkg or .shp [default: iucn_ranges_merged.gpkg]."),
  make_option(c("-p", "--presence"), type = "character", default = "extant,historical",
              help = paste("Presence codes/keywords to KEEP [default: extant,historical",
                           "= 1,2,3,4,5]. e.g. 'extant' or '1,2,3'.")),
  make_option(c("-O", "--origin"), type = "character", default = "native,reintroduced,uncertain",
              help = paste("Origin codes/keywords to KEEP [default: native,reintroduced,",
                           "uncertain = 1,2,5] — drops introduced/vagrant/assisted.")),
  make_option(c("-s", "--seasonal"), type = "character", default = "all",
              help = "Seasonal codes/keywords to KEEP, or 'all' [default: all]."),
  make_option(c("-d", "--dissolve"), type = "character", default = "none",
              help = "Dissolve output: 'none', 'species' (by name), or 'all' [default: none]."),
  make_option("--keep_na_codes", action = "store_true", default = FALSE,
              help = "Keep features whose presence/origin code is missing (default: drop)."),
  make_option("--quiet", action = "store_true", default = FALSE)
)
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Merge & filter IUCN range maps for specimen_map.R."))
say <- function(...) if (!opt$quiet) message(...)
if (is.null(opt$input)) stop("--input is required")

splitcsv <- function(x) trimws(strsplit(x, ",")[[1]])

# map keyword/number tokens -> integer code set
parse_codes <- function(str, kw) {
  if (is.null(str) || tolower(str) == "all") return(NULL)      # NULL = keep all
  toks <- tolower(splitcsv(str))
  out <- integer(0)
  for (t in toks) {
    if (grepl("^[0-9]+$", t)) out <- c(out, as.integer(t))
    else if (!is.null(kw[[t]])) out <- c(out, kw[[t]])
    else warning("Unrecognised code/keyword '", t, "' — ignored")
  }
  sort(unique(out))
}
keep_presence <- parse_codes(opt$presence, PRESENCE_KW)
keep_origin   <- parse_codes(opt$origin,   ORIGIN_KW)
keep_seasonal <- parse_codes(opt$seasonal, SEASONAL_KW)

# case-insensitive column lookup
findcol <- function(df, cands) {
  hit <- match(tolower(cands), tolower(names(df)))
  hit <- hit[!is.na(hit)]
  if (length(hit)) names(df)[hit[1]] else NA_character_
}

# ── Gather input shapefiles (unzip zips to a temp dir) ───────
tmp <- file.path(tempdir(), paste0("iucn_", as.integer(runif(1, 1e5, 9e5))))
dir.create(tmp, showWarnings = FALSE, recursive = TRUE)
items <- splitcsv(opt$input)
shp_paths <- character(0)
for (it in items) {
  if (dir.exists(it)) {
    zips <- list.files(it, pattern = "\\.zip$", full.names = TRUE, recursive = TRUE)
    shps <- list.files(it, pattern = "\\.(shp|gpkg)$", full.names = TRUE, recursive = TRUE)
    for (z in zips) { d <- file.path(tmp, tools::file_path_sans_ext(basename(z)))
                      dir.create(d, showWarnings = FALSE); unzip(z, exdir = d) }
    shp_paths <- c(shp_paths, shps)
  } else if (grepl("\\.zip$", it, ignore.case = TRUE)) {
    d <- file.path(tmp, tools::file_path_sans_ext(basename(it)))
    dir.create(d, showWarnings = FALSE); unzip(it, exdir = d)
  } else if (grepl("\\.(shp|gpkg)$", it, ignore.case = TRUE)) {
    shp_paths <- c(shp_paths, it)
  } else {
    warning("Skipping unrecognised input: ", it)
  }
}
# add everything unzipped into tmp
shp_paths <- unique(c(shp_paths,
                      list.files(tmp, pattern = "\\.(shp|gpkg)$", full.names = TRUE, recursive = TRUE)))
if (!length(shp_paths)) stop("No .shp/.gpkg layers found in --input.")
say("Found ", length(shp_paths), " layer(s).")

# ── Read + standardise each layer ────────────────────────────
std_layer <- function(p) {
  s <- tryCatch(st_read(p, quiet = TRUE), error = function(e) {
    warning("Could not read ", p, ": ", conditionMessage(e)); return(NULL) })
  if (is.null(s) || !nrow(s)) return(NULL)
  s <- st_transform(s, 4326)
  nm_c <- findcol(s, c("sci_name","binomial","SCINAME","sciname","species","SPECIES","tax_comm"))
  pr_c <- findcol(s, c("presence","PRESENCE","presenc"))
  or_c <- findcol(s, c("origin","ORIGIN"))
  se_c <- findcol(s, c("seasonal","SEASONAL","season"))
  data.frame(
    sci_name = if (!is.na(nm_c)) as.character(s[[nm_c]]) else NA_character_,
    presence = if (!is.na(pr_c)) suppressWarnings(as.integer(s[[pr_c]])) else NA_integer_,
    origin   = if (!is.na(or_c)) suppressWarnings(as.integer(s[[or_c]])) else NA_integer_,
    seasonal = if (!is.na(se_c)) suppressWarnings(as.integer(s[[se_c]])) else NA_integer_,
    stringsAsFactors = FALSE
  ) %>% st_sf(geometry = st_geometry(s))
}
layers <- Filter(Negate(is.null), lapply(shp_paths, std_layer))
if (!length(layers)) stop("No readable range features.")
merged <- do.call(rbind, layers)
say("Merged features: ", nrow(merged),
    " from ", length(unique(na.omit(merged$sci_name))), " taxon name(s).")

# ── Filter by presence / origin / seasonal ───────────────────
keeprow <- rep(TRUE, nrow(merged))
codef <- function(col, keep) {
  if (is.null(keep)) return(rep(TRUE, length(col)))
  ifelse(is.na(col), opt$keep_na_codes, col %in% keep)
}
keeprow <- keeprow & codef(merged$presence, keep_presence)
keeprow <- keeprow & codef(merged$origin,   keep_origin)
keeprow <- keeprow & codef(merged$seasonal, keep_seasonal)

say("Filter keep presence={", paste(keep_presence, collapse=","),
    "} origin={", paste(keep_origin, collapse=","),
    "} seasonal={", if (is.null(keep_seasonal)) "all" else paste(keep_seasonal, collapse=","), "}")
say("  kept ", sum(keeprow), " of ", nrow(merged), " features.")
out <- merged[keeprow, ]
if (!nrow(out)) stop("No features left after filtering — loosen --presence/--origin.")

# report what was kept
tab <- as.data.frame(table(presence = out$presence, origin = out$origin))
tab <- tab[tab$Freq > 0, ]
if (!opt$quiet) { cat("\n== kept features by presence x origin ==\n"); print(tab, row.names = FALSE) }

# ── Optional dissolve ────────────────────────────────────────
out <- suppressWarnings(st_make_valid(out))
if (opt$dissolve == "all") {
  out <- st_sf(name = "range", geometry = st_union(st_geometry(out)))
  say("Dissolved to a single polygon.")
} else if (opt$dissolve == "species") {
  out <- out %>% group_by(sci_name) %>%
    summarise(geometry = st_union(geometry), .groups = "drop")
  say("Dissolved by species: ", nrow(out), " polygon(s).")
} else if (opt$dissolve != "none") {
  warning("--dissolve must be none|species|all; leaving undissolved")
}

# ── Write ────────────────────────────────────────────────────
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)
st_write(out, opt$output, delete_dsn = TRUE, quiet = TRUE)
message("Wrote ", opt$output, ": ", nrow(out), " feature(s).")
message("Use it with:  Rscript specimen_map.R --samples <csv> --shapefile ", opt$output)

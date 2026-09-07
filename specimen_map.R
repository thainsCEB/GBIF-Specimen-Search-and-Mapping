#!/usr/bin/env Rscript
# ============================================================
# Specimen locality map — fully flag-driven, cross-project
#
# Single sample sheet required — no separate metadata file.
#
# Minimum required columns in sample sheet:
#   ID, Lat, Lon
#   + 'status' if using --status (values: missing | have | sequenced)
#   + whatever column you set --color_by to
#
# Typical columns: ID, status, Collection, locality, Lat, Lon,
#                  date, sex, species, subspecies, population,
#                  hybrid_status, ... (any extras are fine)
#
# Usage:
#   Rscript specimen_map.R [options]
#
# Examples:
#   # Natural Earth basemap (default, no API key needed)
#   Rscript specimen_map.R \
#     --samples    samples.tsv \
#     --shapefile  range.shp \
#     --show_shapefile \
#     --color_by   subspecies
#
#   # Multiple shapefiles define a combined extent
#   Rscript specimen_map.R \
#     --samples    samples.tsv \
#     --shapefile  "sp1.shp,sp2.shp,sp3.shp" \
#     --show_shapefile \
#     --color_by   subspecies
#
#   # Google Maps, color by Collection, status shapes on
#   Rscript specimen_map.R \
#     --samples   samples.tsv \
#     --shapefile range.shp \
#     --basemap   google --api_key YOUR_KEY \
#     --color_by  Collection \
#     --status
#
#   # Raster basemap, shapefile shown + matched to subspecies
#   Rscript specimen_map.R \
#     --samples    samples.tsv \
#     --shapefile  range.shp \
#     --match \
#     --basemap    /path/to/basemap.tif \
#     --color_by   subspecies
# ============================================================

# install.packages(c("optparse","ggmap","ggplot2","dplyr","readr","sf",
#                    "terra","tidyterra","RColorBrewer","viridis","ggnewscale",
#                    "rnaturalearth","rnaturalearthdata","maptiles"))
#
# The default 'ne' basemap uses maptiles to fetch Esri.WorldPhysical tiles
# (Natural Earth II aesthetic) for the map extent only — fast, small, cached.

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(sf)
  library(RColorBrewer)
  library(viridis)
})

# ── Flags ─────────────────────────────────────────────────────
option_list <- list(
  
  # Input / output
  make_option(c("-s", "--samples"), type = "character",
              help = "Sample sheet TSV/CSV (needs: ID, Lat, Lon) [required]"),
  make_option(c("-o", "--output"), type = "character", default = "specimen_map.png",
              help = "Output PNG path [default: specimen_map.png]"),
  make_option("--sep", type = "character", default = "\t",
              help = "Delimiter: auto-detected from .csv/.tsv extension, or set ',' or '\\t' [default: auto]"),
  
  # Shapefile(s)
  make_option(c("-f", "--shapefile"), type = "character",
              help = paste("Optional: one or more shapefiles (comma-separated) to define the",
                           "map extent and optionally overlay. Multiple bounding boxes are",
                           "unioned. If omitted, the extent comes from the sample points.")),
  make_option(c("-b", "--buffer"), type = "double", default = 5000,
              help = "Padding around the map extent in km (or degrees) [default: 5000]"),
  make_option("--buffer_unit", type = "character", default = "km",
              help = "Buffer unit: 'km' or 'deg' (degrees) [default: km]"),
  make_option(c("-e", "--extent"), type = "character", default = "points",
              help = paste("What the map extent (then buffered) is based on:",
                           "'points' = the sample points [default], 'shapefile' = the",
                           "shapefile bbox, 'both' = union of the two.")),
  make_option(c("-w", "--whole_map"), action = "store_true", default = FALSE,
              help = "Show the entire world instead of framing to the data (ignores --extent/--buffer)."),
  make_option("--show_shapefile", action = "store_true", default = FALSE,
              help = "Overlay shapefile on map as transparent polygon [default: off]"),
  make_option("--shp_alpha", type = "double", default = 0.2,
              help = "Shapefile fill transparency 0-1 [default: 0.2]"),
  make_option(c("-r", "--range"), type = "character", default = NULL,
              help = paste("Range shapefile(s) (comma-separated) to overlay, coloured to MATCH",
                           "the --color_by key (species / subspecies / population). One shared",
                           "colour key covers both points and ranges. Shortcut for",
                           "--shapefile ... --match; use --shp_field if the taxon field isn't",
                           "auto-detected. Does not change the map extent.")),
  make_option("--match", action = "store_true", default = FALSE,
              help = paste("Match shapefile part colors to point colors using --color_by.",
                           "Auto-detects shapefile field from common names (clstr_nm, Name,",
                           "name, clstr_name, fetr_nm, species, subspecies, taxon) or use",
                           "--shp_field to set explicitly. Implies --show_shapefile. [default: off]")),
  make_option("--shp_field", type = "character", default = NULL,
              help = paste("Override auto-detected shapefile field for --match.",
                           "Also used to color shapefile parts independently (without --match).")),
  make_option("--shp_outline", type = "character", default = "white",
              help = "Shapefile polygon outline color [default: white]"),
  
  # Basemap
  make_option(c("-B", "--basemap"), type = "character", default = "cbhypso",
              help = paste("'cbhypso' = Natural Earth 1:50m Cross-blended Hypsometric Tints",
                           "(downloaded & cached on first use) [default];",
                           "'ne' = Esri.WorldPhysical tiles for your extent;",
                           "'google' (needs --api_key); or a raster/vector file path.")),
  make_option("--hypso_url", type = "character",
              default = "https://naciscdn.org/naturalearth/50m/raster/HYP_50M_SR_W.zip",
              help = "Download URL for the Cross-blended Hypsometric Tints raster (for --basemap cbhypso)."),
  make_option("--cache_dir", type = "character", default = path.expand("~/.cache/specimen_map"),
              help = "Where the cbhypso raster is cached [default: ~/.cache/specimen_map]."),
  make_option("--maptype", type = "character", default = "roadmap",
              help = "Google Maps type: roadmap, terrain, satellite, hybrid [default: roadmap]"),
  make_option("--api_key", type = "character",
              default = Sys.getenv("GOOGLE_MAPS_KEY"),
              help = "Google Maps API key or set env var GOOGLE_MAPS_KEY"),
  
  # Country border overlay
  make_option("--borders", type = "character", default = "auto",
              help = paste("Overlay country/state borders on top of raster basemap.",
                           "'auto' = on for raster basemaps (ne + .tif etc.), off for vector.",
                           "'on' = always, 'off' = never. [default: auto]")),
  make_option("--border_color", type = "character", default = "#555555",
              help = "Country border line color [default: #555555]"),
  make_option("--border_size", type = "double", default = 0.3,
              help = "Country border line width [default: 0.3]"),

  # River overlay
  make_option("--rivers", type = "character", default = "auto",
              help = paste("Overlay Natural Earth rivers.",
                           "'auto' = on when --borders is on, 'on', or 'off'. [default: auto]")),
  make_option("--river_color", type = "character", default = "#6baed6",
              help = "River line color [default: #6baed6 (blue)]"),
  make_option("--river_size", type = "double", default = 0.25,
              help = "River line width [default: 0.25]"),
  make_option("--river_scale", type = "integer", default = 50L,
              help = "Natural Earth river scale: 10 (detailed) or 50 (default) [default: 50]"),

  # Geographic label overlays
  make_option(c("-l", "--labels"), type = "character", default = "off",
              help = paste("Overlay geographic labels.",
                           "'off' (default), 'on', or 'auto' (on when --borders is on).",
                           "Controls countries + states + natural features together.",
                           "[default: off]")),
  make_option("--label_countries", action = "store_true", default = TRUE,
              help = "Label country names [default: on]"),
  make_option("--no_label_countries", action = "store_false", dest = "label_countries",
              help = "Suppress country labels"),
  make_option("--label_states", action = "store_true", default = TRUE,
              help = "Label state/province names [default: on]"),
  make_option("--no_label_states", action = "store_false", dest = "label_states",
              help = "Suppress state/province labels"),
  make_option("--label_features", action = "store_true", default = TRUE,
              help = "Label physical features (mountain ranges, deserts, plains…) [default: on]"),
  make_option("--no_label_features", action = "store_false", dest = "label_features",
              help = "Suppress physical feature labels"),
  make_option("--label_size", type = "double", default = 2.5,
              help = "Base label text size (countries scale up, states scale down) [default: 2.5]"),
  make_option("--label_color", type = "character", default = "#333333",
              help = "Label text color [default: #333333]"),
  make_option("--label_scale", type = "integer", default = 50L,
              help = "Natural Earth label scale: 10 (detailed) or 50 [default: 50]"),

  # Taxa filtering
  make_option("--taxa", type = "character", default = NULL,
              help = paste("Comma-separated values to focus on (e.g. 'molinae,frontalis').",
                           "Filtered against --taxa_by column (default: species column).",
                           "Combined with --filter_mode. [default: show all]")),
  make_option("--taxa_by", type = "character", default = NULL,
              help = paste("Column to apply --taxa filter against.",
                           "Defaults to the auto-detected species column, NOT --color_by,",
                           "so you can filter by species while coloring by subspecies.",
                           "[default: auto-detect species column]")),
  make_option("--filter_mode", type = "character", default = "include",
              help = paste("How to apply --taxa: 'include' shows only those taxa,",
                           "'exclude' hides them. [default: include]")),
  make_option("--show_other", action = "store_true", default = FALSE,
              help = paste("When --taxa is set, show non-focal taxa as faded grey",
                           "points in the background. [default: off]")),
  make_option("--other_color", type = "character", default = "#BBBBBB",
              help = "Fill color for background 'other' taxa [default: #BBBBBB]"),
  make_option("--other_alpha", type = "double", default = 0.45,
              help = "Opacity for background 'other' taxa points (0–1) [default: 0.45]"),
  make_option("--filter_shapefile", action = "store_true", default = FALSE,
              help = paste("Also filter shapefile polygons to focal taxa",
                           "(requires --match or a field matching --color_by values).",
                           "[default: off — show all range polygons]")),
  
  # Status shapes (off by default)
  make_option("--status", action = "store_true", default = FALSE,
              help = paste("Color point shapes by 'status' column.",
                           "missing=x  have=filled-circle-pin  sequenced=filled-diamond-pin",
                           "Adds a separate shape legend. [default: off]")),
  make_option("--nudge", type = "double", default = 0.35,
              help = "Degrees to shift pin head above tail [default: 0.35]"),
  
  # Shape by column (alternative to --status)
  make_option(c("-S", "--shape_by"), type = "character", default = NULL,
              help = paste("Column to vary point shape by (e.g. Collection).",
                           "Used alongside --color_by to show two variables simultaneously.",
                           "Uses filled shapes (circle, square, diamond, triangle up/down)",
                           "compatible with the fill color scale.",
                           "Mutually exclusive with --status. [default: off]")),
  make_option("--shape_legend_title", type = "character", default = NULL,
              help = "Shape legend title [default: same as --shape_by]"),
  
  # Color
  make_option(c("-c", "--color_by"), type = "character", default = NULL,
              help = paste("Column to color points by. Auto-detected by taxonomic hierarchy:",
                           "species (if >1 unique value) → subspecies (if only 1 species) →",
                           "collection/institution. Override with any column name.",
                           "[default: auto-detect]")),
  make_option("--na_color", type = "character", default = "#999999",
              help = "Color for specimens where --color_by value is missing/NA [default: #999999]"),
  make_option("--na_label", type = "character", default = "unassigned",
              help = "Legend label for NA/missing specimens [default: unassigned]"),
  make_option(c("-p", "--palette"), type = "character", default = "Dark2",
              help = paste("RColorBrewer palette (Dark2, Set1, Paired...)",
                           "or viridis option (viridis, magma, plasma, turbo)",
                           "[default: Dark2]")),
  make_option("--legend_title", type = "character", default = NULL,
              help = "Color legend title [default: same as --color_by]"),
  
  # Labels
  make_option(c("-t", "--title"), type = "character", default = "",
              help = "Map title [default: none]"),
  make_option("--subtitle", type = "character", default = "",
              help = "Map subtitle [default: none]"),
  
  # Marker
  make_option("--point_size", type = "double", default = 5.0,
              help = "Marker size [default: 5]"),
  
  # Output dimensions
  make_option(c("-W", "--width"),  type = "double",  default = 10.0,
              help = "Output width in inches [default: 10]"),
  make_option(c("-H", "--height"), type = "double",  default = 10.0,
              help = "Output height in inches [default: 10]"),
  make_option(c("-d", "--dpi"),    type = "integer", default = 300L,
              help = "Resolution in DPI [default: 300]")
)

opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Plot specimen localities on a map."))

# ── Validate ──────────────────────────────────────────────────
if (is.null(opt$samples))   stop("--samples is required")
# --range is a shortcut: overlay taxon range polygon(s) coloured to match the
# --color_by key. It adds to --shapefile and turns on --match.
if (!is.null(opt$range) && nzchar(opt$range)) {
  opt$shapefile <- if (!is.null(opt$shapefile) && nzchar(opt$shapefile))
    paste(opt$shapefile, opt$range, sep = ",") else opt$range
  opt$match <- TRUE
}
# --shapefile is optional: without it the map extent comes from the sample
# points and no range polygon is drawn.
have_shp <- !is.null(opt$shapefile) && nzchar(opt$shapefile)
if (!have_shp) { opt$show_shapefile <- FALSE; opt$match <- FALSE }
# --match implies --show_shapefile
if (opt$match) opt$show_shapefile <- TRUE
# Only require API key when google basemap is explicitly chosen
if (opt$basemap == "google" && (is.null(opt$api_key) || opt$api_key == ""))
  stop("--basemap google requires --api_key or env var GOOGLE_MAPS_KEY")
# --status and --shape_by are mutually exclusive
if (opt$status && !is.null(opt$shape_by)) {
  warning("--status and --shape_by are mutually exclusive; --status takes priority")
  opt$shape_by <- NULL
}

# Parse comma-separated shapefile paths (empty when no shapefile given)
shp_paths <- if (have_shp) trimws(strsplit(opt$shapefile, ",")[[1]]) else character(0)

# Resolve the Cross-blended Hypsometric Tints basemap: download once, cache,
# then treat it as a raster path for the rest of the script.
if (identical(opt$basemap, "cbhypso")) {
  dir.create(opt$cache_dir, recursive = TRUE, showWarnings = FALSE)
  tif <- Sys.glob(file.path(opt$cache_dir, "**", "*.tif"))
  tif <- c(tif, Sys.glob(file.path(opt$cache_dir, "*.tif")))
  tif <- tif[grepl("HYP", basename(tif), ignore.case = TRUE)]
  if (!length(tif)) {
    zip <- file.path(opt$cache_dir, basename(opt$hypso_url))
    if (!file.exists(zip)) {
      message("Downloading 1:50m Cross-blended Hypsometric Tints raster (first run) ...")
      old <- options(timeout = 3600); on.exit(options(old), add = TRUE)
      utils::download.file(opt$hypso_url, zip, mode = "wb")
    }
    utils::unzip(zip, exdir = opt$cache_dir)
    tif <- Sys.glob(file.path(opt$cache_dir, "**", "*.tif"))
    tif <- c(tif, Sys.glob(file.path(opt$cache_dir, "*.tif")))
    tif <- tif[grepl("HYP", basename(tif), ignore.case = TRUE)]
    if (!length(tif)) stop("Cross-blended Hypso .tif not found after unzip in ", opt$cache_dir)
  }
  opt$basemap <- tif[1]
  message("Basemap: cross-blended hypso -> ", opt$basemap)
}

# ── Load sample sheet ─────────────────────────────────────────
message("Loading samples...")
# Auto-detect delimiter from file extension if not overridden
if (opt$sep == "\t" && grepl("\\.csv$", opt$samples, ignore.case = TRUE)) {
  opt$sep <- ","
  message("  Auto-detected CSV format from file extension")
}
sep_char <- opt$sep
df <- read_delim(opt$samples, delim = sep_char, show_col_types = FALSE)

# Strip leading/trailing whitespace from all character columns.
# CSV exports often pad values with spaces which silently breaks matching.
df <- df %>% mutate(across(where(is.character), trimws))

# Flexible column name matching (case-insensitive)
flex_col <- function(data, candidates) {
  found <- intersect(tolower(candidates), tolower(colnames(data)))
  if (length(found) == 0) return(NULL)
  colnames(data)[tolower(colnames(data)) == found[1]]
}

lat_col <- flex_col(df, c("Lat","lat","latitude","decimalLatitude"))
lon_col <- flex_col(df, c("Lon","lon","longitude","decimalLongitude"))
id_col  <- flex_col(df, c("ID","id","sampleID","sample_id","SampleID",
                          "voucher_number","voucherNumber","voucher",
                          "catalog_number","catalogNumber"))

if (is.null(lat_col)) stop("No Lat column found. Expected: Lat, lat, latitude, decimalLatitude")
if (is.null(lon_col)) stop("No Lon column found. Expected: Lon, lon, longitude, decimalLongitude")
if (is.null(id_col))  stop("No ID column found. Expected: ID, id, sampleID, sample_id")

# Standardise names internally
df <- df %>% rename(Lat = all_of(lat_col), Lon = all_of(lon_col), ID = all_of(id_col))

# Always detect species / subspecies columns — used by --taxa filter even when
# --color_by is set explicitly, so filtering and coloring can use different columns.
sp_candidates  <- c("species","Species","sp","taxon","Taxon",
                    "sci_name","SCI_NAME","scientific_name","scientificName",
                    "binomial","Binomial")
ssp_candidates <- c("subspecies","Subspecies","subsp","ssp","subsp.",
                    "infraspecificEpithet","race","Race",
                    "population","Population","morph")

auto_sp_col  <- flex_col(df, sp_candidates)    # always set; NULL if no species col found
auto_ssp_col <- flex_col(df, ssp_candidates)   # always set; NULL if no subsp col found

# Auto-detect color_by: prefer species → subspecies → collection
# Skipped entirely when the user sets --color_by explicitly.
if (is.null(opt$color_by)) {
  col_candidates <- c("institutionCode","Collection","collection",
                      "repunit","source_group","Institution","institution")
  col_col <- flex_col(df, col_candidates)

  if (!is.null(auto_sp_col)) {
    n_sp <- length(unique(na.omit(df[[auto_sp_col]])))
    if (n_sp > 1) {
      opt$color_by <- auto_sp_col
      message("  Auto-detected --color_by: '", auto_sp_col,
              "' (", n_sp, " species — species level)")
    } else if (!is.null(auto_ssp_col)) {
      opt$color_by <- auto_ssp_col
      message("  Auto-detected --color_by: '", auto_ssp_col,
              "' (only 1 unique species in '", auto_sp_col,
              "' — falling back to subspecies level)")
    } else {
      opt$color_by <- auto_sp_col
      message("  Auto-detected --color_by: '", auto_sp_col,
              "' (1 species; no subspecies column found)")
    }
  } else if (!is.null(auto_ssp_col)) {
    opt$color_by <- auto_ssp_col
    message("  Auto-detected --color_by: '", auto_ssp_col,
            "' (no species column found — using subspecies)")
  } else if (!is.null(col_col)) {
    opt$color_by <- col_col
    message("  Auto-detected --color_by: '", col_col,
            "' (no taxonomic columns found — using collection)")
  } else {
    stop(paste0(
      "Could not auto-detect a color column.\n",
      "Looked for: species (", paste(sp_candidates[1:4], collapse=", "), "...),\n",
      "            subspecies (", paste(ssp_candidates[1:4], collapse=", "), "...),\n",
      "            collection (institutionCode, Collection, collection...).\n",
      "Please set --color_by explicitly.\n",
      "Available columns: ", paste(colnames(df), collapse = ", ")
    ))
  }
}

if (!opt$color_by %in% colnames(df))
  stop(paste0("--color_by '", opt$color_by, "' not found.\n",
              "Available: ", paste(colnames(df), collapse = ", ")))

if (opt$status) {
  if (!"status" %in% tolower(colnames(df)))
    stop("--status requires a 'status' column in the sample sheet")
  status_col <- colnames(df)[tolower(colnames(df)) == "status"]
  df <- df %>% rename(status = all_of(status_col)) %>%
    mutate(status = factor(status, levels = c("missing","have","sequenced")))
}

df_plot <- df %>% filter(!is.na(Lat), !is.na(Lon))

# ── Taxa filter ───────────────────────────────────────────────
df_other    <- NULL   # non-focal taxa (shown as grey background if --show_other)
focal_taxa  <- NULL   # the resolved set of focal values

if (!is.null(opt[["taxa"]])) {
  taxa_vals <- trimws(strsplit(opt[["taxa"]], ",")[[1]])

  # Resolve which column to filter against.
  # Default: species column (auto_sp_col), NOT color_by — this lets you do
  # --taxa molinae --color_by subspecies without confusion.
  if (!is.null(opt$taxa_by)) {
    taxa_by_col <- opt$taxa_by
    if (!taxa_by_col %in% colnames(df))
      stop("--taxa_by '", taxa_by_col, "' not found. Available: ",
           paste(colnames(df), collapse = ", "))
  } else if (!is.null(auto_sp_col)) {
    taxa_by_col <- auto_sp_col
  } else {
    taxa_by_col <- opt$color_by   # last resort: same as color_by
  }
  message("  Filtering taxa by column: '", taxa_by_col, "'")

  # Check all requested values exist in the data
  all_vals <- unique(na.omit(as.character(df[[taxa_by_col]])))
  bad_vals <- setdiff(taxa_vals, all_vals)
  if (length(bad_vals) > 0)
    warning("--taxa values not found in '", taxa_by_col, "': ",
            paste(bad_vals, collapse = ", "),
            "\nAvailable: ", paste(sort(all_vals), collapse = ", "))

  in_focal <- as.character(df_plot[[taxa_by_col]]) %in% taxa_vals

  if (opt$filter_mode == "include") {
    focal_taxa <- taxa_vals
    df_other   <- df_plot[!in_focal, ]
    df_plot    <- df_plot[ in_focal, ]
  } else if (opt$filter_mode == "exclude") {
    focal_taxa <- setdiff(all_vals, taxa_vals)
    df_other   <- df_plot[ in_focal, ]
    df_plot    <- df_plot[!in_focal, ]
  } else {
    stop("--filter_mode must be 'include' or 'exclude'")
  }

  message("  Taxa filter (", opt$filter_mode, " '", paste(taxa_vals, collapse=","),
          "' in '", taxa_by_col, "'): ",
          nrow(df_plot), " focal  /  ", nrow(df_other), " other")
}

# ── Re-evaluate color_by for focal mode ──────────────────────
# If --taxa narrowed the focal set to a single species and color_by was
# auto-detected at species level, automatically drop to subspecies so the
# focal specimens are distinguished within that one species.
if (!is.null(focal_taxa) && !is.null(auto_sp_col) &&
    opt$color_by == auto_sp_col && !is.null(auto_ssp_col) &&
    auto_ssp_col %in% colnames(df_plot)) {
  n_focal_sp <- length(unique(na.omit(as.character(df_plot[[opt$color_by]]))))
  if (n_focal_sp <= 1) {
    message("  Focal data has ", n_focal_sp, " species — switching color_by from '",
            auto_sp_col, "' to '", auto_ssp_col, "' (subspecies level)")
    opt$color_by <- auto_ssp_col
  }
}

# ── Shapefile(s) — optional (overlay, CRS, optional extent) ──
shp_list <- list(); shp <- NULL; shp_bbox <- NULL
if (have_shp) {
  message("Reading shapefile(s)...")
  # Keep native CRS — coord_sf will use it to drive the display projection
  shp_list <- lapply(shp_paths, function(p) {
    if (!file.exists(p)) stop("Shapefile not found: ", p)
    message("  ", p)
    st_read(p, quiet = TRUE)
  })
  shp_crs <- st_crs(shp_list[[1]])
  message("  Shapefile CRS: ", shp_crs$input)
  shp <- shp_list[[1]]
  all_bbox <- lapply(shp_list, function(s) st_bbox(st_transform(s, 4326)))
  shp_bbox <- c(xmin = min(sapply(all_bbox, `[[`, "xmin")),
                xmax = max(sapply(all_bbox, `[[`, "xmax")),
                ymin = min(sapply(all_bbox, `[[`, "ymin")),
                ymax = max(sapply(all_bbox, `[[`, "ymax")))
} else {
  shp_crs <- sf::st_crs(4326)   # no shapefile -> display in WGS84
}

# Bounding box of the plotted sample points
pts <- df_plot[is.finite(df_plot$Lon) & is.finite(df_plot$Lat), ]
pts_bbox <- if (nrow(pts))
  c(xmin = min(pts$Lon), xmax = max(pts$Lon), ymin = min(pts$Lat), ymax = max(pts$Lat)) else NULL

# ── Choose the raw extent (then buffered below) ──────────────
if (isTRUE(opt$whole_map)) {
  xmin <- -180; xmax <- 180; ymin <- -90; ymax <- 90       # entire world; no buffer
  message("Whole-map mode: showing the entire world.")
} else {
ext_mode <- opt$extent
if (!ext_mode %in% c("points","shapefile","both")) { warning("--extent must be points|shapefile|both; using points"); ext_mode <- "points" }
if (ext_mode == "shapefile" && is.null(shp_bbox)) { warning("--extent shapefile but no shapefile given; using points"); ext_mode <- "points" }
if (ext_mode == "points" && is.null(pts_bbox)) {
  if (!is.null(shp_bbox)) { warning("--extent points but no plottable points; using shapefile"); ext_mode <- "shapefile" }
  else stop("No plottable points and no shapefile — cannot set a map extent.")
}
bb <- switch(ext_mode,
  points    = pts_bbox,
  shapefile = shp_bbox,
  both      = c(xmin = min(c(pts_bbox["xmin"], shp_bbox["xmin"])),
                xmax = max(c(pts_bbox["xmax"], shp_bbox["xmax"])),
                ymin = min(c(pts_bbox["ymin"], shp_bbox["ymin"])),
                ymax = max(c(pts_bbox["ymax"], shp_bbox["ymax"]))))
xmin_raw <- unname(bb["xmin"]); xmax_raw <- unname(bb["xmax"])
ymin_raw <- unname(bb["ymin"]); ymax_raw <- unname(bb["ymax"])
if (xmax_raw == xmin_raw) { xmin_raw <- xmin_raw - 0.5; xmax_raw <- xmax_raw + 0.5 }
if (ymax_raw == ymin_raw) { ymin_raw <- ymin_raw - 0.5; ymax_raw <- ymax_raw + 0.5 }
message("Extent (", ext_mode, ", +", opt$buffer, opt$buffer_unit, " buffer): ",
        round(xmin_raw,2), " to ", round(xmax_raw,2), " lon, ",
        round(ymin_raw,2), " to ", round(ymax_raw,2), " lat")

# Convert buffer to degrees if given in km
# 1 degree lat ≈ 111 km; 1 degree lon ≈ 111 * cos(lat) km
buf_deg <- if (opt$buffer_unit == "km") {
  mid_lat <- mean(c(ymin_raw, ymax_raw))
  buf_lat <- opt$buffer / 111
  buf_lon <- opt$buffer / (111 * cos(mid_lat * pi / 180))
  list(lat = buf_lat, lon = buf_lon)
} else {
  list(lat = opt$buffer, lon = opt$buffer)
}

xmin <- xmin_raw - buf_deg$lon
xmax <- xmax_raw + buf_deg$lon
ymin <- ymin_raw - buf_deg$lat
ymax <- ymax_raw + buf_deg$lat
}   # end !whole_map

# Keep the frame within valid geographic bounds
xmin <- max(xmin, -180); xmax <- min(xmax, 180)
ymin <- max(ymin, -90);  ymax <- min(ymax, 90)

# Shared cache directory (basemap tiles, NE vector layers)
cache_dir <- tools::R_user_dir("specimen_map", "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

# ── Color palette ─────────────────────────────────────────────
# When --taxa is set, build palette from focal specimens only — this keeps
# non-focal taxa' subspecies/values out of the legend.
# Otherwise use ALL rows so every value appears in the legend even if
# some specimens have no coordinates.
viridis_opts <- c("viridis","magma","plasma","inferno","cividis","turbo")
# Auto legend title from the coloured level: Species / Subspecies / Populations.
auto_legend_title <- function(cb) {
  k <- tolower(cb)
  if (k %in% c("species","sp","scientificname","scientific_name","binomial","taxon_species")) "Species"
  else if (k %in% c("subspecies","subsp","ssp","subsp.","infraspecificepithet","race")) "Subspecies"
  else if (k %in% c("population","populations","pop","repunit","source_group")) "Populations"
  else if (k == "taxon") "Species / subspecies"
  else cb
}
legend_ttl   <- if (!is.null(opt$legend_title)) opt$legend_title else auto_legend_title(opt$color_by)

make_palette <- function(values, palette) {
  vals <- sort(unique(na.omit(as.character(values))))
  n    <- length(vals)
  if (palette %in% viridis_opts) {
    cols <- viridis(n, option = palette)
  } else {
    max_n <- tryCatch(brewer.pal.info[palette, "maxcolors"], error = function(e) 8L)
    cols  <- if (n <= max_n) {
      brewer.pal(max(3, n), palette)[seq_len(n)]
    } else {
      colorRampPalette(brewer.pal(8, palette))(n)
    }
  }
  setNames(cols, vals)
}

palette_src <- if (!is.null(focal_taxa)) df_plot else df
point_pal   <- make_palette(palette_src[[opt$color_by]], opt$palette)

# If any focal specimens have NA for color_by, add an explicit "unassigned" entry
# so those points are plotted (as grey) rather than silently dropped.
has_na_focal <- any(is.na(df_plot[[opt$color_by]]))
if (has_na_focal) {
  point_pal <- c(point_pal, setNames(opt$na_color, opt$na_label))
  message("  ", sum(is.na(df_plot[[opt$color_by]])),
          " specimens have no '", opt$color_by, "' value — colored as '",
          opt$na_label, "' (", opt$na_color, ")")
}

# Assign color column — NA values map to the explicit "unassigned" palette entry
df_plot <- df_plot %>%
  mutate(
    .cb_val   = as.character(.data[[opt$color_by]]),
    .cb_val   = ifelse(is.na(.data[[opt$color_by]]), opt$na_label, .cb_val),
    .pt_color = point_pal[.cb_val]
  ) %>%
  select(-.cb_val)

# ── Shapefile coloring ────────────────────────────────────────
shp_has_own_legend <- FALSE
shp_pal            <- NULL
resolved_field     <- NULL

# Helper: case-insensitive palette lookup with trim + suffix word match
# Handles "Pyrrhura devillei" matching palette key "devillei"
safe_color_lookup <- function(pal, values) {
  pal_keys   <- names(pal)
  keys_lower <- tolower(trimws(pal_keys))
  vals_str   <- as.character(values)
  
  out <- vector("character", length(vals_str))
  for (i in seq_along(vals_str)) {
    v <- tolower(trimws(vals_str[i]))
    # 1. Exact match
    idx <- match(v, keys_lower)
    if (!is.na(idx)) { out[i] <- pal[idx]; next }
    # 2. Palette key is a whole word at the end of the shapefile value
    #    e.g. "pyrrhura devillei" ends with "devillei"
    word_hit <- which(endsWith(v, keys_lower) &
                      (nchar(v) == nchar(keys_lower) |
                       substr(v, nchar(v) - nchar(keys_lower), nchar(v) - nchar(keys_lower)) == " "))
    if (length(word_hit) > 0) { out[i] <- pal[word_hit[1]]; next }
    # 3. Any palette key found as a substring
    sub_hit <- which(sapply(keys_lower, function(k) grepl(paste0("\\b", k, "\\b"), v)))
    if (length(sub_hit) > 0) { out[i] <- pal[sub_hit[1]]; next }
    out[i] <- NA_character_
  }
  out
}

color_one_shp <- function(s) {
  geom_col   <- attr(s, "sf_column")
  s_non_geom <- colnames(s)[colnames(s) != geom_col]

  if (!is.null(opt$shp_field)) {
    rf <- opt$shp_field
    if (!rf %in% s_non_geom)
      stop(paste0("--shp_field '", rf, "' not found.\n",
                  "Available columns: ", paste(s_non_geom, collapse = ", ")))
  } else {
    candidates <- c("clstr_nm","Name","name","clstr_name","fetr_nm",
                    "cluster","label","species","sci_name","SCI_NAME",
                    "SCINAME","binomial","BINOMIAL","subspecies","taxon",
                    "id","ID","tax_name","TAX_NAME")
    # Case-insensitive candidate search
    ci_cols  <- tolower(s_non_geom)
    ci_cands <- tolower(candidates)
    hit      <- match(ci_cands, ci_cols)
    rf       <- if (!all(is.na(hit))) s_non_geom[hit[!is.na(hit)][1]] else NULL
    if (is.null(rf)) {
      # Fall back: first non-numeric, non-geometry column
      str_cols <- s_non_geom[!sapply(s_non_geom, function(cc) is.numeric(s[[cc]]))]
      rf <- if (length(str_cols) > 0) str_cols[1] else NULL
    }
  }

  if (opt$match && !is.null(rf)) {
    colors <- safe_color_lookup(point_pal, s[[rf]])
    n_matched <- sum(!is.na(colors))
    message("    Field '", rf, "': ", n_matched, "/", nrow(s),
            " polygons matched to point palette")
    if (n_matched == 0)
      message("    WARNING: no matches — shapefile values: [",
              paste(head(unique(s[[rf]]), 6), collapse=", "),
              "] vs point_pal keys: [",
              paste(head(names(point_pal), 6), collapse=", "), "]")
    s <- s %>% mutate(.shp_color = ifelse(is.na(colors), "#CCCCCC", colors))

  } else if (!is.null(rf)) {
    colors <- safe_color_lookup(shp_pal, s[[rf]])
    s <- s %>% mutate(.shp_color = ifelse(is.na(colors), "#CCCCCC", colors))

  } else {
    message("    WARNING: no usable text field found — using grey fill")
    message("    Shapefile columns: ", paste(s_non_geom, collapse = ", "))
    s <- s %>% mutate(.shp_color = "#AAAAAA")
  }
  s
}

if (opt$show_shapefile) {

  # Detect resolved_field from first shapefile for legend labels
  geom_col_1   <- attr(shp_list[[1]], "sf_column")
  shp_non_geom <- colnames(shp_list[[1]])[colnames(shp_list[[1]]) != geom_col_1]
  message("  Shapefile columns available: ", paste(shp_non_geom, collapse = ", "))

  if (!is.null(opt$shp_field)) {
    resolved_field <- opt$shp_field
  } else {
    candidates <- c("clstr_nm","Name","name","clstr_name","fetr_nm",
                    "cluster","label","species","sci_name","SCI_NAME",
                    "SCINAME","binomial","BINOMIAL","subspecies","taxon",
                    "id","ID","tax_name","TAX_NAME")
    ci_cols    <- tolower(shp_non_geom)
    ci_cands   <- tolower(candidates)
    hit        <- match(ci_cands, ci_cols)
    resolved_field <- if (!all(is.na(hit))) shp_non_geom[hit[!is.na(hit)][1]] else NULL
    if (is.null(resolved_field)) {
      str_cols <- shp_non_geom[!sapply(shp_non_geom, function(cc) is.numeric(shp_list[[1]][[cc]]))]
      resolved_field <- if (length(str_cols) > 0) str_cols[1] else NULL
    }
    if (!is.null(resolved_field))
      message("  Using field '", resolved_field,
              "' for coloring  (override with --shp_field)")
    else
      message("  No text field detected — use --shp_field to specify one")
  }

  if (!opt$match && !is.null(resolved_field)) {
    # Build combined palette across ALL shapefiles
    all_vals        <- unlist(lapply(shp_list, function(s) as.character(s[[resolved_field]])))
    shp_pal         <- make_palette(all_vals, opt$palette)
    shp_has_own_legend <- TRUE
  }

  # Apply coloring to EVERY shapefile
  shp_list <- lapply(shp_list, color_one_shp)
  shp      <- shp_list[[1]]

  # If the shapefile's colors are already captured by the point legend,
  # suppress the shapefile's own legend entry to avoid duplicates.
  # This happens when --match is used OR when the palettes overlap significantly.
  if (shp_has_own_legend) {
    shp_colors   <- na.omit(unique(shp_list[[1]]$.shp_color))
    point_colors <- unname(point_pal)
    overlap_frac <- mean(shp_colors %in% point_colors)
    if (overlap_frac >= 0.5) {
      shp_has_own_legend <- FALSE
      message("  Shapefile colors overlap point palette — using unified legend")
    }
  }

  # Optional: filter shapefile polygons to focal taxa only
  if (opt$filter_shapefile && !is.null(focal_taxa) && !is.null(resolved_field)) {
    shp_list <- lapply(shp_list, function(s) {
      vals_lower <- tolower(trimws(as.character(s[[resolved_field]])))
      keep <- sapply(vals_lower, function(v)
        any(sapply(tolower(trimws(focal_taxa)), function(k)
          v == k | endsWith(v, paste0(" ", k)) | grepl(paste0("\\b", k, "\\b"), v)
        ))
      )
      s[keep, ]
    })
    shp <- shp_list[[1]]
    message("  Shapefile filtered to ", sum(sapply(shp_list, nrow)),
            " polygons matching focal taxa")
  }
}

# ── Basemap ───────────────────────────────────────────────────
message("Building basemap (", opt$basemap, ")...")
is_raster   <- function(p) grepl("\\.(tif|tiff|img|vrt|nc|grd)$", p, ignore.case = TRUE)
is_vector   <- function(p) grepl("\\.(shp|gpkg|geojson|kml)$",   p, ignore.case = TRUE)

if (opt$basemap == "ne") {
  suppressPackageStartupMessages({
    library(maptiles)
    library(terra)
    library(tidyterra)
  })

  # Tile cache — tiles are only downloaded for the visible extent,
  # then stored locally so subsequent runs are instant.
  # (cache_dir is already created in the border/river block above)

  # Build an sf polygon from the buffered extent for tile fetching
  extent_sf <- sf::st_as_sfc(
    sf::st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                crs = sf::st_crs(4326))
  )

  # "Natural Earth II with Shaded Relief and Water" aesthetic.
  # Esri.WorldPhysical tiles are the closest hosted equivalent.
  message("  Fetching Natural Earth basemap tiles (cached after first run)...")
  tiles <- tryCatch(
    get_tiles(
      x        = extent_sf,
      provider = "Esri.WorldPhysical",
      zoom     = NULL,          # auto-select zoom for output resolution
      crop     = TRUE,
      cachedir = cache_dir,
      verbose  = FALSE
    ),
    error = function(e) {
      message("  Tile fetch failed: ", conditionMessage(e),
              "\n  Falling back to vector Natural Earth.")
      NULL
    }
  )

  if (!is.null(tiles)) {
    base_plot <- ggplot() +
      geom_spatraster_rgb(data = tiles, interpolate = TRUE, max_col_value = 255)
  } else {
    # Offline fallback: built-in vector Natural Earth (no download needed)
    suppressPackageStartupMessages({
      library(rnaturalearth); library(rnaturalearthdata)
    })
    world <- ne_countries(scale = "medium", returnclass = "sf")
    base_plot <- ggplot() +
      geom_sf(data = world, fill = "#D8E8D0", colour = "#AAAAAA", linewidth = 0.25)
  }

} else if (opt$basemap == "google") {
  suppressPackageStartupMessages(library(ggmap))
  register_google(key = opt$api_key)
  # Auto-calculate zoom from the buffered extent so the tile fills the frame
  lon_span  <- xmax - xmin
  lat_span  <- ymax - ymin
  zoom_auto <- max(1L, min(18L,
    floor(log2(360 / max(lon_span, lat_span))) + 1L))
  message("  Auto-calculated Google Maps zoom: ", zoom_auto)
  gmap <- get_googlemap(
    center  = c(lon = (xmin + xmax) / 2, lat = (ymin + ymax) / 2),
    zoom    = zoom_auto, maptype = opt$maptype)
  # ggmap() conflicts with coord_sf() (required by geom_sf layers).
  # Solution: extract the tile as a raster and place it with annotation_raster()
  # so the rest of the plot can use coord_sf() normally.
  bb <- attr(gmap, "bb")  # bounding box of the fetched tile
  base_plot <- ggplot() +
    annotation_raster(
      raster = gmap,
      xmin   = bb$ll.lon, xmax = bb$ur.lon,
      ymin   = bb$ll.lat, ymax = bb$ur.lat,
      interpolate = TRUE
    )
  
} else if (is_raster(opt$basemap)) {
  suppressPackageStartupMessages({ library(terra); library(tidyterra) })
  r <- rast(opt$basemap)
  # Crop the raster to the map extent before plotting.
  # Without this, a global or large raster is aggressively downsampled → fuzzy.
  # Build extent in raster's own CRS (transform from WGS84 if needed).
  ext_wgs <- ext(xmin - 1, xmax + 1, ymin - 1, ymax + 1)   # small padding
  r_crs   <- crs(r, describe = TRUE)$code
  if (!is.na(r_crs) && !is.null(r_crs) && r_crs != "4326") {
    crop_ext <- project(ext_wgs, from = "EPSG:4326", to = crs(r))
  } else {
    crop_ext <- ext_wgs
  }
  r_crop <- tryCatch(crop(r, crop_ext), error = function(e) {
    message("  Crop failed (", conditionMessage(e), ") — using full raster")
    r
  })
  base_plot <- ggplot() +
    geom_spatraster_rgb(data = r_crop, interpolate = TRUE,
                        max_col_value = if (max(values(r_crop[[1]]), na.rm=TRUE) > 1) 255 else 1)
  
} else if (is_vector(opt$basemap)) {
  bvec <- st_read(opt$basemap, quiet = TRUE) %>% st_transform(4326)
  base_plot <- ggplot() +
    geom_sf(data = bvec, fill = "#E8E0D5", colour = "#AAAAAA", linewidth = 0.2)
  
} else {
  stop("--basemap must be 'ne', 'google', a raster file, or a vector file")
}

# ── Country border overlay ────────────────────────────────────
# Auto-enable for raster basemaps; controllable with --borders on/off
is_raster_basemap <- opt$basemap == "ne" || is_raster(opt$basemap)
show_borders <- switch(opt$borders,
  "auto" = is_raster_basemap,
  "on"   = TRUE,
  "off"  = FALSE,
  {
    warning("--borders must be 'auto', 'on', or 'off'; using 'auto'")
    is_raster_basemap
  }
)

# Helper: download a Natural Earth vector layer and cache as RDS
ne_cached <- function(scale, type, category) {
  suppressPackageStartupMessages({
    library(rnaturalearth); library(rnaturalearthdata)
  })
  rds <- file.path(cache_dir, paste0("ne_", scale, "_", type, ".rds"))
  if (file.exists(rds)) {
    message("    Using cached: ", basename(rds))
    return(readRDS(rds))
  }
  message("    Downloading: ", type, " (scale ", scale, ")...")
  layer <- tryCatch(
    ne_download(scale = scale, type = type, category = category, returnclass = "sf"),
    error = function(e) {
      message("    Download failed: ", conditionMessage(e)); NULL
    }
  )
  if (!is.null(layer)) saveRDS(layer, rds)
  layer
}

if (show_borders) {
  countries <- ne_cached(50, "admin_0_countries", "cultural")
  if (!is.null(countries)) {
    base_plot <- base_plot +
      geom_sf(data = countries, fill = NA, colour = opt$border_color,
              linewidth = opt$border_size, inherit.aes = FALSE)
    message("  Added country borders")
  }
}

# ── River overlay ─────────────────────────────────────────────
show_rivers <- switch(opt$rivers,
  "auto" = show_borders,
  "on"   = TRUE,
  "off"  = FALSE,
  show_borders
)

if (show_rivers) {
  rivers <- ne_cached(opt$river_scale, "rivers_lake_centerlines", "physical")
  if (!is.null(rivers)) {
    base_plot <- base_plot +
      geom_sf(data = rivers, fill = NA, colour = opt$river_color,
              linewidth = opt$river_size, inherit.aes = FALSE)
    message("  Added rivers (scale ", opt$river_scale, ")")
  }
}

# ── Geographic label layers ───────────────────────────────────
# Labels are collected here and applied AFTER points so they render on top
show_labels <- switch(opt$labels,
  "auto" = show_borders,
  "on"   = TRUE,
  "off"  = FALSE,
  FALSE
)

p_labels <- list()

if (show_labels) {
  suppressPackageStartupMessages({
    library(rnaturalearth); library(rnaturalearthdata)
  })

  # ── Country names ──────────────────────────────────────────
  if (opt$label_countries) {
    ctry <- ne_cached(opt$label_scale, "admin_0_countries", "cultural")
    if (!is.null(ctry)) {
      # Use pre-computed label points (centroids) stored in the dataset
      ctry_pts <- suppressWarnings(sf::st_point_on_surface(ctry))
      p_labels <- c(p_labels, list(
        geom_sf_text(data      = ctry_pts,
                     aes(label = ADMIN),
                     size      = opt$label_size * 1.15,
                     colour    = opt$label_color,
                     fontface  = "bold",
                     check_overlap = TRUE,
                     inherit.aes   = FALSE)
      ))
      message("  Added country labels")
    }
  }

  # ── State / province names ─────────────────────────────────
  if (opt$label_states) {
    states <- tryCatch(
      ne_states(returnclass = "sf"),
      error = function(e) {
        message("  State labels unavailable: ", conditionMessage(e)); NULL
      }
    )
    if (!is.null(states)) {
      state_pts <- suppressWarnings(sf::st_point_on_surface(states))
      p_labels <- c(p_labels, list(
        geom_sf_text(data      = state_pts,
                     aes(label = name),
                     size      = opt$label_size * 0.75,
                     colour    = opt$label_color,
                     fontface  = "plain",
                     check_overlap = TRUE,
                     inherit.aes   = FALSE)
      ))
      message("  Added state/province labels")
    }
  }

  # ── Physical feature names (mountain ranges, deserts, plains…) ──
  if (opt$label_features) {
    # geography_regions_polys: named physical regions (Andes, Amazonia, Patagonia…)
    feat_poly <- ne_cached(opt$label_scale, "geography_regions_polys", "physical")
    if (!is.null(feat_poly)) {
      feat_pts <- suppressWarnings(sf::st_point_on_surface(feat_poly))
      p_labels <- c(p_labels, list(
        geom_sf_text(data      = feat_pts,
                     aes(label = name),
                     size      = opt$label_size * 0.85,
                     colour    = opt$label_color,
                     fontface  = "italic",    # italic = conventional for physical features
                     check_overlap = TRUE,
                     inherit.aes   = FALSE)
      ))
      message("  Added physical feature labels")
    }
    # geography_marine_polys: oceans, seas, gulfs
    marine <- ne_cached(opt$label_scale, "geography_marine_polys", "physical")
    if (!is.null(marine)) {
      marine_pts <- suppressWarnings(sf::st_point_on_surface(marine))
      p_labels <- c(p_labels, list(
        geom_sf_text(data      = marine_pts,
                     aes(label = name),
                     size      = opt$label_size * 0.8,
                     colour    = "#4d9fd6",   # blue for water bodies
                     fontface  = "italic",
                     check_overlap = TRUE,
                     inherit.aes   = FALSE)
      ))
      message("  Added marine feature labels")
    }
  }
}

# ── Shapefile overlay ─────────────────────────────────────────
p <- base_plot

if (opt$show_shapefile) {
  if (shp_has_own_legend) {
    suppressPackageStartupMessages(library(ggnewscale))
  }

  # Render ALL shapefiles with their matched/computed colors
  # The first one carries the scale_fill_identity (with optional legend);
  # subsequent ones reuse identity fill without adding another scale.
  p <- p +
    geom_sf(data = shp_list[[1]], aes(fill = .shp_color),
            colour = opt$shp_outline, alpha = opt$shp_alpha,
            linewidth = 0.5, inherit.aes = FALSE) +
    scale_fill_identity(
      guide  = if (shp_has_own_legend) "legend" else "none",
      name   = if (shp_has_own_legend) resolved_field else NULL,
      breaks = if (shp_has_own_legend) unname(shp_pal) else NULL,
      labels = if (shp_has_own_legend) names(shp_pal)  else NULL
    )
  if (shp_has_own_legend) p <- p + new_scale_fill()

  if (length(shp_list) > 1) {
    for (i in seq(2, length(shp_list))) {
      p <- p +
        geom_sf(data  = shp_list[[i]], aes(fill = .shp_color),
                colour = opt$shp_outline, alpha = opt$shp_alpha,
                linewidth = 0.5, inherit.aes = FALSE) +
        scale_fill_identity()   # identity scale; no extra legend entry
    }
  }
}

# ── Point layers ──────────────────────────────────────────────

# Background layer: non-focal taxa shown as small faded grey points
if (opt$show_other && !is.null(df_other) && nrow(df_other) > 0) {
  p <- p +
    geom_point(
      data   = df_other,
      aes(x = Lon, y = Lat),
      shape  = 21, size = opt$point_size * 0.7,
      fill   = opt$other_color, colour = "white",
      alpha  = opt$other_alpha, stroke = 0.5,
      inherit.aes = FALSE
    )
  message("  Added ", nrow(df_other), " background 'other' taxa points")
}

if (opt$status) {
  
  # Pin tail for have + sequenced
  p <- p +
    geom_point(
      data   = df_plot %>% filter(status %in% c("have","sequenced")),
      aes(x = Lon, y = Lat, fill = .pt_color),
      shape = 25, size = opt$point_size * 0.6,
      colour = "white", stroke = 0.6
    ) +
    # Circle pin head (have)
    geom_point(
      data     = df_plot %>% filter(status == "have"),
      aes(x = Lon, y = Lat, fill = .pt_color),
      shape = 21, size = opt$point_size,
      colour = "white", stroke = 1,
      position = position_nudge(y = opt$nudge)
    ) +
    # Diamond pin head (sequenced)
    geom_point(
      data     = df_plot %>% filter(status == "sequenced"),
      aes(x = Lon, y = Lat, fill = .pt_color),
      shape = 23, size = opt$point_size,
      colour = "white", stroke = 1,
      position = position_nudge(y = opt$nudge)
    ) +
    # Missing: × — white halo first, then colored × on top
    geom_point(
      data   = df_plot %>% filter(status == "missing"),
      aes(x = Lon, y = Lat),
      shape = 4, size = opt$point_size, stroke = 3.5,
      colour = "white"
    ) +
    geom_point(
      data   = df_plot %>% filter(status == "missing"),
      aes(x = Lon, y = Lat, color = .pt_color),
      shape = 4, size = opt$point_size, stroke = 2
    ) +
    # Color scale (identity) — all collections in legend
    scale_fill_identity(
      guide  = "legend",
      name   = legend_ttl,
      breaks = unname(point_pal),
      labels = names(point_pal)
    ) +
    scale_color_identity(guide = "none") +
    # Shape legend as separate manual legend using dummy data
    geom_point(
      data = data.frame(
        x      = rep(-Inf, 3),
        y      = rep(-Inf, 3),
        status = factor(c("missing","have","sequenced"),
                        levels = c("missing","have","sequenced"))
      ),
      aes(x = x, y = y, shape = status),
      fill = "grey50", colour = "white", size = opt$point_size, stroke = 1
    ) +
    scale_shape_manual(
      values = c("missing" = 4, "have" = 21, "sequenced" = 23),
      name   = "Status"
    ) +
    guides(shape = guide_legend(
      override.aes = list(
        fill   = c(NA,       "grey50", "grey50"),  # × has no fill
        colour = c("grey30", "white",  "white")    # × needs dark stroke to be visible
      )
    ))
  
} else if (!is.null(opt$shape_by)) {
  
  # Color by color_by + shape by shape_by — two simultaneous variables.
  # 20-shape pool: 5 filled (shapes 21-25, use fill aesthetic) then
  # 15 stroke/symbol shapes (use colour aesthetic) — gives far more unique options.
  if (!opt$shape_by %in% colnames(df_plot))
    stop(paste0("--shape_by '", opt$shape_by, "' not found in sample sheet.\n",
                "Available: ", paste(colnames(df_plot), collapse = ", ")))
  
  shape_vals <- sort(unique(na.omit(as.character(df[[opt$shape_by]]))))
  n_shp      <- length(shape_vals)
  
  # Ordered pool: filled first (bold, easy to read), then distinctive stroke shapes
  shape_pool <- c(
    21, 22, 23, 24, 25,      # filled: circle, square, diamond, tri-up, tri-down
    3,  4,  8,               # symbols: +  ×  *
    1,  0,  5,  2,  6,       # open: circle, square, diamond, tri-up, tri-down
    7,  9,  10, 12, 13       # compound open shapes
  )
  if (n_shp > length(shape_pool))
    warning("More collections (", n_shp, ") than unique shapes (",
            length(shape_pool), ") — some shapes will repeat")
  
  shape_vals_used <- shape_pool[(seq_len(n_shp) - 1L) %% length(shape_pool) + 1L]
  shape_map  <- setNames(shape_vals_used, shape_vals)
  shp_ttl    <- if (!is.null(opt$shape_legend_title)) opt$shape_legend_title else opt$shape_by
  
  df_plot <- df_plot %>%
    mutate(.pt_shape = as.character(.data[[opt$shape_by]]))
  
  filled_colls <- names(shape_map)[shape_map >= 21]
  stroke_colls <- names(shape_map)[shape_map <  21]
  
  # Layer A — filled shapes (21-25): species color in fill, white border for contrast
  if (length(filled_colls) > 0) {
    df_f <- df_plot %>% filter(.pt_shape %in% filled_colls)
    if (nrow(df_f) > 0)
      p <- p + geom_point(
        data   = df_f,
        aes(x = Lon, y = Lat, fill = .pt_color, shape = .pt_shape),
        size   = opt$point_size, colour = "white", stroke = 1.2,
        inherit.aes = FALSE
      )
  }
  
  # Layer B — stroke/symbol shapes (<21): white halo first, then colored stroke on top.
  # The halo gives contrast against busy raster basemaps.
  if (length(stroke_colls) > 0) {
    df_s <- df_plot %>% filter(.pt_shape %in% stroke_colls)
    if (nrow(df_s) > 0) {
      # Halo: same shape, wider white stroke behind
      p <- p + geom_point(
        data    = df_s,
        aes(x = Lon, y = Lat, shape = .pt_shape),
        size    = opt$point_size, stroke = 4,
        colour  = "white",
        inherit.aes = FALSE
      )
      # Colored stroke on top
      p <- p + geom_point(
        data    = df_s,
        aes(x = Lon, y = Lat, colour = .pt_color, shape = .pt_shape),
        size    = opt$point_size, stroke = 2,
        inherit.aes = FALSE
      )
    }
  }
  
  p <- p +
    scale_fill_identity(
      guide  = "legend",
      name   = legend_ttl,
      breaks = unname(point_pal),
      labels = names(point_pal)
    ) +
    scale_colour_identity(guide = "none") +  # colour used for stroke shapes; no dup legend
    scale_shape_manual(
      values = shape_map,
      name   = shp_ttl,
      breaks = shape_vals,
      labels = shape_vals
    ) +
    guides(
      fill  = guide_legend(order = 1,
                           override.aes = list(shape = 21, size = opt$point_size * 0.8)),
      shape = guide_legend(order = 2,
                           override.aes = list(fill   = "grey60",
                                               colour = "grey20",
                                               size   = opt$point_size * 0.8))
    )

} else {
  
  # No status, no shape_by — pin shape, color legend only
  p <- p +
    geom_point(
      data   = df_plot,
      aes(x = Lon, y = Lat, fill = .pt_color),
      shape = 25, size = opt$point_size * 0.6,
      colour = "white", stroke = 0.6
    ) +
    geom_point(
      data     = df_plot,
      aes(x = Lon, y = Lat, fill = .pt_color),
      shape    = 21, size = opt$point_size,
      colour   = "white", stroke = 1,
      position = position_nudge(y = opt$nudge)
    ) +
    scale_fill_identity(
      guide  = "legend",
      name   = legend_ttl,
      breaks = unname(point_pal),
      labels = names(point_pal)
    )
}

# ── Place labels (on top of everything) ──────────────────────
if (length(p_labels) > 0) {
  for (lyr in p_labels) p <- p + lyr
}

# ── Labels + theme ────────────────────────────────────────────
# xlim/ylim are in WGS84 degrees; crs sets the *display* projection to match
# the shapefile's native CRS. coord_sf reprojects all layers automatically.
p <- p +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax),
           crs = shp_crs, default_crs = sf::st_crs(4326),
           expand = FALSE) +
  labs(
    title    = if (is.null(opt$title) || opt$title == "") NULL else opt$title,
    subtitle = if (opt$subtitle == "") NULL else opt$subtitle,
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", size = 15, hjust = 0.5,
                                   margin = margin(b = 4)),
    plot.subtitle   = element_text(color = "grey40", size = 11, hjust = 0.5,
                                   margin = margin(b = 8)),
    legend.position = "right",
    legend.key.size = unit(1, "lines"),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(10, 10, 10, 10)
  )

# ── Save ──────────────────────────────────────────────────────
message("Saving ", opt$output, "...")
ggsave(filename = opt$output, plot = p,
       width = opt$width, height = opt$height,
       dpi = opt$dpi, bg = "white")
message("Done: ", opt$output)
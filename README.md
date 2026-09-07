# Specimen Mapping Toolkit (R)

R scripts that pull specimen records from GBIF, prepare an IUCN range
polygon, and map them:

1. `gbif_population_sampling.R` — search GBIF → filtered, map-ready sample sheet.
   - `gbif_phylo_sampling.R` — variant: one specimen per taxon for phylogenetics.
2. `merge_iucn_ranges.R` — merge/filter IUCN range downloads → range shapefile.
3. `specimen_map.R` — map the sample sheet (one self-contained mapping script).

## Contents
```
specimen-mapping-toolkit/
├── README.md
├── gbif_population_sampling.R      # GBIF search → sample sheet (+centroid CSV, summary, range shp)
├── gbif_phylo_sampling.R       # one-per-taxon phylo sampler (tissue-first, keep-list, include)
├── merge_iucn_ranges.R         # merge & filter IUCN range zips → shapefile/GeoPackage
├── specimen_map.R              # map a sample sheet (cbhypso basemap, optional shapefile)
└── examples/
    ├── list_template.csv       # --include must-include list template
    └── crocodylus.list         # --keep_taxa_file example (extant Crocodylus)
```

## Requirements
- R (≥ 4.1). Each script is a standalone `Rscript` with `--help`.
- Install the packages once:
```r
install.packages(c(
  "optparse","rgbif","dplyr","readr","stringr",   # search + phylo
  "sf","rnaturalearth","rnaturalearthdata",        # geocoding, IUCN merge, borders
  "ggplot2","RColorBrewer","viridis","ggnewscale", # mapping
  "terra","tidyterra","maptiles"                    # basemaps (cbhypso / ne)
), repos = "https://cloud.r-project.org")
# province-level geocoding also wants: rnaturalearthhires
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev")
```
- Network access to GBIF (search) and, on first use, to Natural Earth
  (basemap/geocoding downloads, cached under `~/.cache/specimen_map`).

## Quick start
```bash
# 1. specimens (one per taxon, tissue-first, priority collections)
Rscript gbif_phylo_sampling.R --taxon Crocodylus \
  --keep_taxa_file examples/crocodylus.list \
  --collections "FMNH,AMNH,UF,LSUMZ,MCZ,USNM" --collection_first \
  --preparations any --country_centroid \
  --output croc_phylo.csv

# 2. (optional) IUCN range polygon for an overlay
Rscript merge_iucn_ranges.R --input iucn_zips/ --dissolve all \
  --output range/croc_range.shp

# 3. map (buffers to the points, cross-blend hypso basemap)
Rscript specimen_map.R --samples croc_phylo_centroid.csv \
  --color_by taxon --shape_by Collection --output croc_map.png
```

## 1. `gbif_population_sampling.R` — search GBIF, write a map-ready sheet

Queries the GBIF occurrence API (via **rgbif**) and writes a CSV whose
columns `specimen_map.R` auto-detects (`ID`, `Lat`, `Lon`, `species`,
`subspecies`, `Collection`, …).

Packages: `optparse`, `rgbif`, `dplyr`, `readr`, `stringr`;
`--geocode`/`--range_shapefile` also need `sf` + `rnaturalearth`
(install `rnaturalearthhires` for province-level centroids).

Full run — ostrich native range, coordinate-less records placed at
province/country centroids, **minimized** to the smallest representative
panel, plus a range shapefile for the map extent:

```bash
Rscript gbif_population_sampling.R \
  --taxon "Struthio camelus,Struthio molybdophanes" \
  --preset ostrich-native \
  --institutions "FMNH,AMNH,NHMUK,MVZ,MSU,MNHN" \
  --preparations "tissue,skin" \
  --geocode state \
  --minimize \
  --priority "FMNH,AMNH,NHMUK,MVZ,MSU,MNHN" \
  --range_shapefile range/ostrich_native_range.shp \
  --output_full ostrich_gbif_full.csv \
  --output ostrich_gbif_min.csv
```
(`--institutions` limits the CSV to those museums, `--preparations` to
tissue/skin, `--preset` to the native-range countries.)

### Minimization (`--minimize`)
Keeps the fewest specimens that still cover the diversity **without dropping
unique localities**: one best record per **stratum**, where a stratum
defaults to `taxon × countryCode × stateProvince × locality` (so every
distinct locality is retained; only exact locality duplicates collapse).
Within each stratum records are ranked by

1. **preparation preference** — `--preparations` order (tissue before skin
   by default, so you get the most sequenceable material);
2. **collection priority** — `--collection_priority` order (so you request
   from collections you can actually access);
3. **coordinate quality** — exact > province > country centroid;
4. **has a date**;
5. smaller coordinate uncertainty.

Tune it:
- `--strata "taxon,countryCode"` — coarser coverage unit (fewer samples),
  or `--strata taxon` for the absolute minimum (one per taxon). The
  default already keeps one per unique locality.
- `--per_stratum 2` — keep N per cell (e.g. replicates).
- `--require_coords` — only choose records that will actually plot.
- `--output_full` — also write the complete, un-minimized set.

### Limiting the CSV (collections, regions, preparations)
The written sheet is restricted to exactly what you filtered for:
- **Regions** — `--countries KE,TZ,…` or `--preset ostrich-native`.
- **Collections** — `--institutions "FMNH,AMNH,NHMUK,MVZ,MSU,MNHN"` keeps
  only those museums (case-insensitive, any number of them).
- **Preparations** — `--preparations "tissue,liquid nitrogen,skin"`
  (default) keeps only records with tissue, liquid-nitrogen, or any skin
  variant (`study skin`, `Mounted skin`, `flat skin & complete skeleton`,
  `frozen tissue`, …), in preference order. The **none / unrecorded**
  bucket is kept by default too (`--drop_unknown_prep` to exclude it);
  eggs, skeletons, feathers, nests, mounts, etc. are dropped.
  `--preparations any` keeps every preparation type (still ranked
  tissue-first). **Non-bird/mammal groups** (reptiles, amphibians, fish) are
  usually preserved as *whole animal / alcohol / skeleton*, not "skin", so
  use `--preparations any` for them — otherwise those specimens are dropped.
- **Institution quirks** — `--exclude "CODE:regex;CODE:regex"` drops
  records from a collection whose **catalogue number** matches a regex.
  Off by default; e.g. `--exclude "NHMUK:E"` removes NHMUK eggs (an `E` in
  the reg number), or `--exclude "NHMUK:E;AMNH:^EGG"` for several rules.

Other flags: `--taxon`/`--taxon_key`, `--basis`,
`--has_coordinate {any,true,false}`, `--no_geospatial_issue`,
`--geocode {none,country,state}`, `--no_jitter`, `--limit_total`.
Records are deduped by `institutionCode` + `catalogNumber`. Every row
carries both `coord_precision` (`exact|state|country`) and a readable
`coord_source` — **GPS (recorded)** for a real specimen coordinate vs
**province centroid** / **country centroid** for a geocoded one (or
**none**) — so you can tell actual GPS from a centroid at a glance.
Taxon-agnostic — works for your parrots too.

**`--country_centroid`** (both GBIF scripts) writes a **separate** CSV
(`<output>_centroid.csv`, or `--centroid_output`) identical to the main
sheet except specimens with no GPS are filled with their country centroid
(tagged `coord_source = country centroid`) so the whole set can be mapped.
The main `--output` is left GPS-only. This is distinct from `--geocode`,
which fills coordinates in place in the main output.

### Must-include specimens (`--include`)
Both GBIF scripts accept `--include list.csv` — a CSV of specimens you want
in the output no matter what. Columns (case-insensitive headers):
- **institution** (`institution`/`collection`) — required
- **catalogue** (`catalogue`/`catalog`/`catalog_number`) — required
- **species** (`species`) — a bare epithet (`niloticus`) or a full binomial
  (`Crocodylus niloticus`); a bare epithet gets the genus from `--taxon`
- **subspecies** (`subspecies`) — optional
- (or a single **taxon**/`scientificName` column instead of species+subspecies)

After the search/selection, any listed specimen not already present is
**looked up on GBIF** by institution + catalogue to fill every column; the
species is used to pick the right record when a catalogue number is reused
across a museum's collections, and if GBIF has no matching record the row is
added from the fields you gave. Rows are flagged `forced_include = TRUE` and,
with `--restrict_taxon` on, are still held to the searched genus.

**With `--include`, the output is the include list by default** — the
automatic genus-wide selection is *not* added. If you also give a species
list (`--keep_taxa`/`--keep_taxa_file`), any target species **not** covered
by an include is **gap-filled** with one auto-picked specimen, chosen from
the priority `--collections` in order, then by preparation. So: your pinned
vouchers, plus one specimen per remaining target species from your
preferred institutions.
- `--include_only` — strictly the include specimens, no gap-fill even with a
  species list.
- `--include_add_all` — also append the full automatic one-per-taxon genus
  selection (the older behavior).

Example `list.csv` (species column, subspecies optional):
```csv
institution,species,subspecies,catalogue
MCZ,siamensis,,R-17574
MCZ,niloticus,,R-54121
FMNH,palustris,,31538
FMNH,niloticus,africanus,00000
```

### Summary counts (always produced)
The script prints counts **per institution**, **per taxon** (species +
subspecies), **per preparation**, and **per coordinate source** (GPS vs
centroid) to stdout, and writes the same as a tidy TSV (`category`,
`value`, `n`) — default `<output>_summary.tsv`, or
set `--summary path.tsv`. Counts describe the written sheet (the minimized
panel when `--minimize`). A specimen with several preparations
(`skin; skeleton`) is counted once per type, so preparation totals are
prep-instances.

### Phylogenetic sampling variant — `gbif_phylo_sampling.R`
A specialised sibling that returns **one specimen per taxon** (species or
subspecies) for building a phylogeny. It prefers **tissue over skin**
(falling back to skin, then anything), searches collections in a
**priority cascade**, excludes eggs, and covers **synonyms** (it resolves
each name to GBIF's accepted key, under which synonym records are indexed).

```bash
Rscript gbif_phylo_sampling.R \
  --taxon "Struthio camelus,Struthio molybdophanes" \
  --collections "FMNH,ANSP" \
  --output ostrich_phylo.csv
```

Per taxon it picks the best specimen ranked by: **collection order**
(`--collections`, always searched in the listed order) → material
(tissue > liquid nitrogen > skin > other) → coordinate quality → has a
date → smaller uncertainty. Flags:
- `--collections "FMNH,ANSP"` — the institution search order, always applied
  as the primary ranking. **By default only these collections are used**
  (`--restrict_collections`, now the default).
- `--collection_first` — search the listed collections in order first, then
  **fall back to other collections** not in the list (ranked after them).
  `--no_restrict_collections` does the same (consider others after the list).
- `--keep_taxa "Crocodylus niloticus,Crocodylus suchus"` (or
  `--keep_taxa_file taxa.txt`) — keep only these taxa and drop the rest.
  Each name is resolved through the GBIF backbone and its **synonyms** are
  pulled, so records published under a synonym are kept too. Handy when you
  search a whole genus but want only certain species. (Naming a subspecies
  keeps just that subspecies; naming a species keeps all its subspecies.)
- `--restrict_taxon` (default **on**) keeps only records whose genus matches
  `--taxon`/`--keep_taxa`, so the output is strictly that group. It also
  guards `--include`: a catalogue number reused across a museum's other
  collections (e.g. `FMNH 31538` is a bat *and* a crocodile) resolves to the
  match whose genus fits the list's taxon, and an off-target hit is replaced
  by a stub from the listed taxon rather than the wrong organism.
  `--no_restrict_taxon` restores full `--include` bypass.
- `--rank_level taxon|species|subspecies` (default `taxon`), `--per_taxon N`.
- `--preparations "tissue,liquid nitrogen,skin"` (preference & filter),
  `--drop_unknown_prep` to forbid no-prep specimens as a last resort.
- Eggs are excluded by default (by egg collectionCode/preparation);
  `--no_exclude_eggs` to keep them. Add `--exclude "NHMUK:E"` for the
  NHMUK egg-catalogue rule (off by default).
- Also takes `--countries`/`--preset`, `--has_coordinate`, `--geocode`.

Writes the same column schema (incl. `coord_source`) plus a summary TSV
with counts per institution / taxon / material / coordinate source. Feed
its CSV straight into `specimen_map.R`.

## 2. `merge_iucn_ranges.R` — merge & filter IUCN range maps

Turns one or more IUCN Red List range downloads into a single range layer
for `--shapefile` (map extent and/or overlay). Unzips `.zip`s, merges all
layers, and keeps only the ranges you want by IUCN `presence`/`origin`.

Packages: `optparse`, `sf`, `dplyr`.

```bash
# Merge several species zips, keep extant + historical native range,
# drop introduced/vagrant, dissolve to one polygon for the map extent:
Rscript merge_iucn_ranges.R \
  --input "Struthio_camelus.zip,Struthio_molybdophanes.zip" \
  --dissolve all \
  --output range/ostrich_iucn_range.shp
```

Defaults: `--presence extant,historical` (codes 1,2,3,4,5) and
`--origin native,reintroduced,uncertain` (codes 1,2,5) — so **introduced
(3), vagrant (4), and assisted-colonisation (6) are dropped**, keeping
extant and historical native range. Override with keywords or raw codes:

- `--presence extant` (1,2,3 only) or `--presence "1,5"`;
- `--origin native` (strictly code 1);
- `--seasonal all` (default) or e.g. `breeding,resident`;
- `--dissolve none|species|all`;
- `--keep_na_codes` to keep features whose code is blank.

IUCN codes — presence: 1 extant · 2 probably extant · 3 possibly extant ·
4 possibly extinct · 5 extinct · 6 uncertain. origin: 1 native ·
2 reintroduced · 3 introduced · 4 vagrant · 5 uncertain · 6 assisted.
Output is `.gpkg` or `.shp` (both work with `specimen_map.R`). Use this in
place of the country-based `--range_shapefile` from script 1 when you have
real IUCN polygons.

## 3. `specimen_map.R` — the mapping script

One self-contained mapping script (flag-driven). Feed it a sample sheet;
a shapefile is optional.

```bash
# simplest — extent from the sample points, cross-blend hypso basemap:
Rscript specimen_map.R --samples ostrich_gbif_min.csv \
  --color_by taxon --shape_by Collection --output ostrich_map.png

# with an IUCN range polygon (fixed extent + overlay):
Rscript specimen_map.R --samples ostrich_gbif_min.csv \
  --shapefile range/ostrich_iucn_range.shp \
  --color_by taxon --shape_by Collection --output ostrich_map.png
```

- **Extent buffers to the points by default** (`--extent points`): the map
  frames to the sample points' bounding box plus `--buffer` (default
  5000 km), so points are never clipped even when a large range polygon is
  overlaid. Use `--extent shapefile` to frame to the shapefile instead, or
  `--extent both` for the union of the two.
- **`--shapefile` is optional** — without it there's no range overlay and
  the extent is the points (same as `--extent points`). Provide one to
  overlay the range (and, with `--extent shapefile`/`both`, to frame by it).
- **`--whole_map`** — show the entire world instead of framing to the data
  (ignores `--extent`/`--buffer`); useful for a global overview.
- **No figure title by default** (`--title ""`); the **colour legend title
  auto-labels** as *Species*, *Subspecies*, *Populations*, or *Species /
  subspecies* (for `--color_by taxon`) based on the level you colour by.
  Override either with `--title "…"` / `--legend_title "…"`.
- **Basemap** defaults to Natural Earth's 1:50m **Cross-blended
  Hypsometric Tints** (`--basemap cbhypso`), downloaded once and cached
  under `~/.cache/specimen_map` (`--cache_dir`; source `--hypso_url`).
  `--basemap ne` uses the lighter tiled basemap; `--basemap /path.tif` any
  raster/vector.
- **Colour key** — the sheet has `species`, `subspecies`, and a combined
  `taxon` column: `--color_by species` (or `subspecies`, `population`,
  `taxon`) makes the point colour key; its title auto-labels as *Species*,
  *Subspecies*, or *Populations*. `--na_color` colours any value-less points.
- **Range shapefiles matched to the key** — `--range ranges.shp` overlays
  range polygon(s) coloured to **match the same key** (one shared legend for
  points and ranges). Give a shapefile with a `species`/`subspecies`/
  taxon field (auto-detected, or set `--shp_field`); comma-separate several.
  It doesn't change the extent. (Equivalent to `--shapefile … --match`.)
- **Shape** — `--shape_by Collection` encodes the museum; omit to disable.
- **Highlight/subset** — `--taxa "FMNH,AMNH,…" --taxa_by Collection`
  (with `--show_other` to add the rest faded); `--status` or `--shape_by`
  for a second variable; `--rivers`, `--labels`, `--title`, `--width/
  --height/--dpi`, etc. Run `--help` for the full list.

Map the **full** set instead of the minimized panel by pointing
`--samples` at `ostrich_gbif_full.csv`.

Requires R with: optparse, sf, ggplot2, dplyr, readr, RColorBrewer,
viridis, ggnewscale; for the basemap: terra, tidyterra (cbhypso / raster),
maptiles (`ne`), rnaturalearth/rnaturalearthdata (borders, rivers, labels).

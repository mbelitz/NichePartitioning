# =============================================================================
# 00_buildDiversityMetrics.R
#
# Builds the site-species matrix and calculates:
#   - Chao1 estimated species richness
#   - Faith's Phylogenetic Diversity (PD)
#   - Abundance-weighted Mean Pairwise Distance (MPD)
#   - Standardized effect size of MPD (SES-MPD)
#
# Outputs:
#   data/derivedData/site_species_matrix.rds
#   data/derivedData/rich.csv          (site metadata + Chao1)
#   data/derivedData/PD_MPD.csv        (PD, SR, MPD per site)
#   data/derivedData/ses_mpd.csv       (SES-MPD null model results)
# =============================================================================

library(tidyverse)
library(vegan)
library(ape)
library(picante)
library(phytools)
library(terra)
library(sf)

# Island sites excluded from all analyses (non-continental dynamics)
EXCLUDE_SITES <- c("LAJA", "GUAN", "PUUM")


# -----------------------------------------------------------------------------
# 1. Load and filter carabid pitfall data
# -----------------------------------------------------------------------------

df <- read.csv("data/neon.df.csv") %>%
  mutate(year = year(collectDate)) %>%
  filter(sampleType %in% c("carabid", "other carabid"))


# -----------------------------------------------------------------------------
# 2. Site coordinates from NEON distributed base plot shapefiles
# -----------------------------------------------------------------------------
# note these shapefiles are ignored in repository due to large file size; 
# must be downloaded for full reproducibility of this script
plots_sf <- vect("data/All_NEON_TOS_Plots_V11/All_NEON_TOS_Plot_Polygons_V11.shp") %>% 
  st_as_sf()

dbp <- plots_sf %>%
  filter(plotType == "distributed", subtype == "basePlot")

dbp_coords <- st_coordinates(st_centroid(dbp))

dbp <- dbp %>%
  mutate(lon = dbp_coords[, 1],
         lat = dbp_coords[, 2])

df <- left_join(df, dbp)


# -----------------------------------------------------------------------------
# 3. Site-level coordinate metadata
# -----------------------------------------------------------------------------

site_meta <- st_drop_geometry(df) %>%
  group_by(siteID) %>%
  summarise(
    meanLat = mean(latitude,  na.rm = TRUE),
    meanLon = mean(longitude, na.rm = TRUE),
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 4. Site-species matrix
# -----------------------------------------------------------------------------
# Retain only species-level identifications and sum total captures per
# site × taxon across all years and sampling events.

species_counts <- df %>%
  filter(!siteID %in% EXCLUDE_SITES,
         taxonRank == "species") %>%
  group_by(siteID, taxonID) %>%
  summarise(totalCount = sum(individualCount, na.rm = TRUE), .groups = "drop")

site_species_matrix <- species_counts %>%
  pivot_wider(names_from  = taxonID,
              values_from = totalCount,
              values_fill = 0) %>%
  column_to_rownames("siteID")

saveRDS(site_species_matrix, "data/derivedData/site_species_matrix.rds")


# -----------------------------------------------------------------------------
# 5. Chao1 estimated richness
# -----------------------------------------------------------------------------
# estimateR returns S.obs, S.chao1, se.chao1, S.ACE, se.ACE as row names;
# transpose so sites are rows.

richness_est <- estimateR(site_species_matrix) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("siteID")

rich <- site_meta %>%
  filter(!siteID %in% EXCLUDE_SITES) %>%
  left_join(richness_est, by = "siteID")

write.csv(rich, "data/derivedData/rich.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 6. Phylogenetic tree preparation
# -----------------------------------------------------------------------------

tree <- read.tree("data/phylogeny/roughTree.DNA.droppedseqs.newick")

# Build taxon lookup: taxonID <-> tree tip label (matched on genus_species)
tax_lookup <- df %>%
  filter(taxonRank == "species") %>%
  select(taxonID, scientificName) %>%
  distinct() %>%
  mutate(genus_species = paste(word(scientificName, 1),
                               word(scientificName, 2),
                               sep = "_"))

tree_taxa <- tibble(tip = tree$tip.label) %>%
  mutate(genus_species = str_extract(tip, "^[^_]+_[^_]+"))

tax_lookup <- tax_lookup %>%
  inner_join(tree_taxa, by = "genus_species")

# Subset SSM to taxa present in the tree and rename columns to tip labels
ssm_tree <- site_species_matrix %>%
  select(any_of(tax_lookup$taxonID))

colnames(ssm_tree) <- tax_lookup$tip[match(colnames(ssm_tree), tax_lookup$taxonID)]

# Prune tree to taxa represented in the SSM, then midpoint-root
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, colnames(ssm_tree)))
tree_rooted <- midpoint.root(tree_pruned)

stopifnot(all(colnames(ssm_tree) %in% tree_rooted$tip.label))


# -----------------------------------------------------------------------------
# 7. Faith's PD, abundance-weighted MPD, and SES-MPD
# -----------------------------------------------------------------------------

ssm_pa <- (ssm_tree > 0) * 1  # presence/absence matrix for PD

pd_res  <- pd(ssm_pa, tree_rooted, include.root = TRUE)
mpd_res <- mpd(ssm_tree, cophenetic(tree_pruned), abundance.weighted = TRUE)

ses_mpd_res <- ses.mpd(
  ssm_pa,
  cophenetic(tree_rooted),
  null.model = "taxa.labels",
  runs       = 999
)

pd_df <- data.frame(
  siteID = rownames(ssm_pa),
  PD     = pd_res$PD,
  SR     = pd_res$SR,
  MPD    = mpd_res
)

write.csv(pd_df,       "data/derivedData/PD_MPD.csv",  row.names = FALSE)
write.csv(ses_mpd_res, "data/derivedData/ses_mpd.csv", row.names = FALSE)
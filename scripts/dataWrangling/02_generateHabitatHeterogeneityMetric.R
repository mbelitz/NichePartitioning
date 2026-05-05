# =============================================================================
# 01_buildHabitatHeterogeneity.R
#
# Calculates two site-level plant community metrics, both restricted to
# distributed base plots where beetles were also sampled:
#
#   (a) Within-site plant beta diversity (H_plant_bc): mean pairwise Jaccard
#       dissimilarity between plots — used as a covariate for habitat
#       heterogeneity in the turnover and co-occurrence models (scripts 04-05).
#
#   (b) Between-site habitat position (habitat_PC1): first axis of a PCoA on
#       between-site Jaccard distances, using site-level pooled presence/absence.
#       Used in script 06 as the habitat gradient along which each species'
#       niche position is calculated (abundance-weighted mean PC1 across sites
#       of capture). PC1 sign is arbitrary; its ecological interpretation
#       depends on which NEON sites load at each extreme.
#
# Outputs:
#   data/derivedData/plant_heterogeneity.csv      (siteID, H_plant_bc, n_plots)
#   data/derivedData/site_habitat_position.csv    (siteID, habitat_PC1)
# =============================================================================

library(neonDivData)
library(vegan)
library(ape)
library(tidyverse)

# Island sites excluded from all analyses (non-continental dynamics)
EXCLUDE_SITES <- c("LAJA", "GUAN", "PUUM")


# -----------------------------------------------------------------------------
# 1. Get plot IDs where beetles were sampled
# -----------------------------------------------------------------------------

df <- read.csv("data/neon.df.csv") %>%
  mutate(year = year(collectDate)) %>%
  filter(sampleType %in% c("carabid", "other carabid"),
         !siteID %in% EXCLUDE_SITES)

beetle_plots <- df %>%
  distinct(siteID, plotID)


# -----------------------------------------------------------------------------
# 2. Load plant data and restrict to beetle-sampled plots
# -----------------------------------------------------------------------------

plantData <- neonDivData::data_plant

plant_filtered <- plantData %>%
  filter(presence_absence == 1,
         taxon_rank == "species",
         plotID %in% beetle_plots$plotID)


# -----------------------------------------------------------------------------
# 3. Calculate mean pairwise Jaccard dissimilarity per site
# -----------------------------------------------------------------------------

plant_hetero <- plant_filtered %>%
  group_by(siteID, plotID, taxon_id) %>%
  summarise(present = 1, .groups = "drop") %>%
  pivot_wider(names_from  = taxon_id,
              values_from = present,
              values_fill = 0) %>%
  group_by(siteID) %>%
  group_split() %>%
  map_dfr(function(x) {
    
    site <- unique(x$siteID)
    
    mat <- x %>%
      select(-siteID) %>%
      column_to_rownames("plotID")
    
    # Need at least 2 plots to compute dissimilarity
    if (nrow(mat) < 2) {
      return(tibble(siteID = site, H_plant_bc = NA_real_, n_plots = nrow(mat)))
    }
    
    mat <- mat[, colSums(mat) > 0, drop = FALSE]
    bc  <- vegdist(mat, method = "jaccard")
    
    tibble(
      siteID     = site,
      H_plant_bc = mean(bc, na.rm = TRUE),
      n_plots    = nrow(mat)
    )
  })


# -----------------------------------------------------------------------------
# 4. Diagnostic: how many sites lost plots relative to the full plant dataset?
# -----------------------------------------------------------------------------

plant_hetero_all <- plantData %>%
  filter(presence_absence == 1, taxon_rank == "species") %>%
  group_by(siteID, plotID, taxon_id) %>%
  summarise(present = 1, .groups = "drop") %>%
  pivot_wider(names_from  = taxon_id,
              values_from = present,
              values_fill = 0) %>%
  group_by(siteID) %>%
  summarise(n_plots_all = n(), .groups = "drop")

comparison <- plant_hetero %>%
  left_join(plant_hetero_all, by = "siteID") %>%
  mutate(plots_dropped = n_plots_all - n_plots)

message("Sites with fewer plots after beetle-plot filter:")
print(filter(comparison, plots_dropped > 0))


# -----------------------------------------------------------------------------
# 5. Save H_plant_bc
# -----------------------------------------------------------------------------

write.csv(plant_hetero,
          "data/derivedData/plant_heterogeneity.csv",
          row.names = FALSE)

message("H_plant_bc written to data/derivedData/plant_heterogeneity.csv")


# -----------------------------------------------------------------------------
# 6. Between-site PCoA: continuous habitat position gradient
# -----------------------------------------------------------------------------
# Pools plant presence/absence to the site level (a species is present at a
# site if it occurred in any beetle-sampled plot), then runs PCoA on
# between-site Jaccard distances. PC1 captures the dominant axis of
# compositional variation across NEON continental sites (typically a
# forest <-> open-ground or mesic <-> xeric gradient).
#
# NOTE: PC1 sign is arbitrary. Check which sites load at each extreme to
# interpret direction before reporting.

site_plant_mat <- plant_filtered %>%
  group_by(siteID, taxon_id) %>%
  summarise(present = 1L, .groups = "drop") %>%
  pivot_wider(names_from  = taxon_id,
              values_from = present,
              values_fill = 0L) %>%
  column_to_rownames("siteID")

# Drop any all-zero species columns (shouldn't occur given presence_absence == 1 filter)
site_plant_mat <- site_plant_mat[, colSums(site_plant_mat) > 0, drop = FALSE]

if (nrow(site_plant_mat) < 3) {
  stop("Fewer than 3 sites with plant data — cannot compute between-site PCoA.")
}

site_jac  <- vegdist(site_plant_mat, method = "jaccard")
pcoa_res  <- pcoa(site_jac)

pct_pc1 <- round(pcoa_res$values$Relative_eig[1] * 100, 1)
pct_pc2 <- round(pcoa_res$values$Relative_eig[2] * 100, 1)
message(sprintf("Between-site plant PCoA: PC1 = %.1f%%, PC2 = %.1f%% of variance",
                pct_pc1, pct_pc2))

# Show which sites anchor each end of PC1 (aids interpretation)
site_hab_pos <- tibble(
  siteID      = rownames(pcoa_res$vectors),
  habitat_PC1 = pcoa_res$vectors[, 1]
)

message("PC1 extremes:")
print(arrange(site_hab_pos, habitat_PC1) %>% slice(c(1:3, (n()-2):n())))

write.csv(site_hab_pos,
          "data/derivedData/site_habitat_position.csv",
          row.names = FALSE)

message("Done. Habitat position written to data/derivedData/site_habitat_position.csv")
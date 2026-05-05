# =============================================================================
# 02_buildSpatialTurnover.R
#
# Calculates mean pairwise Bray-Curtis dissimilarity between distributed
# plots within each NEON site as a measure of spatial community turnover
# (i.e., spatial niche partitioning). Counts are pooled across all years
# and sampling events before computing dissimilarity.
#
# Prerequisites (must be run first):
#   00_buildDiversityMetrics.R  → data/derivedData/rich.csv, PD_MPD.csv
#   01_buildHabitatHeterogeneity.R → data/derivedData/plant_heterogeneity.csv
#
# NOTE: data/derivedData/ must exist before running. Create it with:
#   dir.create("data/derivedData", showWarnings = FALSE)
#
# Outputs:
#   data/derivedData/spatial_turnover.csv
# =============================================================================

library(vegan)
library(tidyverse)
library(sf)
library(terra)

# Island sites excluded from all analyses
EXCLUDE_SITES <- c("LAJA", "GUAN", "PUUM")


# -----------------------------------------------------------------------------
# 1. Load and filter carabid pitfall data
# -----------------------------------------------------------------------------

df <- read.csv("data/neon.df.csv") %>%
  mutate(year = year(collectDate)) %>%
  filter(sampleType %in% c("carabid", "other carabid"),
         !siteID %in% EXCLUDE_SITES)


# -----------------------------------------------------------------------------
# 2. Attach plot-level coordinates from NEON distributed base plot shapefiles
# -----------------------------------------------------------------------------

# NOTE: shapefiles are not tracked in the repository due to file size;
# must be downloaded separately for full reproducibility.
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
# 3. Aggregate counts by site × plot × taxon (pooled across years and weeks)
# -----------------------------------------------------------------------------

plot_counts <- df %>%
  filter(taxonRank == "species") %>%
  group_by(siteID, plotID, taxonID) %>%
  summarise(totalCount = sum(individualCount, na.rm = TRUE), .groups = "drop")


# -----------------------------------------------------------------------------
# 4. Spatial beta diversity function
# -----------------------------------------------------------------------------

calc_spatial_beta <- function(site_df) {
  
  if (is.null(site_df) || nrow(site_df) == 0) return(NULL)
  
  # Build plot × species matrix
  plot_matrix <- site_df %>%
    select(plotID, taxonID, totalCount) %>%
    pivot_wider(names_from  = taxonID,
                values_from = totalCount,
                values_fill = 0) %>%
    arrange(plotID) %>%
    column_to_rownames("plotID")
  
  # Require at least 3 plots
  if (nrow(plot_matrix) < 3) return(NULL)
  
  # Drop species absent at all plots in this site
  plot_matrix <- plot_matrix[, colSums(plot_matrix) > 0, drop = FALSE]
  if (ncol(plot_matrix) == 0) return(NULL)
  
  # Hellinger transformation then Bray-Curtis dissimilarity
  hel_matrix <- decostand(plot_matrix, method = "hellinger")
  bc_dist    <- vegdist(hel_matrix, method = "bray")
  
  all_pairs <- as.vector(bc_dist)
  
  tibble(
    mean_spatial_turnover = mean(all_pairs, na.rm = TRUE),
    sd_spatial_turnover   = sd(all_pairs,   na.rm = TRUE),
    n_plots               = nrow(plot_matrix)
  )
}


# -----------------------------------------------------------------------------
# 5. Apply across sites
# -----------------------------------------------------------------------------

spatial_beta <- plot_counts %>%
  group_by(siteID) %>%
  group_split() %>%
  map_dfr(function(x) {
    result <- calc_spatial_beta(x)
    if (!is.null(result)) mutate(result, siteID = unique(x$siteID))
  })

message(sprintf("Spatial turnover calculated for %d of %d sites.",
                nrow(spatial_beta),
                n_distinct(plot_counts$siteID)))


# -----------------------------------------------------------------------------
# 6. Join diversity metrics, habitat heterogeneity, and season length
# -----------------------------------------------------------------------------

rich    <- read.csv("data/derivedData/rich.csv")
pd      <- read.csv("data/derivedData/PD_MPD.csv")
plant_h <- read.csv("data/derivedData/plant_heterogeneity.csv")

sl <- read.csv("data/spatial.partioning.nov.csv") %>% 
  distinct(siteID, seasonLength)

mdf <- rich %>%
  left_join(pd,           by = "siteID") %>%
  left_join(spatial_beta, by = "siteID") %>%
  left_join(plant_h,      by = "siteID") %>%
  left_join(sl,           by = "siteID")


# -----------------------------------------------------------------------------
# 7. Save output
# -----------------------------------------------------------------------------

write.csv(mdf, "data/derivedData/spatial_turnover.csv", row.names = FALSE)

message("Done. Output written to data/derivedData/spatial_turnover.csv")

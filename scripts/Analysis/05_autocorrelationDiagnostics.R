# =============================================================================
# 08_diagnosticsAutocorrelation.R
#
# Tests for spatial autocorrelation in residuals of the four main brms models:
#   fit_temporal_turnover   (script 04)
#   fit_spatial_turnover    (script 04)
#   fit_diversity_richness  (script 07)
#   fit_diversity_pd        (script 07)
#
# Method: Moran's I on posterior mean residuals using inverse geographic
# distance weights between site centroids. Site coordinates are derived from
# mean plot-level lat/lon in the carabid data (neon.df.csv).
#
# Temporal autocorrelation is not tested separately: turnover and diversity
# models use site-level means averaged across years (temporal structure
# collapsed before modelling), and the co-occurrence models (script 05)
# include a (1|year) random effect.
#
# Phylogenetic autocorrelation is not tested for site-level models because
# the unit of analysis is sites, not species.
#
# If Moran's I is significant for any model, the recommended fix is to add
# latitude as a covariate (or a thin-plate spline of lat/lon) and refit.
# With n ≈ 44 sites, a full spatial covariance model (GP/CAR) is likely
# overparameterized.
#
# Prerequisites:
#   04_analysisTurnoverModels.R → modelOutputs/fit_temporal_turnover.rds
#                               → modelOutputs/fit_spatial_turnover.rds
#   07_linkingTurnoverToDiversity.R → modelOutputs/fit_diversity_richness.rds
#                                   → modelOutputs/fit_diversity_pd.rds
#
# Outputs:
#   data/derivedData/morans_i_results.csv
#   figures/residuals_vs_latitude.png
# =============================================================================

library(tidyverse)
library(ape)
library(brms)
library(terra)
library(sf)

EXCLUDE_SITES <- c("LAJA", "GUAN", "PUUM")


# -----------------------------------------------------------------------------
# 1. Site coordinates from carabid data
# -----------------------------------------------------------------------------
# Average plot-level lat/lon within each site to get site centroids.
# This avoids dependence on the NEON shapefile (not tracked in repo).

plots_sf <- vect("data/All_NEON_TOS_Plots_V11/All_NEON_TOS_Plot_Polygons_V11.shp") %>% 
  st_as_sf()

dbp <- plots_sf %>%
  filter(plotType == "distributed", subtype == "basePlot")

dbp_coords <- st_coordinates(st_centroid(dbp))

dbp <- dbp %>%
  mutate(lon = dbp_coords[, 1],
         lat = dbp_coords[, 2])

df <- read.csv("data/neon.df.csv")

df <- left_join(df, dbp)

# Site coordinates from carabid data (averaged across plots within site)
site_coords <- df %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  group_by(siteID) %>%
  summarise(
    lat = mean(lat,  na.rm = TRUE),
    lon = mean(lon, na.rm = TRUE),
    .groups = "drop"
  )

message(sprintf("Site coordinates available for %d sites", nrow(site_coords)))

# -----------------------------------------------------------------------------
# 2. Reconstruct modelling data frames (mirrors scripts 04 and 07)
# -----------------------------------------------------------------------------

scale_predictors <- function(df) {
  df %>%
    mutate(
      S.chao1      = as.numeric(scale(S.chao1)),
      MPD          = as.numeric(scale(MPD)),
      seasonLength = as.numeric(scale(seasonLength)),
      H_habitat    = as.numeric(scale(H_plant_bc))
    )
}

mdf_temporal <- read.csv("data/derivedData/temporal_turnover.csv") %>%
  scale_predictors() %>%
  filter(!is.na(mean_consecutive_turnover), !is.na(seasonLength),
         !is.na(H_habitat), !is.na(MPD))

mdf_spatial <- read.csv("data/derivedData/spatial_turnover.csv") %>%
  scale_predictors() %>%
  filter(!is.na(mean_spatial_turnover), !is.na(seasonLength),
         !is.na(H_habitat), !is.na(MPD))

# Script 07 data: raw turnover files, only predictor scaling
temporal_raw <- read.csv("data/derivedData/temporal_turnover.csv")
spatial_raw  <- read.csv("data/derivedData/spatial_turnover.csv") %>%
  select(siteID, mean_spatial_turnover)

mdf_diversity <- temporal_raw %>%
  left_join(spatial_raw, by = "siteID") %>%
  filter(!is.na(S.chao1), !is.na(PD),
         !is.na(mean_consecutive_turnover), !is.na(mean_spatial_turnover),
         S.chao1 > 0, PD > 0) %>%
  mutate(
    temporal_sc = as.numeric(scale(mean_consecutive_turnover)),
    spatial_sc  = as.numeric(scale(mean_spatial_turnover))
  )


# -----------------------------------------------------------------------------
# 3. Moran's I helper
# -----------------------------------------------------------------------------
# Builds an inverse-distance weight matrix from site coordinates, then runs
# ape::Moran.I() on the supplied residual vector.
# Returns a one-row tibble with the test statistics.

run_morans <- function(resid_vec, site_ids, site_coords, model_label) {
  
  coords <- site_coords %>%
    filter(siteID %in% site_ids) %>%
    arrange(match(siteID, site_ids))
  
  if (nrow(coords) != length(resid_vec)) {
    warning(sprintf("%s: coordinate/residual length mismatch — skipping", model_label))
    return(NULL)
  }
  
  coord_mat <- as.matrix(coords[, c("lon", "lat")])
  dist_mat  <- as.matrix(dist(coord_mat))
  inv_dist  <- 1 / dist_mat
  diag(inv_dist) <- 0
  
  mi <- Moran.I(resid_vec, inv_dist)
  
  tibble(
    model      = model_label,
    n_sites    = length(resid_vec),
    observed_I = round(mi$observed, 4),
    expected_I = round(mi$expected, 4),
    sd_I       = round(mi$sd,       4),
    p_value    = round(mi$p.value,  4)
  )
}


# -----------------------------------------------------------------------------
# 4. Extract residuals and run Moran's I for each model
# -----------------------------------------------------------------------------

fit_temp  <- readRDS("modelOutputs/fit_temporal_turnover.rds")
fit_space <- readRDS("modelOutputs/fit_spatial_turnover.rds")
fit_rich  <- readRDS("modelOutputs/fit_diversity_richness.rds")
fit_pd    <- readRDS("modelOutputs/fit_diversity_pd.rds")

# Posterior mean residuals (same row order as modelling data)
resid_temp  <- residuals(fit_temp,  summary = TRUE)[, "Estimate"]
resid_space <- residuals(fit_space, summary = TRUE)[, "Estimate"]
resid_rich  <- residuals(fit_rich,  summary = TRUE)[, "Estimate"]
resid_pd    <- residuals(fit_pd,    summary = TRUE)[, "Estimate"]

morans_results <- bind_rows(
  run_morans(resid_temp,  mdf_temporal$siteID,  site_coords, "Temporal Turnover"),
  run_morans(resid_space, mdf_spatial$siteID,   site_coords, "Spatial Turnover"),
  run_morans(resid_rich,  mdf_diversity$siteID, site_coords, "Species Richness"),
  run_morans(resid_pd,    mdf_diversity$siteID, site_coords, "Phylogenetic Diversity")
)

message("\n--- Moran's I results ---")
print(morans_results)

write.csv(morans_results, "data/derivedData/morans_i_results.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 5. Residuals vs. latitude plot
# -----------------------------------------------------------------------------
# Visual check for latitudinal trend in residuals. A systematic trend
# suggests spatial structure not captured by the model.

resid_long <- bind_rows(
  tibble(siteID = mdf_temporal$siteID,  resid = resid_temp,  model = "Temporal Turnover"),
  tibble(siteID = mdf_spatial$siteID,   resid = resid_space, model = "Spatial Turnover"),
  tibble(siteID = mdf_diversity$siteID, resid = resid_rich,  model = "Species Richness"),
  tibble(siteID = mdf_diversity$siteID, resid = resid_pd,    model = "Phylogenetic Diversity")
) %>%
  left_join(site_coords, by = "siteID")

ggplot(resid_long, aes(x = lat, y = resid)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.6, size = 2, color = "grey30") +
  geom_smooth(method = "loess", span = 1, color = "#2E86AB",
              fill = "#2E86AB", alpha = 0.15, linewidth = 0.8) +
  facet_wrap(~ model, scales = "free_y", ncol = 2) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold")) +
  labs(x = "Latitude", y = "Posterior mean residual",
       title = "Residuals vs. latitude — spatial autocorrelation check")

ggsave("figures/residuals_vs_latitude.png", width = 7, height = 5)

message("\nDone. Results written to data/derivedData/ and figures/")

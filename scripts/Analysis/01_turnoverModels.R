# =============================================================================
# Bayesian regression models testing how season length, habitat heterogeneity
# (plant beta diversity), and MPD predict temporal and spatial community
# turnover and the partitioning ratio (spatial / temporal turnover).
#
# Spatial random effect: Moran's I tests (script 08) detected significant
# spatial autocorrelation in residuals for both turnover models (temporal:
# I = 0.220, p = 0.003; spatial: I = 0.145, p = 0.039). A 2D isotropic
# Gaussian Process term gp(lon_sc, lat_sc) is included in all turnover and
# partitioning ratio models to absorb residual spatial covariance. The GP
# models spatial covariance as a continuous function of geographic distance
# (squared exponential kernel) without assuming a directional gradient.
# Coordinates are scaled to mean 0, SD 1 before fitting for numerical
# stability. If GP length-scale posteriors are poorly identified (very wide
# CIs), consider falling back to MEM eigenvectors (adespatial::mem()) as
# an alternative spatial filter.
#
# Likelihood families:
#   Temporal & spatial turnover [0,1]: Beta
#     - Theoretically appropriate for bounded dissimilarity values.
#     - LOO comparison against Gaussian showed negligible differences
#       (temporal: ΔELPD = 0.8 ± 0.3 SE in favour of Beta;
#        spatial:  ΔELPD = 0.6 ± 0.5 SE, within noise). Beta is used
#       on principled grounds and for consistency across both turnover models.
#     - Observed ranges (spatial: 0.29–0.81; temporal: 0.15–0.70) are well
#       within (0, 1) with no boundary values, so no transformation is needed.
#   Partitioning ratio (positive, unbounded): log-normal
#     - LOO strongly preferred log-normal over Gaussian
#       (ΔELPD = 10.7 ± 2.3 SE).
#
# Prerequisites:
#   02_buildSpatialTurnover.R  → data/derivedData/spatial_turnover.csv
#   03_buildTemporalTurnover.R → data/derivedData/temporal_turnover.csv
#
# Outputs:
#   modelOutputs/fit_temporal_turnover.rds
#   modelOutputs/fit_temporal_turnover_rich.rds
#   modelOutputs/fit_spatial_turnover.rds
#   modelOutputs/fit_spatial_turnover_rich.rds
#   modelOutputs/fit_partitioning_ratio.rds
#   modelOutputs/fit_partitioning_ratio_mpd.rds
#   data/derivedData/partitioning_ratio.csv
#   figures/temporalTurnover.png
#   figures/spatialTurnover.png
#   figures/partitioningRatio.png
#   figures/turnoverModels_allPanel.png
# =============================================================================

library(tidyverse)
library(brms)
library(tidybayes)
library(sf)
library(terra)
library(adespatial)
library(spdep)

# -----------------------------------------------------------------------------
# 1. Load and prepare data
# -----------------------------------------------------------------------------

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

mdf_temporal <- read.csv("data/derivedData/temporal_turnover.csv")
mdf_spatial  <- read.csv("data/derivedData/spatial_turnover.csv")

scale_predictors <- function(df) {
  df %>%
    left_join(site_coords, by = "siteID") %>%
    mutate(
      S.chao1      = as.numeric(scale(S.chao1)),
      MPD          = as.numeric(scale(MPD)),
      seasonLength = as.numeric(scale(seasonLength)),
      H_habitat    = as.numeric(scale(H_plant_bc)),
      lat_sc       = as.numeric(scale(lat)),
      lon_sc       = as.numeric(scale(lon))
    )
}

mdf_temporal <- scale_predictors(mdf_temporal)
mdf_spatial  <- scale_predictors(mdf_spatial)

# Partitioning ratio
mdf_tot <- mdf_temporal %>%
  left_join(select(mdf_spatial, siteID, mean_spatial_turnover), by = "siteID") %>%
  mutate(partitioning_ratio = mean_spatial_turnover / mean_consecutive_turnover) %>%
  filter(is.finite(partitioning_ratio))

write.csv(mdf_tot, "data/derivedData/partitioning_ratio.csv", row.names = FALSE)



# -----------------------------------------------------------------------------
# 2. MEM spatial filtering
# -----------------------------------------------------------------------------
# Forward-selects Moran's Eigenvector Maps that significantly reduce spatial
# autocorrelation in OLS residuals of the baseline (non-spatial) model.
# Returns: list with $mem_names (character) and $mem_df (data frame with
# siteID + selected MEM columns for joining to the modelling data).

select_mems <- function(lm_residuals, site_ids, coord_df,
                        alpha = 0.05, nperm = 999) {
  
  coords <- coord_df %>%
    filter(siteID %in% site_ids) %>%
    arrange(match(siteID, site_ids))
  
  # Inverse-distance spatial weights (all sites connected)
  dist_mat  <- as.matrix(dist(coords[, c("lon_sc", "lat_sc")]))
  inv_dist  <- 1 / dist_mat
  diag(inv_dist) <- 0
  listw_obj <- mat2listw(inv_dist, style = "W")
  
  # All positive MEMs (positive eigenvalues → positive Moran's I by construction)
  mem_all    <- mem(listw_obj)
  mem_matrix <- as.data.frame(mem_all)
  
  # Forward selection against OLS residuals
  fs <- tryCatch(
    forward.sel(Y     = as.data.frame(lm_residuals),
                X     = mem_matrix,
                alpha = alpha,
                nperm = nperm),
    error = function(e) {
      message("  forward.sel() error: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(fs) || nrow(fs) == 0) {
    message("  → 0 MEMs selected")
    return(list(
      mem_names = character(0),
      mem_df    = tibble(siteID = coords$siteID)
    ))
  }
  
  sel_names <- fs$variables
  message(sprintf("  → %d MEM(s) selected: %s",
                  length(sel_names), paste(sel_names, collapse = ", ")))
  
  mem_out        <- mem_matrix[, sel_names, drop = FALSE]
  mem_out$siteID <- coords$siteID
  
  list(mem_names = sel_names,
       mem_df    = select(mem_out, siteID, everything()))
}

# Coordinate data frame for MEM computation (scaled, one row per site)
coords_df <- mdf_temporal %>% select(siteID, lon_sc, lat_sc) %>% distinct()

# --- Temporal turnover MEMs ---
message("Selecting MEMs — temporal turnover:")
lm_temp_base <- lm(mean_consecutive_turnover ~ seasonLength + H_habitat + MPD,
                   data = mdf_temporal)
mems_temp  <- select_mems(residuals(lm_temp_base), mdf_temporal$siteID, coords_df)
mdf_temporal <- left_join(mdf_temporal, mems_temp$mem_df, by = "siteID")

# --- Spatial turnover MEMs ---
message("Selecting MEMs — spatial turnover:")
lm_space_base <- lm(mean_spatial_turnover ~ seasonLength + H_habitat + MPD,
                    data = mdf_spatial)
mems_space <- select_mems(residuals(lm_space_base), mdf_spatial$siteID, coords_df)
mdf_spatial <- left_join(mdf_spatial, mems_space$mem_df, by = "siteID")

# --- Partitioning ratio MEMs ---
# Update mdf_tot with any newly added MEM columns from mdf_temporal join
mdf_tot <- mdf_temporal %>%
  left_join(select(mdf_spatial, siteID, mean_spatial_turnover), by = "siteID") %>%
  mutate(partitioning_ratio = mean_spatial_turnover / mean_consecutive_turnover) %>%
  filter(is.finite(partitioning_ratio))

message("Selecting MEMs — partitioning ratio:")
lm_pr_base <- lm(log(partitioning_ratio) ~ seasonLength + H_habitat,
                 data = mdf_tot)
mems_pr <- select_mems(residuals(lm_pr_base), mdf_tot$siteID, coords_df)
# mdf_tot inherits MEM columns from mdf_temporal; only join MEMs not already present
new_pr_mems <- setdiff(mems_pr$mem_names, names(mdf_tot))
if (length(new_pr_mems) > 0) {
  mdf_tot <- left_join(mdf_tot,
                       select(mems_pr$mem_df, siteID, all_of(new_pr_mems)),
                       by = "siteID")
}


# -----------------------------------------------------------------------------
# 3. Shared priors and MCMC settings
# -----------------------------------------------------------------------------

priors <- c(
  prior(normal(0, 1),   class = Intercept),
  prior(normal(0, 0.5), class = b)
)

mcmc_args <- list(
  chains  = 4,
  iter    = 2000,
  warmup  = 1000,
  control = list(adapt_delta = 0.99),
  cores   = 4,
  seed    = 1234,
  backend = "cmdstanr"
)


# -----------------------------------------------------------------------------
# 4. Fit models
# -----------------------------------------------------------------------------
# Formulas built dynamically so the number of MEM terms adjusts automatically
# to however many were selected per response.

build_formula <- function(response, fixed_terms, mem_names) {
  rhs <- paste(c(fixed_terms, mem_names), collapse = " + ")
  as.formula(paste(response, "~", rhs))
}

# --- Temporal turnover (Beta) ---
fit_temp <- do.call(brm, c(
  list(formula = bf(build_formula("mean_consecutive_turnover",
                                  c("seasonLength", "H_habitat", "MPD"),
                                  mems_temp$mem_names)),
       family   = Beta(),
       data     = mdf_temporal,
       prior    = priors),
  mcmc_args
))

fit_temp_rich <- do.call(brm, c(
  list(formula = bf(build_formula("mean_consecutive_turnover",
                                  c("seasonLength", "H_habitat", "MPD", "S.chao1"),
                                  mems_temp$mem_names)),
       family   = Beta(),
       data     = mdf_temporal,
       prior    = priors),
  mcmc_args
))

saveRDS(fit_temp,      "modelOutputs/fit_temporal_turnover.rds")
saveRDS(fit_temp_rich, "modelOutputs/fit_temporal_turnover_rich.rds")

summary(fit_temp)
summary(fit_temp_rich)

# --- Spatial turnover (Beta) ---
fit_space <- do.call(brm, c(
  list(formula = bf(build_formula("mean_spatial_turnover",
                                  c("seasonLength", "H_habitat", "MPD"),
                                  mems_space$mem_names)),
       family   = Beta(),
       data     = mdf_spatial,
       prior    = priors),
  mcmc_args
))

fit_space_rich <- do.call(brm, c(
  list(formula = bf(build_formula("mean_spatial_turnover",
                                  c("seasonLength", "H_habitat", "MPD", "S.chao1"),
                                  mems_space$mem_names)),
       family   = Beta(),
       data     = mdf_spatial,
       prior    = priors),
  mcmc_args
))

saveRDS(fit_space,      "modelOutputs/fit_spatial_turnover.rds")
saveRDS(fit_space_rich, "modelOutputs/fit_spatial_turnover_rich.rds")

summary(fit_space)
summary(fit_space_rich)

# --- Partitioning ratio (log-normal) ---
fit_pr <- do.call(brm, c(
  list(formula = bf(build_formula("partitioning_ratio",
                                  c("seasonLength", "H_habitat"),
                                  mems_pr$mem_names)),
       family   = lognormal(),
       data     = mdf_tot,
       prior    = priors),
  mcmc_args
))

fit_pr_mpd <- do.call(brm, c(
  list(formula = bf(build_formula("partitioning_ratio",
                                  c("seasonLength", "H_habitat", "MPD"),
                                  mems_pr$mem_names)),
       family   = lognormal(),
       data     = mdf_tot,
       prior    = priors),
  mcmc_args
))

saveRDS(fit_pr,     "modelOutputs/fit_partitioning_ratio.rds")
saveRDS(fit_pr_mpd, "modelOutputs/fit_partitioning_ratio_mpd.rds")

summary(fit_pr)
summary(fit_pr_mpd)


# -----------------------------------------------------------------------------
# 5. Posterior visualization helpers
# -----------------------------------------------------------------------------
# MEM columns are excluded from figures — they are spatial nuisance terms.

tidy_brms <- function(fit, label) {
  as_draws_df(fit) %>%
    select(starts_with("b_"), -contains("Intercept"), -matches("b_MEM\\d+")) %>%
    pivot_longer(everything(), names_to = "term", values_to = "value") %>%
    mutate(
      term  = str_remove(term, "^b_"),
      term  = recode(term,
                     "seasonLength" = "Season Length",
                     "H_habitat"    = "Habitat Heterogeneity",
                     "MPD"          = "MPD",
                     "S.chao1"      = "Species Richness"),
      model = label
    )
}

add_significance <- function(draws_df) {
  draws_df %>%
    group_by(model, term) %>%
    mutate(significant = !(quantile(value, 0.025) < 0 &
                             quantile(value, 0.975) > 0)) %>%
    ungroup()
}

coef_plot <- function(draws_df) {
  draws_df %>%
    add_significance() %>%
    ggplot(aes(x = value, y = term, fill = significant)) +
    stat_halfeye(
      .width         = c(0.5, 0.95),
      point_interval = median_qi,
      alpha          = 0.8
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "grey70"),
                      guide  = "none") +
    facet_wrap(~ model, scales = "free_x") +
    theme_classic() +
    theme(strip.background = element_blank(),
          strip.text       = element_text(face = "bold")) +
    labs(x = "Posterior Estimate", y = NULL)
}


# -----------------------------------------------------------------------------
# 6. Figures
# -----------------------------------------------------------------------------

tt <- tidy_brms(fit_temp,   "Temporal Turnover")  %>% coef_plot()
ggsave("figures/temporalTurnover.png", width = 4, height = 3)

st <- tidy_brms(fit_space,  "Spatial Turnover")   %>% coef_plot()
ggsave("figures/spatialTurnover.png",  width = 4, height = 3)

tidy_brms(fit_pr_mpd, "Partitioning Ratio") %>% coef_plot()
ggsave("figures/partitioningRatio.png", width = 4, height = 3)

bind_rows(
  tidy_brms(fit_temp,       "Temporal Turnover"),
  tidy_brms(fit_temp_rich,  "Temporal Turnover\n(+ Richness)"),
  tidy_brms(fit_space,      "Spatial Turnover"),
  tidy_brms(fit_space_rich, "Spatial Turnover\n(+ Richness)"),
  tidy_brms(fit_pr,         "Partitioning Ratio"),
  tidy_brms(fit_pr_mpd,     "Partitioning Ratio\n(+ MPD)")
) %>%
  coef_plot()
ggsave("figures/turnoverModels_allPanel.png", width = 12, height = 4)

## two panel of temporal turnover and spatial turnover
cowplot::plot_grid(tt, st, labels = c("A", "B"))
ggsave("figures/turnoverModels_spatialTemporal.png", width = 8, height = 4)

## now plus richness
ttr <- tidy_brms(fit_temp_rich,  "Temporal Turnover\n(+ Richness)") %>% coef_plot()
str <- tidy_brms(fit_space_rich, "Spatial Turnover\n(+ Richness)")%>% coef_plot()
cowplot::plot_grid(ttr, str, labels = c("A", "B"))
ggsave("figures/turnoverModels_spatialTemporal_plusRichness.png", width = 8, height = 4)

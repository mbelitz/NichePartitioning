# =============================================================================
# Tests whether temporal and spatial community turnover predict site-level
# diversity. The core hypothesis: if niche partitioning facilitates
# coexistence, sites with greater turnover (more partitioning) should
# support higher species richness and phylogenetic diversity.
#
# Spatial random effect: Moran's I tests (script 08) detected strong spatial
# autocorrelation in residuals for both diversity models (Richness: I = 0.362,
# p ≈ 0; PD: I = 0.336, p ≈ 0). A 2D isotropic Gaussian Process term
# gp(lon_sc, lat_sc) is included in both models to absorb residual spatial
# covariance without assuming a directional (latitudinal) gradient. Coordinates
# are scaled to mean 0, SD 1 for numerical stability.
#
# Response variables:
#   S.chao1  — Chao1 species richness estimate
#   PD       — Faith's phylogenetic diversity
#
# NOTE: MPD is excluded here because it was used as a predictor of turnover
# in scripts 04–05. Using it as a response to turnover would create
# circular inference.
#
# Predictors (scaled):
#   mean_consecutive_turnover  — temporal niche partitioning
#   mean_spatial_turnover      — spatial niche partitioning
#
# Likelihood: lognormal for both responses (positive, right-skewed).
#
# Prerequisites:
#   03_buildTemporalTurnover.R → data/derivedData/temporal_turnover.csv
#   02_buildSpatialTurnover.R  → data/derivedData/spatial_turnover.csv
#
# Outputs:
#   modelOutputs/fit_diversity_richness.rds
#   modelOutputs/fit_diversity_pd.rds
#   figures/linkingTurnoverToDiversity.png
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
# Load raw (unscaled) turnover files so that S.chao1 and PD are on their
# natural scales as response variables. Only the predictors are scaled.


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

temporal <- read.csv("data/derivedData/temporal_turnover.csv")
spatial  <- read.csv("data/derivedData/spatial_turnover.csv") %>%
  select(siteID, mean_spatial_turnover)

mdf <- temporal %>%
  left_join(spatial,      by = "siteID") %>%
  left_join(site_coords,  by = "siteID") %>%
  filter(
    !is.na(S.chao1),
    !is.na(PD),
    !is.na(mean_consecutive_turnover),
    !is.na(mean_spatial_turnover),
    !is.na(lat), !is.na(lon),
    S.chao1 > 0,
    PD      > 0
  ) %>%
  mutate(
    temporal_sc = as.numeric(scale(mean_consecutive_turnover)),
    spatial_sc  = as.numeric(scale(mean_spatial_turnover)),
    lat_sc      = as.numeric(scale(lat)),
    lon_sc      = as.numeric(scale(lon))
  )

message(sprintf("Linking diversity model dataset: %d sites", nrow(mdf)))


# -----------------------------------------------------------------------------
# 2. MEM spatial filtering
# -----------------------------------------------------------------------------
# Computes the full MEM matrix once from inverse-distance weights, then
# forward-selects MEMs separately for each response against OLS residuals.
# The union of selected MEMs is joined to mdf in a single operation to avoid
# duplicate column names (which arise when both responses select the same MEM
# and sequential left_joins create .x/.y suffixes that break brms formulas).

coords_df <- mdf %>% select(siteID, lon_sc, lat_sc) %>% distinct()

# Build full MEM matrix (positive MEMs only; positive eigenvalues → positive
# spatial autocorrelation by construction)
dist_mat_mem <- as.matrix(dist(coords_df[, c("lon_sc", "lat_sc")]))
inv_dist_mem <- 1 / dist_mat_mem
diag(inv_dist_mem) <- 0
listw_mem    <- mat2listw(inv_dist_mem, style = "W")
mem_all      <- mem(listw_mem)
mem_matrix   <- as.data.frame(mem_all)  # rows aligned to coords_df row order

# Helper: returns names of forward-selected MEMs (no data frame, no join here)
select_mem_names <- function(lm_residuals, mem_matrix,
                             alpha = 0.05, nperm = 999) {
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
    return(character(0))
  }
  sel <- fs$variables
  message(sprintf("  → %d MEM(s) selected: %s",
                  length(sel), paste(sel, collapse = ", ")))
  sel
}

# --- Species richness MEMs ---
message("Selecting MEMs — species richness:")
lm_rich_base <- lm(S.chao1 ~ temporal_sc + spatial_sc, data = mdf)
mems_rich    <- select_mem_names(residuals(lm_rich_base), mem_matrix)

# --- PD MEMs ---
message("Selecting MEMs — phylogenetic diversity:")
lm_pd_base <- lm(PD ~ temporal_sc + spatial_sc, data = mdf)
mems_pd    <- select_mem_names(residuals(lm_pd_base), mem_matrix)

# Single join of the union of all selected MEMs — no duplicate columns
all_mem_names <- union(mems_rich, mems_pd)
if (length(all_mem_names) > 0) {
  mem_join_df <- bind_cols(
    tibble(siteID = coords_df$siteID),
    mem_matrix[, all_mem_names, drop = FALSE]
  )
  mdf <- left_join(mdf, mem_join_df, by = "siteID")
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

build_formula <- function(response, fixed_terms, mem_names) {
  rhs <- paste(c(fixed_terms, mem_names), collapse = " + ")
  as.formula(paste(response, "~", rhs))
}

fit_rich <- do.call(brm, c(
  list(
    formula = bf(build_formula("S.chao1",
                               c("temporal_sc", "spatial_sc"),
                               mems_rich)),
    family  = lognormal(),
    data    = mdf,
    prior   = priors
  ),
  mcmc_args
))

summary(fit_rich)
saveRDS(fit_rich, "modelOutputs/fit_diversity_richness.rds")


fit_pd <- do.call(brm, c(
  list(
    formula = bf(build_formula("PD",
                               c("temporal_sc", "spatial_sc"),
                               mems_pd)),
    family  = lognormal(),
    data    = mdf,
    prior   = priors
  ),
  mcmc_args
))

summary(fit_pd)
saveRDS(fit_pd, "modelOutputs/fit_diversity_pd.rds")


# -----------------------------------------------------------------------------
# 5. Coefficient plot
# -----------------------------------------------------------------------------
# MEM terms are excluded — only the ecological fixed effects are shown.

tidy_brms <- function(fit, label) {
  as_draws_df(fit) %>%
    select(b_temporal_sc, b_spatial_sc) %>%
    pivot_longer(everything(), names_to = "term", values_to = "value") %>%
    mutate(
      term  = recode(term,
                     "b_temporal_sc" = "Temporal Turnover",
                     "b_spatial_sc"  = "Spatial Turnover"),
      model = label
    )
}

draws <- bind_rows(
  tidy_brms(fit_rich, "Species Richness (Chao1)"),
  tidy_brms(fit_pd,   "Phylogenetic Diversity")
) %>%
  group_by(model, term) %>%
  mutate(significant = !(quantile(value, 0.025) < 0 &
                           quantile(value, 0.975) > 0)) %>%
  ungroup()

ggplot(draws, aes(x = value, y = term, fill = significant)) +
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
  labs(x = "Posterior Estimate (log scale)", y = NULL)

ggsave("figures/linkingTurnoverToDiversity.png", width = 6, height = 3)

message("\nDone. Results written to modelOutputs/ and figures/")
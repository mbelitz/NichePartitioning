# =============================================================================
# 05_analysisCoOccurrenceModels.R
#
# Bayesian multilevel models testing how season length, pairwise phylogenetic
# distance, and habitat heterogeneity (plant beta diversity) predict pairwise
# species co-occurrence effects across NEON sites.
#
# Three models are fit:
#   fit_tp    — temporal co-occurrence effects
#   fit_sp    — spatial co-occurrence effects
#   fit_joint — joint model with process (temporal/spatial) × predictor
#               interactions, enabling direct comparison across axes
#
# Random intercepts account for non-independence of pairwise observations:
#   (1|siteID), (1|sp1_clean), (1|pair_key), (1|year)
#
# Prerequisites:
#   01_buildHabitatHeterogeneity.R → data/derivedData/plant_heterogeneity.csv
#   Input co-occurrence files are pre-computed pairwise correlation effects:
#     data/time.partioning.nov.csv    (temporal co-occurrence)
#     data/spatial.partioning.nov.csv (spatial co-occurrence)
#
# Outputs:
#   modelOutputs/fit_temporal_cooccurrence.rds
#   modelOutputs/fit_spatial_cooccurrence.rds
#   modelOutputs/fit_joint_cooccurrence.rds
#   figures/temporal_coOccurrenceResults.png
#   figures/spatial_coOccurrenceResults.png
#   figures/joint_coOccurrenceResults.png
# =============================================================================

library(tidyverse)
library(brms)
library(tidybayes)


# -----------------------------------------------------------------------------
# 1. Load co-occurrence data and join habitat heterogeneity
# -----------------------------------------------------------------------------

plant_h <- read.csv("data/derivedData/plant_heterogeneity.csv")  # H_plant_bc

tp_raw <- read.csv("data/time.partioning.nov.csv")
sp_raw <- read.csv("data/spatial.partioning.nov.csv")

tp_raw <- left_join(tp_raw, plant_h, by = "siteID")
sp_raw <- left_join(sp_raw, plant_h, by = "siteID")


# -----------------------------------------------------------------------------
# 2. Prepare modeling data frames
# -----------------------------------------------------------------------------
# Scales seasonLength, pairDist, and H_plant_bc; removes rows with missing
# values in any modelling variable.

prepare_mdf <- function(df) {
  df %>%
    filter(
      !is.na(effects),
      !is.na(seasonLength),
      !is.na(pairDist),
      !is.na(siteID),
      !is.na(sp1_clean),
      !is.na(H_plant_bc)
    ) %>%
    mutate(
      seasonLength_sc = as.numeric(scale(seasonLength)),
      pairDist_sc     = as.numeric(scale(pairDist)),
      H_habitat_sc    = as.numeric(scale(H_plant_bc))
    )
}

tp_mdf <- prepare_mdf(tp_raw)
sp_mdf <- prepare_mdf(sp_raw)

message(sprintf("Temporal dataset: %d pairwise observations across %d sites",
                nrow(tp_mdf), n_distinct(tp_mdf$siteID)))
message(sprintf("Spatial dataset:  %d pairwise observations across %d sites",
                nrow(sp_mdf), n_distinct(sp_mdf$siteID)))


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
  threads = threading(4),
  backend = "cmdstanr"
)


# -----------------------------------------------------------------------------
# 4. Temporal co-occurrence model
# -----------------------------------------------------------------------------

fit_tp <- do.call(brm, c(
  list(
    formula = bf(effects ~ seasonLength_sc + pairDist_sc + H_habitat_sc +
                   (1 | siteID) + (1 | sp1_clean) + (1 | pair_key) + (1 | year)),
    family  = gaussian(),
    data    = tp_mdf,
    prior   = priors
  ),
  mcmc_args
))

summary(fit_tp)
saveRDS(fit_tp, "modelOutputs/fit_temporal_cooccurrence.rds")


# -----------------------------------------------------------------------------
# 5. Spatial co-occurrence model
# -----------------------------------------------------------------------------

fit_sp <- do.call(brm, c(
  list(
    formula = bf(effects ~ seasonLength_sc + pairDist_sc + H_habitat_sc +
                   (1 | siteID) + (1 | sp1_clean) + (1 | pair_key) + (1 | year)),
    family  = gaussian(),
    data    = sp_mdf,
    prior   = priors
  ),
  mcmc_args
))

summary(fit_sp)
saveRDS(fit_sp, "modelOutputs/fit_spatial_cooccurrence.rds")


# -----------------------------------------------------------------------------
# 6. Joint model (process × predictor interactions)
# -----------------------------------------------------------------------------
# Stacks temporal and spatial data with a process indicator. Interaction terms
# (process × predictor) capture whether effects differ between axes. The
# reference level is spatial; temporal effects are the reference + interaction.

joint_long <- bind_rows(
  tp_mdf %>%
    mutate(process = "temporal") %>%
    select(effects, seasonLength_sc, pairDist_sc, H_habitat_sc,
           siteID, sp1_clean, year, pair_key, process),
  sp_mdf %>%
    mutate(process = "spatial") %>%
    select(effects, seasonLength_sc, pairDist_sc, H_habitat_sc,
           siteID, sp1_clean, year, pair_key, process)
) %>%
  mutate(process = factor(process, levels = c("spatial", "temporal")))

fit_joint <- do.call(brm, c(
  list(
    formula = bf(effects ~ process * (seasonLength_sc + pairDist_sc + H_habitat_sc) +
                   (1 | siteID) + (1 | sp1_clean) + (1 | year) + (1 | pair_key)),
    family  = gaussian(),
    data    = joint_long,
    prior   = priors
  ),
  mcmc_args
))

summary(fit_joint)
saveRDS(fit_joint, "modelOutputs/fit_joint_cooccurrence.rds")


# -----------------------------------------------------------------------------
# 7. Individual model figures (temporal and spatial)
# -----------------------------------------------------------------------------
# Reusable helper: extract fixed effect posteriors and plot as halfeye.

plot_cooccurrence <- function(fit, title) {
  as_draws_df(fit) %>%
    select(b_seasonLength_sc, b_pairDist_sc, b_H_habitat_sc) %>%
    pivot_longer(everything(), names_to = "term", values_to = "value") %>%
    mutate(term = recode(term,
                         "b_seasonLength_sc" = "Season Length",
                         "b_pairDist_sc"     = "Phylogenetic Distance",
                         "b_H_habitat_sc"    = "Habitat Heterogeneity")) %>%
    group_by(term) %>%
    mutate(significant = !(quantile(value, 0.025) < 0 &
                             quantile(value, 0.975) > 0)) %>%
    ungroup() %>%
    ggplot(aes(x = value, y = term, fill = significant)) +
    stat_halfeye(
      .width         = c(0.5, 0.95),
      point_interval = median_qi,
      alpha          = 0.8
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "grey70"),
                      guide  = "none") +
    theme_classic() +
    labs(x = "Posterior Estimate", y = NULL, title = title)
}

plot_cooccurrence(fit_tp, "Predictors of Temporal Co-occurrence")
ggsave("figures/temporal_coOccurrenceResults.png", width = 5, height = 4)

plot_cooccurrence(fit_sp, "Predictors of Spatial Co-occurrence")
ggsave("figures/spatial_coOccurrenceResults.png", width = 5, height = 4)


# -----------------------------------------------------------------------------
# 8. Joint model figure (derived process-specific effects)
# -----------------------------------------------------------------------------
# The joint model uses spatial as reference. Temporal effects are derived as:
#   temporal effect = b_predictor + b_processtemporal:predictor

as_draws_df(fit_joint) %>%
  transmute(
    spatial_seasonLength  = b_seasonLength_sc,
    temporal_seasonLength = b_seasonLength_sc  + `b_processtemporal:seasonLength_sc`,
    spatial_pairDist      = b_pairDist_sc,
    temporal_pairDist     = b_pairDist_sc      + `b_processtemporal:pairDist_sc`,
    spatial_habitat       = b_H_habitat_sc,
    temporal_habitat      = b_H_habitat_sc     + `b_processtemporal:H_habitat_sc`
  ) %>%
  pivot_longer(
    everything(),
    names_to  = c("process", "variable"),
    names_sep = "_"
  ) %>%
  mutate(variable = recode(variable,
                           "seasonLength" = "Season Length",
                           "pairDist"     = "Phylogenetic Distance",
                           "habitat"      = "Habitat Heterogeneity")) %>%
  ggplot(aes(x = value, y = variable, fill = process)) +
  stat_halfeye(alpha = 0.7, position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_fill_manual(values = c("spatial" = "#2E86AB", "temporal" = "#E07B39"),
                    labels  = c("Spatial", "Temporal")) +
  theme_classic() +
  labs(x = "Effect on co-occurrence", y = NULL, fill = "Process")

ggsave("figures/joint_coOccurrenceResults.png", width = 6, height = 4)
# =============================================================================
# 03_buildTemporalTurnover.R
#
# Calculates mean consecutive-week Bray-Curtis dissimilarity within each
# NEON site-year as a measure of temporal community turnover (i.e., temporal
# niche partitioning). Site-year values are then averaged across years to
# produce a single site-level estimate.
#
# Also retains mean_all_pairs_bc (all pairwise weeks, not just consecutive)
# and season_length_weeks (observed sampling weeks per site-year) for
# use in sensitivity analyses.
#
# Prerequisites (must be run first):
#   00_buildDiversityMetrics.R  → data/derivedData/rich.csv, PD_MPD.csv
#   01_buildHabitatHeterogeneity.R → data/derivedData/plant_heterogeneity.csv
#
# Outputs:
#   data/derivedData/temporal_turnover.csv
#
# FLAGS:
#
#   [1] season_length_weeks (computed here) is the observed number of weekly
#       sampling events per site-year, averaged across years. This differs
#       from seasonLength in spatial.partioning.nov.csv, which is a modeled
#       climate-based estimate. Both are retained in the output.
# =============================================================================

library(vegan)
library(tidyverse)
library(lubridate)

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
# 2. Build weekly counts per site-year
# -----------------------------------------------------------------------------

weekly_counts <- df %>%
  filter(taxonRank == "species") %>%
  mutate(week = floor_date(as.Date(collectDate), unit = "week")) %>%
  group_by(siteID, year, week, taxonID) %>%
  summarise(totalCount = sum(individualCount, na.rm = TRUE), .groups = "drop")


# -----------------------------------------------------------------------------
# 3. Temporal beta diversity function (one site-year)
# -----------------------------------------------------------------------------

calc_temporal_beta <- function(site_year_df) {
  
  if (is.null(site_year_df) || nrow(site_year_df) == 0) return(NULL)
  
  # Pivot to week × species matrix
  week_matrix <- site_year_df %>%
    select(week, taxonID, totalCount) %>%
    pivot_wider(names_from  = taxonID,
                values_from = totalCount,
                values_fill = 0) %>%
    arrange(week) %>%
    column_to_rownames("week")
  
  # Require at least 3 sampling weeks and at least 1 species
  if (nrow(week_matrix) < 3 || ncol(week_matrix) == 0) return(NULL)
  
  # Drop species absent this site-year
  week_matrix <- week_matrix[, colSums(week_matrix) > 0, drop = FALSE]
  if (ncol(week_matrix) == 0) return(NULL)
  
  # Hellinger transform then Bray-Curtis dissimilarity
  hel_matrix <- decostand(week_matrix, method = "hellinger")
  bc_dist    <- vegdist(hel_matrix, method = "bray")
  bc_matrix  <- as.matrix(bc_dist)
  
  n <- nrow(bc_matrix)
  
  # Consecutive week pairs (week t → week t+1)
  consecutive_bc <- sapply(1:(n - 1), function(i) bc_matrix[i, i + 1])
  
  # All pairwise weeks
  all_pairs_bc <- bc_matrix[lower.tri(bc_matrix)]
  
  tibble(
    mean_consecutive_turnover = mean(consecutive_bc, na.rm = TRUE),
    sd_consecutive_turnover   = sd(consecutive_bc,   na.rm = TRUE),
    mean_all_pairs_bc         = mean(all_pairs_bc,   na.rm = TRUE),
    season_length_weeks       = n
  )
}


# -----------------------------------------------------------------------------
# 4. Apply across all site-years
# -----------------------------------------------------------------------------

temporal_beta <- weekly_counts %>%
  group_by(siteID, year) %>%
  group_split() %>%
  map_dfr(function(x) {
    result <- calc_temporal_beta(x)
    if (!is.null(result)) {
      mutate(result, siteID = unique(x$siteID), year = unique(x$year))
    }
  })

message(sprintf("Temporal turnover calculated for %d site-years across %d sites.",
                nrow(temporal_beta),
                n_distinct(temporal_beta$siteID)))


# -----------------------------------------------------------------------------
# 5. Average across years per site
# -----------------------------------------------------------------------------

temporal_beta_site <- temporal_beta %>%
  group_by(siteID) %>%
  summarise(
    across(
      c(mean_consecutive_turnover, sd_consecutive_turnover,
        mean_all_pairs_bc, season_length_weeks),
      ~ mean(.x, na.rm = TRUE)
    ),
    n_years = n(),
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# 6. Join diversity metrics, habitat heterogeneity, and season length
# -----------------------------------------------------------------------------

rich    <- read.csv("data/derivedData/rich.csv")
pd      <- read.csv("data/derivedData/PD_MPD.csv")
plant_h <- read.csv("data/derivedData/plant_heterogeneity.csv")

sl <- read.csv("data/spatial.partioning.nov.csv") %>%  
  distinct(siteID, seasonLength)

mdf <- rich %>%
  left_join(pd,                  by = "siteID") %>%
  left_join(temporal_beta_site,  by = "siteID") %>%
  left_join(plant_h,             by = "siteID") %>%
  left_join(sl,                  by = "siteID")


# -----------------------------------------------------------------------------
# 7. Save output
# -----------------------------------------------------------------------------

write.csv(mdf, "data/derivedData/temporal_turnover.csv", row.names = FALSE)

message("Done. Output written to data/derivedData/temporal_turnover.csv")
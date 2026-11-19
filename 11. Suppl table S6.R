
# ============================================================
# SUPPLEMENTARY TABLE 6
# Single-moderator and pairwise comparison analyses
# Overall, RCTs, and QEDs
# ============================================================

# 1. LIBRARIES
# ============================================================

library(metafor)
library(emmeans)
library(dplyr)
library(purrr)
library(tidyr)
library(writexl)


# 2. LOAD AND CLEAN DATA
# ============================================================

mydata_cog <- read.csv(
  "D:/PhD Nurul/Script submission/mydata_cog.csv",
  stringsAsFactors = FALSE
)

mydata_cog <- mydata_cog %>%
  mutate(
    final_d_effect = as.numeric(final_d_effect),
    variance_final_d_effect = as.numeric(variance_final_d_effect)
  )

mydata_cog[] <- lapply(mydata_cog, function(x) {
  if (is.character(x)) {
    x[trimws(x) == ""] <- NA
  }
  x
})


# 3. DEFINE MODERATORS AND SAMPLES
# ============================================================

moderators <- c(
  "outcome_content",
  "hours_cat",
  "weeks_duration_cat_1"
)

samples <- c("All", "RCTs", "QEDs")


# 4. FIT SINGLE-MODERATOR MODELS
# ============================================================

fit_moderator_model <- function(data, moderator, sample) {
  
  dat <- data %>%
    filter(
      !is.na(.data[[moderator]]),
      !is.na(final_d_effect),
      !is.na(variance_final_d_effect),
      !is.na(study_id),
      !is.na(effect_id)
    )
  
  if (sample == "RCTs") {
    dat <- dat %>% filter(study_type == "rct")
  }
  
  if (sample == "QEDs") {
    dat <- dat %>% filter(study_type == "quasi")
  }
  
  if (n_distinct(dat[[moderator]]) < 2) {
    return(NULL)
  }
  
  V <- vcalc(
    vi = variance_final_d_effect,
    cluster = study_id,
    obs = effect_id,
    data = dat,
    rho = 0.5
  )
  
  model <- rma.mv(
    yi = final_d_effect,
    V = V,
    mods = as.formula(paste("~", moderator)),
    random = ~ 1 | study_id/effect_id,
    method = "REML",
    data = dat
  )
  
  list(
    model = model,
    data = dat,
    moderator = moderator,
    sample = sample
  )
}

model_grid <- expand_grid(
  moderator = moderators,
  sample = samples
)

all_models <- pmap(
  model_grid,
  function(moderator, sample) {
    fit_moderator_model(
      data = mydata_cog,
      moderator = moderator,
      sample = sample
    )
  }
)

names(all_models) <- paste(
  model_grid$moderator,
  model_grid$sample,
  sep = "_"
)

all_models <- compact(all_models)


# 5. OMNIBUS MODERATOR TESTS
# ============================================================

omnibus_results <- map_dfr(all_models, function(x) {
  
  m <- x$model
  
  tibble(
    moderator = x$moderator,
    sample = x$sample,
    k = m$k,
    n_studies = n_distinct(x$data$study_id),
    QM = m$QM,
    df = m$m,
    p_value = m$QMp
  )
})


# 6. POOLED EFFECT SIZES BY CATEGORY
# ============================================================

get_pooled_results <- function(x) {
  
  emm <- emmeans(
    emmprep(x$model),
    as.formula(paste("~", x$moderator))
  )
  
  cld_results <- multcomp::cld(
    emm,
    adjust = "sidak",
    Letters = letters
  ) %>%
    as.data.frame()
  
  counts <- x$data %>%
    group_by(category = .data[[x$moderator]]) %>%
    summarise(
      k = n(),
      n_studies = n_distinct(study_id),
      .groups = "drop"
    ) %>%
    mutate(category = as.character(category))
  
  cld_results %>%
    mutate(
      category = as.character(.data[[x$moderator]]),
      moderator = x$moderator,
      sample = x$sample,
      .group = gsub(" ", "", .group)
    ) %>%
    left_join(counts, by = "category") %>%
    transmute(
      moderator,
      sample,
      category,
      k,
      n_studies,
      cohens_d = emmean,
      SE,
      CI_lower = asymp.LCL,
      CI_upper = asymp.UCL,
      letter = .group
    )
}

pooled_results <- map_dfr(
  all_models,
  get_pooled_results
)


# 7. PAIRWISE COMPARISONS (SIDAK ADJUSTMENT)
# ============================================================

get_pairwise_results <- function(x) {
  
  emm <- emmeans(
    emmprep(x$model),
    as.formula(paste("~", x$moderator))
  )
  
  pairs(emm, adjust = "sidak") %>%
    as.data.frame() %>%
    transmute(
      moderator = x$moderator,
      sample = x$sample,
      contrast,
      estimate,
      SE,
      z_ratio = z.ratio,
      p_value = p.value
    )
}

pairwise_results <- map_dfr(
  all_models,
  get_pairwise_results
)


# 8. DISPLAY RESULTS
# ============================================================

print(omnibus_results, n = Inf)
print(pooled_results, n = Inf)
print(pairwise_results, n = Inf)



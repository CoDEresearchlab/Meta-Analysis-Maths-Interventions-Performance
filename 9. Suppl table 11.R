library(metafor)
library(emmeans)
library(multcomp)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(writexl)


# Loading the data
#######
mydata_cog <- read.csv(
  "D:/PhD Nurul/Script submission/mydata_cog.csv",
  stringsAsFactors = FALSE
)

# Ensure numeric columns
mydata_cog <- mydata_cog %>%
  mutate(
    final_d_effect = as.numeric(final_d_effect),
    variance_final_d_effect = as.numeric(variance_final_d_effect)
  )

# Replace blank cells with NA across all columns
mydata_cog[] <- lapply(mydata_cog, function(x) {
  if (is.character(x)) {
    x[trimws(x) == ""] <- NA
  }
  x
})
#######

## 1. Filter MD data
dat_MD <- mydata_cog %>%
  filter(risk_status2 == "maths difficulties (MD)")

## 2. Moderators and labels
moderators <- c(
  "outcome_content",
  "weeks_duration_cat_1",
  "hours_cat"
)

moderator_labels <- c(
  outcome_content = "Content domain",
  weeks_duration_cat_1 = "Program duration",
  hours_cat = "Session length"
)

sample_grid <- tibble(
  sample = c("All", "RCTs", "QEDs"),
  study_type_filter = c(NA, "rct", "quasi")
)

## 3. Fit one single-moderator model
fit_single_moderator_model <- function(data, moderator, sample_name, study_type_filter = NA, rho = 0.5) {
  
  dat <- data %>%
    filter(!is.na(.data[[moderator]]))
  
  if (!is.na(study_type_filter)) {
    dat <- dat %>%
      filter(study_type == study_type_filter)
  }
  
  ## Skip if there are too few rows or fewer than 2 categories
  if (nrow(dat) < 2 || length(unique(dat[[moderator]])) < 2) {
    return(NULL)
  }
  
  V <- vcalc(
    vi      = variance_final_d_effect,
    cluster = study_id,
    obs     = effect_id,
    data    = dat,
    rho     = rho
  )
  
  model <- rma.mv(
    yi     = final_d_effect,
    V      = V,
    mods   = as.formula(paste("~", moderator)),
    random = ~ 1 | study_id/effect_id,
    method = "REML",
    data   = dat
  )
  
  list(
    model = model,
    data = dat,
    V = V,
    moderator = moderator,
    sample = sample_name
  )
}

## 4. Fit all MD models: moderator × sample
model_grid_MD <- expand_grid(
  moderator = moderators,
  sample_grid
)

all_models_MD <- pmap(
  model_grid_MD,
  function(moderator, sample, study_type_filter) {
    fit_single_moderator_model(
      data = dat_MD,
      moderator = moderator,
      sample_name = sample,
      study_type_filter = study_type_filter,
      rho = 0.5
    )
  }
)

names(all_models_MD) <- paste(
  model_grid_MD$moderator,
  model_grid_MD$sample,
  sep = "_"
)

all_models_MD <- compact(all_models_MD)

## 5. Extract CLD results
get_cld_results_MD <- function(model_object, adjust = "sidak") {
  
  model <- model_object$model
  moderator <- model_object$moderator
  sample_name <- model_object$sample
  
  emm <- emmeans(
    emmprep(model),
    as.formula(paste("~", moderator))
  )
  
  cld_obj <- cld(
    emm,
    adjust = adjust,
    Letters = letters
  )
  
  df <- as.data.frame(cld_obj)
  
  df$moderator <- moderator
  df$moderator_label <- moderator_labels[[moderator]]
  df$sample <- sample_name
  df$category <- as.character(df[[moderator]])
  df$.group <- gsub(" ", "", df$.group)
  
  df$category_label <- dplyr::case_when(
    moderator == "outcome_content" &
      df$category == "algebra and prealgebra" ~ "AP",
    
    moderator == "outcome_content" &
      df$category == "geometry" ~ "Geom.",
    
    moderator == "outcome_content" &
      df$category == "measurement and data and statistics" ~ "MDS",
    
    moderator == "outcome_content" &
      df$category == "mixed" ~ "Mixed",
    
    moderator == "outcome_content" &
      df$category == "number sense and arithmetic" ~ "NsA",
    
    moderator == "outcome_content" &
      df$category == "rational numbers and fraction" ~ "RnF",
    
    moderator == "weeks_duration_cat_1" &
      df$category == "Long" ~ ">9 weeks",
    
    moderator == "weeks_duration_cat_1" &
      df$category == "Middle" ~ "5-8 weeks",
    
    moderator == "weeks_duration_cat_1" &
      df$category == "Short" ~ "1-4 weeks",
    
    moderator == "hours_cat" ~ gsub("minutes", "min.", df$category),
    
    TRUE ~ df$category
  )
  
  df
}

combined_cld_results_MD <- map_dfr(
  all_models_MD,
  get_cld_results_MD,
  adjust = "sidak"
)

## 6. Count k and number of studies
get_counts_MD <- function(data, moderator, sample_name, study_type_filter = NA) {
  
  dat <- data %>%
    filter(!is.na(.data[[moderator]]))
  
  if (!is.na(study_type_filter)) {
    dat <- dat %>%
      filter(study_type == study_type_filter)
  }
  
  dat %>%
    group_by(
      category = .data[[moderator]]
    ) %>%
    summarise(
      k = n(),
      n_studi = n_distinct(study_id),
      .groups = "drop"
    ) %>%
    mutate(
      moderator = moderator,
      sample = sample_name,
      category = as.character(category),
      count_label = paste0(k, "(", n_studi, ")")
    )
}

count_grid_MD <- expand_grid(
  moderator = moderators,
  sample_grid
)

count_results_MD <- pmap_dfr(
  count_grid_MD,
  function(moderator, sample, study_type_filter) {
    get_counts_MD(
      data = dat_MD,
      moderator = moderator,
      sample_name = sample,
      study_type_filter = study_type_filter
    )
  }
)

## 7. Join counts to CLD results
combined_cld_results_MD <- combined_cld_results_MD %>%
  dplyr::select(-any_of(c("k", "n_studi", "count_label"))) %>%
  mutate(
    category = as.character(category),
    sample = as.character(sample)
  ) %>%
  left_join(
    count_results_MD,
    by = c("moderator", "sample", "category")
  )


## 11. Save estimated effect size table
effect_size_results_MD <- combined_cld_results_MD %>%
  dplyr::select(
    moderator,
    moderator_label,
    sample,
    category,
    category_label,
    emmean,
    SE,
    asymp.LCL,
    asymp.UCL,
    .group,
    k,
    n_studi
  ) %>%
  arrange(
    moderator_label,
    sample,
    category
  )


# get pairwise
get_pairwise_pvalues_MD <- function(x, adjust = "sidak") {
  
  emm_grid <- metafor::emmprep(x$model)
  
  emm <- emmeans::emmeans(
    emm_grid,
    specs = x$moderator
  )
  
  emmeans::contrast(
    emm,
    method = "pairwise",
    adjust = adjust
  ) %>%
    as.data.frame() %>%
    mutate(
      moderator = x$moderator,
      sample = x$sample,
      .before = 1
    )
}

combined_pairwise_results_MD <- purrr::map_dfr(
  all_models_MD,
  get_pairwise_pvalues_MD,
  adjust = "sidak",
  .id = "model_name"
)

pairwise_pvalue_results_MD <- combined_pairwise_results_MD %>%
  dplyr::mutate(
    p_value_label = dplyr::case_when(
      is.na(p.value) ~ NA_character_,
      p.value < .001 ~ "< .001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  dplyr::select(
    model_name,
    moderator,
    sample,
    contrast,
    estimate,
    SE,
    df,
    z.ratio,
    p.value,
    p_value_label
  ) %>%
  dplyr::arrange(
    moderator,
    sample,
    contrast
  )


pairwise_pvalue_results_MD




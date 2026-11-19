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

## 8. Prepare plot variables
combined_cld_results_MD <- combined_cld_results_MD %>%
  mutate(
    sample = factor(
      sample,
      levels = c("All", "QEDs", "RCTs"),
      labels = c("All (MD only)", "QEDs (MD only)", "RCTs (MD only)")
    ),
    moderator_label = factor(
      moderator_label,
      levels = c(
        "Content domain",
        "Program duration",
        "Session length"
      )
    ),
    category_axis = paste(
      moderator_label,
      category_label,
      sep = "___"
    )
  )

category_axis_labels <- function(x) {
  sub(".*___", "", x)
}

common_x_limits_MD <- range(
  c(
    combined_cld_results_MD$asymp.LCL,
    combined_cld_results_MD$asymp.UCL
  ),
  na.rm = TRUE
)

x_padding_MD <- 0.80 * diff(common_x_limits_MD)

common_x_limits_MD <- c(
  common_x_limits_MD[1] - 0.10 * diff(common_x_limits_MD),
  common_x_limits_MD[2] + x_padding_MD
)

combined_cld_results_MD <- combined_cld_results_MD %>%
  mutate(
    letter_x = asymp.UCL + 0.04 * diff(common_x_limits_MD),
    count_x  = common_x_limits_MD[2] - 0.02 * diff(common_x_limits_MD)
  )

## 9. Horizontal MD plot
combined_moderator_plot_horizontal_MD <- ggplot(
  combined_cld_results_MD,
  aes(
    x = emmean,
    y = category_axis,
    color = sample
  )
) +
  geom_vline(
    data = distinct(combined_cld_results_MD, sample),
    aes(
      xintercept = 0,
      color = sample
    ),
    linetype = "dashed",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  geom_errorbarh(
    aes(
      xmin = asymp.LCL,
      xmax = asymp.UCL
    ),
    height = 0.15, , width = 0.08
  ) +
  geom_point(
    size = 3
  ) +
  geom_text(
    aes(
      x = count_x,
      label = count_label
    ),
    size = 3,
    hjust = 1
  ) +
  facet_grid(
    moderator_label ~ sample,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_y_discrete(
    labels = category_axis_labels
  ) +
  coord_cartesian(
    xlim = common_x_limits_MD,
    clip = "off"
  ) +
  
  scale_color_manual(
    values = c(
      "All (MD only)"  = "#a6cee3",
      "QEDs (MD only)" = "dodgerblue",
      "RCTs (MD only)" = "dodgerblue4"
    )
  ) +
  labs(
    x = "Effect size (Cohen's d)",
    y = NULL
  ) +
  theme_test() +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    strip.text.x = element_text(
      face = "bold",
      size = 10
    ),
    strip.text.y.left = element_text(
      face = "bold",
      size = 10,
      angle = 90
    ),
    strip.placement = "outside",
    legend.position = "none",
    axis.text.y = element_text(
      size = 8
    ),
    axis.text.x = element_text(
      size = 9
    ),
    plot.margin = margin(
      t = 5,
      r = 40,
      b = 5,
      l = 5
    )
  )


combined_moderator_plot_horizontal_MD


# Supplementary figure 9
# Fit the models
##############
## Moderators
moderators <- c(
  outcome_content = "outcome_content",
  duration = "weeks_duration_cat_1",
  session = "hours_cat"
)

## Study samples
samples <- tibble(
  sample = c("all", "rct", "quasi"),
  study_type_filter = c(NA, "rct", "quasi")
)

fit_moderator_model <- function(data, moderator, sample, study_type_filter = NA, rho = 0.5) {
  
  dat <- data %>%
    filter(!is.na(.data[[moderator]]))
  
  if (!is.na(study_type_filter)) {
    dat <- dat %>%
      filter(study_type == study_type_filter)
  }
  
  ## Optional safety check
  if (nrow(dat) == 0 || length(unique(dat[[moderator]])) < 2) {
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
    sample = sample
  )
}

model_grid <- tidyr::expand_grid(
  moderator_name = names(moderators),
  moderator = unname(moderators),
  samples
)

moderator_models <- pmap(
  model_grid,
  function(moderator_name, moderator, sample, study_type_filter) {
    fit_moderator_model(
      data = mydata_cog,
      moderator = moderator,
      sample = sample,
      study_type_filter = study_type_filter,
      rho = 0.5
    )
  }
)

names(moderator_models) <- paste(
  model_grid$moderator_name,
  model_grid$sample,
  sep = "_"
)

## Remove skipped NULL models, if any
moderator_models <- purrr::compact(moderator_models)

mod_content_all   <- moderator_models$content_all$model
mod_content_rct   <- moderator_models$content_rct$model
mod_content_quasi <- moderator_models$content_quasi$model

mod_dur_all       <- moderator_models$duration_all$model
mod_dur_rct       <- moderator_models$duration_rct$model
mod_dur_quasi     <- moderator_models$duration_quasi$model

mod_sess_all      <- moderator_models$session_all$model
mod_sess_rct      <- moderator_models$session_rct$model
mod_sess_quasi    <- moderator_models$session_quasi$model
##############

# Supplementary figure 9
# PLOTTING, ADDING LETTER TO SIGNIFICANT POINTS
##############

## Moderators to analyse
moderators <- c(
  content = "outcome_content",
  duration = "weeks_duration_cat_1",
  session = "hours_cat"
)

moderator_labels <- c(
  outcome_content = "Content domain",
  weeks_duration_cat_1 = "Program duration",
  hours_cat = "Session length"
)

fit_single_moderator_model <- function(data, moderator, sample_name) {
  
  dat <- data %>%
    filter(!is.na(.data[[moderator]]))
  
  if (sample_name == "RCTs") {
    dat <- dat %>% filter(study_type == "rct")
  }
  
  if (sample_name == "QEDs") {
    dat <- dat %>% filter(study_type == "quasi")
  }
  
  dat <- dat %>%
    filter(!is.na(.data[[moderator]]))
  
  ## Skip model if fewer than 2 categories are available
  if (length(unique(dat[[moderator]])) < 2) {
    return(NULL)
  }
  
  V <- vcalc(
    vi      = variance_final_d_effect,
    cluster = study_id,
    obs     = effect_id,
    data    = dat,
    rho     = 0.5
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
    moderator = moderator,
    sample = sample_name
  )
}

## These labels must match the count_results and facet labels
sample_names <- c(
  "All",
  "QEDs",
  "RCTs"
)

model_grid <- expand.grid(
  moderator = moderators,
  sample = sample_names,
  stringsAsFactors = FALSE
)

all_models <- pmap(
  list(
    moderator = model_grid$moderator,
    sample = model_grid$sample
  ),
  ~ fit_single_moderator_model(
    data = mydata_cog,
    moderator = ..1,
    sample_name = ..2
  )
)

names(all_models) <- paste(
  model_grid$moderator,
  model_grid$sample,
  sep = "_"
)

## Remove failed/skipped models
all_models <- compact(all_models)

get_cld_results <- function(model_object, adjust = "sidak") {
  
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
  
  ## Create moderator-specific x-axis labels
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

## 1. Extract CLD results
combined_cld_results <- map_dfr(
  all_models,
  get_cld_results,
  adjust = "sidak"
)

## 2. Count k and number of studies for each moderator/category/sample
get_counts <- function(data, moderator, sample_name, study_type_filter = NA) {
  
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

count_grid <- expand_grid(
  moderator = c(
    "outcome_content",
    "weeks_duration_cat_1",
    "hours_cat"
  ),
  sample = c("All", "QEDs", "RCTs")
) %>%
  mutate(
    study_type_filter = case_when(
      sample == "All"  ~ NA_character_,
      sample == "RCTs" ~ "rct",
      sample == "QEDs" ~ "quasi"
    )
  )

count_results <- pmap_dfr(
  count_grid,
  function(moderator, sample, study_type_filter) {
    get_counts(
      data = mydata_cog,
      moderator = moderator,
      sample_name = sample,
      study_type_filter = study_type_filter
    )
  }
)

## 3. Join counts to CLD results
combined_cld_results <- combined_cld_results %>%
  dplyr::select(-any_of(c("k", "n_studi", "count_label"))) %>%
  mutate(
    category = as.character(category),
    sample = as.character(sample)
  ) %>%
  left_join(
    count_results,
    by = c("moderator", "sample", "category")
  )

## 4. Prepare plotting variables
combined_cld_results <- combined_cld_results %>%
  mutate(
    sample = factor(
      sample,
      levels = c("All", "QEDs", "RCTs")
    ),
    moderator_label = factor(
      moderator_label,
      levels = c(
        "Content domain",
        "Program duration",
        "Session length"
      )
    ),
    category_label = case_when(
      moderator == "outcome_content" ~ factor(
        category_label,
        levels = c("AP", "Geom.", "MDS", "Mixed", "NsA", "RnF")
      ) %>% as.character(),
      
      moderator == "weeks_duration_cat_1" ~ factor(
        category_label,
        levels = c(">9 weeks", "1-4 weeks", "5-8 weeks")
      ) %>% as.character(),
      
      moderator == "hours_cat" ~ factor(
        category_label,
        levels = c(">60 min.", "≤30 min.", "31–60 min.")
      ) %>% as.character(),
      
      TRUE ~ category_label
    ),
    category_axis = paste(
      moderator_label,
      category_label,
      sep = "___"
    )
  )

## Manual y-axis order to match expected plot
y_levels <- c(
  paste("Content domain", c("AP", "Geom.", "MDS", "Mixed", "NsA", "RnF"), sep = "___"),
  paste("Program duration", c(">9 weeks", "1-4 weeks", "5-8 weeks"), sep = "___"),
  paste("Session length", c(">60 min.", "≤30 min.", "31–60 min."), sep = "___")
)

combined_cld_results <- combined_cld_results %>%
  mutate(
    category_axis = factor(category_axis, levels = y_levels)
  )

category_axis_labels <- function(x) {
  sub(".*___", "", x)
}

## 5. Define x-axis limits and text positions
common_x_limits <- range(
  c(
    combined_cld_results$asymp.LCL,
    combined_cld_results$asymp.UCL
  ),
  na.rm = TRUE
)

x_padding <- 0.80 * diff(common_x_limits)

common_x_limits <- c(
  common_x_limits[1] - 0.10 * diff(common_x_limits),
  common_x_limits[2] + x_padding
)

combined_cld_results <- combined_cld_results %>%
  mutate(
    letter_x = asymp.UCL + 0.04 * diff(common_x_limits),
    count_x  = common_x_limits[2] - 0.02 * diff(common_x_limits)
  )

## 6. Final horizontal plot
combined_moderator_plot_horizontal <- ggplot(
  combined_cld_results,
  aes(
    x = emmean,
    y = category_axis,
    color = sample
  )
) +
  geom_vline(
    data = distinct(combined_cld_results, sample),
    aes(
      xintercept = 0,
      color = sample
    ),
    linetype = "dashed",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  geom_errorbarh(
    aes(
      xmin = asymp.LCL,
      xmax = asymp.UCL
    ),
    height = 0.15
  ) +
  geom_point(
    size = 2
  ) +
  geom_text(
    aes(
      x = letter_x,
      label = .group
    ),
    size = 3.5,
    hjust = 0
  ) +
  geom_text(
    aes(
      x = count_x,
      label = count_label
    ),
    size = 3,
    hjust = 1
  ) +
  facet_grid(
    moderator_label ~ sample,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_y_discrete(
    labels = category_axis_labels
  ) +
  coord_cartesian(
    xlim = common_x_limits,
    clip = "off"
  ) +
  scale_color_manual(
    values = c(
      "All"  = "#A80000",
      "QEDs" = "#00ACFC",
      "RCTs" = "#262A76"
    )
  ) +
  labs(
    x = "Effect size (Cohen's d)",
    y = NULL
  ) +
  theme_test() +
  theme(
    strip.background = element_rect(
      fill = "grey85",
      color = "grey30"
    ),
    strip.text.x = element_text(
      face = "bold",
      size = 10
    ),
    strip.text.y.left = element_text(
      face = "bold",
      size = 10,
      angle = 90
    ),
    strip.placement = "outside",
    legend.position = "none",
    axis.text.y = element_text(
      size = 8
    ),
    axis.text.x = element_text(
      size = 9
    ),
    axis.title.x = element_text(
      size = 12
    ),
    plot.margin = margin(
      t = 5,
      r = 35,
      b = 5,
      l = 5
    )
  )

combined_moderator_plot_horizontal
##############

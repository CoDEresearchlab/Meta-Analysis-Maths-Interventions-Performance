library(metafor)
library(dplyr)
library(tibble)

# Loading the data
#######
mydata_cog <- read.csv(
  "D:/PhD Nurul/Script submission/mydata_cog_updated.csv",
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

V_cog <- vcalc(
  vi = variance_final_d_effect,
  cluster = study_id,
  obs = effect_id,
  data = mydata_cog,
  rho = 0.5
)
#########


# Grade level
########
# Supplementary table 1 - omnibus test
mod_RQ2_main_grade <- rma.mv(
  yi = final_d_effect,
  V  = V_cog,
  mods = ~factor(grade_cat),            
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(mod_RQ2_main_grade)
count_int3 <- mydata_cog %>%
  group_by(grade_cat) %>%
  summarise(
    k = n(),                              
    n_studi = n_distinct(study_id)        
  )
count_int3
# Supplementary table 1 - estimate per category
mod_RQ2_main_grade_category <- rma.mv(
  yi = final_d_effect,
  V  = V_cog,
  mods = ~0 + factor(grade_cat),            
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(mod_RQ2_main_grade_category)
#########


# Single-grade level
##########
# Supplementary table 3
dat_grade <- mydata_cog %>%
  mutate(
    row_id = row_number(),
    grade = suppressWarnings(as.numeric(as.character(grade_level)))
  ) %>%
  filter(
    !is.na(grade),
    grade %in% 1:6
  )
# Subset variance-covariance matrix to the same observations
V_grade <- V_cog[dat_grade$row_id, dat_grade$row_id]
# Treat Grade 1–6 as categories
dat_grade$grade <- factor(
  dat_grade$grade,
  levels = 1:6
)
# Meta-regression
mod_grade <- rma.mv(
  yi = final_d_effect,
  V = V_grade,
  mods = ~ grade,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = dat_grade
)
summary(mod_grade)

# Count k and n for single-grade level
count_int3_single <- dat_grade %>%
  group_by(grade) %>%
  summarise(
    k = n(),
    n_studi = n_distinct(study_id),
    .groups = "drop"
  )

print(count_int3_single, n = Inf)

# Fit model
mod_grade_category <- rma.mv(
  yi = final_d_effect,
  V = V_grade,
  mods = ~ 0 + grade,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = dat_grade
)
summary(mod_grade_category)
# Extract results
grade_results <- tibble(
  grade_level = levels(factor(dat_grade$grade)),
  cohens_d = as.numeric(mod_grade$b),
  SE = mod_grade$se,
  CI_lower = mod_grade$ci.lb,
  CI_upper = mod_grade$ci.ub,
  p = mod_grade$pval
)

# Calculate k and n for each grade
grade_counts <- dat_grade %>%
  group_by(grade) %>%
  summarise(
    k = n(),
    n = n_distinct(study_id),
    .groups = "drop"
  ) %>%
  rename(grade_level = grade) %>%
  mutate(grade_level = as.character(grade_level))

# Combine and format
grade_results <- grade_results %>%
  mutate(grade_level = as.character(grade_level)) %>%
  left_join(grade_counts, by = "grade_level") %>%
  mutate(
    cohens_d = round(cohens_d, 3),
    SE = round(SE, 3),
    `95% CI` = paste0(
      "[", round(CI_lower, 3),
      ", ", round(CI_upper, 3), "]"
    ),
    p = ifelse(p < .001, "<.001", sprintf("%.3f", p))
  ) %>%
  select(
    grade_level,
    k,
    n,
    cohens_d,
    SE,
    `95% CI`,
    p
  )
grade_results
##########


# Risk status
###########
# Supplementary table 1 - omnibus test
mod_RQ2_main_status <- rma.mv(
  yi = final_d_effect,
  V  = V_cog,
  mods = ~risk_status2,            
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(mod_RQ2_main_status)
count_int15 <- mydata_cog %>%
  group_by(risk_status2) %>%
  summarise(
    k       = n(),                   
    n_studi = n_distinct(study_id)   
  )
count_int15
# Supplementary table 1 - estimate per category
mod_RQ2_main_status_category <- rma.mv(
  yi = final_d_effect,
  V  = V_cog,
  mods = ~0 + risk_status2,            
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(mod_RQ2_main_status_category)
###########


# Gender
##########
# Supplementary table 1
mydata_gender <- mydata_cog %>%
  mutate(
    distribution_male = as.numeric(distribution_male),
    distribution_female = as.numeric(distribution_female),
    total_gender = distribution_male + distribution_female,
    prop_female = distribution_female / total_gender
  ) %>%
  filter(
    !is.na(prop_female),
    is.finite(prop_female),
    total_gender > 0,
    prop_female >= 0,
    prop_female <= 1)
V_gender <- vcalc(
  vi = variance_final_d_effect,
  cluster = study_id,
  obs = effect_id,
  data = mydata_gender,
  rho = 0.5)
# Omnibus test
meta_mod_gender <- rma.mv(
  yi = final_d_effect,
  V = V_gender,
  mods = ~ prop_female,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_gender)
summary(meta_mod_gender)
# Count effect sizes (k) and studies (n) for gender
count_gender <- mydata_gender %>%
  summarise(
    k = n(),
    n_studi = n_distinct(study_id)
  )

print(count_gender)
# Estimate per category
meta_mod_gender_category <- rma.mv(
  yi = final_d_effect,
  V = V_gender,
  mods = ~ 0 + prop_female,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_gender)
summary(meta_mod_gender_category)
##########


# Intervention type
########
meta_mod_intervention <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ factor(intervention_type2),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_intervention)
count_int <- mydata_cog %>%
  group_by(intervention_type2) %>%
  summarise(
    k = n(),                              
    n_studi = n_distinct(study_id)        
  )
count_int
# estimate per category
meta_mod_intervention_category <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ 0 + factor(intervention_type2),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_intervention_category)
########


# session length
########
mydata_cog <- mydata_cog %>%
  mutate(
    hours_cat = case_when(
      is.na(hours_each_duration) ~ NA_character_,
      hours_each_duration <= 0.5 ~ "≤30 minutes",
      hours_each_duration <= 1 ~ "31–60 minutes",
      hours_each_duration > 1  ~ ">60 minutes"
    ),
    hours_cat = factor(hours_cat, levels = c("≤30 minutes", "31–60 minutes", ">60 minutes"))
  )

meta_mod_each_hour1 <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ factor (hours_cat),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_each_hour1)
count_int_SL <- mydata_cog %>%
  group_by(hours_cat) %>%
  summarise(
    k       = n(),                   
    n_studi = n_distinct(study_id)   
  )
count_int_SL
# estimate per category
meta_mod_each_hour1_category <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ 0 + factor (hours_cat),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_each_hour1_category)
########


# program duration
#########
meta_mod_duration2 <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ factor (weeks_duration_cat_1),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_duration2)
count_int4b <- mydata_cog %>%
  group_by(weeks_duration_cat_1) %>%
  summarise(
    k = n(),                             
    n_studi = n_distinct(study_id)        
  )
count_int4b
# estimate per category
meta_mod_duration2_category <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ 0 + factor (weeks_duration_cat_1),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_duration2_category)
#########


# format
#########
meta_mod_format <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ factor(format),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_format)
count_int9 <- mydata_cog %>%
  group_by(format) %>%
  summarise(
    k       = n(),                   
    n_studi = n_distinct(study_id)  
  )
count_int9
meta_mod_format_category <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ 0 + factor(format),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_format_category)
#########


# interventionist
#########
meta_mod_interventionist <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ factor(interventionist_cat),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_interventionist)
count_int10 <- mydata_cog %>%
  group_by(interventionist_cat) %>%
  summarise(
    k       = n(),                   
    n_studi = n_distinct(study_id)   
  )
count_int10
meta_mod_interventionist_category <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ 0 + factor(interventionist_cat),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_interventionist_category)
#########


# outcome
#########
meta_mod_outcome_content <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ factor(outcome_content),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_outcome_content)
count_int11 <- mydata_cog %>%
  group_by(outcome_content) %>%
  summarise(
    k       = n(),                   
    n_studi = n_distinct(study_id)  
  )
count_int11
meta_mod_outcome_content_category <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ 0 + factor(outcome_content),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_outcome_content_category)
#########


library(tibble)
library(metafor)
library(dplyr)
library(tibble)
library(rnaturalearth)
library(sf)
library(ggplot2)
library(scales)

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

V_cog <- vcalc(
  vi = variance_final_d_effect,
  cluster = study_id,
  obs = effect_id,
  data = mydata_cog,
  rho = 0.5
)
#########


# Results of effectiveness by features
########
mod_RQ3_main<- rma.mv(
  yi   = final_d_effect,
  V    = V_cog,
  mods = ~ intervention_type2 +
    hours_cat +
    weeks_duration_cat_1 +
    interventionist_cat +
    format +
    outcome_content,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data  = mydata_cog
)
# Supplementary table 7
metafor::vif(mod_RQ3_main)

mod_RQ3_main1<- rma.mv(
  yi   = final_d_effect,
  V    = V_cog,
  mods = ~ intervention_type2 +
    hours_cat +
    weeks_duration_cat_1 +
    format +
    outcome_content,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data  = mydata_cog
)
summary(mod_RQ3_main1)
# Session length
anova(mod_RQ3_main1, btt = 5:6)
# Programme duration
anova(mod_RQ3_main1, btt = 7:8)
# Content domain
anova(mod_RQ3_main1, btt = 11:15)
#######


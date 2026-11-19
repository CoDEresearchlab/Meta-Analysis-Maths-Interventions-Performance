library(ggplot2)
library(RColorBrewer)
library(emmeans)
library(metafor)
library(dplyr)
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

# Preparing data for Figure 3a
########
studypercountry <- mydata %>%
  group_by(country1) %>%
  summarise(num_of_study = n_distinct(study_id))%>%
  arrange(desc(num_of_study))
# Add continent info to studypercountry by matching country names
studypercountry <- studypercountry %>%
  mutate(continent = countrycode(sourcevar = country1,
                                 origin = "country.name",
                                 destination = "continent"))
# Group by continent and sum number of studies
studypercontinent <- studypercountry %>%
  group_by(continent) %>%
  summarise(num_of_study = sum(num_of_study, na.rm = TRUE)) %>%
  arrange(desc(num_of_study))


# Get world map data
world_map <- map_data("world")
countries <- sort(unique(world_map$region))
# 1. Normalization helper (as you already have)
norm <- function(x) {
  x %>%
    str_trim() %>%
    str_squish() %>%
    str_to_lower()
}

# 2. Build lookup table from map regions
countries_std <- tibble(
  country_std      = countries,           # original region names from map_data("world")
  country_std_norm = norm(countries)      # normalized version
)

# 3. Match your country1 to map regions via normalized names
study_fixed <- studypercountry %>%
  mutate(country1_norm = norm(country1)) %>%
  left_join(
    countries_std,
    by = c("country1_norm" = "country_std_norm")
  ) %>%
  mutate(
    # Use matched map name when available; fall back to original if NA
    country1_map = coalesce(country_std, country1)
  )

# 4. Inspect any that still did not match exactly to map regions
unmatched <- study_fixed %>% filter(is.na(country_std))
unmatched %>% select(country1, num_of_study, continent)

studypercountry <- studypercountry %>%
  mutate(
    country1_map = recode(
      country1,
      "United States"       = "USA",
      "United Kingdom"      = "UK",
      "Trinidad and Tobago" = "Trinidad",
      "Northern Cyprus"     = "Cyprus",
      "Hong Kong"           = "China",
      .default              = country1  # keep others as-is
    )
  )

# ensure one row per country1_map for mapping
studypercountry_map <- studypercountry %>%
  group_by(country1_map, continent) %>%
  summarise(num_of_study = sum(num_of_study), .groups = "drop")

map_joined <- world_map %>%
  left_join(studypercountry_map, by = c("region" = "country1_map"))
#########

# Plot (Figure 3a)
#######
p1 <- ggplot(map_joined,
             aes(x = long, y = lat, group = group, fill = num_of_study)) +
  geom_polygon(color = "black", linewidth = 0.1) +
  coord_fixed(
    ratio = 1.3,
    xlim = c(-180, 180),
    ylim = c(-55, 90),
    expand = TRUE
  ) +
  scale_fill_gradient2(
    low = "#2166AC", 
    mid = "lightblue", 
    high = "blue",
    na.value = "grey80"  # countries with zero studies
  ) +
  labs(
    fill = "Number of\nstudies",
    x = NULL, y = NULL
  ) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),       # legend text size  
        legend.key.size = unit(0.4, "cm"))
#######

# Preparing data for figure 3b
#########
mydata_cog <- mydata_cog %>%
  mutate(
    country1_map = dplyr::recode(
      country1_fix,
      "China, Hong Kong Special Administrative Region"        = "China",
      "Iran (Islamic Republic of)"                            = "Iran", 
      "Russian Federation"                                    = "Russia",
      "State of Palestine"                                    = "Palestine",
      "Trinidad and Tobago"                                   = "Trinidad",
      "United Kingdom of Great Britain and Northern Ireland"  = "UK",
      "United Republic of Tanzania"                           = "Tanzania",
      "United States of America"                              = "USA"
    )
  )
meta_mod_country <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ 0 + factor(country1_map),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_country)

# 1. Extract country estimates from your model
country_estimates <- data.frame(
  country1_map = rownames(summary(meta_mod_country)$b),
  estimate = as.numeric(summary(meta_mod_country)$b),
  stringsAsFactors = FALSE
)
country_estimates <- country_estimates %>%
  mutate(
    country1_map = str_remove(country1_map, "^factor\\([^)]+\\)")
  )
# 2. Verify the names match (should be high now)
countries <- sort(unique(world_map$region))
# 3. Join estimates to world map
map_model <- map_joined %>%
  left_join(country_estimates, by = c("region" = "country1_map"))
# Pooled cohen's d ALL
p3 <- ggplot(map_model, 
             aes(x = long, y = lat, group = group, fill = estimate)) +
  geom_polygon(color = "black", linewidth = 0.1) +
  coord_fixed(
    ratio = 1.3,
    xlim = c(-180, 180),
    ylim = c(-55, 90),
    expand = TRUE
  ) +
  scale_fill_gradient2(
    low = "#2166AC",     # negative effects (blue)
    mid = "#F7F7F7",     # neutral (light grey)
    high = "#D7301F",    # positive effects (red)
    na.value = "grey80",
    name = "Pooled\neffect size\n(Cohen's d)"
  ) +
  labs(
    x = NULL, y = NULL
  ) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),       # legend text size  
        legend.key.size = unit(0.4, "cm"))
##########

# Prepating data for Supplementary table 4 and 5
##########
# Run model for continent
mydata_cog$continent <- factor(mydata_cog$continent, ordered = FALSE)
mydata_cog$continent <- relevel(mydata_cog$continent, ref = "Americas")

mod_RQ2_main7 <- rma.mv(
  yi = final_d_effect,
  V  = V_cog,
  mods = ~ continent,            
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(mod_RQ2_main7)

# Estimated marginal means by continent
emm_continent <- emmeans(
  emmprep(mod_RQ2_main7),
  ~ continent)
# Pairwise comparisons
pairwise_continent <- pairs(
  emm_continent,
  adjust = "sidak")
pairwise_continent_df <- as.data.frame(pairwise_continent)
pairwise_continent_df
# Supplementary table 5
count_int_cont <- mydata_cog %>%
  group_by(continent) %>%
  summarise(
    k       = n(),                   
    n_studi = n_distinct(study_id))
count_int_cont
#######

# Supplementary figure 8
########
mydata_cog$continent <- factor(mydata_cog$continent, ordered = FALSE)
mydata_cog$continent <- relevel(mydata_cog$continent, ref = "Americas")
# Fit the model
mod_RQ2_main7 <- rma.mv(
  yi = final_d_effect,
  V  = V_cog,
  mods = ~ continent,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
cont_levels <- levels(mydata_cog$continent)   
newdat <- data.frame(continent = factor(cont_levels, levels = cont_levels))
X <- model.matrix(~ continent, data = newdat)
b <- coef(mod_RQ2_main7)
cont_est <- tibble(
  continent = cont_levels,
  pred_es   = as.numeric(X %*% b)
)
world <- ne_countries(scale = "medium", returnclass = "sf")
world_map <- world %>%
  mutate(continent5 = case_when(
    continent %in% c("North America", "South America") ~ "Americas",
    continent %in% c("Africa", "Asia", "Europe", "Oceania") ~ continent,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(continent5)) %>%
  left_join(cont_est, by = c("continent5" = "continent"))
mx <- max(world_map$pred_es, na.rm = TRUE)
ggplot(world_map) +
  geom_sf(aes(fill = pred_es), color = "white", linewidth = 0.2) +
  scale_fill_gradient(
    limits = c(0, mx),
    oob = squish,          
    name = "Cohen's d"
  ) +
  coord_sf(crs = "EPSG:4326", expand = FALSE) +
  labs(title = "Pooled Effect Sizes by Continent") +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
  )
#########

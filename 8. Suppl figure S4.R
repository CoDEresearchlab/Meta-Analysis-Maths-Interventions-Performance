library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(stringr)

# Load data
########
mydata <- read.csv(
  "D:/PhD Nurul/Script submission/mydata_processed.csv",
  stringsAsFactors = FALSE
)

mydata1 <- mydata %>%
  mutate(country1 = recode(country1,
                           "US" = "United States",
                           "USA" = "United States",
                           "UK" = "United Kingdom",
                           "England" = "United Kingdom",
                           "Kosovo" = "Serbia",
                           "German" = "Germany",
                           "NIgeria" = "Nigeria"))
########

# Preparing data for suppl figure 4a
#######
rct_data <- mydata1 %>%
  filter(
    study_type == "rct",
    outcome_cat == "cognitive",
    notes_in_notes == "Include"
  ) %>%
  select(study_id, rob2_d1:rob2_d5, final_rob)

rct_long <- rct_data %>%
  pivot_longer(
    cols = rob2_d1:rob2_d5,
    names_to  = "domain",
    values_to = "rating"
  ) %>%
  distinct(study_id, domain, rating)

domain_labels <- c(
  rob2_d1 = "D1-Bias arising from the randomisation process",
  rob2_d2 = "D2-Bias due to deviations from intended interventions",
  rob2_d3 = "D3-Bias due to missing outcome data",
  rob2_d4 = "D4-Bias in measurement of the outcome",
  rob2_d5 = "D5-Bias in selection of the reported result"
)

rct_long <- rct_long %>%
  mutate(
    domain = recode(domain, !!!domain_labels),
    domain = factor(domain, levels = rev(domain_labels))
  )

rct_long <- rct_long %>%
  mutate(
    rating = str_to_lower(str_trim(rating)),
    rating = case_when(
      rating %in% c("low", "low risk", "low risk of bias") ~ "Low risk of bias",
      rating %in% c("some concerns", "some concern") ~ "Some concerns",
      rating %in% c("high", "high risk", "high risk of bias") ~ "High risk of bias",
      rating %in% c("no information", "ni", "not reported") ~ "No information",
      TRUE ~ "No information"
    ),
    rating = factor(
      rating,
      levels = c(
        "Low risk of bias",
        "Some concerns",
        "High risk of bias",
        "No information"
      )
    )
  )

rb_cols <- c(
  "Low risk of bias"  = "dodgerblue4",
  "Some concerns"     = "lightblue",
  "High risk of bias" = "blue",
  "No information"    = "navy"
)
########

# Plot (Supplementary figure 4a)
######
plot_rct <- ggplot(rct_long, aes(x = domain, fill = rating)) +
  geom_bar(position = "fill", colour = "black", width = 0.7) +
  coord_flip() +
  scale_y_continuous(
    labels = percent,
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.025))
  ) +
  scale_fill_manual(values = rb_cols, name = "") +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 9),
    legend.position = "bottom",
    panel.border = element_rect(colour = "black", fill = NA)
  )
#######

# preparing data for figure 4b
#######
quasi_data <- mydata1 %>%
  filter(
    study_type == "quasi",
    notes_in_notes == "Include",
    outcome_cat == "cognitive"
  ) %>%
  select(study_id, robin_d1:robin_d7)

quasi_long <- quasi_data %>%
  pivot_longer(
    cols      = robin_d1:robin_d7,
    names_to  = "domain",
    values_to = "rating"
  ) %>%
  distinct(study_id, domain, rating)

domain_labels_quasi <- c(
  robin_d1 = "D1 – Bias due to confounding",
  robin_d2 = "D2 – Bias in selection of participants into the study",
  robin_d3 = "D3 – Bias in classification of interventions",
  robin_d4 = "D4 – Bias due to deviations from intended interventions",
  robin_d5 = "D5 – Bias due to missing data",
  robin_d6 = "D6 – Bias in measurement of outcomes",
  robin_d7 = "D7 – Bias in selection of the reported result"
)

quasi_long <- quasi_long %>%
  mutate(
    domain = recode(domain, !!!domain_labels_quasi),
    domain = factor(domain, levels = rev(domain_labels_quasi))
  )

quasi_long <- quasi_long %>%
  mutate(
    rating = str_to_lower(str_trim(rating)),
    rating = case_when(
      rating %in% c("low", "low risk", "low risk of bias") ~ "Low risk of bias",
      rating %in% c("moderate", "moderate risk", "moderate risk of bias") ~ "Moderate risk of bias",
      rating %in% c("serious", "serious risk", "serious risk of bias") ~ "Serious risk of bias",
      rating %in% c("critical", "critical risk", "critical risk of bias") ~ "Critical risk of bias",
      TRUE ~ NA_character_
    ),
    rating = factor(
      rating,
      levels = c(
        "Low risk of bias",
        "Moderate risk of bias",
        "Serious risk of bias",
        "Critical risk of bias"
      )
    )
  ) %>%
  filter(!is.na(rating))

rb_cols_quasi <- c(
  "Low risk of bias"      = "dodgerblue4",
  "Moderate risk of bias" = "lightblue",
  "Serious risk of bias"  = "mediumblue",
  "Critical risk of bias" = "midnightblue"
)
##########

# plot (supplementary figure 4b)
#######
plot_quasi <- ggplot(quasi_long, aes(x = domain, fill = rating)) +
  geom_bar(position = "fill", colour = "black", width = 0.7) +
  coord_flip() +
  scale_y_continuous(
    labels = percent,
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.025))
  ) +
  scale_fill_manual(values = rb_cols_quasi, name = "") +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 9),
    legend.position = "bottom",
    panel.border = element_rect(colour = "black", fill = NA)
  )
#######

# combining plots
#########
library(patchwork)
combined_plot <- plot_rct / plot_quasi +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 14, face = "bold")
  )
########
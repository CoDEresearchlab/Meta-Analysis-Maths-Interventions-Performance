library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(emmeans)
library(metafor)
library(dplyr)
library(tibble)

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



mod_RQ3_main1 <- rma.mv(
  yi   = final_d_effect,
  V    = V_cog,
  mods = ~ 0 + intervention_type2 +
    hours_cat +
    weeks_duration_cat_1 +
    format +
    outcome_content,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data  = mydata_cog
)

# Session length
anova(mod_RQ3_main1, btt = 5:6)
# Programme duration
anova(mod_RQ3_main1, btt = 7:8)
# Content domain
anova(mod_RQ3_main1, btt = 11:15)

# Supplementary table 8
summary(mod_RQ3_main1)

emm_mod <- emmprep(mod_RQ3_main1)

moderators <- c("intervention_type2", "hours_cat", "weeks_duration_cat_1", "format", "outcome_content")
moderator_labels <- c("Intervention type", "Session length", "Program duration", "Format", "Content domain")

x_axis_labels <- list(
  intervention_type2 = c(
    "Other" = "Other",
    "Schema-Based Instruction (SBI)" = "SBI",
    "Strategy Instruction (SI)" = "SI",
    "Technology-Based Instruction (TBI)" = "TBI"
  ),
  hours_cat = function(x) gsub("minutes", "min.", x),
  format = c(
    "class" = "Class",
    "individual (one-to-one)" = "Indiv.",
    "small group" = "Small group"
  ),
  weeks_duration_cat_1 = c(
    "Long" = ">9 weeks",
    "Middle" = "5-8 weeks",
    "Short" = "1-4 weeks"
  ),
  outcome_content = c(
    "algebra and prealgebra" = "AP",
    "geometry" = "Geom.",
    "measurement and data and statistics" = "MDS",
    "mixed" = "Mixed",
    "number sense and arithmetic" = "NsA",
    "rational numbers and fraction" = "RnF"
  )
)

get_plot_df <- function(moderator, adjust = "sidak") {
  emm <- emmeans(emm_mod, as.formula(paste("~", moderator)))
  as.data.frame(emm)
}

# Supplementary table 9
pairwise_results <- lapply(moderators, function(m) {
  emm <- emmeans(emm_mod, as.formula(paste("~", m)))
  as.data.frame(pairs(emm, adjust = "sidak"))
})
names(pairwise_results) <- moderators

get_sig_df <- function(moderator, plot_df, pairwise_results) {
  pw <- pairwise_results[[moderator]] %>%
    filter(p.value < 0.05) %>%
    separate(contrast, into = c("group1", "group2"), sep = " - ")
  
  if (nrow(pw) == 0) return(NULL)
  
  rng <- range(c(plot_df$asymp.LCL, plot_df$asymp.UCL), na.rm = TRUE)
  step <- 0.10 * diff(rng)
  top <- max(plot_df$asymp.UCL, na.rm = TRUE)
  
  pw %>%
    mutate(
      xmin = match(group1, plot_df[[moderator]]),
      xmax = match(group2, plot_df[[moderator]]),
      stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01 ~ "**",
        TRUE ~ "*"
      )
    ) %>%
    arrange(xmin, xmax) %>%
    mutate(
      y = top + step * seq_len(n()),
      y_text = y + 0.02 * diff(rng)
    )
}

all_plot_dfs <- lapply(moderators, get_plot_df)
names(all_plot_dfs) <- moderators

common_y_limits <- range(
  unlist(lapply(all_plot_dfs, function(x) c(x$asymp.LCL, x$asymp.UCL))),
  na.rm = TRUE
)

max_n_sig <- max(sapply(moderators, function(m) {
  sum(pairwise_results[[m]]$p.value < 0.05, na.rm = TRUE)
}))

base_diff <- diff(common_y_limits)
common_y_limits <- c(
  common_y_limits[1] - 0.15 * base_diff,
  common_y_limits[2] + (0.15 + 0.12 * max_n_sig) * base_diff
)

plot_moderator <- function(moderator, x_label = NULL, adjust = "sidak",
                           x_axis_labels = NULL, y_limits = NULL) {
  plot_df <- get_plot_df(moderator = moderator, adjust = adjust)
  sig_df <- get_sig_df(moderator, plot_df, pairwise_results)
  if (is.null(x_label)) x_label <- moderator
  
  p <- ggplot(plot_df, aes(x = .data[[moderator]], y = emmean)) +
    #geom_hline(yintercept = pretty(y_limits), color = "lightgrey", alpha = 0.5) +
    geom_point(color = "dodgerblue", size = 2) +
    geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), color = "dodgerblue", width = 0.08) +
    labs(title = x_label, x = NULL, y = NULL) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    theme_test() +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "none",
      axis.text.x = element_text(size = 7),
      plot.margin = margin(5.5, 12, 5.5, 5.5)
    )
  
  if (!is.null(x_axis_labels[[moderator]])) {
    p <- p + scale_x_discrete(labels = x_axis_labels[[moderator]])
  }
  
  if (!is.null(sig_df)) {
    rng <- diff(range(c(plot_df$asymp.LCL, plot_df$asymp.UCL), na.rm = TRUE))
    tick <- 0.03 * rng
    
    p <- p +
      geom_segment(data = sig_df,
                   aes(x = xmin, xend = xmax, y = y, yend = y),
                   inherit.aes = FALSE, color = "black") +
      geom_segment(data = sig_df,
                   aes(x = xmin, xend = xmin, y = y - tick, yend = y),
                   inherit.aes = FALSE, color = "black") +
      geom_segment(data = sig_df,
                   aes(x = xmax, xend = xmax, y = y - tick, yend = y),
                   inherit.aes = FALSE, color = "black") +
      geom_text(data = sig_df,
                aes(x = (xmin + xmax) / 2, y = y_text, label = stars),
                inherit.aes = FALSE, color = "black", size = 3)
  }
  
  p
}

plots <- Map(
  f = plot_moderator,
  moderator = moderators,
  x_label = moderator_labels,
  MoreArgs = list(x_axis_labels = x_axis_labels, y_limits = common_y_limits)
)

plots <- lapply(seq_along(plots), function(i) {
  if (i %% 2 == 1) plots[[i]] + labs(y = "ES") else plots[[i]] + labs(y = NULL)
})

# Figure 4
combined_plot <- plot_grid(
  plotlist = plots,
  ncol = 2,
  align = "none"
)


# Supplementary table 9
pairwise_df <- bind_rows(
  lapply(names(pairwise_results), function(x) {
    as.data.frame(pairwise_results[[x]]) %>%
      mutate(predictor = x)
  })
)
pairwise_df

library(writexl)

write_xlsx(
  pairwise_df,
  "D:/PhD Nurul/Supplementary table 9.xlsx"
)

library(metafor)
library(dplyr)



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
#########

# Running REML model for creating Caterpilar and funnel plots
#######
# Build V matrix
V_cog <- vcalc(
  vi = variance_final_d_effect,
  cluster = study_id,
  obs = effect_id,
  data = mydata_cog,
  rho = 0.5
)

# Run multilevel model (REML)
meta_model_cog <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_model_cog)
# Prediction interval
pred_overall <- predict(meta_model_cog)
pred_overall
#######

# Caterpilar plot (figure 2a)
########
# Ringkasan model untuk caption
sum_m <- summary(meta_model_cog)

pooled_est <- as.numeric(sum_m$beta[1])
ci_lb      <- as.numeric(sum_m$ci.lb[1])
ci_ub      <- as.numeric(sum_m$ci.ub[1])

k_eff      <- nrow(mydata_cog)
k_studies  <- length(unique(mydata_cog$study_id))

# >>> gabung subtitle ke caption + REML tetap ada
caption_text <- sprintf(
  "Vertical lines = 95%% CI; solid line = pooled effect; dashed line = 0\nREML pooled effect: d = %.2f [%.2f; %.2f], k = %d effects from %d studies",
  pooled_est, ci_lb, ci_ub, k_eff, k_studies
)

# Siapkan data untuk plot
plot_dat <- mydata_cog %>%
  mutate(
    se      = sqrt(variance_final_d_effect),
    ci_low  = final_d_effect - 1.96 * se,
    ci_high = final_d_effect + 1.96 * se
  ) %>%
  arrange(final_d_effect) %>%
  mutate(effect_order = row_number())

# Plot
caterpilarplot <- ggplot(plot_dat, aes(x = effect_order, y = final_d_effect)) +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high),
                 color = "#a6cee3", linewidth = 0.7) +
  geom_point(
    color = "#1f78b4",
    size = 1, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = pooled_est, linetype = "solid",
             color = "black", linewidth = 0.8) +
  labs(
    x = "Effect (ordered by size)",
    y = "Cohen's d"
  ) +
  theme_test(base_size = 14) +
  theme(
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank()
  )
########

# Funnel plot (figure 2b)
#######
# Pooled effect
pooled_effect <- as.numeric(coef(meta_model_cog)[1])
# Map colors to study types
col_map <- c(
  quasi = "dodgerblue",
  rct   = "dodgerblue4")
legend_labels <- c(
  quasi = "QEDs",
  rct   = "RCTs")
# Data for points
funnel_df <- data.frame(
  yi = meta_model_cog$yi,
  sei = sqrt(meta_model_cog$vi),
  study_type = meta_model_cog$data$study_type)
# Maximum SE
max_se <- max(funnel_df$sei, na.rm = TRUE)
# Sequence of SE values for contour lines/shading
sei_seq <- seq(0, max_se, length.out = 500)
# z-values for 90%, 95%, 99% contours
z90 <- qnorm(0.95)
z95 <- qnorm(0.975)
z99 <- qnorm(0.995)
# Build contour boundaries
shade_df <- data.frame(
  sei = sei_seq,
  c90_l = pooled_effect - z90 * sei_seq,
  c90_u = pooled_effect + z90 * sei_seq,
  c95_l = pooled_effect - z95 * sei_seq,
  c95_u = pooled_effect + z95 * sei_seq,
  c99_l = pooled_effect - z99 * sei_seq,
  c99_u = pooled_effect + z99 * sei_seq)
# Plot
funnel_gg <- ggplot() +
  # Outer shaded band: 95% to 99%
  geom_ribbon(data = shade_df, aes(y = sei, xmin = c99_l, xmax = c95_l), fill = "grey75") +
  geom_ribbon(data = shade_df, aes(y = sei, xmin = c95_u, xmax = c99_u), fill = "grey75") +
  # Middle shaded band: 90% to 95%
  geom_ribbon(data = shade_df, aes(y = sei, xmin = c95_l, xmax = c90_l), fill = "grey85") +
  geom_ribbon(data = shade_df, aes(y = sei, xmin = c90_u, xmax = c95_u), fill = "grey85") +
  # Inner shaded band: center to 90%
  geom_ribbon(data = shade_df, aes(y = sei, xmin = c90_l, xmax = c90_u), fill = "grey93") +
  # Contour lines
  geom_line(data = shade_df, aes(x = c90_l, y = sei), linetype = "dashed", colour = "grey50") +
  geom_line(data = shade_df, aes(x = c90_u, y = sei), linetype = "dashed", colour = "grey50") +
  geom_line(data = shade_df, aes(x = c95_l, y = sei), linetype = "dashed", colour = "grey40") +
  geom_line(data = shade_df, aes(x = c95_u, y = sei), linetype = "dashed", colour = "grey40") +
  geom_line(data = shade_df, aes(x = c99_l, y = sei), linetype = "dashed", colour = "grey30") +
  geom_line(data = shade_df, aes(x = c99_u, y = sei), linetype = "dashed", colour = "grey30") +
  # Vertical pooled-effect line
  geom_vline(xintercept = pooled_effect, colour = "black", linewidth = 0.8) +
  # Study points
  geom_point(data = funnel_df, aes(x = yi, y = sei, colour = study_type), size = 2) +
  scale_colour_manual(name = "Study design", values = col_map, breaks = names(col_map),
                      labels = unname(legend_labels[names(col_map)])) +
  scale_shape_manual(name = "Study design", values = c(quasi = 15, rct = 17),
                     breaks = names(col_map), labels = unname(legend_labels[names(col_map)])) +
  # Reverse y-axis for funnel plot
  scale_y_reverse(limits = c(max_se, 0), expand = c(0, 0)) +
  labs(x = "Cohen's d", y = "Standard error") +
  theme_test(base_size = 14) +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 11),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = 0.1),
        legend.key = element_blank())
funnel_gg
#########

# Pooled effect size (figure 2c)
# We will be running REML models with all study designs, QED only, and RCT only
# In each of them we will be running three different models...
## 1. with all data available
## 2. excluding outliers from dataset
## 3. including high-quality data only
########
# Model for all data
# V matrix
V_overall <- vcalc(
  vi      = variance_final_d_effect,
  cluster = study_id,
  obs     = effect_id,
  data    = mydata_cog,
  rho     = 0.5
)
meta_overall_full <- rma.mv(
  yi     = final_d_effect,
  V      = V_overall,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = mydata_cog
)
summary(meta_overall_full)
# Model for all data excluding outliers
mydata_cog$stud_resid_overall <- rstandard(meta_overall_full)$z
mydata_cog$abs_z_overall <- abs(mydata_cog$stud_resid_overall)
out_overall <- mydata_cog %>%
  arrange(desc(abs_z_overall)) %>%
  select(study_id, effect_id, final_d_effect, stud_resid_overall, abs_z_overall)
out_overall_ids <- mydata_cog %>%
  filter(abs_z_overall > 2.5) %>%
  pull(effect_id)
mydata_cog_no_out <- mydata_cog %>%
  filter(!effect_id %in% out_overall_ids)
V_overall_no_out <- vcalc(
  vi      = variance_final_d_effect,
  cluster = study_id,
  obs     = effect_id,
  data    = mydata_cog_no_out,
  rho     = 0.5
)
meta_overall_no_out <- rma.mv(
  yi     = final_d_effect,
  V      = V_overall_no_out,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = mydata_cog_no_out
)
# Model for all data with high-quality data only
mydata_cog_hq <- mydata_cog %>%
  filter(final_rob_cat %in% c("Low risk", "Moderate"))
V_overall_hq <- vcalc(
  vi      = variance_final_d_effect,
  cluster = study_id,
  obs     = effect_id,
  data    = mydata_cog_hq,
  rho     = 0.5
)
meta_overall_hq <- rma.mv(
  yi     = final_d_effect,
  V      = V_overall_hq,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = mydata_cog_hq
)

# Model for RCT
mydata_rct <- mydata_cog %>%
  filter(study_type == "rct")
V_rct <- vcalc(
  vi      = variance_final_d_effect,
  cluster = study_id,
  obs     = effect_id,
  data    = mydata_rct,
  rho     = 0.5
)

meta_rct_full <- rma.mv(
  yi     = final_d_effect,
  V      = V_rct,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = mydata_rct
)

# Model for RCT excluding outliers
mydata_rct$stud_resid_rct <- rstandard(meta_rct_full)$z
mydata_rct$abs_z_rct <- abs(mydata_rct$stud_resid_rct)
out_rct <- mydata_rct %>%
  arrange(desc(abs_z_rct)) %>%
  select(study_id, effect_id, final_d_effect, stud_resid_rct, abs_z_rct)
out_rct_ids <- mydata_rct %>%
  filter(abs_z_rct > 2.5) %>%
  pull(effect_id)
mydata_rct_no_out <- mydata_rct %>%
  filter(!effect_id %in% out_rct_ids)
V_rct_no_out <- vcalc(
  vi      = variance_final_d_effect,
  cluster = study_id,
  obs     = effect_id,
  data    = mydata_rct_no_out,
  rho     = 0.5
)
meta_rct_no_out <- rma.mv(
  yi     = final_d_effect,
  V      = V_rct_no_out,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = mydata_rct_no_out
)

# Model for RCT with high-quality data only
mydata_rct_hq <- mydata_rct %>%
  filter(final_rob_cat %in% c("Low risk", "Moderate"))
V_rct_hq <- vcalc(
  vi      = variance_final_d_effect,
  cluster = study_id,
  obs     = effect_id,
  data    = mydata_rct_hq,
  rho     = 0.5
)
meta_rct_hq <- rma.mv(
  yi     = final_d_effect,
  V      = V_rct_hq,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = mydata_rct_hq
)

# Model for QUASI all dataset
mydata_quasi <- mydata_cog %>%
  filter(study_type == "quasi")
V_quasi <- vcalc(
  vi      = variance_final_d_effect,
  cluster = study_id,
  obs     = effect_id,
  data    = mydata_quasi,
  rho     = 0.5
)
meta_quasi_full <- rma.mv(
  yi     = final_d_effect,
  V      = V_quasi,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = mydata_quasi
)

# Model for QUASI excluding outliers
mydata_quasi$stud_resid_quasi <- rstandard(meta_quasi_full)$z
mydata_quasi$abs_z_quasi <- abs(mydata_quasi$stud_resid_quasi)
out_quasi <- mydata_quasi %>%
  arrange(desc(abs_z_quasi)) %>%
  select(study_id, effect_id, final_d_effect, stud_resid_quasi, abs_z_quasi)
out_quasi_ids <- mydata_quasi %>%
  filter(abs_z_quasi > 2.5) %>%
  pull(effect_id)
mydata_quasi_no_out <- mydata_quasi %>%
  filter(!effect_id %in% out_quasi_ids)
V_quasi_no_out <- vcalc(
  vi      = variance_final_d_effect,
  cluster = study_id,
  obs     = effect_id,
  data    = mydata_quasi_no_out,
  rho     = 0.5
)
meta_quasi_no_out <- rma.mv(
  yi     = final_d_effect,
  V      = V_quasi_no_out,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = mydata_quasi_no_out
)

# Model for QUASI with high-quality data only
mydata_quasi_hq <- mydata_quasi %>%
  filter(final_rob_cat %in% c("Low risk", "Moderate"))
V_quasi_hq <- vcalc(
  vi      = variance_final_d_effect,
  cluster = study_id,
  obs     = effect_id,
  data    = mydata_quasi_hq,
  rho     = 0.5
)
meta_quasi_hq <- rma.mv(
  yi     = final_d_effect,
  V      = V_quasi_hq,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = mydata_quasi_hq
)

# Compile all models
get_es <- function(model, group_label, variant_label) {
  s <- summary(model)
  tibble(
    Group   = group_label,       # All / RCT / Quasi
    Variant = variant_label,     # Overall / No outlier / High quality
    est     = as.numeric(s$beta[1]),
    ci_lb   = as.numeric(s$ci.lb[1]),
    ci_ub   = as.numeric(s$ci.ub[1])
  )
}
plot_df <- bind_rows(
  # ALL
  get_es(meta_overall_full,   "All",  "Overall"),
  get_es(meta_overall_no_out, "All",  "No outlier"),
  get_es(meta_overall_hq,     "All",  "High quality"),
  
  # RCT
  get_es(meta_rct_full,       "RCT",  "Overall"),
  get_es(meta_rct_no_out,     "RCT",  "No outlier"),
  get_es(meta_rct_hq,         "RCT",  "High quality"),
  
  # QUASI
  get_es(meta_quasi_full,     "Quasi","Overall"),
  get_es(meta_quasi_no_out,   "Quasi","No outlier"),
  get_es(meta_quasi_hq,       "Quasi","High quality")
) %>%
  mutate(
    Group = factor(Group, levels = c("All", "RCT", "Quasi")),
    Variant = factor(Variant, levels = c("Overall", "No outlier", "High quality"))
  )
group_cols <- c("All" = "#a6cee3", "RCT" = "dodgerblue4", "Quasi" = "dodgerblue")
variant_shapes <- c("Overall" = 16, "No outlier" = 17, "High quality" = 15)

pos <- position_dodge(width = 0.6)

plot_df <- plot_df %>%
  dplyr::mutate(
    Group = factor(Group, levels = c("Quasi", "RCT", "All")),
    Variant = factor(Variant, levels = c("High quality", "No outlier", "Overall"))
  )

p <- ggplot(plot_df, aes(x = est, y = Group, colour = Group, shape = Variant)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lb, xmax = ci_ub), height = 0.15, position = pos, linewidth = 0.8) +
  geom_point(position = pos, size = 3) +
  scale_colour_manual(values = group_cols, labels = c("Quasi" = "QEDs", "RCT" = "RCTs", "All" = "All"), guide = "none") +
  scale_shape_manual(values = variant_shapes, breaks = c("Overall", "No outlier", "High quality")) +
  scale_y_discrete(labels = c("Quasi" = "QEDs", "RCT" = "RCTs", "All" = "All")) +
  labs(x = "Pooled effect size (Cohen's d)", y = NULL, shape = NULL) +
  theme_test(base_size = 14) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

p
########

# Forest plot (Supplementary figure S1a-l)
#######
# Overall forest plot (optional: for on-screen viewing)
forest(
  meta_model_cog,
  slab = paste("Study", mydata_cog$author_year),
  xlab = "Effect size (Cohen's d)",
  mlab = "Pooled effect (Cognitive)"
)

# --- Split the data into several chunks ---
chunk_size <- 48
n_effects <- length(meta_model_cog$yi)
chunk_indices <- split(seq_len(n_effects), ceiling(seq_len(n_effects) / chunk_size))

# --- Loop through each chunk ---
for (i in seq_along(chunk_indices)) {
  inds <- chunk_indices[[i]]
  
  effect_sizes <- meta_model_cog$yi[inds]
  variances <- diag(meta_model_cog$V)[inds]
  slabs <- paste(mydata_cog$author_year[inds])
  
  # Specify the file path for each chunk
  file_name <- paste0("D:/9. PhD/1. PhD/2nd Year/2. Nurul_META ANALYSIS/1. data extraction/analysis_2/forest_plot_chunk_study", i, ".png")
  
  # Save the plot as a PNG file
  png(file = file_name, width = 800, height = 600)
  
  forest(
    x = effect_sizes,
    vi = variances,
    slab = slabs,
    xlab = "Effect size (Cohen's d)",
    main = paste("Pooled effect (Cognitive) - chunk", i),
    cex.main = 0.7
  )
  
  dev.off()  # Close the graphics device to save the file
  
  cat("✅ Saved:", file_name, "\n")
  
  # Optional: pause between chunks
  # readline(prompt = "Press [Enter] to continue to the next chunk...")
}

# --- (Optional) Additional plot: overall pooled effect ---
png(file = "D:/9. PhD/1. PhD/2nd Year/2. Nurul_META ANALYSIS/1. data extraction/analysis_2/forest_plot_overall.png",
    width = 800, height = 600)

forest(
  meta_model_cog,
  slab = paste("Study", mydata_cog$study_id),
  xlab = "Effect size (Cohen's d)",
  mlab = "Overall pooled effect (Cognitive)",
  main = "Overall pooled effect (Cognitive)"
)

dev.off()
########

# Supplementary figure 2
########
meta_model_cog_no_out <- rma.mv(
  yi = final_d_effect,
  V = V_cog_no_out,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog_no_out
)
summary(meta_model_cog_no_out)

library(metafor)
library(orchaRd)
library(gridExtra)
library(ggplot2)
orch1 <- orchard_plot(
  meta_model_cog, 
  mod = "1",           # intercept only
  group = "study_id",
  xlab = "Cohens' d"
) + ggtitle("Full Model")

orch2 <- orchard_plot(
  meta_model_cog_no_out, 
  mod = "1",
  group = "study_id",
  xlab = "Cohens' d"
) + ggtitle("No Outliers")

grid.arrange(orch1, orch2, ncol = 2)
########









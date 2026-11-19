

# Load data
##########
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


# PUB BIAS: OVERALL - NO OUTLIERS (Suppl. figure S3)
########
library(clubSandwich)
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

mydata_cog$stud_resid_overall <- rstandard(meta_overall_full)$z
mydata_cog$abs_z_overall <- abs(mydata_cog$stud_resid_overall)

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

dat <- mydata_cog_no_out
Vmat <- V_overall_no_out   # ganti jika nama V beda (mis. V_cog_no_out)

dat <- dat %>%
  mutate(
    sei       = sqrt(diag(Vmat)),
    precision = 1 / sei,
    SND       = final_d_effect / sei
  )

# Egger-type MV + CR2
egger_overall_noout <- rma.mv(
  yi     = SND,
  V      = 1,
  mods   = ~ precision,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = dat
)

egger_overall_noout_CR2 <- coef_test(
  egger_overall_noout,
  vcov    = "CR2",
  cluster = dat$study_id,
  test    = "Satterthwaite"
)
print(egger_overall_noout_CR2)

# Funnel plots (pakai model mv no-outliers supaya funnel() bisa langsung)
meta_overall_noout <- rma.mv(
  yi     = final_d_effect,
  V      = Vmat,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = dat
)
summary(meta_overall_noout)

col_map <- c(quasi="dodgerblue", rct="dodgerblue4")
idx <- meta_overall_noout$not.na
stype <- tolower(trimws(dat$study_type[idx]))
pt_col <- unname(col_map[stype]); pt_col[is.na(pt_col)] <- "grey20"

funnel(meta_overall_noout, yaxis="sei", level=c(90,95,99),
       xlab="Effect size (Cohen's d)", main="Funnel - Overall NO OUTLIERS",
       pch=19, col=pt_col)
legend("topright", legend=c("QEDs","RCTs"),
       col=unname(col_map), pch=19, bty="n", title="Study Type", inset=0.02)

funnel(meta_overall_noout, yaxis="sei", level=c(90,95,99),
       contour=c(0.90,0.95,0.99), col.contour=c("grey75","grey60","grey40"),
       xlab="Effect size (Cohen's d)", main="Contour Funnel - Overall NO OUTLIERS",
       pch=19, col=pt_col)
legend("topright", legend=c("QEDs","RCTs"),
       col=unname(col_map), 
       pch=19, 
       bty="n", 
       title="Study Type", 
       x.intersp = 0.4,   # jarak simbol-ke-teks
       y.intersp = 0.48,   # jarak antar baris
       text.width = 0.5)  # lebar kolom teks 



# Trim-and-fill on aggregated
agg_overall_noout <- dat %>%
  group_by(study_id) %>%
  summarise(
    yi = weighted.mean(final_d_effect, w = 1/variance_final_d_effect, na.rm = TRUE),
    vi = 1 / sum(1/variance_final_d_effect, na.rm = TRUE),
    .groups = "drop"
  )

res_agg_overall_noout <- rma(yi, vi, data=agg_overall_noout, method="REML")
tf_overall_noout <- trimfill(res_agg_overall_noout)
summary(tf_overall_noout)
##########


# PUB BIAS: OVERALL - HIGH QUALITY (Suppl. figure S5)
##########

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

dat <- mydata_cog_hq
Vmat <- V_overall_hq   # ganti jika nama V beda (mis. V_cog_lm)

dat <- dat %>%
  mutate(
    sei       = sqrt(diag(Vmat)),
    precision = 1 / sei,
    SND       = final_d_effect / sei
  )

egger_overall_hq <- rma.mv(
  yi     = SND,
  V      = 1,
  mods   = ~ precision,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = dat
)

egger_overall_hq_CR2 <- coef_test(
  egger_overall_hq,
  vcov    = "CR2",
  cluster = dat$study_id,
  test    = "Satterthwaite"
)
print(egger_overall_hq_CR2)

meta_overall_hq_mv <- rma.mv(
  yi     = final_d_effect,
  V      = Vmat,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = dat
)
summary(meta_overall_hq_mv)


col_map <- c(quasi="dodgerblue", rct="dodgerblue4")
idx <- meta_overall_hq_mv$not.na
stype <- tolower(trimws(dat$study_type[idx]))
pt_col <- unname(col_map[stype]); pt_col[is.na(pt_col)] <- "grey20"

funnel(meta_overall_hq_mv, yaxis="sei", level=c(90,95,99),
       xlab="Effect size (Cohen's d)", main="Funnel - Overall HIGH QUALITY",
       pch=19, col=pt_col)
legend("topright", legend=c("QEDs","RCTs"),
       col=unname(col_map), pch=19, bty="n", title="Study Type", inset=0.02)

funnel(meta_overall_hq_mv, yaxis="sei", level=c(90,95,99),
       contour=c(0.90,0.95,0.99), col.contour=c("grey75","grey60","grey40"),
       xlab="Effect size (Cohen's d)", main="Contour Funnel - Overall HIGH QUALITY",
       pch=19, col=pt_col)
legend("topright", legend=c("QEDs","RCTs"),
       col=unname(col_map), pch=19, bty="n", title="Study Type",
       x.intersp=0.4, y.intersp=0.48, text.width=0.5)

agg_overall_hq <- dat %>%
  group_by(study_id) %>%
  summarise(
    yi = weighted.mean(final_d_effect, w = 1/variance_final_d_effect, na.rm = TRUE),
    vi = 1 / sum(1/variance_final_d_effect, na.rm = TRUE),
    .groups = "drop"
  )

res_agg_overall_hq <- rma(yi, vi, data=agg_overall_hq, method="REML")
tf_overall_hq <- trimfill(res_agg_overall_hq)
summary(tf_overall_hq)
##########


# summary of moderator by study design
#######
meta_mod_design <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ factor(study_type),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_design)

count_int13 <- mydata_cog %>%
  group_by(study_type) %>%
  summarise(
    k       = n(),                  
    n_studi = n_distinct(study_id)   
  )

meta_mod_design_category <- rma.mv(
  yi = final_d_effect,
  V = V_cog,
  mods = ~ 0 + factor(study_type),
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data = mydata_cog
)
summary(meta_mod_design_category)
##########


# PUB BIAS: RCT - FULL (Suppl. figure S6)
##########
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

summary(meta_rct_full)
predict(meta_rct_full)


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

summary(meta_rct_no_out)
predict(meta_rct_no_out)


dat <- mydata_rct
Vmat <- V_rct

dat <- dat %>%
  mutate(
    sei       = sqrt(diag(Vmat)),
    precision = 1 / sei,
    SND       = final_d_effect / sei
  )

egger_rct_full <- rma.mv(
  yi     = SND,
  V      = 1,
  mods   = ~ precision,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = dat
)

egger_rct_full_CR2 <- coef_test(
  egger_rct_full,
  vcov    = "CR2",
  cluster = dat$study_id,
  test    = "Satterthwaite"
)
print(egger_rct_full_CR2)

meta_rct_full_mv <- rma.mv(
  yi     = final_d_effect,
  V      = Vmat,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = dat
)
summary(meta_rct_full_mv)


col_map <- c(rct="dodgerblue4")
idx <- meta_rct_full_mv$not.na
stype <- rep("rct", sum(idx))              # karena ini model RCT saja
pt_col <- unname(col_map[stype])

funnel(meta_rct_full_mv, yaxis="sei", level=c(90,95,99),
       xlab="Effect size (Cohen's d)", main="Funnel - RCT FULL",
       pch=19, col=pt_col)
legend("topright", legend="RCTs", col=col_map["rct"],
       pch=19, bty="n", title="Study Type", inset=0.02)

funnel(meta_rct_full_mv, yaxis="sei", level=c(90,95,99),
       contour=c(0.90,0.95,0.99), col.contour=c("grey75","grey60","grey40"),
       xlab="Effect size (Cohen's d)", main="Contour Funnel - RCT FULL",
       pch=19, col=pt_col)
legend("topright", legend="RCTs", col=col_map["rct"],
       pch=19, bty="n", title="Study Type",
       x.intersp=0.4, y.intersp=0.48, text.width=0.5)

agg_rct_full <- dat %>%
  group_by(study_id) %>%
  summarise(
    yi = weighted.mean(final_d_effect, w = 1/variance_final_d_effect, na.rm = TRUE),
    vi = 1 / sum(1/variance_final_d_effect, na.rm = TRUE),
    .groups = "drop"
  )

res_agg_rct_full <- rma(yi, vi, data=agg_rct_full, method="REML")
tf_rct_full <- trimfill(res_agg_rct_full)
summary(tf_rct_full)
###########


# PUB BIAS: QUASI - FULL (Suppl. figure S7)
########

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

summary(meta_quasi_full)
predict(meta_quasi_full)

dat <- mydata_quasi
Vmat <- V_quasi

dat <- dat %>%
  mutate(
    sei       = sqrt(diag(Vmat)),
    precision = 1 / sei,
    SND       = final_d_effect / sei
  )

egger_quasi_full <- rma.mv(
  yi     = SND,
  V      = 1,
  mods   = ~ precision,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = dat
)

egger_quasi_full_CR2 <- coef_test(
  egger_quasi_full,
  vcov    = "CR2",
  cluster = dat$study_id,
  test    = "Satterthwaite"
)
print(egger_quasi_full_CR2)

meta_quasi_full_mv <- rma.mv(
  yi     = final_d_effect,
  V      = Vmat,
  random = ~ 1 | study_id/effect_id,
  method = "REML",
  data   = dat
)
summary(meta_quasi_full_mv)

col_map <- c(quasi="dodgerblue")
idx <- meta_quasi_full_mv$not.na
pt_col <- rep(unname(col_map["quasi"]), sum(idx))

funnel(meta_quasi_full_mv, yaxis="sei", level=c(90,95,99),
       xlab="Effect size (Cohen's d)", main="Funnel - QUASI FULL",
       pch=19, col=pt_col)
legend("topright", legend="QEDs", col=col_map["quasi"],
       pch=19, bty="n", title="Study Type", inset=0.02)

funnel(meta_quasi_full_mv, yaxis="sei", level=c(90,95,99),
       contour=c(0.90,0.95,0.99), col.contour=c("grey75","grey60","grey40"),
       xlab="Effect size (Cohen's d)", main="Contour Funnel - QUASI FULL",
       pch=19, col=pt_col)
legend("topright", legend="QEDs", col=col_map["quasi"],
       pch=19, bty="n", title="Study Type",
       x.intersp=0.4, y.intersp=0.48, text.width=0.5)

agg_quasi_full <- dat %>%
  group_by(study_id) %>%
  summarise(
    yi = weighted.mean(final_d_effect, w = 1/variance_final_d_effect, na.rm = TRUE),
    vi = 1 / sum(1/variance_final_d_effect, na.rm = TRUE),
    .groups = "drop"
  )

res_agg_quasi_full <- rma(yi, vi, data=agg_quasi_full, method="REML")
tf_quasi_full <- trimfill(res_agg_quasi_full)
summary(tf_quasi_full)
#########
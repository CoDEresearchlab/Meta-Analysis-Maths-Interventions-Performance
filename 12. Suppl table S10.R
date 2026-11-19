library(dplyr)

# ============================================================
# 1. LOAD DATA
# ============================================================

mydata_cog <- read.csv(
  "D:/PhD Nurul/Script submission/mydata_cog.csv",
  stringsAsFactors = FALSE
)

# Replace blank cells with NA
mydata_cog[] <- lapply(mydata_cog, function(x) {
  if (is.character(x)) {
    x[trimws(x) == ""] <- NA
  }
  x
})

# ============================================================
# 2. COUNT STUDIES (n) AND EFFECT SIZES (k)
# ============================================================

count_risk_format <- mydata_cog %>%
  filter(!is.na(risk_status2), !is.na(format)) %>%
  group_by(risk_status2, format) %>%
  summarise(
    k = n(),
    n = n_distinct(study_id),
    .groups = "drop"
  )

print(count_risk_format, n = Inf)

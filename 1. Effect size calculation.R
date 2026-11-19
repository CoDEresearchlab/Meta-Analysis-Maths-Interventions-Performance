library(readxl)
library(dplyr)
library(tidyr)
library(dplyr)
library(stringr)
library(countrycode)
library(zoo)


# LOAD DATA
#######
# Load as csv
mydata <- read.csv(
  "D:/PhD Nurul/Script submission/mydata.csv",
  stringsAsFactors = FALSE
)

# Load as xlsx
mydata <- read_excel(
  "D:/PhD Nurul/Script submission/mydata.xlsx"
)
#######


# CONVERT HEDGE'S G TO COHEN'S D
########
mydata <- mydata %>%
  mutate(
    J = ifelse(
      effect_size_calculation_method == "g" & notes_in_notes == "Include",
      1 - (3 / (4 * (n_exp_post + n_control_post) - 9)),
      NA_real_ 
    )
  )
#######

# CREATE FINAL D EFFECT SIZE
########
mydata <- mydata %>%
  mutate(
    final_d_effect = case_when(
      
      # 1️⃣ Use existing d values (prefer pre-post first, then post-only)
      notes_in_notes == "Include" & !is.na(d_effect_pre_post) ~ d_effect_pre_post,
      notes_in_notes == "Include" & is.na(d_effect_pre_post) & !is.na(d_effect_post) ~ d_effect_post,
      
      # 2️⃣ Convert from Hedges' g → Cohen's d
      # For pre-post comparison
      notes_in_notes == "Include" & !is.na(es_comparison_between_pre_post) &
        effect_size_calculation_method == "g" & !is.na(J) ~ es_comparison_between_pre_post / J,
      
      # For post-only experimental vs control comparison
      notes_in_notes == "Include" & !is.na(es_post_only_exp_con) &
        effect_size_calculation_method_1 == "g" &
        !is.na(n_exp_post) & !is.na(n_control_post) ~
        es_post_only_exp_con / (1 - (3 / (4 * (n_exp_post + n_control_post) - 9))),
      
      # 3️⃣ Convert r(from U via Z) → Cohen's d
      notes_in_notes == "Include" & !is.na(es_comparison_between_pre_post) &
        effect_size_calculation_method == "r" ~
        2 * es_comparison_between_pre_post / sqrt(pmax(0, 1 - es_comparison_between_pre_post^2)),
      
      notes_in_notes == "Include" & !is.na(es_post_only_exp_con) &
        effect_size_calculation_method_1 == "r" ~
        2 * es_post_only_exp_con / sqrt(pmax(0, 1 - es_post_only_exp_con^2)),
      
      # 4️⃣ Convert from standardized regression coefficient (β) → d
      notes_in_notes == "Include" & !is.na(es_comparison_between_pre_post) &
        effect_size_calculation_method == "Standardized regression coefficient β" ~
        2 * es_comparison_between_pre_post / sqrt(pmax(0, 1 - (es_comparison_between_pre_post)^2)),
      
      notes_in_notes == "Include" & !is.na(es_post_only_exp_con) &
        effect_size_calculation_method_1 == "Standardized regression coefficient β" ~
        2 * es_post_only_exp_con / sqrt(pmax(0, 1 - (es_post_only_exp_con)^2)),
      
      # 5️⃣ Convert from eta squared / partial eta squared → d
      notes_in_notes == "Include" & !is.na(es_comparison_between_pre_post) &
        effect_size_calculation_method %in% c("eta squared", "partial eta squared") ~
        2 * sqrt(pmax(0, es_comparison_between_pre_post / (1 - es_comparison_between_pre_post))),
      
      notes_in_notes == "Include" & !is.na(es_post_only_exp_con) &
        effect_size_calculation_method_1 %in% c("eta squared", "partial eta squared") ~
        2 * sqrt(pmax(0, es_post_only_exp_con / (1 - es_post_only_exp_con))),
      
      # 6️⃣ Default → NA if not matched or not included
      TRUE ~ NA_real_
    )
  )
#########

# Step 1. Fill missing SDs from SE (standard error)
#######
mydata <- mydata %>%
  mutate(
    pretest_sd_exp = if_else(
      notes_in_notes == "Include" & is.na(pretest_sd_exp) & !is.na(se) & !is.na(n_exp),
      se * sqrt(n_exp),
      pretest_sd_exp
    ),
    pretest_sd_control = if_else(
      notes_in_notes == "Include" & is.na(pretest_sd_control) & !is.na(se_1) & !is.na(n_control),
      se_1 * sqrt(n_control),
      pretest_sd_control
    ),
    post_sd_exp = if_else(
      notes_in_notes == "Include" & is.na(post_sd_exp) & !is.na(se_2) & !is.na(n_exp_post),
      se_2 * sqrt(n_exp_post),
      post_sd_exp
    ),
    post_sd_control = if_else(
      notes_in_notes == "Include" & is.na(post_sd_control) & !is.na(se_3) & !is.na(n_control_post),
      se_3 * sqrt(n_control_post),
      post_sd_control
    )
  )
########

# Step 2. Fill missing SDs from SS (sum of squares)
###########
mydata <- mydata %>%
  mutate(
    pretest_sd_exp = if_else(
      notes_in_notes == "Include" & is.na(pretest_sd_exp) & !is.na(ss) & !is.na(n_exp),
      sqrt(ss / (n_exp - 1)),
      pretest_sd_exp
    ),
    pretest_sd_control = if_else(
      notes_in_notes == "Include" & is.na(pretest_sd_control) & !is.na(ss_1) & !is.na(n_control),
      sqrt(ss_1 / (n_control - 1)),
      pretest_sd_control
    ),
    post_sd_exp = if_else(
      notes_in_notes == "Include" & is.na(post_sd_exp) & !is.na(ss_2) & !is.na(n_exp_post),
      sqrt(ss_2 / (n_exp_post - 1)),
      post_sd_exp
    ),
    post_sd_control = if_else(
      notes_in_notes == "Include" & is.na(post_sd_control) & !is.na(ss_3) & !is.na(n_control_post),
      sqrt(ss_3 / (n_control_post - 1)),
      post_sd_control
    )
  )
############


# STEP 3: Compute pooled SDs for pre–post and post-only
########
mydata <- mydata %>%
  mutate(
    sd_pooled_pre = ifelse(
      notes_in_notes == "Include" &
        !is.na(pretest_sd_exp) & !is.na(pretest_sd_control),
      sqrt(
        ((n_exp_post - 1) * pretest_sd_exp^2 +
           (n_control_post - 1) * pretest_sd_control^2) /
          (n_exp_post + n_control_post - 2)
      ),
      NA_real_
    ),
    s_within_post = ifelse(
      notes_in_notes == "Include" &
        !is.na(post_sd_exp) & !is.na(post_sd_control),
      sqrt(
        ((n_exp_post - 1) * post_sd_exp^2 +
           (n_control_post - 1) * post_sd_control^2) /
          (n_exp_post + n_control_post - 2)
      ),
      NA_real_
    )
  )
###########

# STEP 4: Compute manual pre–post d (priority 1, guarded)
#########
mydata <- mydata %>%
  mutate(
    d_manual_pre_post = ifelse(
      notes_in_notes == "Include" &                    # guard first
        is.na(final_d_effect) &
        !is.na(pretest_mean_exp) & !is.na(pretest_mean_control) &
        !is.na(post_mean_exp) & !is.na(post_mean_control) &
        !is.na(sd_pooled_pre) & sd_pooled_pre != 0,
      ((post_mean_exp - pretest_mean_exp) -            # safe math: only runs if Include
         (post_mean_control - pretest_mean_control)) / sd_pooled_pre,
      NA_real_
    )
  )
##########

## STEP 5: Compute manual post-only d (priority 2)
#######
mydata <- mydata %>%
  mutate(
    d_manual_post = ifelse(
      notes_in_notes == "Include" &
        is.na(final_d_effect) & is.na(d_manual_pre_post) &  # ⬅ strict guard
        is.na(pretest_mean_exp) & is.na(pretest_mean_control) &
        !is.na(post_mean_exp) & !is.na(post_mean_control) &
        !is.na(s_within_post) & s_within_post != 0,
      (post_mean_exp - post_mean_control) / s_within_post,
      NA_real_
    )
  )
########

# STEP 6: Compute t-based d (priority 3)
##########
mydata <- mydata %>%
  mutate(
    d_from_t = ifelse(
      notes_in_notes == "Include" &
        is.na(final_d_effect) & is.na(d_manual_pre_post) & is.na(d_manual_post) &  # ⬅ strict guard
        !is.na(p_postexp_con) & !is.na(n_exp_post) & !is.na(n_control_post),
      qt(1 - as.numeric(p_postexp_con) / 2, n_control_post + n_exp_post - 2) *
        sqrt((1 / n_exp_post) + (1 / n_control_post)),
      NA_real_
    )
  )
###########

# STEP 6b: Compute d from known t-value (directly, not from p)
#######
mydata <- mydata %>%
  mutate(
    d_from_t_direct = ifelse(
      notes_in_notes == "Include" &
        is.na(final_d_effect) &
        is.na(d_manual_pre_post) & is.na(d_manual_post) &
        !is.na(t_postexp_con) &
        !is.na(n_exp_post) & !is.na(n_control_post),
      t_postexp_con * sqrt((1 / n_exp_post) + (1 / n_control_post)),
      NA_real_
    )
  )
###########

# STEP 6c: Compute d from known F-value (postexp_con)
#######
mydata <- mydata %>%
  mutate(
    d_from_f = ifelse(
      notes_in_notes == "Include" &
        is.na(final_d_effect) &
        is.na(d_manual_pre_post) & is.na(d_manual_post) &
        is.na(d_from_t_direct) &
        !is.na(f_postexp_con) &
        !is.na(n_exp_post) & !is.na(n_control_post),
      sqrt(f_postexp_con * ((n_exp_post + n_control_post) /
                              (n_exp_post * n_control_post))),
      NA_real_
    )
  )
#########

# STEP 6d: Compute d from F (time × group interaction)
##########
mydata <- mydata %>%
  mutate(
    d_from_f_time = ifelse(
      notes_in_notes == "Include" &
        is.na(final_d_effect) &
        is.na(d_manual_pre_post) & is.na(d_manual_post) &
        is.na(d_from_t_direct) & is.na(d_from_f) &
        !is.na(ftimexgroup) &
        !is.na(n_exp_post) & !is.na(n_control_post),
      sqrt(ftimexgroup * ((n_exp_post + n_control_post) /
                            (n_exp_post * n_control_post))),
      NA_real_
    )
  )
##########

## STEP 6e: Calculate t-value from p-value (as last fallback before coalesce)
#########
mydata <- mydata %>%
  mutate(
    p_postexp_con = ifelse(
      notes_in_notes == "Include" & is.na(final_d_effect),
      as.numeric(p_postexp_con),
      p_postexp_con
    ),
    p_value = coalesce(as.numeric(p_postexp_con), as.numeric(ptimexgroup)),
    df = n_control_post + n_exp_post - 2,
    t_value_from_p = ifelse(
      notes_in_notes == "Include" &
        is.na(final_d_effect) &
        is.na(d_manual_pre_post) & is.na(d_manual_post) &
        is.na(d_from_t_direct) & is.na(d_from_f) & is.na(d_from_f_time) &
        !is.na(p_value),
      qt(1 - p_value / 2, df),
      NA_real_
    ),
    d_from_t_p = ifelse(
      notes_in_notes == "Include" &
        is.na(final_d_effect) &
        !is.na(t_value_from_p) &
        !is.na(n_exp_post) & !is.na(n_control_post),
      t_value_from_p * sqrt((1 / n_exp_post) + (1 / n_control_post)),
      NA_real_
    )
  )
##########

# STEP 7: Combine into final_d_effect (complete with all sources)
#########
mydata <- mydata %>%
  mutate(
    final_d_effect = coalesce(
      final_d_effect,     # existing or converted
      d_manual_pre_post,  # pre–post (best)
      d_manual_post,      # post-only
      d_from_t,           # t from p_postexp_con (old)
      d_from_t_direct,    # known t-value
      d_from_f,           # known F-value
      d_from_f_time,      # F time×group
      d_from_t_p          # computed t from p-value (fallback)
    )
  )
#######

# STEP 8: Adjust direction if "low is good"
#######
mydata <- mydata %>%
  mutate(
    final_d_effect = if_else(
      !is.na(notes) & str_detect(notes, "low is good") &
        !is.na(final_d_effect) & notes_in_notes == "Include",
      abs(final_d_effect),
      final_d_effect
    )
  )
########

# STEP 9: Track source of final_d_effect (complete)
##########
mydata <- mydata %>%
  mutate(
    final_d_effect_source = case_when(
      !is.na(d_manual_pre_post) ~ "manual_pre_post",
      !is.na(d_manual_post) ~ "manual_post",
      !is.na(d_from_t_direct) ~ "from_t_value_direct",
      !is.na(d_from_f) ~ "from_f_value",
      !is.na(d_from_f_time) ~ "from_f_time_group",
      !is.na(d_from_t_p) ~ "from_t_value_from_p",
      !is.na(d_from_t) ~ "from_t_value_from_p_postexp_con",
      !is.na(final_d_effect) ~ "existing_or_converted",
      TRUE ~ NA_character_
    )
  )
##########

# CALCULATING VARIANS
#######
mydata <- mydata %>%
  mutate(
    variance_final_d_effect = ifelse(notes_in_notes == "Include",
                                     ((n_exp_post + n_control_post) / (n_exp_post * n_control_post)) +
                                       (final_d_effect^2 / (2 * (n_exp_post + n_control_post))),
                                     NA_real_
    )
  )
#########

# ADDING ID TO OUTCOME AND EFFECT SIZE
###########
# create outcome_id
mydata <- mydata %>%
  group_by(study_id) %>%
  mutate(
    outcome_num = as.integer(factor(x_1, levels = unique(x_1))),
    outcome_id = ifelse(is.na(outcome_num), NA, paste0("o", outcome_num))
  ) %>%
  ungroup() %>%
  select(-outcome_num)



# create effect_id only if final_d_effect is not NA
mydata <- mydata %>% 
  mutate(
    effect_id = ifelse(
      !is.na(final_d_effect),
      paste(study_id, x, outcome_id, sep = "."),
      NA_character_
    )
  )
#########

# COUNTRY RENAMING
#######
mydata <- mydata %>%
  # Fill author_year missing values with last observation carried forward
  mutate(author_year = zoo::na.locf(author_year, na.rm = FALSE)) %>%
  
  # Extract Year from author_year and propagate non-NA year per study_id
  mutate(year = as.integer(str_extract(author_year, "\\d{4}"))) %>%
  group_by(study_id) %>%
  mutate(year = first(year[!is.na(year)])) %>%
  ungroup() %>%
  
  # Recode country names in country1 to correct names
  mutate(country1 = recode(country1,
                           "US" = "United States",
                           "USA" = "United States",
                           "UK" = "United Kingdom",
                           "England" = "United Kingdom",
                           "Kosovo" = "Serbia",
                           "German" = "Germany",
                           "NIgeria" = "Nigeria"
  )) %>%
  
  
  # colomn continent
  mutate(continent = countrycode(sourcevar = country1,
                                 origin = "country.name",
                                 destination = "continent"))
#########

# MATCHING COUNTRY NAMES
#########
country_name <- read_xlsx("D:/PhD Nurul/script/country_name.xlsx")
# make both lowercase (or uppercase)
m_c1  <- tolower(mydata$country1)
m_cn  <- tolower(country_name$country)

# how many of the unique country1 are present in country_name$country?
sum(unique(m_c1) %in% unique(m_cn))
countries_in    <- unique(m_c1)[unique(m_c1) %in% unique(m_cn)]
countries_notin <- unique(m_c1)[!unique(m_c1) %in% unique(m_cn)]
length(countries_in)      # count present
length(countries_notin)   # count absent
countries_notin           # inspect mismatches (spelling issues etc.)

# inspect country names from UN
grep("united", country_name$country, ignore.case = TRUE, value = TRUE)
# Fix for united kingdom and states
# working copy to edit
mydata$country1_fix <- mydata$country1
# recode to match your country_name$country values
mydata$country1_fix[grepl("^united states$", mydata$country1_fix,
                          ignore.case = TRUE)] <-
  "United States of America"
mydata$country1_fix[grepl("^united kingdom$", mydata$country1_fix,
                          ignore.case = TRUE)] <-
  "United Kingdom of Great Britain and Northern Ireland"
mydata$country1_fix[grepl("^tanzania$", mydata$country1_fix,
                          ignore.case = TRUE)] <-
  "United Republic of Tanzania"
# check remaining non-matches
m_c1  <- tolower(mydata$country1_fix)
m_cn  <- tolower(country_name$country)
countries_notin <- unique(m_c1)[!unique(m_c1) %in% unique(m_cn)]
countries_notin
#check for russia
grep("russia", country_name$country, ignore.case = TRUE, value = TRUE)
grep("federation", country_name$country, ignore.case = TRUE, value = TRUE)
# continue using the working copy
mydata$country1_fix <- mydata$country1_fix
mydata$country1_fix[grepl("^russia$", mydata$country1_fix,
                          ignore.case = TRUE)] <- "Russian Federation"
# now recompute remaining non-matches
m_c1  <- tolower(mydata$country1_fix)
m_cn  <- tolower(country_name$country)
countries_notin <- unique(m_c1)[!unique(m_c1) %in% unique(m_cn)]
countries_notin
# Check for taiwan and hongkong
grep("taiwan", country_name$country, ignore.case = TRUE, value = TRUE)
grep("hong",   country_name$country, ignore.case = TRUE, value = TRUE)
grep("korea",  country_name$country, ignore.case = TRUE, value = TRUE)
# continue from your working copy
mydata$country1_fix <- mydata$country1_fix
# Hong Kong -> ISO-style name
mydata$country1_fix[grepl("^hong kong$", mydata$country1_fix,
                          ignore.case = TRUE)] <-
  "China, Hong Kong Special Administrative Region"
# South Korea -> Republic of Korea
mydata$country1_fix[grepl("^south korea$", mydata$country1_fix,
                          ignore.case = TRUE)] <-
  "Republic of Korea"
# recompute remaining non-matches
m_c1  <- tolower(mydata$country1_fix)
m_cn  <- tolower(country_name$country)
countries_notin <- unique(m_c1)[!unique(m_c1) %in% unique(m_cn)]
countries_notin
# CHeck for iran, palestine and cyprus
grep("iran",      country_name$country, ignore.case = TRUE, value = TRUE)
grep("palest",    country_name$country, ignore.case = TRUE, value = TRUE)
grep("cyprus",    country_name$country, ignore.case = TRUE, value = TRUE)
# continue from working copy
mydata$country1_fix <- mydata$country1_fix
# Iran
mydata$country1_fix[grepl("^iran$", mydata$country1_fix,
                          ignore.case = TRUE)] <-
  "Iran (Islamic Republic of)"
# Palestine
mydata$country1_fix[grepl("^palestine$", mydata$country1_fix,
                          ignore.case = TRUE)] <-
  "State of Palestine"
# Northern Cyprus -> map to Cyprus
mydata$country1_fix[grepl("^northern cyprus$", mydata$country1_fix,
                          ignore.case = TRUE)] <-
  "Cyprus"
# recompute remaining non-matches
m_c1  <- tolower(mydata$country1_fix)
m_cn  <- tolower(country_name$country)
countries_notin <- unique(m_c1)[!unique(m_c1) %in% unique(m_cn)]
countries_notin
grep("china",      country_name$country, ignore.case = TRUE, value = TRUE)
grep("china",      mydata$country1_fix, ignore.case = TRUE, value = TRUE)
# Change Taiwan with China
mydata$country1_fix <- mydata$country1_fix
# merge Taiwan into China
mydata$country1_fix[grepl("^taiwan$", mydata$country1_fix,
                          ignore.case = TRUE)] <- "China"
# final non-match check
m_c1  <- tolower(mydata$country1_fix)
m_cn  <- tolower(country_name$country)
countries_notin <- unique(m_c1)[!unique(m_c1) %in% unique(m_cn)]
countries_notin

# add country status
mydata_with_status <- mydata %>%
  # join on standardized names
  left_join(country_name,
            by = c("country1_fix" = "country")) %>%
  # rename status column
  rename(country_status = status)
#########

# Save dataset
write.csv(mydata,
          "D:/PhD Nurul/Script submission/mydata_processed_updated.csv",
          row.names = FALSE)
write_xlsx(
  mydata,
  "D:/PhD Nurul/Script submission/mydata_processed_updated.xlsx"
)

# Clean data for running the model
#########
mydata_cog <- mydata %>%
  filter(
    !is.na(final_d_effect),
    !is.na(variance_final_d_effect),
    is.finite(final_d_effect),
    is.finite(variance_final_d_effect),
    variance_final_d_effect > 0
  )
########

# Save filtered data
write.csv(mydata_cog,
          "D:/PhD Nurul/Script submission/mydata_cog_updated.csv",
          row.names = FALSE)
write_xlsx(
  mydata,
  "D:/PhD Nurul/Script submission/mydata_cog_updated.xlsx"
)



# ---- Compare the immune responses between randomization groups using GAM(M) ----

# Strain-specific immune responses are compared using GAM models.
# Overall time-of-day effect across vaccine strains is estimated using a mixed-
#   effects framework (GAMM), taking within-participant correlation into account.



# ---- antibody response ----

# 1. Compare strain-specific antibody responses between the vaccination time groups
# A/H1N1 (adjusted and unadjusted)
GAM_Ah1n1 <- gam( 
  fold_change_Ah1n1 ~ time_group + 
    # s(T0_titer_Ah1n1, k=3) + # baseline titers (continuous)
    # s(age) + 
    # vaccination_history + # number of previous influenza vaccinations (categorical) 
    # sex, 
  data = data, method = "REML")
# Get a different time group as the reference level
# data$time_group <- relevel(data$time_group, ref = "09:00-11:40")
# data$time_group <- relevel(data$time_group, ref = "11:40-14:20")

# get p-values from the model
summary(GAM_Ah1n1)

# Use emmeans to get the means per time group, including 95% CI
Ah1n1_emm <- emmeans(GAM_Ah1n1, specs = "time_group", type = "response")
confint(Ah1n1_emm, method = "Wald") # get 95% CI for the estimated means
test(Ah1n1_emm) # get p-value (H0: mean = 0)
confint(pairs(Ah1n1_emm, adjust="none")) # get 95% confidence interval for the pairwise comparisons
pairs(Ah1n1_emm, adjust="tukey") # Tukey correction for pairwise comparison

# Plot the results
Ah1n1_emm_df <- as.data.frame(Ah1n1_emm)
ggplot(Ah1n1_emm_df, aes(x = time_group, y = emmean)) +
  geom_point(size = 4, color = "#B71C1C") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 0.5, color = "#B71C1C") +
  geom_line(group = 1, color = "#B71C1C", linewidth = 0.75) +
  labs(x = "Vaccination time group", y = "Fold change in antibody titers", title = "Influenza A/H1N1") +
  theme_bw(base_size = 14) 

# Using the fold change (FC) in antibody titers or the post vaccination titers as the 
#   outcome results in the same output (for the effect of vaccination timing)
# Therefore, we pick the FC as the outcome, as it is the primary outcome of our trial



### Repeat this for the other 2 strains: A/H3N2 and B/victoria ###
  




# 2. Test for an overall time-of-day effect across vaccine strains using a GAMM
# Put the data in long format
data_long <- data %>%
  select(participant_id, time_group, vac_time, time_MSFsc, 
         age, vaccination_history, 
         sex, fold_change_Ah1n1, fold_change_Ah3n2, fold_change_Bvict, 
         T0_titer_Ah1n1, T0_titer_Ah3n2, T0_titer_Bvict) %>%
  pivot_longer(
    cols = c(fold_change_Ah1n1, fold_change_Ah3n2, fold_change_Bvict),
    names_to = "strain",
    values_to = "FC_titer"
    ) %>%
  mutate(
    strain = recode(strain,
                    "fold_change_Ah1n1" = "A/H1N1",
                    "fold_change_Ah3n2" = "A/H3N2",
                    "fold_change_Bvict" = "B/Vict"))
data_long <- data_long %>%
  mutate(
    T0_gmean = case_when(
      strain == "A/H1N1"  ~ T0_titer_Ah1n1,
      strain == "A/H3N2"   ~ T0_titer_Ah3n2,
      strain == "B/Vict"   ~ T0_titer_Bvict,
      TRUE ~ NA_real_)
    )
data_long$strain <- as.factor(data_long$strain)

# Fit a GAM with mixed-effects structure, as the 3 outcomes are correlated within participants
# A participant-level random effect was included to account for correlation between 
#   strain-specific antibody responses within individuals.

GAMM_interaction <- gam(
  FC_titer ~ # fold change in antibody titers
    time_group * strain +
    # s(T0_gmean) +
    # s(age) +
    # vaccination_history +
    # sex +
    s(participant_id, bs = "re"), # random effect
  data = data_long, method = "REML")

summary(GAMM_interaction) 
# interaction is not significant, suggesting a similar time-of-day effect across strains

GAMM <- gam(
  FC_titer ~
    time_group +
    strain +
    # s(T0_gmean) +
    # s(age) +
    # vaccination_history +
    # sex +
    s(participant_id, bs = "re"),
  data = data_long, method = "REML")
# Get a different time group as the reference level
# data_long$time_group <- relevel(data_long$time_group, ref = "09:00-11:40")
# data_long$time_group <- relevel(data_long$time_group, ref = "11:40-14:20")



# Get p-values form the model
summary(GAMM)
# Use emmeans to get the means per time group, including 95% CI  
GAMM_emm <- emmeans(GAMM, specs = "time_group", type = "response")
confint(GAMM_emm, method = "Wald") 
test(GAMM_emm) 
confint(pairs(GAMM_emm, adjust="none")) 
pairs(GAMM_emm, adjust="tukey") 


# check normality of residuals
hist(resid(GAMM), main="Residuals histogram")
qqnorm(resid(GAMM)); qqline(resid(GAMM)) 
plot(fitted(GAMM), resid(GAMM), main="Residuals vs fitted")







# ---- T-cell response ----

# Note: this is similar to the method used for the antibody response

# 1. Compare strain-specific T-cell responses between the vaccination time groups
# A/H1N1 (adjusted and unadjusted)
GAM_Ah1n1_Tcell <- gam( 
  Tcell_response_Ah1n1 ~ time_group + 
    # s(age) + 
    # vaccination_history + # number of previous influenza vaccinations (categorical) 
    # sex, 
    data = data, method = "REML")
# Get a different time group as the reference level
# data$time_group <- relevel(data$time_group, ref = "09:00-11:40")
# data$time_group <- relevel(data$time_group, ref = "11:40-14:20")

# get p-values from the model
summary(GAM_Ah1n1_Tcell)

# Use emmeans to get the means per time group, including 95% CI
Ah1n1_emm_Tcell <- emmeans(GAM_Ah1n1_Tcell, specs = "time_group", type = "response")
confint(Ah1n1_emm_Tcell, method = "Wald") 
test(Ah1n1_emm_Tcell) 
confint(pairs(Ah1n1_emm_Tcell, adjust="none"))
pairs(Ah1n1_emm_Tcell, adjust="tukey") 


### Repeat this for the other 2 strains: A/H3N2 and B/victoria ###






# 2. Test for an overall time-of-day effect across vaccine strains using a GAMM
# Put the data in long format
data_long_T <- data %>%
  select(participant_id, time_group, vac_time, time_MSFsc, 
         age, vaccination_history, sex, 
         Tcell_response_Ah1n1, Tcell_response_Ah3n2, Tcell_response_Bvict) %>%
  pivot_longer(
    cols = c(Tcell_response_Ah1n1, Tcell_response_Ah3n2, Tcell_response_Bvict),
    names_to = "strain",
    values_to = "Tcell_response"
  ) %>%
  mutate(
    strain = recode(strain,
                    "Tcell_response_Ah1n1" = "A/H1N1",
                    "Tcell_response_Ah3n2" = "A/H3N2",
                    "Tcell_response_Bvict" = "B/Vict"))
data_long_T$strain <- as.factor(data_long_T$strain)

# Fit a GAM with mixed-effects structure, as the 3 FC are correlated within participants
# A participant-level random effect was included to account for correlation between 
#   strain-specific antibody responses within individuals.

Tcell_GAMM_interaction <- gam(
  Tcell_response ~ 
    time_group * strain +
    # s(age) +
    # vaccination_history +
    # sex +
    s(participant_id, bs = "re"), # random effect
  data = data_long_T, method = "REML")

summary(Tcell_GAMM_interaction) 
# interaction is not significant, suggesting a similar time-of-day effect across strains

Tcell_GAMM <- gam(
  Tcell_response ~
    time_group +
    strain +
    # s(age) +
    # vaccination_history +
    # sex +
    s(participant_id, bs = "re"),
  data = data_long_T, method = "REML")
# Get a different time group as the reference level
# data_long_T$time_group <- relevel(data_long_T$time_group, ref = "09:00-11:40")
# data_long_T$time_group <- relevel(data_long_T$time_group, ref = "11:40-14:20")

# Get p-values form the model
summary(Tcell_GAMM)
# Use emmeans to get the means per time group, including 95% CI  
Tcell_GAMM_emm <- emmeans(Tcell_GAMM, specs = "time_group", type = "response")
confint(Tcell_GAMM_emm, method = "Wald") 
test(Tcell_GAMM_emm) 
confint(pairs(Tcell_GAMM_emm, adjust="none")) 
pairs(Tcell_GAMM_emm, adjust="tukey") 


# check normality of residuals
hist(resid(Tcell_GAMM), main="Residuals histogram")
qqnorm(resid(Tcell_GAMM)); qqline(resid(Tcell_GAMM)) 
plot(fitted(Tcell_GAMM), resid(Tcell_GAMM), main="Residuals vs fitted")



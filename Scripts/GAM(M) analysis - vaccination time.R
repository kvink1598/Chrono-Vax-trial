# ---- Compare immune responses from 09:00 to 17:00 using GAM(M) ----

# Overall time-of-day effect across vaccine strains is estimated using a mixed-
#   effects framework (GAMM), taking within-participant correlation into account.

# Instead of time_group (= randomization group), we will use vac_time (continuous
#   vaccination time), which is defined as the number of minutes after 09:00.



# ---- antibody response ----

# Fit a GAM with mixed-effects structure, as the 3 outcomes are correlated within participants
# A participant-level random effect was included to account for correlation between 
#   strain-specific antibody responses within individuals.

overall_interaction <- gam(
  FC_titer ~ # fold change in antibody titers
    s(vac_time, by = strain, bs = "ps") + # penalized spline
    strain +
    # s(T0_gmean) +
    # s(age) +
    # vaccination_history +
    # sex +
    s(participant_id, bs = "re"), # random effect
  data = data_long, method = "REML")

summary(overall_interaction) 
# interaction is not significant, suggesting a similar time-of-day effect across strains

overall_GAMM <- gam(
  FC_titer ~
    s(vac_time, bs = "ps") +
    strain +
    # s(T0_gmean) +
    # s(age) +
    # vaccination_history +
    # sex +
    s(participant_id, bs = "re"),
  data = data_long, method = "REML")

# Check significance of the spline for vac_time 
summary(overall_GAMM) 

# Get the estimated mariginal means using emmeans
overall_GAMM_emm <- emmeans(
  object = overall_GAMM,
  specs = ~ vac_time,
  # For the entire time range
  at = list(vac_time = seq(from = 9, to = 17, by = 0.2))) # 09:00-17:00

# Use as_tibble for printing
overall_GAMM_emm_df <- as_tibble(overall_GAMM_emm)

# Visualise the difference in antibody responses over time.
ggplot(
  data = overall_GAMM_emm_df,
  mapping = aes(x = vac_time, y = emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_ribbon(fill = "#2C2C2C", alpha = 0.2) +
  geom_line(color = "#2C2C2C", linewidth = 1) +
  labs(x="Time of vaccine administration", y="Antibody response") +
  theme_bw() +
  scale_y_continuous(
    breaks = seq(0.5, 2, by = 0.5), 
    limits = c(0.4, 2.1)
  ) +
  scale_x_continuous(
    limits = c(9, 17),
    breaks = c(9, 11, 13, 15, 17),
    labels = c("09:00", "11:00", "13:00", "15:00", "17:00"))


# check normality of residuals
hist(resid(overall_GAMM), main="Residuals histogram")
qqnorm(resid(overall_GAMM)); qqline(resid(overall_GAMM)) 
plot(fitted(overall_GAMM), resid(overall_GAMM), main="Residuals vs fitted")





# ---- T-cell response ----

overall_interaction_Tcell <- gam(
  Tcell_response ~ 
    s(vac_time, by = strain, bs = "ps") +
    strain +
    # s(age) +
    # vaccination_history +
    # sex +
    s(participant_id, bs = "re"), 
  data = data_long_T, method = "REML")

summary(overall_interaction_Tcell) 
# interaction is not significant, suggesting a similar time-of-day effect across strains

overall_GAMM_Tcell <- gam(
  Tcell_response ~
    s(vac_time, bs = "ps") +
    strain +
    # s(age) +
    # vaccination_history +
    # sex +
    s(participant_id, bs = "re"),
  data = data_long_T, method = "REML")

# Check significance of the spline for vac_time 
summary(overall_GAMM_Tcell) 

# Get the estimated mariginal means using emmeans
overall_GAMM_emm_Tcell <- emmeans(
  object = overall_GAMM_Tcell,
  specs = ~ vac_time,
  # For the entire time range
  at = list(vac_time = seq(from = 9, to = 17, by = 0.2))) # 09:00-17:00

# Use as_tibble for printing
overall_GAMM_emm_Tcell_df <- as_tibble(overall_GAMM_emm_Tcell)

# Visualise the difference in T-cell responses over time.
ggplot(
  data = overall_GAMM_emm_Tcell_df,
  mapping = aes(x = vac_time, y = emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_ribbon(fill = "#2C2C2C", alpha = 0.2) +
  geom_line(color = "#2C2C2C", linewidth = 1) +
  labs(x="Time of vaccine administration", y="T-cell response") +
  theme_bw() +
  scale_y_continuous(
    breaks = seq(0.5, 2, by = 0.5), 
    limits = c(0.4, 2.1)
  ) +
  scale_x_continuous(
    limits = c(9, 17),
    breaks = c(9, 11, 13, 15, 17),
    labels = c("09:00", "11:00", "13:00", "15:00", "17:00"))






# ---- Compare immune responses over time of vaccination relative to MSFsc ----

# Chronotype can be quantified by the sleep-corrected midpoint of sleep (MSFsc), 
#  using the micro-MCTQ questionnaire.  
# Vaccination time relative to the MSFsc is another way of treating vaccination 
#  timing. This maps everyone on the same scale, relative to their own MSFsc and, 
#  therefore, takes variation in chronotype between participants into account.

# The MSFsc can be calculated using the following formula:
# MSFsc = MSF - (SDf - SDw)/2
# with SDf = sleepduration on free days (no alarm in the morning);
#  SDw = sleep duration on workdays (or days you have to put an alarm in the morning)
#  MSF = sleep onset on free days + (SDf/2)

# Time of vaccination since MSFsc (time_MSFsc) is the amount of hours that vaccination
#  occured after MSFsc


# These GAMM analyses are similar to the analyses in which we treated time as 
#  a continuous variable.



# ---- antibody response ----

chrono_interaction <- gam(
  FC_titer ~ # fold change in antibody titers
    s(time_MSFsc, by = strain, bs = "ps") + # penalized spline
    strain +
  # s(T0_gmean) +
  # s(age) +
  # vaccination_history +
  # sex +
  s(participant_id, bs = "re"), # random effect
  data = data_long, method = "REML")

summary(chrono_interaction) 
# interaction is not significant, suggesting a similar timing effect across strains

chrono_GAMM <- gam(
  FC_titer ~
    s(time_MSFsc, bs = "ps") +
    strain +
    # s(T0_gmean) +
    # s(age) +
    # vaccination_history +
    # sex +
    s(participant_id, bs = "re"),
  data = data_long, method = "REML")

# Check significance of the spline for time_MSFsc 
summary(chrono_GAMM) 

# Get the estimated mariginal means using emmeans
max_time_MSFsc <- max(data_long$time_MSFsc)
min_time_MSFsc <- min(data_long$time_MSFsc)
chrono_GAMM_emm <- emmeans(
  object = chrono_GAMM,
  specs = ~ time_MSFsc,
  # For the entire time range
  at = list(time_MSFsc = seq(from = min_time_MSFsc, to = max_time_MSFsc, by = 0.5)))

# Use as_tibble for printing
chrono_GAMM_emm_df <- as_tibble(chrono_GAMM_emm)

# Visualise the difference in antibody responses over time relative to MSFsc.
ggplot(
  data = chrono_GAMM_emm_df,
  mapping = aes(x = time_MSFsc, y = emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_ribbon(fill = "#2C2C2C", alpha = 0.2) +
  geom_line(color = "#2C2C2C", linewidth = 1) +
  labs(x="Time of vaccination relative to MSFsc", y="Antibody response") +
  theme_bw() 


# check normality of residuals
hist(resid(chrono_GAMM), main="Residuals histogram")
qqnorm(resid(chrono_GAMM)); qqline(resid(chrono_GAMM)) 
plot(fitted(chrono_GAMM), resid(chrono_GAMM), main="Residuals vs fitted")




# ---- T-cell response ----

chrono_interaction_Tcell <- gam(
  Tcell_response ~ 
    s(time_MSFsc, by = strain, bs = "ps") + 
    strain +
    # s(age) +
    # vaccination_history +
    # sex +
    s(participant_id, bs = "re"), # random effect
  data = data_long_T, method = "REML")

summary(chrono_interaction_Tcell) 
# interaction is not significant, suggesting a similar timing effect across strains

chrono_GAMM_Tcell <- gam(
  Tcell_response ~
    s(time_MSFsc, bs = "ps") +
    strain +
    # s(age) +
    # vaccination_history +
    # sex +
    s(participant_id, bs = "re"),
  data = data_long_T, method = "REML")

# Check significance of the spline for time_MSFsc 
summary(chrono_GAMM_Tcell) 

# Get the estimated mariginal means using emmeans
max_time_MSFsc_T <- max(data_long_T$time_MSFsc)
min_time_MSFsc_T <- min(data_long_T$time_MSFsc)
chrono_GAMM_emm_Tcell <- emmeans(
  object = chrono_GAMM_Tcell,
  specs = ~ time_MSFsc,
  # For the entire time range
  at = list(time_MSFsc = seq(from = min_time_MSFsc_T, to = max_time_MSFsc_T, by = 0.5)))

# Use as_tibble for printing
chrono_GAMM_emm_Tcell_df <- as_tibble(chrono_GAMM_emm_Tcell)

# Visualise the difference in T-cell responses over time relative to MSFsc.
ggplot(
  data = chrono_GAMM_emm_Tcell_df,
  mapping = aes(x = time_MSFsc, y = emmean, ymin = lower.CL, ymax = upper.CL)) +
  geom_ribbon(fill = "#2C2C2C", alpha = 0.2) +
  geom_line(color = "#2C2C2C", linewidth = 1) +
  labs(x="Time of vaccination relative to MSFsc", y="T-cell response") +
  theme_bw() 
  







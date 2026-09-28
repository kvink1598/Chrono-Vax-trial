# Chrono-Vax

The Chrono-Vax trial is a randomized controlled trial investigating the effect
of inlfuenza vaccination timing (time of day) on both antibody and T-cell 
responses in older adults (aged 60-58 years).

Participants received influenza vaccination at a randomized time:
- 09:00-11:40
- 11:40-14:20
- 14:20-17:00

The statistical analysis consists of two parts:
1. the primary analysis: permutation-based test
2. exploratory analysis: GAM(M) analysis

The primary analysis was conducted in accordance with the trial protocol.
However, upon inspection of the data, significant variation in baseline antibody
titers was observed. Those vaccinated in the late afternoon had higher baseline
titers, which affects the fold change. As the permutation test is unable to
adjust for baseline titers, a generalized additive (mixed-effects) model (GAM(M))
was used to get unbiased estimates of the effect of vaccination timing on
vaccine-induced immune responses.


# Scripts:
There are .. scripts in the 'Scripts' folder:
 - Load packages
    This script shows which packages need to be installed and loaded for the analysis  
 - Permutation test 
    The primary analysis
 - GAM(M) analysis - vaccination time group
    GAM(M) analyses comparing the ranomdization groups
 - GAM(M) analysis - vaccination time 
    GAM(M) analyses treating time as a continuous variable
 - GAM(M) analysis - vaccination time relative to MSFsc
   
 
 
 
 
 
 
 
 
 
 
 
 
# Chrono-Vax

The Chrono-Vax trial is a randomized controlled trial investigating the effect
of inlfuenza vaccination timing (time of day) on both antibody and T-cell 
responses in older adults (aged 60-58 years).

Participants received influenza vaccination at a randomized time:
- 09:00-11:40
- 11:40-14:20
- 14:20-17:00

# Statistical analysis

Primary analysis:

- Comparing randomization groups using a permutation-based test (primary analysis)
- Comparing vaccine-induced immune responses from 09:00 to 17:00 using GAM(M)s

Secondary analysis:

- Comparing randomization groups using generalized additive (mixed-effect) models
    (GAM(M)) 
    
Exploratory analysis:

- Assess the relationship between vaccine-induced immune responses and 
    an individual's internal time (circadian phase).


The primary analysis was conducted in accordance with the trial protocol.
However, upon inspection of the data, significant variation in baseline antibody
titers was observed. Those vaccinated in the late afternoon had higher baseline
titers, which affects the fold change. As the permutation test is unable to
adjust for baseline titers, a generalized additive (mixed-effects) model (GAM(M))
was used to reduce bias in the comparison of immune responses between the 
randomization groups.


# Scripts:
There are 5 scripts in the 'Scripts' folder:
 - Load packages
    This script shows which packages need to be installed and loaded for the analyses  
 - Permutation test 
    The primary analysis
 - GAM(M) analysis - group comparison
    GAM(M) analyses comparing the randomization groups
 - GAM(M) analysis - vaccination time 
    GAM(M) analyses treating time as a continuous variable
 - GAMM analysis - vaccination time relative to MSFsc
    GAMM analyses comparing immune responses over vaccination time relative
    to the sleep-corrected midpoint of sleep (MSFsc), which is another way of
    of treating vaccination timing. This maps everyone on the same scale, relative
    to their own MSFsc and therefore takes variation in chronotype between 
    participants into account.
    

   
# Sensitivity analyses 
 
The scripts can be reused to perform the sensitivity analyses in which different
values for titers below the lower limit of detection (<10) were used. 
 

 
 
 
 
 
 
 
 
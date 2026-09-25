
### Permutation test to compare immune responses between vaccination time groups ###

# Data required:  time of vaccination
#                 Antibody response (fold change in anitbody titers)
#                 T-cell response


## --- Permutation-based omnibus test 

# --- Antibody response

# Test for significant differences in fold change between the randomization groups
# Randomization groups = time_group

# A/H1N1
set.seed(123)
prm_Ah1n1 <- independence_test(
  fold_change_Ah1n1 ~ time_group, 
  data=data, 
  teststat="quadratic", # testat=quadratic for multiple group comparisons 
  distribution=approximate(nresample = 10000)) #  10,000 random permutations to 
                                # approximate the null distribution and p-value
prm_Ah1n1

# A/H3N2
set.seed(123)
prm_Ah3n2 <- independence_test(
  fold_change_Ah3n2 ~ time_group, 
  data=data, 
  teststat="quadratic",  
  distribution=approximate(nresample = 10000))
prm_Ah3n2

# B/Victoria
set.seed(123)
prm_Bvict <- independence_test(
  fold_change_Bvict ~ time_group, 
  data=data, 
  teststat="quadratic",  
  distribution=approximate(nresample = 10000))
prm_Ah3n2



# Optional:
# permutation test of B_fold vs time_group, with permutations restricted within 
# age_group × sex strata (controls for age and sex)
independence_test(
  fold_change_Ah1n1 ~ time_group | interaction(age_group, sex), 
  data=data, 
  teststat="quadratic", 
  distribution=approximate(nresample = 10000))


# Optional:
# Compare antibody responses between vaccinations before and after 1 pm 
# (= 1pm_split)
independence_test(
  fold_change_Ah1n1 ~ time_group | interaction(age_group, sex), 
  data=data, 
  teststat="scalar", # scalar: used for comparing two groups 
  distribution=approximate(nresample = 10000))


# Optional:
# The permutation test tests for differences among the 3 vaccination time groups
# Perform pairwise group comparisons to see which groups differ
# Get all pairwise combinations of the time groups
levels <- levels(d$time_group)
combinations <- combn(levels, 2, simplify = FALSE)
# Perform tests for each pair
pvals <- sapply(combinations, function(pair) {
  subset_data <- subset(d, time_group %in% pair)
  
  test_Ah1n1 <- independence_test(fold_change_Ah1n1 ~ time_group, 
                                  data = subset_data, 
                                  teststat="scalar", 
                                  distribution = approximate(nresample = 10000))  
  pvalue(test_Ah1n1)
})
# Adjust p-values for multiple comparisons using Holm's method
adjusted_pvals <- p.adjust(pvals, method = "holm")
# Combine results
results <- data.frame(
  Comparison = sapply(combinations, paste, collapse = " vs "),
  Raw_p_value = pvals,
  Adj_p_value = adjusted_pvals
)
print(results)





# --- T-cell response

# Test for significant differences in change in T-cell response between the 
# randomization groups (= time_group)

# A/H1N1
set.seed(123)
independence_test(
  Tcell_response_Ah1n1 ~ time_group, 
  data=data, 
  teststat="quadratic", 
  distribution=approximate(nresample = 10000)) 

# A/H3N2
set.seed(123)
independence_test(
  Tcell_response_Ah3n2 ~ time_group, 
  data=data, 
  teststat="quadratic", 
  distribution=approximate(nresample = 10000)) 

# B/Victoria
set.seed(123)
independence_test(
  Tcell_response_Bvict ~ time_group, 
  data=data, 
  teststat="quadratic", 
  distribution=approximate(nresample = 10000)) 






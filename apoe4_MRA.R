library(metafor)

MRA <- function(beta_full, se_full, beta_sub, se_sub){
  wi = 1 / (se_sub^2)  
  SE2 = se_full^2      
  beta_red = beta_full + (SE2 * (beta_full * wi - beta_sub * wi) / (1 - SE2 * wi))
  SE_red = se_full / (sqrt(1 - SE2 * wi))
  return(list(beta_red = beta_red, SE_red = SE_red))
}

studies <- data.frame(
  study = c("aa", "eas", "hi", "nhw"),
  beta = c(0.978326, 1.597365, 0.815365, 1.163151),
  se = c(0.05071, 0.041235, 0.046304, 0.019909)
)

full_meta_results <- rma(yi = studies$beta, sei = studies$se, method = "FE")
beta_full <- coef(full_meta_results)
se_full <- sqrt(vcov(full_meta_results)[1, 1])

mra_results <- matrix(ncol = 5, nrow = nrow(studies))


for (i in 1:nrow(studies)) {
  current_study <- studies[i, ]
  other_studies <- studies[-i, ]
  
  reduced_meta_results <- rma(yi = other_studies$beta, sei = other_studies$se, method = "FE")
  beta_reduced_full <- coef(reduced_meta_results)
  se_reduced_full <- sqrt(vcov(reduced_meta_results)[1, 1])
  
  mra_result <- MRA(beta_full, se_full, current_study$beta, current_study$se)
  mra_results[i,1] = current_study$study
  mra_results[i,2] = mra_result$beta_red
  mra_results[i,3] = mra_result$SE_red
  mra_results[i,4] = beta_reduced_full
  mra_results[i,5] = se_reduced_full
  
}
colnames(mra_results) = c("study removed","MRA_beta","MRA_se","Traditional_beta", "Traditional_se")
mra_results = data.frame(mra_results)

print(mra_results)

############################################################
library(metafor)

MRA <- function(beta_full, se_full, beta_sub, se_sub){
  wi = 1 / (se_sub^2)  
  SE2 = se_full^2     
  beta_red = beta_full - (SE2 * (beta_full * wi - beta_sub * wi) / (1 + SE2 * wi))
  SE_red = se_full / (sqrt(1 + SE2 * wi))
  return(list(beta_red = beta_red, SE_red = SE_red))
}


studies <- data.frame(
  study = c("aa", "eas", "hi", "nhw"),
  beta = c(0.978326, 1.597365, 0.815365, 1.163151),
  se = c(0.05071, 0.041235, 0.046304, 0.019909)
)


full_meta_results <- rma(yi = studies$beta, sei = studies$se, method = "FE")
print("Full meta-analysis results (including all studies):")
print(full_meta_results)

studies_excl_aa <- studies[-1,]
excl_aa_meta_results <- rma(yi = studies_excl_aa$beta, sei = studies_excl_aa$se, method = "FE")
print("Meta-analysis results excluding 'aa':")
print(excl_aa_meta_results)

aa_study <- studies[1,]  
mra_results <- MRA(coef(full_meta_results), sqrt(vcov(full_meta_results)[1, 1]),
                   aa_study$beta, aa_study$se)
print("MRA adjusted results (excluding 'aa' using MRA function):")
print(mra_results)


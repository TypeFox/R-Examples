getVar <-
function(covs, tIndex){

  ## create treatment vector
  treat = numeric(nrow(covs))
  treat[tIndex] = 1
  
  ## model propensity score, calculate variance
  psModel = glm(treat~., family=binomial('logit'), data=covs)
  psVar = var(psModel$fitted)
  return(psVar)
}

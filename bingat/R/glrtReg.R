glrtReg <-
function(data, type, groups){
regCoeffHA <- glmReg(data, type, groups)
loglikHA <- regCoeffHA$loglik

regCoeffH0 <- glmReg(data, type, rep(0, ncol(data)))
loglikH0 <- regCoeffH0$loglik

glrt <- 2 * (loglikHA-loglikH0)
return(glrt)
}

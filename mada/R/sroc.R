sroc <- function(fit, ...) UseMethod("sroc")

### calculate naive sroc curves
calc.sroc <- function(fpr, alpha.sens, alpha.fpr, mu1, mu2, sigma2, sigma){
  theta <- sigma/sigma2
  return(inv.trafo(alpha.sens, (mu1 - theta*mu2) + theta*trafo(alpha.fpr,fpr)))
}

sroc2 <- function(fit, fpr = 1:99/100){
  if(!attr(fit$logLik,"df") == 5){stop("AUC can not be calculated for meta-regression")}
  estimate <- fit$coefficients
  alpha.sens <- fit$alphasens
  alpha.fpr <- fit$alphafpr
  mu1 <- estimate[1]
  mu2 <- estimate[2]
  sigma2 <- fit$Psi[2,2]
  sigma <- fit$Psi[1,2]  
  return(cbind(fpr, sens = calc.sroc(fpr, alpha.sens, alpha.fpr, mu1, mu2, sigma2, sigma)))
}

mcsroc <- function(fit, ...) UseMethod("mcsroc")

mcsroc<- function(fit, replications = 10000, lambda = 100){
  if(!attr(fit$logLik,"df") == 5){stop("SROC can not be calculated for meta-regression")}
  estimate <- fit$coefficients
  alpha.sens <- fit$alphasens
  alpha.fpr <- fit$alphafpr
  mu <- estimate[1:2]
  Sigma <- fit$Psi
  stud.pars <- rmvnorm(replications, mu, Sigma)
  sens <- inv.trafo(alpha.sens, stud.pars[,1])
  fpr <- inv.trafo(alpha.fpr, stud.pars[,2])
  N.sens <- rpois(replications, lambda)
  N.fpr <- rpois(replications, lambda)
  TN <- rbinom(replications, N.sens, sens)
  FP <- rbinom(replications, N.fpr, fpr)
  obs.sens <- TN/N.sens
  obs.fpr <- FP/N.fpr
  return(lowess(obs.fpr, obs.sens))
}

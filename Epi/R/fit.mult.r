fit.mult <-
function(y, rates.frame, cov.frame, start)
{
  if (missing(start)) {
    ## Fit model without covariates to get initial rates estimates
    glm.out.rates <- fit.baseline(y, rates.frame)
    ## Initial values for iterative fitting
    lambda <- coef(glm.out.rates)
    beta <- rep(0, ncol(cov.frame))
  }
  else {
    lambda <- start[1:ncol(rates.frame)]
    beta <- start[ncol(rates.frame) + 1:ncol(cov.frame)]
  }

  niter <- 1
  cy <- 1 - y
  while(TRUE) {
    ## covariates model
    off <- log(-as.matrix(rates.frame) %*% lambda)
    glm.out.cov <- glm(cy ~ -1 + offset(off) + .,
                       family=binomial(link=cloglog),
                       data=cov.frame, start=beta, maxit=100)
    beta <- coef(glm.out.cov)
    
    ## rates model
    wgt <- exp(as.matrix(cov.frame) %*% beta)
    temp.rates.frame <- wgt * rates.frame
    glm.out.rates <- glm(y ~ -1 + ., family=binomial(link=log),
                         data=temp.rates.frame, start=lambda, maxit=100)
    lambda <- coef(glm.out.rates)

    ## Check convergence 
    ## Convergence <==> deviances are equal
    TOL <- max(glm.out.cov$control$epsilon, glm.out.rates$control$epsilon)
    dev1 <- glm.out.cov$deviance
    dev2 <- glm.out.rates$deviance
    if (abs(dev1 - dev2)/(0.1 + abs(dev1)) < TOL) 
      break
    else
      niter <- niter + 1
  }

  return(list( rates=glm.out.rates, cov=glm.out.cov, niter=niter))
}

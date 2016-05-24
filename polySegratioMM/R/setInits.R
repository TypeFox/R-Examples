`setInits` <-
function(model, priors, seed=1)
{
  ## at moment only one chain - think about it
  ## from JAGS 0.97 seed is set in inits as vector with one seed for
  ##                each chain
  ## previously JAGS 0.90 was latest version for Windows so
  ## didn't set seed here but in "jags.cmd" file
  ## outdated now that JAGS 1.0 available
  
  if (class(model) != "modelSegratioMM")
    stop("'model' must be of class 'modelSegratioMM'")

  if (class(priors) != "priorsSegratioMM")
    stop("'priors' must be of class 'priorsSegratioMM'")
  
  n <- model$n.components

  mu <- c(0, rep(NA,n-1))
  P <- c(0.7, 0.6*c((n-1):1)/(n-1)/n)  # nb: denom=sum arithmetic prog.
  if (model$equal.variances) {
    tau <- priors$params$logit.prec[1]
  } else {
    tau <- priors$params$logit.prec[1:n]
  }
  theta <- diff(gtools::logit(model$E.segRatio$ratio[1:n]))

  res <- list(mu=mu, P=P, tau=tau, theta=theta)

  ## superseded now that JAGS Version 1.0 required
  ##if (.Platform$OS.type != "windows"){# to be fixed after JAGS 0.90 superseded
  ##  res$seed <- seed
  ##}
  
  if (model$random.effect) {
    res$taub <-  priors$params$logit.prec[1]
    names(res$taub) <- NULL
  }
  return(res)
}


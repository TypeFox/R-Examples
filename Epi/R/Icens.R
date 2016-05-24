Icens <- function(first.well, last.well, first.ill, formula,
                  model.type=c("MRR","AER"), breaks,
                  boot=FALSE, alpha=0.05, keep.sample=FALSE, data)
{
  ## Create follow-up matrix containing three event times
  fu.expression <- substitute(cbind(first.well, last.well, first.ill))
  fu <- if (missing(data)) {
    eval(fu.expression)
  }
  else {
    eval(fu.expression, data)
  }
  ## Check consistency of arguments
  missing.f1 <- is.na(fu[,1])
  missing.f2 <- is.na(fu[,2])
  if (any(missing.f1 & missing.f2)) {
    stop("You must supply at least one of \"first.well\" and \"last.well\"")
  }
  if (any(fu[,1] > fu[,2], na.rm=TRUE) | any(fu[,2] > fu[,3], na.rm=TRUE)) {
    stop("Some units do not meet: first.well < last.well < first.ill" )
  }
  ## Fill in any gaps
  fu[,1][missing.f1] <- fu[,2][missing.f1]
  fu[,2][missing.f2] <- fu[,1][missing.f2]
  ## Recensor cases that fall after the last break point
  is.censored <- fu[,3] > max(breaks)
  is.censored[is.na(is.censored)] <- FALSE
  fu[is.censored,3] <- NA
  
  exp.dat <- expand.data(fu, formula, breaks, data)

  model.type <- match.arg(model.type)
  if (missing(formula)) {
    fit.out <- with(exp.dat, fit.baseline(y, rates.frame))
    lambda <- coef(fit.out)
  }
  else {
    fit.out <- switch(model.type,
                      "MRR"=with(exp.dat, fit.mult(y, rates.frame, cov.frame)),
                      "AER"=with(exp.dat, fit.add(y, rates.frame, cov.frame)))
    lambda <- coef(fit.out$rates)
  }
  
  beta <- if (is.null(fit.out$cov)) numeric(0) else coef(fit.out$cov)
  params <- c(lambda,beta)
  if (boot) {
    nboot <- ifelse (is.numeric(boot), boot, 100)
    boot.coef <- matrix(NA, nrow=nboot, ncol=length(lambda) + length(beta))
    colnames(boot.coef) <- names(params)
    
    for (i in 1:nboot) {
      subsample <- sample(nrow(fu), replace=TRUE)
      exp.dat <- expand.data(fu[subsample,], formula, breaks, data[subsample,])
      if (missing(formula)) {
        sim.out <- with(exp.dat, fit.baseline(y, rates.frame, params))
        boot.coef[i,] <- coef(sim.out)
      }
      else {
        sim.out <- switch(model.type,
                          "MRR"=with(exp.dat,
                            fit.mult(y, rates.frame, cov.frame, params)),
                          "AER"=with(exp.dat,
                            fit.add(y, rates.frame, cov.frame, params)))
        boot.coef[i,] <- switch(model.type,
                                "MRR"=c(coef(sim.out[[1]]),
                                  coef(sim.out[[2]])),
                                "AER"=coef(sim.out[[1]]))
      }
    }
    ci.quantiles=c(0.5, alpha/2, 1 - alpha/2)
    boot.ci <- t(apply(boot.coef,2,quantile,ci.quantiles))

    lower.ci.lab <-
      paste("lower ", formatC(100 * alpha/2, format="f", digits=1),"%", sep="")
    upper.ci.lab <- 
      paste("upper ", formatC(100 * (1-alpha/2), format="f", digits=1),"%",
            sep="")
    
    colnames(boot.ci) <- c("median", lower.ci.lab, upper.ci.lab)
    fit.out$boot.ci <- boot.ci
    if (keep.sample) {
      fit.out$sample <- boot.coef
    }
  }
  
  class( fit.out ) <- "Icens"
  attr( fit.out, "model" ) <- model.type
  return( fit.out )
}

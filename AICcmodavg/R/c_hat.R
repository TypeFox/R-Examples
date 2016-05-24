##create generic c_hat
c_hat <- function(mod, method = "pearson", ...) {
  ##format list according to model class
  UseMethod("c_hat", mod)
}

##default to indicate when object class not supported
c_hat.default <- function(mod, method = "pearson", ...) {
  stop("\nFunction not yet defined for this object class\n")
}



##function to compute c-hat from Poisson or binomial GLM with success/total syntax
c_hat.glm <- function(mod, method = "pearson", ...){

  ##determine family of model
  fam <- family(mod)$family

  ##if binomial, check if n > 1 for each case
  if(fam == "binomial") {
    ##extract number of trials
    n.trials <- mod$prior.weights
    if(identical(unique(n.trials), 1)) {
      stop("\nWith a binomial distribution, the number of successes must be summarized for valid computation of c-hat\n")
    }
  }

  ##Poisson or binomial
  if(!any(fam == c("poisson", "binomial"))) {
    stop("\nEstimation of c-hat only valid for Poisson or binomial GLM's\n")
  }
  
  ##Pearson chi-square
  chisq <- sum(residuals(mod, type = "pearson")^2)

  ##return estimate based on Pearson chi-square
  if(method == "pearson") {
    c_hat.est <- chisq/mod$df.residual
    attr(c_hat.est, "method") <- "pearson estimator"
  }

  
  ##return estimate based on deviance estimator
  if(method == "deviance") {
    ##estimate deviance
    mod.deviance <- sum(residuals(mod, type = "deviance")^2)
    
    c_hat.est <- mod.deviance/mod$df.residual
    attr(c_hat.est, "method") <- "deviance estimator"
  }

  ##extract raw residuals
  raw.res <- residuals(mod, type = "response")
  ##extract fitted values
  fit.vals <- fitted(mod)
      
  ##estimate s.bar for Poisson
  if(fam == "poisson") {
    si <- 1/fit.vals * raw.res
    s.bar <- mean(si)
  }

  ##estimate s.bar for binomial
  if(fam == "binomial") {
    si <- (1 - 2 * fit.vals)/((n.trials * fit.vals) * (1 - fit.vals))
    s.bar <- mean(si)
  }

  
  ##return estimate based on Farrington estimator
  if(method == "farrington") {
    c_hat.est <- (chisq - sum(si))/mod$df.residual
    attr(c_hat.est, "method") <- "farrington estimator"
  }
  
  
  ##return estimate based on Fletcher estimator
  if(method == "fletcher") {
    c_hat.est <- (chisq/mod$df.residual)/(1 + s.bar)
    attr(c_hat.est, "method") <- "fletcher estimator"
  }

  class(c_hat.est) <- "c_hat"
  return(c_hat.est)
}



##function to compute c-hat from Poisson or binomial GLM with success/total syntax
c_hat.vglm <- function(mod, method = "pearson", ...){

  ##determine family of model
  fam <- mod@family@vfamily
  if(length(fam) > 1) fam <- fam[1]

  ##if binomial, check if n > 1 for each case
  if(fam == "binomialff") {
    ##extract number of trials
    n.trials <- mod@prior.weights
    if(identical(nrow(mod@prior.weights), 0)) {
      stop("\nWith a binomial distribution, the number of successes must be summarized for valid computation of c-hat\n")
    }
  }


  ##Poisson or binomial
  if(!any(fam == c("poissonff", "binomialff"))) {
    stop("\nEstimation of c-hat only valid for Poisson or binomial GLM's\n")
  }
  
  ##Pearson chi-square
  chisq <- sum(residuals(mod, type = "pearson")^2)

  ##return estimate based on Pearson chi-square
  if(method == "pearson") {
    c_hat.est <- chisq/mod@df.residual
    attr(c_hat.est, "method") <- "pearson estimator"
  }

  
  ##return estimate based on deviance estimator
  if(method == "deviance") {
    ##estimate deviance
    mod.deviance <- sum(residuals(mod, type = "deviance")^2)
    
    c_hat.est <- mod.deviance/mod@df.residual
    attr(c_hat.est, "method") <- "deviance estimator"
  }

  ##extract raw residuals
  raw.res <- residuals(mod, type = "response")
  ##extract fitted values
  fit.vals <- fitted(mod)
      
  ##estimate s.bar for Poisson
  if(fam == "poisson") {
    si <- 1/fit.vals * raw.res
    s.bar <- mean(si)
  }

  ##estimate s.bar for binomial
  if(fam == "binomialff") {
    si <- (1 - 2 * fit.vals)/((n.trials * fit.vals) * (1 - fit.vals))
    s.bar <- mean(si)
  }

  
  ##return estimate based on Farrington estimator
  if(method == "farrington") {
    c_hat.est <- (chisq - sum(si))/mod@df.residual
    attr(c_hat.est, "method") <- "farrington estimator"
  }
  
  
  ##return estimate based on Fletcher estimator
  if(method == "fletcher") {
    c_hat.est <- (chisq/mod@df.residual)/(1 + s.bar)
    attr(c_hat.est, "method") <- "fletcher estimator"
  }

  class(c_hat.est) <- "c_hat"
  return(c_hat.est)
}



##method for GLMM from lme4
c_hat.merMod <- function(mod, method = "pearson", ...) {

  #determine family of model
  fam <- family(mod)$family

  ##if binomial, check if n > 1 for each case
  if(fam == "binomial") {
    if(identical(unique(mod@resp$n), 1)) {
      stop("\nWith a binomial distribution, the number of successes must be summarized for valid computation of c-hat\n")
    }
  }
      
  ##Poisson or binomial
  if(!any(fam == c("poisson", "binomial"))) {
    stop("\nEstimation of c-hat only valid for Poisson or binomial GLMM's\n")
  }
    
  ##number of parameters estimated
  n.parms <- attr(logLik(mod), "df")
  
  ##total number of observations
  n.obs <- nrow(model.frame(mod))

  ##residual df
  res.df <- n.obs - n.parms

  if(method == "pearson") {
    chisq <- sum(residuals(mod, type = "pearson")^2)
    c_hat.est <- chisq/res.df
    attr(c_hat.est, "method") <- "pearson estimator"
    
  } else {stop("\nOnly Pearson estimator is currently supported for GLMM's\n")}

  class(c_hat.est) <- "c_hat"
  return(c_hat.est)
}


##print method
print.c_hat <- function(x, digits = 2, ...) {
  cat("'c-hat' ", paste(format(c(x), digits = digits), collapse = ", "),
      " (method: ", format(attr(x, "method")), ")\n", sep = "")
}



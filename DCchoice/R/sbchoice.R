# a binary choice model for Single-bounded data. a simple logit or probit model
sbchoice <- function(formula, data, subset, dist = "log-logistic", ...){
  if (!inherits(formula, "Formula"))
  formula <- Formula(formula)

  # evaluating the formula and stops if the formula is not defined correctly
  if (!inherits(formula, "formula")) stop("invalid formula")

  cl <- match.call()            
#  if(missing(data)) data <- environment(formula)
  if(missing(data)) stop("the name of the data frame object must be supplied in the 'data' argument")
  
#  data <- eval(data, parent.frame())
  mf <- match.call(expand.dots = TRUE)
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  original.data <- data
  data <- mf
  mm.data <- model.matrix(formula, data = data, rhs = 1:2)

  # removing observations with missing values
  na.num <- max(sum(as.numeric(is.na(data))))
  if(na.num != 0){ 
    d1 <- nrow(data)
    data <- na.omit(data)
    d2 <- nrow(data)
    warning(paste("Missing values detected.", d1 - d2, "rows are removed.", sep = " "))
  }

  # defining the dependent variable 
  y1 <- model.part(formula, data = data, lhs = 1)[[1]]  # yes/no to the bids

  nobs <- length(y1)
  
  BID <- model.frame(formula, data = data, lhs = 0, rhs = 2)[[1]]  # the suggested bid
  
  # handling the data matrix
  f2 <- formula(formula, lhs = 0, rhs = 1)
  X <- model.frame(f2, data)
  mmX <- model.matrix(f2, X)

  form <- formula(terms(formula))

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ##
if(dist == "logistic" | dist == "log-logistic"){  # logistic or log-logistic error distribution
    glm.out <- glm(form, family = binomial(link = "logit"), data = data)  # unrestricted model
    glm.null <- update(glm.out, .~ 1)  # restricted (only the intercept) model
    
    npar <- length(glm.out$coefficients) # the number of parameters
    if (substr(dist, 1, 4) == "log-") {  # changing the name label if the model is log-logistic
        names(glm.out$coefficients)[npar] <- "log(bid)"
    } else {
        names(glm.out$coefficients)[npar] <- "BID"
    }
    names(glm.out$coefficients)[1] <- "(Intercept)"
    estimates <- glm.out$coefficients  # taking out the estimates
    
} else if(dist == "normal" | dist == "log-normal") {  # normal or log-normal error distribution
    glm.out <- glm(form, family = binomial(link = "probit"), data = data)
    glm.null <- update(glm.out, .~ 1)
    
    npar <- length(glm.out$coefficients)
    if (substr(dist, 1, 4) == "log-") {
        names(glm.out$coefficients)[npar] <- "log(bid)"
    } else {
        names(glm.out$coefficients)[npar] <- "BID"
    }
    names(glm.out$coefficients)[1] <- "(Intercept)"
    estimates <- glm.out$coefficients
} else if(dist == "weibull"){
    # likelihood function
    sbLL <- function(param, dvar, ivar){
          y1 <- dvar
          X <- ivar
          ll <- 
            sum(pweibull(exp(-X[y1==1, , drop=FALSE]%*%param), shape = 1, lower.tail = FALSE, log.p=TRUE)) + 
            sum(pweibull(exp(-X[y1==0, , drop=FALSE]%*%param), shape = 1, lower.tail = TRUE, log.p=TRUE))
          ifelse(is.finite(ll), return(-ll), NaN) 
      }
    # initial parameter values
    ini <- glm(form, family = binomial(link = "probit"), data = data)
    ini.par <- ini$coefficients
    ini.par.null <- update(ini, . ~ 1)$coefficients

    
    # ML estimation
    suppressWarnings( # "glm." is nothing to do with GLM. The naming is merely because of compatibility
#        glm.out <- optim(ini.par, fn = sbLL, method="BFGS", hessian = TRUE, dvar = y1, ivar = cbind(X, BID), control=list(abstol=10^(-30)))
        glm.out <- optim(ini.par, fn = sbLL, method="BFGS", hessian = TRUE, dvar = y1, ivar = mm.data)
    )
    suppressWarnings(
#        glm.null <- optim(ini.par.null, fn = sbLL, method = "BFGS", hessian = TRUE, dvar = y1, ivar = matrix(1, nobs, 1), control=list(abstol=10^(-30)))
        glm.null <- optim(ini.par.null, fn = sbLL, method = "BFGS", hessian = TRUE, dvar = y1, ivar = matrix(1, nobs, 1))
    )
    names(glm.out$par)[1] <- "(Intercept)"
    if(dist == "weibull") {
      npar <- length(glm.out$par)     # the number of parameters
      names(glm.out$par)[npar] <- "log(bid)"
      colnames(glm.out$hessian)[npar] <- "log(bid)"
      rownames(glm.out$hessian)[npar] <- "log(bid)"
    }
    estimates <- glm.out$par
    # compatibility with other distributions
    glm.out$converged <- ifelse(glm.out$convergence == 0, TRUE, FALSE)  # convergence status of the ML estimation
    glm.out$iter <- glm.out$counts
    glm.out$coefficients <- glm.out$par
    glm.null$coefficients <- glm.null$par
} else {
    stop("dist must be logistic, normal or weibull")
}

  terms <- terms(formula)
  fac <- which(attr(attr(X, "terms"), "dataClasses") == "factor")
  xlevels <- as.list(fac)
  j <- 0
  for (i in fac) {
    j <- j + 1
    xlevels[[j]] <- levels(X[[i]])
  }
  contrasts <- attr(mmX, "contrasts")


  # arranging outcomes into a single list variable
   output <- list(
      coefficients = estimates, # the coefficient estimates
      call = cl,            # the function call
      formula = formula,    # the defined model formula
      glm.out = glm.out,    # the outcome of the unrestricted model
      glm.null = glm.null,  # the outcome of the null model
      distribution = dist,  # the specified error distribution
      covariates = mmX,       # a matrix of the covariates
      bid = BID,            # the suggested bid
      nobs = nobs,          # the number of observations
      yn = y1,              # the acceptance/rejection variable
      data.name = data,     # the data matrix
      terms = terms,
      contrasts = contrasts,
      xlevels = xlevels
      )

  class(output) <- "sbchoice"   # setting the object class
  return(output)
}

# summarizing the output
summary.sbchoice <- function(object, ...){
  dist <- object$distribution
  coef <- object$coefficients
  npar <- length(coef)
  X <- object$covariates
  bid <- object$bid

    # function for obrtaining AIC and BIC
    akaike <- function(loglik, npar, k ){
      -2*loglik + k*npar
    }

  if(dist == "weibull"){
    # creating a table for coefficients, se, etc. 
    se <- sqrt(diag(solve(object$glm.out$hessian)))
    zstat <- coef/se  # z statistics
    pval <- round(2*pnorm(-abs(zstat)), 6)  # p-value
    coefmat <- cbind(coef, se, zstat, pval)
    colnames(coefmat) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    object$coefmat <- coefmat
    
    object$loglik <- -object$glm.out$value
    object$loglik.null <- -object$glm.null$value
  } else {
    # obtaining necessary components from the object
    object$glm.summary <- summary.glm(object$glm.out)
    object$glm.null.summary <- summary.glm(object$glm.null)
    object$loglik <- logLik(object$glm.out)[1]
    object$loglik.null <- logLik(object$glm.null)[1]
  } 
    # computing various mean estimates for different error distributions
    WTPs <- wtp(object = X, b = coef, bid = bid, dist = dist)
    object$medianWTP <- WTPs$medianWTP
    object$meanWTP <- WTPs$meanWTP
    object$trunc.meanWTP <- WTPs$trunc.meanWTP
    object$adj.trunc.meanWTP <- WTPs$adj.trunc.meanWTP

    # computing pseudo-R^2
    object$psdR2 <- 1 - object$loglik/object$loglik.null
    names(object$psdR2) <- "pseudo-R^2 measure"
    object$adjpsdR2 <- 1 - (object$loglik - npar)/object$loglik.null
    names(object$adjpsdR2) <- "adjusted pseudo-R^2 measure"

    # computing Likelihood Ratio Statistic and its p-value
    LR <- 2*(object$loglik - object$loglik.null)    # the LR statistic
    d.f <- length(object$glm.out$coefficients) - length(object$glm.null$coefficients) # the degrees of freedom
    pvalLR <- pchisq(LR, df = d.f, lower.tail = FALSE)  # p-value
    object$LR.test <- c(LR, d.f, pvalLR)

    # computing AIC and BIC by the function AKAIKE
    object$AIC <- akaike(object$loglik, npar, k = c(2, log(object$nobs)))
    names(object$AIC) <- c("AIC", "BIC")  

    class(object) <- "summary.sbchoice"
    return(object)
}


print.sbchoice <- function(x, digits = max(3, getOption("digits") - 1), ...){
  cat("\nDistribution:", x$distribution, "\n", sep = " ")
#  print(x$glm.out$coefficients)
  print(x$coefficients)
  invisible(x)
}


print.summary.sbchoice <- function(x, digits = max(3, getOption("digits") - 1), ...){
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  cat("Formula:", deparse(x$formula, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  
  if(!x$glm.out$converged)  cat("The optimization did not converge\n")
  cat("Coefficients:", "\n")
  if(x$distribution == "weibull" | x$distribution == "log-weibull") {
    printCoefmat(x$coefmat, digits = 4, dig.tst = 3)
  } else {
    printCoefmat(x$glm.summary$coefficients, digits = 4, dig.tst = 3)    # printCoefmat() is defined in stats package
  }
  cat("\nDistribution:", x$distribution, "", sep = " ")
  cat("\nNumber of Obs.:", formatC(x$nobs, digits = 0), " ")
  cat("\nlog-likelihood:", x$loglik, "\n", sep = " ")
  cat("pseudo-R^2:", formatC(x$psdR2, format="f", digits = 4), 
      ", adjusted pseudo-R^2:", formatC(x$adjpsdR2, format="f", digits = 4), "")
  cat("\nLR statistic:", round(x$LR.test[1], 3), "on", formatC(x$LR.test[2], digits = 0), 
    "DF, p-value:", formatC(x$LR.test[3], format="f", digits = 3), "\n")
  cat("AIC:", formatC(x$AIC[1], format="f", digits = digits), ", BIC:", formatC(x$AIC[2], format="f", digits = digits), "\n")
  cat("\nIterations:", x$glm.out$iter, " ")
  cat("\nConvergence:", x$glm.out$converged, "\n")
  
  cat("\nWTP estimates:")
  if(is.finite(x$meanWTP)){
    cat("\n Mean :", x$meanWTP, "", sep = " ")
  } else {
    cat("\n Mean :", x$meanWTP, "(because of |beta_bid| < 1)", sep = " ")
  }
  cat("\n Mean :", x$trunc.meanWTP, "(truncated at the maximum bid)", sep = " ")
  cat("\n Mean :", x$adj.trunc.meanWTP, "(truncated at the maximum bid with adjustment)", sep = " ")
  cat("\n Median :", x$medianWTP, "\n", sep = " ")

}



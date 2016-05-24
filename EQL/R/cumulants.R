cumulants <- function(saddlef, cgf=NULL, kappa2f=NULL, rho3f=NULL, rho4f=NULL,
                      cgf.deriv=NULL, domain=interval(-Inf, Inf), ...) {
  if (!is.function(saddlef)) {
    stop("'saddlef' is not a function!")
  }

  if (class(domain) != "interval") {
    warning("'domain' must be an object of class 'interval'!\n",
            "'domain' will bet set to (-Inf, Inf)!")
    domain <- interval(-Inf, Inf)
  }

  mu.inv <- function(x) {
    if (any(!liesWithin(x, domain))) {
      stop("Variable 'x' outside the domain!")
    }
    saddlef(x, ...)
  }
  
  isMissing <- FALSE
  
  # Derive cumulant functions from cgf.deriv
  if (!is.null(cgf.deriv)) {
    if (!is.function(cgf.deriv)) {
      stop("The generic derivative of the cgf must be given as a function!")
    }
    if (length(formals(cgf.deriv)) < 2) {
      stop("The generic derivative function must take at least two arguments!")
    }
    if (any(!is.null(c(cgf, kappa2f, rho3f, rho4f)))) {
      warning("The generic derivative of the cgf is given!\n  ",
              "Thus, any explicit given cumulants will be overwritten!")
    }
    type <- "implicit"
    K      <- function(x) cgf.deriv(0, x, ...)
    kappa2 <- function(x) cgf.deriv(2, x, ...)
    rho3   <- function(x) cgf.deriv(3, x, ...) / (kappa2(x) ^ (3 / 2)) 
    rho4   <- function(x) cgf.deriv(4, x, ...) / (kappa2(x) ^ 2)
  } else { # Funtions are given explicitly
    cumulantf <- list(cgf=cgf, kappa2f=kappa2f)
    isFunction <- sapply(cumulantf, is.function)
    if (any(!isFunction)) {
      badFunctions <- paste("'", names(cumulantf)[!isFunction], "'", sep="")
      if (length(badFunctions) == 1) {
        msg <- " is not a function!"
      } else {
        msg <- " are not functions!"
      }
      stop(paste(badFunctions, collapse=", "), msg)
    }
    check <- .missingFormals(cgf, list(...), 1)
    if (!is.null(check)) {
      check <- paste("'", check, "'", sep="")
      if (length(check) == 1) { 
        msg <- paste("Argument", check, "needed in 'cgf' but not supplied!")
      } else {
        msg <- paste("Arguments", check, "needed in 'cgf' but not supplied!",
                     collapse="")
      }
      warning(msg)
    }
    check <- .missingFormals(kappa2f, list(...), 1)
    if (!is.null(check)) {
      check <- paste("'", check, "'", sep="")
      if (length(check) == 1) { 
        msg <- paste("Argument", check, "needed in 'kappa2f' but not supplied!")
      } else {
        msg <- paste("Arguments", check, "needed in 'kappa2f' but not supplied!",
                     collapse="")
      }
      warning(msg)
    }
    type <- "explicit"
    K      <- function(x) cgf(x, ...)
    kappa2 <- function(x) kappa2f(x, ...)
    if (is.null(rho3f)) {
      isMissing <- TRUE
      rho3 <- function(x) rep(0, length(x))
      warning("Function 'rho3f' not supplied!")
    } else {
       check <- .missingFormals(rho3f, list(...), 1)
       if (!is.null(check)) {
         check <- paste("'", check, "'", sep="")
         if (length(check) == 1) { 
           msg <- paste("Argument", check, "needed in 'rho3f' but not supplied!")
         } else {
           msg <- paste("Arguments", check, "needed in 'rho3f' but not supplied!",
                        collapse="")
         }
         warning(msg)
       }
      rho3 <- function(x) rho3f(x, ...)
    }
    if (is.null(rho4f)) {
      isMissing <- TRUE
      rho4 <- function(x) rep(0, length(x))
      warning("Function 'rho4f' not supplied!")
    } else {
       check <- .missingFormals(rho4f, list(...), 1)
       if (!is.null(check)) {
         check <- paste("'", check, "'", sep="")
         if (length(check) == 1) { 
           msg <- paste("Argument", check, "needed in 'rho4f' but not supplied!")
         } else {
           msg <- paste("Arguments", check, "needed in 'rho4f' but not supplied!",
                        collapse="")
         }
         warning(msg)
       }
      rho4 <- function(x) rho4f(x, ...)
    }
  }
  cumulantsList <- list(K=K, mu.inv=mu.inv, kappa2=kappa2, rho3=rho3,
                        rho4=rho4, domain=domain, extra.params=c(...),
                        type=type, missing=isMissing)
  class(cumulantsList) <- "cumulants"
  return(cumulantsList)
}

gammaCumulants <- function(shape, scale) {
  K.deriv <- function(n, x, shape, scale) {
    if (n == 0) {
      return(-shape * log(1 - scale * x))
    } else {
      return(factorial(n - 1) * shape * scale ^ n / (1 - scale * x) ^ n)
    }
  }
  theta.hat <- function(x, shape, scale) {
    return((x - shape * scale) / (x * scale))
  }
  return(cumulants(theta.hat, cgf.deriv=K.deriv, domain=interval(0, Inf),
                   shape=shape, scale=scale))
}

gaussianCumulants <- function(mu, sigma2) {
  K.deriv <- function(n, x, mu, sigma2) {
    if (n <= 2) {
      switch(n + 1,
             return(mu * x + sigma2 * x ^ 2 / 2), # n == 0
             return(mu + sigma2 * x),             # n == 1
             return(rep(sigma2, length(x))))      # n == 2
    } else {
      return(rep(0, length(x)))                   # n >= 3
    }
  }
  theta.hat <- function(x, mu, sigma2) {
    return((x - mu) / sigma2)
  }
  return(cumulants(theta.hat, cgf.deriv=K.deriv, mu=mu, sigma2=sigma2))
}

inverseGaussianCumulants <- function(lambda, nu) {
  K.deriv <- function(n, x, lambda, nu) {
    if (n == 0) {
      return((lambda / nu) * (1 - sqrt(1 - 2 * nu ^ 2 * x / lambda)))
    } else {
      coef <- factorial(2 * n) / (factorial(n) * 2 ^ n)
      return(coef * nu ^ (2 * n - 1) /
             (lambda ^ (n - 1) * (1 - 2 * nu ^ 2 * x / lambda) ^ (n - 1 / 2)))
    }
  }
  theta.hat <- function(x, lambda, nu) {
    return(lambda * (x ^ 2 - nu ^ 2) / (2 * nu ^ 2 * x ^ 2))
  }
  return(cumulants(theta.hat, cgf.deriv=K.deriv, domain=interval(0, Inf),
                   lambda=lambda, nu=nu))
}

check.cumulants <- function(object, ...) {
  n <- 10
  isCumulant <- TRUE
  functNames <- c("K", "mu.inv", "kappa2", "rho3", "rho4")
  extraNames <- c("domain", "extra.params", "type", "missing")
  standardNames <- c(functNames, extraNames)
  
  haveSameElements <- identical(names(object), standardNames)
  if (!haveSameElements) {
    isCumulant <- FALSE
    warning("The fields of the object do not equate the fields of an ",
            "'cumulants' object!", call.=FALSE)
    return(isCumulant)
  }

  x <- runif(n, min=0.1, max=0.9)
  valueList <- try(suppressWarnings(list(object$K(x), object$mu.inv(x), object$kappa2(x),
                    object$rho3(x), object$rho4(x))), silent=TRUE)
  if (class(valueList) == "try-error") {
    isCumulant <- FALSE
    warning("A function call did not work properly:\n\t",
            valueList, "Does the function accept extra arguments?", call.=FALSE)
    return(isCumulant)
  }
  names(valueList) <- functNames
  isVectorized <- sapply(functNames, function(name) length(valueList[[name]]) == n)

  if (any(!isVectorized)) {
    isCumulant <- FALSE
    badNames <- paste("'", functNames[!isVectorized], "'", sep="")
    if (sum(!isVectorized) > 1) {
      msg <- paste("The functions ", badNames, " are not properly verctorized!",
                   collapse=", ")
    } else {
      msg <- paste("The function", badNames, "is not properly vectorized!")
    }
    warning(msg, call.=FALSE)
    return(isCumulant)
  }
  return(isCumulant)
}
    
  


  
                  

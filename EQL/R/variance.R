varianceFamily <- function(varf, devf=NULL, link="log", initf=NULL, 
                           validmuf=NULL, name="default") {
  if (!is.function(varf)) {
    stop("'varf' is not a function!")
  }
  if (is.null(validmuf)) {
    warning("'validmuf' is not given! Should be provided ",
            "to check for valid mus!")
    validmuf <- function(mu, ...) TRUE
  }
  if (is.null(initf)) {
    warning("No initialize sequence is given! Should be provided to set up ", 
            "starting values correctly!")
    initf <- function(...) {
      return(expression({
        n <- rep.int(1, nobs)
        mustart <- y
      }))
    }
  }
  # Getting the family parameters, first parameter of V is mu, so skip it
  params <- formals(varf)
  params[1] <- NULL
  params["..."] <- NULL

  if (is.null(devf)) {
    warning("Deviance function is determined numerically!\n",
            "For more accurate results supply an explicit function!", call.=FALSE)
    devf <- function(y, mu, ...) {
      mapply(function(y, mu, ...) {
        integrand <- function(u, y, ...)
          return((y - u) / varf(u, ...))
        int <- try(suppressWarnings(integrate(integrand, lower=y, upper=mu, y=y, ...)),
                   silent=TRUE)
        if (class(int) == "try-error") {
          stop("Deviance could not be evaluated!\n", int)
        } else {
          return(-2 * int$value)
        }
      }, y, mu, ...)
    }
    type <- "numerical"
  } else if (!is.function(devf)) {
    stop("'devf' is not a function!")
  } else {
    # Function given explicitly - check if parameters are passed properly
    # Remove first argument and "..." to get a list of all 'nonlinear' params
    check <- .missingFormals(devf, params, 2)
    if (!is.null(check)) {
      check <- paste("'", check, "'", sep="")
      if (length(check) == 1) { 
        msg <- paste("Argument", check, "supplied in 'varf' but not in 'devf'!")
      } else {
        msg <- paste("Arguments", check, "supplied in 'varf' but not in 'devf'!",
                     collapse="")
      }
      warning(msg)
    }
    type <- "explicit"
  }

  # Get link, code snippet from stats::glm
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  }
  stats <- try(make.link(linktemp), silent=T)
  if (class(stats) == "try-error") {
    if (is.character(link)) {
      stats <- make.link(link)
      linktemp <- link
    } else {
      if (inherits(link, "link-glm")) {
        stats <- link
        if (!is.null(stats$name)) 
          linktemp <- stats$name
      }
      else {
        stop("Link '", link, "' not recognized!") 
      }
    }
  }
  # end snippet
  family <- function(...) {
    # Check if params specified correctly
    familyParams <- list(...)
    if (length(familyParams) < length(params)) {
      stop("Too few arguments specified!\nNeeded arguments:",
           paste("'", names(params), "'", sep="", collapse=", "), "!")
    } else if (length(familyParams) > length(params)) {
      stop("Too many arguments specified!\nNeeded arguments:",
           paste("'", names(params), "'", sep="", collapse=", "),
           "\nPassed arguments:",
           paste("'", familyParams, "'", sep="", collapse=", "), "!")
    }
    variance <- function(mu) varf(mu, ...)
    validmu <- function(mu) validmuf(mu, ...)
    dev.resids <- function(y, mu, wt) wt * devf(y, mu, ...)
    aic <- function(y, n, mu, wt, dev) NA
    initialize <- initf(...)
    if (is.null(names(familyParams))) {
      names(familyParams) <- names(params)
    }
    return(structure(list(family = name, link = linktemp, 
                          linkfun = stats$linkfun, linkinv = stats$linkinv,
                          variance = variance, dev.resids = dev.resids,
                          aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                          validmu = validmu, valideta = stats$valideta,
                          params=familyParams), 
                     class = c("extFamily", "family")))
  }
  value <- list(family=family, name=name, params=params, type.dev=type)
  class(value) <- "varianceFamily"
  return(value)
}

powerVarianceFamily <- function(link="log") {
  init <- function(theta,...) {
    if (theta == 0) {
      return(expression({
        mustart <- y
        n <- rep.int(1, nobs)}))
    } else if (theta == 1) {
      return(expression({
        if (any(y < 0)) {
          stop("'y' must not be negative for theta=1!")
        }
        mustart <- y + 0.1
        n <- rep.int(1, nobs)}))
    } else {
      return(expression({
        if (any(y <= 0)) {
          stop("'y' must be postive!")
        }
        mustart <- y
        n <- rep.int(1, nobs)}))
    }
  }
  validmuf <- function(mu, theta) {
    if (theta == 0) {
      return(TRUE)
    } else if (theta == 1 | theta == 2  | theta == 3) {
      return(all(mu > 0))
    } else {
      return(all(mu > 0))
    }
  }
  varf <- function(y, theta) y^theta
  devf <- function(y, mu, theta) {
    if(theta == 1) {
      return(2 * (y * log(ifelse(y == 0, 1, y / mu)) - (y - mu)))
    } else if(theta == 2) {
      return(2 * (y / mu - log(ifelse(y == 0, 1, y / mu)) - 1))
    } else {
      return(1/ ((1 - theta) * (2 - theta)) * 2 *
             (y ^ (2 - theta) - (2 - theta) * y * mu ^ (1 - theta) +
            (1 - theta) * mu ^ (2 - theta)))
    }
  }
  return(varianceFamily(varf=varf, devf=devf, link=link, initf=init,
                        validmuf=validmuf, name="Power-Family"))
}

extBinomialVarianceFamily <- function(link="logit") {
  # Exact integral not available
  init <- function(k, l) {
    return(expression({
      mustart <- (weights * y + 0.5)/(weights + 1)
      n <- rep.int(1, nobs)}))
  }
  validmuf <- function(mu, k, l) {
    return(all(mu > 0) && all(mu < 1))
  }
  varf <- function(y, k, l)  y ^ k * (1 - y) ^ l
  return(varianceFamily(varf=varf, link=link, initf=init, 
                        validmuf=validmuf, name="Extended-Binomial-Family"))
}

print.extFamily <- function(x, ...) {
  NextMethod()
  cat("Family-Parameters:\n\t")
  z <- x$param
  paramNames <- names(z)
  cat(paste(lapply(names(z),
                     function(name) paste(name, z[[name]], sep=" = ")) ,
        collapse="\n\t"), "\n")
  invisible(x)
}

dim.varianceFamily <- function(x) {
  return(length(x$params))
}

print.varianceFamily <- function(x, ...) {
  cat("Variance-Family:", x$name, "\n")
  cat("Deviance       :", x$type.dev , "\n")
  cat("Parameter      :", names(x$params), "\n")
  invisible(x)
}


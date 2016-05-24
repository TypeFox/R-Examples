############################################################################################### Definitions

# my_dx = E[Y | D, X]; md_x = E[D | X]; mz_x = E[Z | X]; md_zx = E[D |
# Z, X]; my_zx = E[Y | Z, X]; md = E[D]; mz = E[Z]; myd1_zx = E[YD | Z,
# X]; myd0_zx = E[Y(1-D) | Z, X];


#' Functions for estimation of treatment effects
#'
#' This class of functions estimates the average treatment effect (ATE), the ATE of the tretated (ATET), the local average treatment effects (LATE) and the LATE of
#' the tretated (LATET). The estimation methods rely on immunized / orthogonal moment
#' conditions which guarantee valid post-selection inference in a high-dimensional setting. Further details can be found in Belloni et al. (2014).
#'
#' Details can be found in Belloni et al. (2014).
#'
#' @aliases late latet ate atet LATE LATET ATE ATET
#' @param y outcome variable / dependent variable
#' @param d treatment variable (binary)
#' @param x exogenous variables
#' @param z instrumental variables (binary)
#' @param bootstrap boostrap method which should be employed: 'none', 'Bayes',
#' 'normal', 'wild'
#' @param nRep number of replications for the bootstrap
#' @param ... arguments passed, e.g. \code{intercept} and \code{post}
#' @return Functions return an object of class \code{rlassoTE} with estimated effects, standard errors and
#' individual effects in the form of a \code{list}.
#' @references A. Belloni, V. Chernozhukov, I. Fernandez-Val, and C. Hansen
#' (2014). Program evaluation with high-dimensional data. Working Paper.
#' @keywords local average treatment effect late latt local average structural
#' effect lasf lasft
#' @rdname TE
#' @export

rlassoATE <- function(x, d, y, bootstrap = "none", nRep = 500, ...) {
  z <- d
  res <- rlassoLATE(x, d, y, z, bootstrap = bootstrap, nRep = nRep, ...)
  res$type <- "ATE"
  return(res)
}

#' @export
#' @rdname TE
rlassoATET <- function(x, d, y, bootstrap = "none", nRep = 500, ...) {
  z <- d
  res <- rlassoLATET(x, d, y, z, bootstrap = bootstrap, nRep = nRep, 
                     ...)
  res$type <- "ATET"
  return(res)
}
#' @export
#' @param post logical. If \code{TRUE}, post-lasso estimation is conducted.
#' @param intercept logical. If \code{TRUE}, intercept is included which is not
#' penalized.
#' @rdname TE

rlassoLATE <- function(x, d, y, z, bootstrap = "none", nRep = 500, post = TRUE, 
                       intercept = TRUE) {
  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  checkmate::checkChoice(bootstrap, c("none", "Bayes", "normal", "wild"))
  lambda <- 2.2 * sqrt(n) * qnorm(1 - (1/log(n))/(2 * (2 * p)))
  control <- list(numIter = 15, tol = 10^-5)
  # penalty <- list(method = 'none', lambda.start = rep(lambda, p), c =
  # 1.1, gamma = 0.1)
  penalty <- list(homoscedastic = "none", lambda.start = rep(lambda, 
                                                             p), c = 1.1, gamma = 0.1)
  indz1 <- (z == 1)
  indz0 <- (z == 0)
  # E[Y|Z = 1,X] = my_z1x
  b_y_z1xL <- rlasso(y[indz1] ~ x[indz1, , drop = FALSE], post = post, 
                     intercept = intercept, control = control, penalty = penalty)
  my_z1x <- predict(b_y_z1xL, newdata = x)
  # E[Y|Z = 0,X] = my_z0x
  b_y_z0xL <- rlasso(y[indz0] ~ x[indz0, , drop = FALSE], post = post, 
                     intercept = intercept, control = control, penalty = penalty)
  my_z0x <- predict(b_y_z0xL, newdata = x)
  # E[D|Z = 1,X] = md_z1x
  lambda <- 2.2 * sqrt(n) * qnorm(1 - (1/log(n))/(2 * (2 * p)))
  penalty <- list(lambda.start = lambda, c = 1.1, gamma = 0.1)
  if (sum(d - z) != 0) {
    b_d_z1xL <- rlassologit(d[indz1] ~ x[indz1, , drop = FALSE], post = post, 
                            intercept = intercept, penalty = penalty)
    md_z1x <- predict(b_d_z1xL, newdata = x)
  } else {
    md_z1x <- rep(1, n)
  }
  # E[D|Z = 0,X] = md_z0x
  md_z0x <- rep(0, n)
  
  # E[Z|X] = mz_x
  b_z_xL <- rlassologit(z ~ x, post = post, intercept = intercept)
  mz_x <- predict(b_z_xL, newdata = x)
  mz_x <- mz_x * (mz_x > 1e-12 & mz_x < 1 - 1e-12) + (1 - 1e-12) * (mz_x > 
                                                                      1 - 1e-12) + 1e-12 * (mz_x < 1e-12)
  
  eff <- (z * (y - my_z1x)/mz_x - ((1 - z) * (y - my_z0x)/(1 - mz_x)) + 
            my_z1x - my_z0x)/mean(z * (d - md_z1x)/mz_x - ((1 - z) * (d - md_z0x)/(1 - 
                                                                                     mz_x)) + md_z1x - md_z0x)
  
  se <- sqrt(var(eff))/sqrt(n)
  late <- mean(eff)
  individual <- eff
  
  object <- list(se = se, te = late, individual = individual, type = "LATE", 
                 call = match.call(), samplesize = n)
  
  if (bootstrap != "none") {
    boot <- rep(NA, nRep)
    for (i in 1:nRep) {
      if (bootstrap == "Bayes") {
        weights <- rexp(n, rate = 1) - 1
      }
      if (bootstrap == "normal") {
        weights <- rnorm(n)
      }
      if (bootstrap == "wild") {
        weights <- rnorm(n)/sqrt(2) + (rnorm(n)^2 - 1)/2
      }
      weights <- weights + 1
      boot[i] <- mean(weights * (z * (y - my_z1x)/mz_x - ((1 - z) * 
                                                            (y - my_z0x)/(1 - mz_x)) + my_z1x - my_z0x))/mean(weights * 
                                                                                                                (z * (d - md_z1x)/mz_x - ((1 - z) * (d - md_z0x)/(1 - mz_x)) + 
                                                                                                                   md_z1x - md_z0x))
    }
    object$boot.se <- sqrt(var(boot))
    object$type_boot <- bootstrap
  }
  object$type <- "LATE"
  class(object) <- "rlassoTE"
  return(object)
}

#' @export
#' @rdname TE
rlassoLATET <- function(x, d, y, z, bootstrap = "none", nRep = 500, post = TRUE, 
                        intercept = TRUE) {
  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  checkmate::checkChoice(bootstrap, c("none", "Bayes", "normal", "wild"))
  lambda <- 2.2 * sqrt(n) * qnorm(1 - (1/log(n))/(2 * (2 * p)))
  control <- list(numIter = 15, tol = 10^-5)
  # penalty <- list(method = 'none', lambda.start = rep(lambda, p), c =
  # 1.1, gamma = 0.1)
  penalty <- list(homoscedastic = "none", lambda.start = rep(lambda, 
                                                             p), c = 1.1, gamma = 0.1)
  indz1 <- (z == 1)
  indz0 <- (z == 0)
  # E[Y|Z = 0,X] = my_z0x
  b_y_z0xL <- rlasso(y[indz0] ~ x[indz0, ], post = post, intercept = intercept, 
                     control = control, penalty = penalty)
  my_z0x <- predict(b_y_z0xL, newdata = x)
  # E[D|Z = 0,X] = md_z0x
  md_z0x <- rep(0, n)
  # E[Z|X] = mz_x
  lambdaP <- 2.2 * sqrt(n) * qnorm(1 - (1/log(n))/(2 * p))
  # penalty <- list(lambda.start = lambdaP, c = 1.1, gamma = 0.1)
  penalty <- list(homoscedastic = "none", lambda.start = p, c = 1.1, 
                  gamma = 0.1)
  b_z_xL <- rlassologit(z ~ x, post = post, intercept = intercept, penalty = penalty)
  mz_x <- predict(b_z_xL, newdata = x)
  mz_x <- mz_x * (mz_x > 1e-12 & mz_x < 1 - 1e-12) + (1 - 1e-12) * (mz_x > 
                                                                      1 - 1e-12) + 1e-12 * (mz_x < 1e-12)
  
  
  ## 
  eff <- ((y - my_z0x) - (1 - z) * (y - my_z0x)/(1 - mz_x))/mean((d - 
                                                                    md_z0x) - (1 - z) * (d - md_z0x)/(1 - mz_x))
  
  se <- sqrt(var(eff))/sqrt(n)
  latet <- mean(eff)
  individual <- eff
  
  object <- list(se = se, te = latet, individual = individual, type = "LATET", 
                 call = match.call(), samplesize = n)
  
  if (bootstrap != "none") {
    boot <- rep(NA, nRep)
    for (i in 1:nRep) {
      if (bootstrap == "Bayes") {
        weights <- rexp(n, rate = 1) - 1
      }
      if (bootstrap == "normal") {
        weights <- rnorm(n)
      }
      if (bootstrap == "wild") {
        weights <- rnorm(n)/sqrt(2) + (rnorm(n)^2 - 1)/2
      }
      weights <- weights + 1
      # boot[i] <- mean(weights*(z*(y-my_z1x)/mz_x -
      # ((1-z)*(y-my_z0x)/(1-mz_x)) + my_z1x - my_z0x))/
      # mean(weights*(z*(d-md_z1x)/mz_x - ((1-z)*(d-md_z0x)/(1-mz_x)) +
      # md_z1x - md_z0x))
      boot[i] <- mean(weights * ((y - my_z0x) - (1 - z) * (y - my_z0x)/(1 - 
                                                                          mz_x)))/mean(weights * ((d - md_z0x) - (1 - z) * (d - md_z0x)/(1 - 
                                                                                                                                           mz_x)))
    }
    object$boot.se <- sqrt(var(boot))
    object$type_boot <- bootstrap
  }
  object$type <- "LATET"
  class(object) <- "rlassoTE"
  return(object)
}


################# Methods for rlassoTE

#' Methods for S3 object \code{rlassoTE}
#'
#' Objects of class \code{rlassoTE} are constructed by  \code{rlassoATE},  \code{rlassoATET}, \code{rlassoLATE},  \code{rlassoLATET}.
#' \code{print.rlassoTE} prints and displays some information about fitted \code{rlassoTE} objects.
#' \code{summary.rlassoTE} summarizes information of a fitted \code{rlassoTE} object.
#' \code{confint.rlassoTE} extracts the confidence intervals.
#' @param object an object of class \code{rlassoTE}
#' @param x an object of class \code{rlassoTE}
#' @param digits number of significant digits in printout
#' @param ... arguments passed to the print function and other methods
#' @param parm a specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector of names. If missing, all parameters are considered.
#' @param level confidence level required.
#' @keywords methods rlassoTE
#' @rdname methods.rlassoTE
#' @aliases methods.rlassoTE print.rlassoTE summary.rlassoTE
#' @export

print.rlassoTE <- function(x, digits = max(3L, getOption("digits") - 3L), 
                           ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  if (length(x$te)) {
    cat("Treatment Effect\n")
    cat(paste("Type:", x$type), "\n")
    cat("Value:\n")
    print.default(format(x$te, digits = digits), print.gap = 2L, quote = FALSE)
  } else cat("No treatment effect\n")
  cat("\n")
  invisible(x$te)
}

#' @rdname methods.rlassoTE
#' @export

summary.rlassoTE <- function(object, digits = max(3L, getOption("digits") - 
                                                    3L), ...) {
  if (length(object$te)) {
    table <- matrix(NA, ncol = 4, nrow = 1)
    rownames(table) <- "TE"
    colnames(table) <- c("coeff.", "se.", "t-value", "p-value")
    table[, 1] <- object$te
    if (is.null(object$type_boot)) {
      table[, 2] <- object$se
    } else {
      table[, 2] <- object$boot.se
    }
    table[, 3] <- table[, 1]/table[, 2]
    table[, 4] <- 2 * pnorm(-abs(table[, 3]))
    cat("Estimation and significance testing of the treatment effect", 
        "\n")
    cat(paste("Type:", object$type), "\n")
    cat(paste("Bootstrap:", ifelse(is.null(object$type_boot), "not applicable", 
                                   object$type_boot)), "\n")
    printCoefmat(table, digits = digits, P.values = TRUE, has.Pvalue = TRUE)
    cat("\n")
  } else {
    cat("No coefficients\n")
  }
  cat("\n")
  invisible(table)
}

#' @rdname methods.rlassoTE
#' @export

confint.rlassoTE <- function(object, parm, level = 0.95, ...) {
  n <- object$samplesize
  k <- 1
  cf <- object$te
  pnames <- "TE"
  a <- (1 - level)/2
  a <- c(a, 1 - a)
  fac <- qt(a, n - k)
  pct <- format.perc(a, 3)
  ci <- array(NA, dim = c(length(pnames), 2L), dimnames = list(pnames, 
                                                               pct))
  if (is.null(object$type_boot)) {
    ses <- object$se
  } else {
    ses <- object$boot.se
  }
  ses <- object$se[parm]
  ci[] <- cf[parm] + ses %o% fac
  print(ci)
  invisible(ci)
}

#' Calibration for Linear and Nonlinear Regression Models.
#' 
#' The function \code{invest} computes the inverse estimate and a condfidence 
#' interval for the unknown predictor value that corresponds to an observed 
#' value of the response (or vector thereof) or specified value of the mean 
#' response. See the references listed below for more details. 
#' 
#' @rdname invest
#' 
#' @importFrom nlme fixef
#' @importFrom stats coef family fitted formula na.omit numericDeriv 
#' @importFrom stats model.matrix predict qnorm qt quantile resid residuals 
#' @importFrom stats rnorm runif sd simulate uniroot update var vcov
#' @importFrom utils setTxtProgressBar txtProgressBar
#' 
#' @export
#' 
#' @param object An object that inherits from class \code{"lm"}, \code{"glm"},
#'               \code{"nls"}, or \code{"lme"}.
#' @param y0 The value of the observed response(s) or specified value of the 
#'           mean response. For \code{"glm"} objects, \code{y0} should be on 
#'           scale of the response variable.
#' @param interval The type of interval required.
#' @param level A numeric scalar between 0 and 1 giving the confidence level 
#'              for the interval to be calculated. 
#' @param mean.response Logical indicating whether confidence intervals should
#'                      correspond to an individual response (\code{FALSE}) or 
#'                      a mean response (\code{TRUE}). For \code{glm} objects,
#'                      this is always \code{TRUE}.
#' @param x0.name For multiple linear regression, a character string giving the
#'   the name of the predictor variable of interest.
#' @param newdata For multiple linear regression, a \code{data.frame} giving the
#'   values of interest for all other predictor variables (i.e., those other 
#'   than \code{x0.name}).
#' @param data An optional data frame. This is required if \code{object$data} 
#'             is \code{NULL}.
#' @param nsim Positive integer specifying the number of bootstrap simulations; 
#'             the bootstrap B (or R).
#' @param boot.type Character string specifying the type of bootstrap to use 
#'                  when \code{interval = "percentile"}. Options are 
#'                  \code{"parametric"} and \code{"nonparametric"}.
#' @param seed Optional argument to \code{set.seed}.
#' @param progress Logical indicating whether to display a text-based progress
#'                 bar during the bootstrap simulation.
#' @param lower The lower endpoint of the interval to be searched.
#' @param upper The upper endpoint of the interval to be searched.
#' @param extendInt Character string specifying if the interval 
#'   \code{c(lower, upper)} should be extended or directly produce an error when
#'   the inverse of the prediction function does not have differing signs at the
#'   endpoints. The default, \code{"no"}, keeps the search interval and hence 
#'   produces an error. Can be abbreviated. See the documentation for the 
#'   \code{base} R function \code{uniroot} for details.
#' @param q1 Optional lower cutoff to be used in forming confidence intervals. 
#'           Only used when \code{object} inherits from class \code{"lme"}. 
#'           Defaults to \code{qnorm((1+level)/2)}.
#' @param q2 Optional upper cutoff to be used in forming confidence intervals. 
#'           Only used when \code{object} inherits from class \code{"lme"}. 
#'           Defaults to \code{qnorm((1-level)/2)}.
#' @param tol The desired accuracy passed on to \code{uniroot}. Recommend a 
#'            minimum of 1e-10.
#' @param maxiter The maximum number of iterations passed on to \code{uniroot}. 
#' @param adjust A logical value indicating if an adjustment should be made to
#'               the critical value used in calculating the confidence interval.
#'               This is useful for when the calibration curve is to be used 
#'               multiple, say k, times.
#' @param k The number times the calibration curve is to be used for computing 
#'          a confidence interval. Only needed when 
#'          \code{adjust = "Bonferroni"}.
#' @param ... Additional optional arguments. At present, no optional arguments 
#'            are used.
#' @return \code{invest} returns an object of class \code{"invest"} or, if
#'         \code{interval = "percentile"}, of class 
#'         \code{c("invest", "bootCal")}. The generic function \code{plot} 
#'         can be used to plot the output of the bootstrap simulation when 
#'         \code{interval = "percentile"}.
#'         
#'         An object of class \code{"invest"} contains the following 
#'         components:
#'   \itemize{
#'     \item \code{estimate} The estimate of x0.
#'     \item \code{lwr} The lower confidence limit for x0.
#'     \item \code{upr} The upper confidence limit for x0.
#'     \item \code{se} An estimate of the standard error (Wald and percentile 
#'                     intervals only).
#'     \item \code{bias} The bootstrap estimate of bias (percentile interval 
#'                       only).
#'     \item \code{bootreps} Vector of bootstrap replicates (percentile 
#'                           interval only).
#'     \item \code{nsim} The number of bootstrap replicates (percentile 
#'                       interval only).
#'     \item \code{interval} The method used for calculating \code{lower} and 
#'           \code{upper} (only used by \code{print} method).
#'   }
#'
#' @references
#' Greenwell, B. M., and Schubert Kabban, C. M. (2014). investr: An R Package 
#' for Inverse Estimation. \emph{The R Journal}, \bold{6}(1), 90--100. 
#' URL http://journal.r-project.org/archive/2014-1/greenwell-kabban.pdf.
#'
#' Graybill, F. A., and Iyer, H. K. (1994).
#' \emph{Regression analysis: Concepts and Applications}. Duxbury Press. 
#'
#' Huet, S., Bouvier, A., Poursat, M-A., and Jolivet, E.  (2004)
#' \emph{Statistical Tools for Nonlinear Regression: A Practical Guide with S-PLUS and R Examples}. Springer. 
#' 
#' Norman, D. R., and Smith H. (2014).
#' \emph{Applied Regression Analysis}. John Wiley & Sons.
#' 
#' Oman, Samuel D. (1998).
#' Calibration with Random Slopes.
#' \emph{Biometrics} \bold{85}(2): 439--449.
#' doi:10.1093/biomet/85.2.439.
#' 
#' Seber, G. A. F., and Wild, C. J. (1989)
#' \emph{Nonlinear regression}. Wiley.
#' 
#' @examples
#' #
#' # Dobson's beetle data (generalized linear model)
#' #
#' 
#' # Complementary log-log model
#' mod <- glm(cbind(y, n-y) ~ ldose, data = beetle, 
#'            family = binomial(link = "cloglog"))
#' plotFit(mod, pch = 19, cex = 1.2, lwd = 2, 
#'         xlab = "Log dose of carbon disulphide",
#'         interval = "confidence", shade = TRUE, 
#'         col.conf = "lightskyblue")
#' 
#' # Approximate 95% confidence intervals and standard error for LD50
#' invest(mod, y0 = 0.5)
#' invest(mod, y0 = 0.5, interval = "Wald")
#' 
#' #
#' # Nasturtium example (nonlinear least-squares with replication)
#' #
#' 
#' # Log-logistic model
#' mod <- nls(weight ~ theta1/(1 + exp(theta2 + theta3 * log(conc))),
#'            start = list(theta1 = 1000, theta2 = -1, theta3 = 1),
#'            data = nasturtium)
#' plotFit(mod, lwd.fit = 2)
#'            
#' # Compute approximate 95% calibration intervals
#' invest(mod, y0 = c(309, 296, 419), interval = "inversion")
#' invest(mod, y0 = c(309, 296, 419), interval = "Wald")  
#' 
#' # Bootstrap calibration intervals. In general, nsim should be as large as 
#' # reasonably possible (say, nsim = 9999).
#' boo <- invest(mod, y0 = c(309, 296, 419), interval = "percentile", 
#'               nsim = 999, seed = 101)
#' boo  # print bootstrap summary
#' plot(boo)  # plot results
invest <- function(object, ...) {
  UseMethod("invest")
} 


#' @rdname invest
#' @export
#' @method invest lm
invest.lm <- function(object, y0, 
                      interval = c("inversion", "Wald", "percentile", "none"), 
                      level = 0.95, mean.response = FALSE, 
                      x0.name, newdata, data, 
                      boot.type = c("parametric", "nonparametric"), nsim = 999, 
                      seed = NULL, progress = FALSE, lower, upper, 
                      extendInt = "no", tol = .Machine$double.eps^0.25, 
                      maxiter = 1000, adjust = c("none", "Bonferroni"), 
                      k,  ...) {
  
  # Extract data, variable names, etc.
  .data  <- if (!missing(data)) data else eval(object$call$data, 
                                               envir = parent.frame())
  yname <- all.vars(formula(object)[[2]])
  
  # Predictor variable(s)
  multi <- FALSE
  xnames <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  if (length(xnames) != 1) {
    multi <- TRUE
    if (missing(x0.name)) {
      stop("'x0.name' is missing, please select a valid predictor variable")
    }
    xnames <- xnames[xnames != x0.name]
    if (missing(newdata)) {
      stop("'newdata' must be supplied when multiple predictor variables exist!")
    }
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data frame")
    }
    if (nrow(newdata) != 1) {
      stop("'newdata' must have a single row")
    }
    if (ncol(newdata) != length(xnames) || !all(xnames %in% names(newdata))) {
      stop(paste0("'newdata' must contain a column for each predictor variable", 
                  " used by ", deparse(substitute(object)),
                  " (except ", x0.name, ")"))
    }
  } else {
    x0.name <- xnames
  }

  # End-points for 'uniroot'  
  if (missing(lower)) lower <- min(.data[, x0.name])  # lower limit default
  if (missing(upper)) upper <- max(.data[, x0.name])  # upper limit default
  
  # Define constants
  m <- length(y0)  # number of unknowns 
  if(mean.response && m > 1) stop("Only one mean response value allowed.")
  eta <- mean(y0)  # mean unknown
  n <- length(resid(object))  # in case of missing values
  p <- length(coef(object))   # number of regression coefficients
  df1 <- n - p  # stage 1 degrees of freedom
  df2 <- m - 1  # stage 2 degrees of freedom
  var1 <- Sigma(object)^2  # stage 1 variance estimate
  var2 <- if (m == 1) 0 else var(y0)  # stage 2 variance estimate
  var_pooled <- (df1*var1 + df2*var2) / (df1 + df2)  # pooled estimate
  rat <- var_pooled / var1  # right variance?
  
  # Calculate point estimate by inverting fitted model
  x0_est <- try(uniroot(function(x) {
    nd <- if (multi) {
            cbind(newdata, makeData(x, x0.name))  # append newdata
          } else {
            makeData(x, x0.name)
          }
    predict(object, newdata = nd) - eta  #  solve yhat(x0) - eta = 0 for x0
  }, interval = c(lower, upper), extendInt = extendInt,
  tol = tol, maxiter = maxiter)$root,
  silent = TRUE)

  # Provide (informative) error message if point estimate is not found
  if (inherits(x0_est, "try-error")) {
    stop(paste("Point estimate not found in the search interval (", lower,
               ", ", upper, "). ",
               "Try tweaking the values of lower and upper. ",
               "Use plotFit for guidance.", sep = ""),
         call. = FALSE)
  }

  # Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") {  # || multi) {
    return(setNames(x0_est, x0.name))
  }
  
  # Bootstrap intervals -------------------------------------------------------
  if (interval == "percentile") {
    
    # Sanity check
    stopifnot((nsim <- as.integer(nsim[1])) > 0)
    
    # Set up progress bar (if requested)
    if (progress) { 
      pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    }
    
    # Initialize random number generator
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
      runif(1)
    if (is.null(seed)) 
      RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
      R.seed <- get(".Random.seed", envir = .GlobalEnv)
      set.seed(seed)
      RNGstate <- structure(seed, kind = as.list(RNGkind()))
      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    
    # Simulate new response vectors
    ftd <- fitted(object)  # fitted values
    res <- residuals(object) # redisuals
    boot.type <- match.arg(boot.type)
    if (boot.type == "parametric") {  
      ss <- simulate(object, nsim = nsim)
    } else {
      ss <- replicate(nsim, ftd + sample(res, replace = TRUE), simplify = FALSE)
    }
    
    # Bootstrap function
    x0Fun <- function(i) {
      
      # Update model using simulated response data
      boot_data <- eval(object$call$data)  # copy data
      boot_data[, yname] <- ss[[i]]  # simulated response vector
      boot_object <- tryCatch(update(object, data = boot_data),
                              error = function(e) NULL)
      
      # If updating the model fails, then return value is NA
      if (is.null(boot_object)) {
        ret <- NA
      } else {
        
        # Simulate new response (different from simulated response vector)
        if (mean.response) {  # regulation
          y0_star <- y0  # hold constant in bootstrap replications
        } else {  # calibration
          if (boot.type == "parametric") {
            y0_star <- y0 + rnorm(length(y0), sd = Sigma(object))
          } else {
            y0_star <- y0 + sample(res, size = length(y0), replace = TRUE)
          }
        }
        
        # Calculate point estimate
        ret <- tryCatch(uniroot(function(x) {
          predict(boot_object, newdata = makeData(x, x0.name)) - mean(y0_star)
        }, interval = c(lower, upper), extendInt = extendInt, 
        tol = tol, maxiter = maxiter)$root, 
        error = function(e) NA)
      }
      
      # Update progress bar
      if (progress) { 
        setTxtProgressBar(pb, i) 
      }
      
      # Return estimate
      ret
      
    }
    
    # Calculate bootstrap replicates
    x0_star <- sapply(seq_len(nsim), x0Fun)
    
    # Check for errors and return the runs that did not fail
    if (AnyNA(x0_star)) {
      num_fail <- sum(is.na(x0_star))
      warning("some bootstrap runs failed (", num_fail, "/", nsim, 
              ")")
      x0_star <- na.omit(x0_star)  # remove runs that failed
      attributes(x0_star) <- NULL  # remove attributes
    } else {
      num_fail <- NULL
    }
    
    # Percentile interval
    perc <- unname(quantile(x0_star, probs = c((1-level)/2, (1+level)/2)))
        
    # Create and return a bootCal object (essentially a list that can later be 
    # converted to an object of class boot)
    res <- list("estimate" = x0_est,   # original estimate
                "lower"    = perc[1],  # lower percentile
                "upper"    = perc[2],  # upper percentile
                "se"       = sd(x0_star),  # standard error
                "bias"     = mean(x0_star) - x0_est,  # estimated bias
                "bootreps" = x0_star,  # bootstrap replicates
                "nsim"     = nsim,     # number of simulations
                "level"    = level,    # desired confidence level
                "interval" = interval)  # type of interval requested
    class(res) = c("invest", "bootCal")
    attr(res, "bootFail") <- num_fail
    return(res)
     
  }
  
  # Critical value for confidence interval computations
  adjust <- match.arg(adjust)
  crit <- if (adjust == "Bonferroni" && m == 1) {
            qt((level + 2*k - 1) / (2*k), n+m-p-1)  # Bonferroni adjustment
          } else {
            qt((level+1) / 2, n+m-p-1)  # no adjustment
          }
  
  # inversion interval --------------------------------------------------------
  if (interval == "inversion") { 
    
    # Inversion function
    inversionFun <- function(x) {
      nd <- if (multi) {
        cbind(newdata, makeData(x, x0.name))  # append newdata
      } else {
        makeData(x, x0.name)
      }
      pred <- predict(object, newdata = nd, se.fit = TRUE)
      denom <- if (mean.response) pred$se.fit^2 else var_pooled/m + 
        rat*pred$se.fit^2
      (eta - pred$fit)^2/denom - crit^2
    }
    
    # Compute lower and upper confidence limits (i.e., the roots of the 
    # inversion function)
    lwr <- try(uniroot(inversionFun, interval = c(lower, x0_est), 
                       extendInt = extendInt, tol = tol, 
                       maxiter = maxiter)$root, silent = TRUE)
    upr <- try(uniroot(inversionFun, interval = c(x0_est, upper), 
                       extendInt = extendInt, tol = tol, 
                       maxiter = maxiter)$root, silent = TRUE)
    
    # Provide (informative) error message if confidence limits not found
    if (inherits(lwr, "try-error")) {
      stop(paste("Lower confidence limit not found in the search interval (", 
                 lower, ", ", upper, 
                 "). ", "Try tweaking the values of lower and upper. ", 
                 "Use plotFit for guidance.", sep = ""), 
           call. = FALSE)
    }
    if (inherits(upr, "try-error")) {
      stop(paste("Upper confidence limit not found in the search interval (", 
                 lower, ", ", upper, 
                 "). ", "Try tweaking the values of lower and upper. ", 
                 "Use plotFit for guidance.", sep = ""), 
           call. = FALSE)
    }
    
    # Store results in a list
    res <- list("estimate" = x0_est, 
                "lower" = lwr, 
                "upper" = upr, 
                "interval" = interval)
  } 
  
  # Wald interval --------------------------------------------------------------
  if (interval == "Wald") { 
    
    # Function of parameters whose gradient is required
    object_copy <- object # FIXME: Is a copy really needed?
    dmFun <- function(params) {
      if (mean.response) {
        object_copy$coefficients <- params
        z <- eta
      } else {
        object_copy$coefficients <- params[-length(params)]
        z <- params[length(params)]
      }
      uniroot(function(x) { 
        nd <- if (multi) {
          cbind(newdata, makeData(x, x0.name))  # append newdata
        } else {
          makeData(x, x0.name)
        }
        predict(object_copy, newdata = nd) - z
      }, interval = c(lower, upper), extendInt = extendInt, 
      tol = tol, maxiter = maxiter)$root
    }
    
    # Variance-covariane matrix
    if (mean.response) {
      params <- coef(object)
      covmat <- vcov(object)
    } else {
      params <- c(coef(object), eta)
      covmat <- diag(p + 1)
      covmat[p + 1, p + 1] <- var_pooled / m
      covmat[1:p, 1:p] <- vcov(object)
    }
    
    # Calculate gradient, and return standard error
    gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
    se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
    
    # Store results in a list
    res <- list("estimate" = x0_est, 
                "lower" = x0_est - crit * se, 
                "upper" = x0_est + crit * se, 
                "se" = se,
                "interval" = interval)
  }
  
  # Assign class label and return results
  class(res) <- "invest"
  res
  
}


#' @rdname invest
#' @export
#' @method invest glm
invest.glm <- function(object, y0, 
                       interval = c("inversion", "Wald", "percentile", "none"),  
                       level = 0.95, lower, upper, 
                       x0.name, newdata, data,
                       tol = .Machine$double.eps^0.25, maxiter = 1000, ...) {
  
  # NOTE: Currently, this function only works for the case 
  #       mean.response = TRUE. 
  
  # Extract data, variable names, etc.
  .data  <- if (!missing(data)) data else eval(object$call$data, 
                                               envir = parent.frame())
  
  # End-points for 'uniroot'  
  if (missing(lower)) lower <- min(.data[, x0.name])  # lower limit default
  if (missing(upper)) upper <- max(.data[, x0.name])  # upper limit default
  
  # Calculations should be done on the "link" scale!
  eta <- family(object)$linkfun(y0)
  
  # Predictor variable(s)
  multi <- FALSE
  xnames <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  if (length(xnames) != 1) {
    multi <- TRUE
    if (missing(x0.name)) {
      stop("'x0.name' is missing, please select a valid predictor variable")
    }
    xnames <- xnames[xnames != x0.name]
    if (missing(newdata)) {
      stop("'newdata' must be supplied when multiple predictor variables exist!")
    }
    if (!is.data.frame(newdata)) {
      stop("'newdata' must be a data frame")
    }
    if (nrow(newdata) != 1) {
      stop("'newdata' must have a single row")
    }
    if (ncol(newdata) != length(xnames) || !all(xnames %in% names(newdata))) {
      stop(paste0("'newdata' must contain a column for each predictor variable", 
                  " used by ", deparse(substitute(object)),
                  " (except ", x0.name, ")"))
    }
  } else {
    x0.name <- xnames
  }
  
  # Calculate point estimate by inverting fitted model
  x0_est <- try(uniroot(function(x) {
    nd <- if (multi) {
      cbind(newdata, makeData(x, x0.name))  # append newdata
    } else {
      makeData(x, x0.name)
    }
    predict(object, newdata = nd, type = "link") - eta  #  solve yhat(x0) - eta = 0 for x0
  }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root, 
  silent = TRUE)
  
  # Provide (informative) error message if point estimate is not found
  if (inherits(x0_est, "try-error")) {
    stop(paste("Point estimate not found in the search interval (", lower, 
               ", ", upper, "). ", 
               "Try tweaking the values of lower and upper. ",
               "Use plotFit for guidance.", sep = ""), 
         call. = FALSE)
  }
  
  # Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") {
    return(setNames(x0_est, x0.name))
  }
  
  # Stop and print error if user requests bootstrap intervals
  if (interval == "percentile") {
    stop("Bootstrap intervals not available for 'glm' objects.", call. = FALSE)
  }
  
  # Critical value
  crit <- qnorm((level + 1) / 2)  # quantile from standard normal
  
  # Inversion interval
  #
  # Based on exercise 5.31 on pg. 207 of Categorical Data Analysis (2nd ed.) by 
  # Alan Agresti.
  if (interval == "inversion") { 

    # Inversion function
    inversionFun <- function(x) {
      nd <- if (multi) {
        cbind(newdata, makeData(x, x0.name))  # append newdata
      } else {
        makeData(x, x0.name)
      }
      pred <- predict(object, newdata = nd, se.fit = TRUE, type = "link")
      ((eta - pred$fit) ^ 2) / (pred$se.fit ^ 2) - crit^2
    }
    
    # Compute lower and upper confidence limits (i.e., the roots of the 
    # inversion function)
    lwr <- try(uniroot(inversionFun, interval = c(lower, x0_est), tol = tol, 
                       maxiter = maxiter)$root, silent = TRUE)
    upr <- try(uniroot(inversionFun, interval = c(x0_est, upper), tol = tol, 
                       maxiter = maxiter)$root, silent = TRUE)
    
    # Provide (informative) error message if confidence limits not found
    if (inherits(lwr, "try-error")) {
      stop(paste("Lower confidence limit not found in the search interval (", 
                 lower, ", ", upper, 
                 "). ", "Try tweaking the values of lower and upper. ", 
                 "Use plotFit for guidance.", sep = ""), 
           call. = FALSE)
    }
    if (inherits(upr, "try-error")) {
      stop(paste("Upper confidence limit not found in the search interval (", 
                 lower, ", ", upper, 
                 "). ", "Try tweaking the values of lower and upper. ", 
                 "Use plotFit for guidance.", sep = ""), 
           call. = FALSE)
    }
    
    # Store results in a list
    res <- list("estimate" = x0_est, 
                "lower" = lwr, 
                "upper" = upr, 
                "interval" = interval)
    
  }
  
  # Wald interval -------------------------------------------------------------
  if (interval == "Wald") { 
    
    # Function of parameters whose gradient is required
    object_copy <- object # FIXME: Is a copy really needed?
    dmFun <- function(params) {
      object_copy$coefficients <- params
      z <- eta
      uniroot(function(x) { 
        nd <- if (multi) {
          cbind(newdata, makeData(x, x0.name))  # append newdata
        } else {
          makeData(x, x0.name)
        }
        predict(object_copy, newdata = nd, type = "link") - z
      }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root
    }
    
    # Variance-covariane matrix
    params <- coef(object)
    covmat <- vcov(object)
    
    # Calculate gradient, and return standard error
    gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
    se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
    
    # Store results in a list
    res <- list("estimate" = x0_est, 
                "lower" = x0_est - crit * se, 
                "upper" = x0_est + crit * se, 
                "se" = se,  # FIXME: How should this get transformed?
                "interval" = interval)
  }
  
  # Assign class label and return results
  class(res) <- "invest"
  res
  
}


#' @rdname invest
#' @export
#' @method invest nls
invest.nls <- function(object, y0, 
                       interval = c("inversion", "Wald", "percentile", "none"),  
                       level = 0.95, mean.response = FALSE, data, 
                       boot.type = c("parametric", "nonparametric"), nsim = 1, 
                       seed = NULL, progress = FALSE, lower, upper, 
                       tol = .Machine$double.eps^0.25, maxiter = 1000, 
                       adjust = c("none", "Bonferroni"), k, ...) {
  
  # No support for the Golub-Pereyra algorithm for partially linear 
  # least-squares models
  if (object$call$algorithm == "plinear") {
    stop(paste("The Golub-Pereyra algorithm for partially linear least-squares 
               models is currently not supported."))
  }
  
  # Extract data, variable names, etc.
  .data  <- if (!missing(data)) data else eval(object$call$data, 
                                               envir = parent.frame())
  x0.name <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(x0.name) != 1) stop("Only one independent variable allowed.")
  if (missing(lower)) lower <- min(.data[, x0.name])  # lower limit default
  if (missing(upper)) upper <- max(.data[, x0.name])  # upper limit default
  
  # Set up for inverse estimation
  m <- length(y0)  # number of unknowns 
  if(mean.response && m > 1) stop("Only one mean response value allowed.")
  eta <- mean(y0)  # mean response
  n <- length(resid(object))  # sample size
  p <- length(coef(object))  # number of parameters
  var_pooled <- Sigma(object)^2  # residual variance
  
  # Calculate point estimate by inverting fitted model
  x0_est <- try(uniroot(function(x) {
    predict(object, newdata = makeData(x, x0.name)) - eta
  }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root, 
  silent = TRUE)
  
  # Provide (informative) error message if point estimate is not found
  if (inherits(x0_est, "try-error")) {
    stop(paste("Point estimate not found in the search interval (", lower, 
               ", ", upper, "). ", 
               "Try tweaking the values of lower and upper. ",
               "Use plotFit for guidance.", sep = ""), 
         call. = FALSE)
  }
  
  # Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") return(x0_est)
  
  # Bootstrap -----------------------------------------------------------------
  if (interval == "percentile") {
    
    # Sanity check
    stopifnot((nsim <- as.integer(nsim[1])) > 0)
    
    # Set up progress bar (if requested)
    if (progress) { 
      pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    }
    
    # Initialize random number generator
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
      runif(1)
    if (is.null(seed)) 
      RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
      R.seed <- get(".Random.seed", envir = .GlobalEnv)
      set.seed(seed)
      RNGstate <- structure(seed, kind = as.list(RNGkind()))
      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    
    # Simulate new response vectors
    ftd <- fitted(object)  # fitted values
    res <- residuals(object) # redisuals
    boot.type <- match.arg(boot.type)
    if (boot.type == "parametric") {  
      ss <- simulate(object, nsim = nsim)
    } else {
      ss <- replicate(nsim, ftd + sample(res, replace = TRUE), simplify = FALSE)
    }
    
    # Function to calculate inverse estimate
    x0Fun <- function(i) {
      
      # Update model using simulated response data
      boot_data <- eval(object$call$data)  # copy data
      boot_data[, yname] <- ss[[i]]  # simulated response vector
      boot_object <- tryCatch(update(object, data = boot_data),
                              error = function(e) NULL)
      
      # If updating the model fails, then return value is NA
      if (is.null(boot_object)) {
        ret <- NA
      } else {
        
        # Simulate new response (different from simulated response vector)
        if (mean.response) {  # regulation
          y0_star <- y0  # hold constant in bootstrap replications
        } else {  # calibration
          if (boot.type == "parametric") {
            y0_star <- y0 + rnorm(length(y0), sd = Sigma(object))
          } else {
            y0_star <- y0 + sample(res, size = length(y0), replace = TRUE)
          }
        }
        
        # Calculate point estimate
        ret <- tryCatch(uniroot(function(x) {
            predict(boot_object, newdata = makeData(x, x0.name)) - mean(y0_star)
          }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root, 
          error = function(e) NA)
      }
      
      # Update progress bar
      if (progress) { 
        setTxtProgressBar(pb, i) 
      }
      
      # Return estimate
      ret
      
    }
    
    # Calculate bootstrap replicates
    x0_star <- sapply(seq_len(nsim), x0Fun)
    
    # Check for errors and return the runs that did not fail
    if (AnyNA(x0_star)) {
      num_fail <- sum(is.na(x0_star))
      warning("some bootstrap runs failed (", num_fail, "/", nsim, 
              ")")
      x0_star <- na.omit(x0_star)  # remove runs that failed
      attributes(x0_star) <- NULL  # remove attributes
    } else {
      num_fail <- NULL
    }

    # Percentile interval
    perc <- unname(quantile(x0_star, probs = c((1-level)/2, (1+level)/2)))
    
    # Create and return a bootCal object (essentially a list that can later be 
    # converted to an object of class boot)
    res <- list("estimate" = x0_est,   # original estimate
                "lower"    = perc[1],  # lower percentile
                "upper"    = perc[2],  # upper percentile
                "se"       = sd(x0_star),  # standard error
                "bias"     = mean(x0_star) - x0_est,  # estimated bias
                "bootreps" = x0_star,  # bootstrap replicates
                "nsim"     = nsim,     # number of simulations
                "level"    = level,    # desired confidence level
                "interval" = interval)  # type of interval requested
    class(res) = c("invest", "bootCal")
    attr(res, "bootFail") <- num_fail
    return(res)

  }
  
  # Critical value for confidence interval computations
  adjust <- match.arg(adjust)
  crit <- if (adjust == "Bonferroni" && m == 1) {
    qt((level + 2*k - 1) / (2*k), n+m-p-1)  # Bonferroni adjustment
  } else {
    qt((level + 1) / 2, n+m-p-1)  # no adjustment
  }

  # Inversion interval --------------------------------------------------------
  if (interval == "inversion") {
    
    # Inversion function
    inversionFun <- function(x) {
      pred <- predFit(object, newdata = makeData(x, x0.name), se.fit = TRUE) 
      denom <- if (mean.response) {
                 pred$se.fit^2 
               } else {
                 var_pooled/m + pred$se.fit^2
               }
      (eta - pred$fit)^2/denom - crit^2
    }
    
    # Compute lower and upper confidence limits (i.e., the roots of the 
    # inversion function)
    lwr <- try(uniroot(inversionFun, interval = c(lower, x0_est), tol = tol, 
                       maxiter = maxiter)$root, silent = TRUE)
    upr <- try(uniroot(inversionFun, interval = c(x0_est, upper), tol = tol, 
                       maxiter = maxiter)$root, silent = TRUE)
    
    # Provide (informative) error message if confidence limits not found
    if (inherits(lwr, "try-error")) {
      stop(paste("Lower confidence limit not found in the search interval (", 
                 lower, ", ", upper, 
                 "). ", "Try tweaking the values of lower and upper. ", 
                 "Use plotFit for guidance.", sep = ""), 
           call. = FALSE)
    }
    if (inherits(upr, "try-error")) {
      stop(paste("Upper confidence limit not found in the search interval (", 
                 lower, ", ", upper, 
                 "). ", "Try tweaking the values of lower and upper. ", 
                 "Use plotFit for guidance.", sep = ""), 
           call. = FALSE)
    }
    
    # Store results in a list
    res <- list("estimate" = x0_est, 
                "lower" = lwr, 
                "upper" = upr, 
                "interval" = interval)
  } 
  
  # Wald interval --------------------------------------------------------------
  if (interval == "Wald") { 
    
    # Function of parameters whose gradient is required
    object_copy <- object # FIXME: Is a copy really needed?
    dmFun <- function(params) {
      if (mean.response) {
        object_copy$m$setPars(params)
        z <- eta
      } else {
        object_copy$m$setPars(params[-length(params)])
        z <- params[length(params)]
      }
      uniroot(function(x) { 
        predict(object_copy, newdata = makeData(x, x0.name)) - z
      }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root
    }
    
    # Variance-covariance matrix
    if (mean.response) {
      params <- coef(object)
      covmat <- vcov(object)
    } else {
      params <- c(coef(object), eta)
      covmat <- diag(p + 1)
      covmat[p + 1, p + 1] <- var_pooled/m
      covmat[1:p, 1:p] <- vcov(object)
    }
    
    # Calculate gradient, and return standard error
    gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
    se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
    
    # Store results in a list
    res <- list("estimate" = x0_est, 
                "lower" = x0_est - crit * se, 
                "upper" = x0_est + crit * se, 
                "se" = se,
                "interval" = interval)
  }
  
  # Assign class label and return results
  class(res) <- "invest"
  res
  
}


#' @rdname invest
#' @export
#' @method invest lme
invest.lme <- function(object, y0, 
                       interval = c("inversion", "Wald", "percentile", "none"),  
                       level = 0.95, mean.response = FALSE, data, lower, upper, 
                       q1, q2, tol = .Machine$double.eps^0.25, maxiter = 1000, 
                       ...) 
{
  
  # Extract data, variable names, etc.
  .data  <- if (!missing(data)) data else object$data
  x0.name <- intersect(all.vars(formula(object)[[3]]), colnames(.data)) 
  yname <- all.vars(formula(object)[[2]])
  if (length(x0.name) != 1) stop("Only one independent variable allowed.")
  if (missing(lower)) lower <- min(.data[, x0.name])  # lower limit default
  if (missing(upper)) upper <- max(.data[, x0.name])  # upper limit default
  
  # Set up for inverse estimation
#   if(m > 1) stop("Only one response value allowed.")
  m <- length(y0)
  if(mean.response && m > 1) stop("Only one mean response value allowed.")
  eta <- mean(y0)
  if (m != 1) stop('Only a single unknown allowed for objects of class "lme".')
  N <- length(resid(object)) 
  p <- length(fixef(object))
#   res.var <- Sigma(object)^2  # residual variance
  
  # Critical value. Oman (1998. pg. 445) suggests a t(1-alpha/2, N-1) dist.
  if (missing(q1)) q1 <- qnorm((1-level) / 2)
  if (missing(q2)) q2 <- qnorm((1+level) / 2)
  
  # Calculate point estimate by inverting fitted model
  x0_est <- try(uniroot(function(x) {
    predict(object, newdata = makeData(x, x0.name), level = 0) - eta
  }, interval = c(lower, upper), tol = tol, maxiter = maxiter)$root, 
  silent = TRUE)
  
  # Provide (informative) error message if point estimate is not found
  if (inherits(x0_est, "try-error")) {
    stop(paste("Point estimate not found in the search interval (", lower, 
               ", ", upper, "). ", 
               "Try tweaking the values of lower and upper.", 
               sep = ""), 
         call. = FALSE)
  }
  
  # Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") return(x0_est)
  
  # Stop and print error if user requested bootstrap intervals
  if (interval == "percentile") {
    stop("Bootstrap intervals not available for 'lme' objects.", call. = FALSE)
  }

  # Estimate variance of new response
  if (!mean.response) var.y0 <- varY(object, newdata = makeData(x0_est, x0.name))
  
  # Inversion interval --------------------------------------------------------
  if (interval == "inversion") { 
    
    # Inversion function
    inversionFun <- function(x, bound = c("lower", "upper")) {
      pred <- predFit(object, newdata = makeData(x, x0.name), se.fit = TRUE)
      denom <- if (mean.response) {
                 pred[, "se.fit"] 
               } else {
                 sqrt(var.y0 + pred[, "se.fit"]^2)
               }
      bound <- match.arg(bound)
      if (bound == "upper") {
        (eta - pred[, "fit"])/denom - q1
      } else {
        (eta - pred[, "fit"])/denom - q2
      }
    }
    
    # Compute lower and upper confidence limits (i.e., the roots of the 
    # inversion function)
    lwr <- try(uniroot(inversionFun, interval = c(lower, x0_est), 
                       bound = "lower", tol = tol, maxiter = maxiter)$root, 
               silent = TRUE)
    upr <- try(uniroot(inversionFun, interval = c(x0_est, upper), 
                       bound = "upper", tol = tol, maxiter = maxiter)$root, 
               silent = TRUE)
    
    # Provide (informative) error message if confidence limits not found
    if (inherits(lwr, "try-error")) {
      stop(paste("Lower confidence limit not found in the search interval (", 
                 lower, ", ", upper, 
                 "). ", "Try tweaking the values of lower and upper.", 
                 sep = ""), 
           call. = FALSE)
    }
    if (inherits(upr, "try-error")) {
      stop(paste("Upper confidence limit not found in the search interval (", 
                 lower, ", ", upper, 
                 "). ", "Try tweaking the values of lower and upper.",
                 sep = ""), 
           call. = FALSE)
    }
    
    # Store results in a list
    res <- list("estimate" = x0_est, 
                "lower" = lwr, 
                "upper" = upr, 
                "interval" = interval)
    
  }
  
  # Wald interval -------------------------------------------------------------
  if (interval == "Wald") { 
    
    # Function of parameters whose gradient is required
    dmFun <- function(params) {
      fun <- function(x) {
        X <- model.matrix(eval(object$call$fixed)[-2], 
                          data = makeData(x, x0.name))
        if (mean.response) {
          X %*% params - eta
        } else {
          X %*% params[-length(params)] - params[length(params)]
        }
      }
      uniroot(fun, lower = lower, upper = upper, tol = tol, 
              maxiter = maxiter)$root
    }
    
    # Variance-covariance matrix
    if (mean.response) {
      params <- fixef(object)
      covmat <- vcov(object)
    } else {
      params <- c(fixef(object), eta)
      covmat <- diag(p + 1)
      covmat[p + 1, p + 1] <- var.y0
      covmat[1:p, 1:p] <- vcov(object)
    }
    
    # Calculate gradient, and return standard error
    gv <- attr(numericDeriv(quote(dmFun(params)), "params"), "gradient")
    se <- as.numeric(sqrt(gv %*% covmat %*% t(gv)))
    
    # Store results in a list
    res <- list("estimate" = x0_est, 
                "lower" = x0_est - q2 * se, 
                "upper" = x0_est + q2 * se,
                "se" = se,
                "interval" = interval)
    
  } 
  
  # Assign class label and return results
  class(res) <- "invest"
  res
  
}


#' @keywords internal
#' @export
print.invest <- function(x, digits = getOption("digits"), ...) {
  if (x$interval == "inversion") print(round(unlist(x[1:3]), digits))
  if (x$interval == "Wald") print(round(unlist(x[1:4]), digits))
  if (x$interval == "percentile") print(round(unlist(x[1:5]), digits))
  invisible(x)
} 


#' Plots of the Output of a Bootstrap Calibration Simulation
#' 
#' This takes a bootstrap calibration object and produces plots for the 
#' bootstrap replicates of the inverse estimate.
#' 
#' @rdname plot.bootCal
#' @importFrom graphics hist par plot
#' @importFrom stats qqline qqnorm
#' @export
#' @method plot bootCal
#' 
#' @param x An object that inherits from class \code{"bootCal"}.
#' @param ... Additional optional arguments. At present, no optional arguments 
#'            are used.
#' @export
plot.bootCal <- function(x, ...) {
  
  t <- x$bootreps  # bootstrap replicates
  t0 <- x$estimate  # original estimate
  
  # Calculate number of histogram breaks (better than default)
  if (!is.null(t0)) {
    nclass <- min(max(ceiling(length(t)/25), 10), 100)
    rg <- range(t)
    if (t0 < rg[1L]) 
      rg[1L] <- t0
    else if (t0 > rg[2L]) 
      rg[2L] <- t0
    rg <- rg + 0.05 * c(-1, 1) * diff(rg)
    lc <- diff(rg)/(nclass - 2)
    n1 <- ceiling((t0 - rg[1L])/lc)
    n2 <- ceiling((rg[2L] - t0)/lc)
    bks <- t0 + (-n1:n2) * lc
  }
  
  # Plots
  par(mfrow = c(1, 2))
  hist(t, breaks = bks, probability = TRUE, main = "",
       xlab = "Bootstrap replicate")
  qqnorm(t, xlab = "Standard normal quantile", ylab = "Bootstrap quantile",
         main = "")
  qqline(t)
  par(mfrow = c(1, 1))
  invisible(x)
  
}

#' Calibration for the simple linear regression model.
#' 
#' The function \code{calibrate} computes the maximum likelihood estimate and a
#' condfidence interval for the unknown predictor value that corresponds to an 
#' observed value of the response (or vector thereof) or specified value of the 
#' mean response. See the reference listed below for more details.
#'  
#' @param object An object that inherits from class \code{"lm"}, a matrix, a 
#' list, or a data frame.
#' @param formula A formula of the form \code{y ~ x}.
#' @param data an optional data frame, list or environment (or object coercible 
#' by \code{as.data.frame} to a data frame) containing the variables in the 
#' model. If not found in data, the variables are taken from 
#' \code{environment(formula)}, typically the environment from which \code{lm}
#' is called. 
#' @param subset An optional vector specifying a subset of observations to be 
#' used in the fitting process.
#' @param na.action a function which indicates what should happen when the data 
#' contain \code{NA}s. 
#' @param y0 The value of the observed response(s) or specified value of the
#'           mean response.
#' @param interval The method to use for forming a confidence interval.
#' @param level A numeric scalar between 0 and 1 giving the confidence level for 
#'              the interval to be calculated. 
#' @param mean.response Logicial indicating whether confidence intervals should 
#' correspond to an observed response(s) (\code{FALSE}) or a specified value of 
#' the mean response (\code{TRUE}). Default is \code{FALSE}.
#' @param adjust A logical value indicating if an adjustment should be made to
#'               the critical value used in calculating the confidence interval.
#'               This useful for when the calibration curve is to be used 
#'               multiple, say k, times.
#' @param k The number times the calibration curve is to be used for computing a 
#'          confidence interval. Only needed when \code{adjust = TRUE}.
#' @param ... Additional optional arguments. At present, no optional arguments 
#'            are used.
#'            
#' @return An object of class \code{"invest"} containing the following 
#'         components:
#' \itemize{
#'   \item \code{estimate} The estimate of x0.
#'   \item \code{lwr} The lower confidence limit for x0.
#'   \item \code{upr} The upper confidence limit for x0.
#'   \item \code{se} An estimate of the standard error (Wald interval only).
#'   \item \code{interval} The method used for calculating \code{lower} and 
#'                   \code{upper} (only used by \code{print} method).
#' }
#' 
#' @references 
#' Graybill, F. A., and Iyer, H. K. (1994)
#' \emph{Regression analysis: Concepts and Applications}. Duxbury Press.
#' 
#' Miller, R. G. (1981)
#' \emph{Simultaneous Statistical Inference}. Springer-Verlag.
#' 
#' @rdname calibrate
#' 
#' @aliases print.calibrate
#' 
#' @importFrom stats coef formula lm lm.fit model.extract model.frame 
#' @importFrom stats model.matrix model.response na.fail qf qt resid terms var
#' @export
#'
#' @note The function \code{invest} is more general, but based on numerical
#' techniques to find the solution. When the underlying model is that of the 
#' simple linear regression model with normal errors, closed-form expressions
#' exist which are utilized by the function \code{calibrate}.
#' 
#' @examples
#' #
#' # Arsenic example (simple linear regression with replication)
#' #
#' 
#' # Inverting a prediction interval for an individual response
#' arsenic.lm <- lm(measured ~ actual, data = arsenic)
#' plotFit(arsenic.lm, interval = "prediction", shade = TRUE, 
#'         col.pred = "lightblue")
#' (cal <- calibrate(arsenic.lm, y0 = 3, interval = "inversion"))
#' abline(h = 3)
#' segments(cal$estimate, 3, cal$estimate, par()$usr[3])
#' arrows(cal$lower, 3, cal$lower, par()$usr[3])
#' arrows(cal$upper, 3, cal$upper, par()$usr[3])
#' 
#' #
#' # Crystal weight example (simple linear regression)
#' #
#' 
#' # Inverting a confidence interval for the mean response
#' crystal.lm <- lm(weight ~ time, data = crystal)
#' plotFit(crystal.lm, interval = "confidence", shade = TRUE,
#'         col.conf = "lightblue")
#' (cal <- calibrate(crystal.lm, y0 = 8, interval = "inversion", 
#'                   mean.response = TRUE))
#' abline(h = 8)
#' segments(cal$estimate, 8, cal$estimate, par()$usr[3])
#' arrows(cal$lower, 8, cal$lower, par()$usr[3])
#' arrows(cal$upper, 8, cal$upper, par()$usr[3])
#'
#' # Wald interval and approximate standard error based on the delta method
#' calibrate(crystal.lm, y0 = 8, interval = "Wald", mean.response = TRUE)
calibrate <- function(object, ...) {
  UseMethod("calibrate")
}


#' @rdname calibrate
#' @export
#' @method calibrate default
calibrate.default <- function(object, y0, 
                              interval = c("inversion", "Wald", "none"), 
                              level = 0.95, mean.response = FALSE, 
                              adjust = c("none", "Bonferroni", "Scheffe"), k, 
                              ...) {

  # Extract needed components from fitted model
  if (inherits(object, "matrix")) {
    x <- object[, 1]
    y <- object[, 2]
  } else if (inherits(object, "data.frame")) {
    object <- data.matrix(object)
    x <- object[, 1]
    y <- object[, 2]
  } else if (inherits(object, "list")) {
    x <- object[[1]]
    y <- object[[2]]
    if (length(x) != length(y)) {
      stop(paste("Components of '", deparse(substitute(object)), 
                 "' not of same length.", sep = ""), call. = FALSE)
    }
  } else {
    stop("'", paste(deparse(substitute(object)), 
                    "' is not a valid matrix, list, or data frame.", sep = ""),
         call. = FALSE)
  }
  eta <- mean(y0)  # mean of new observations
  m <- length(y0)  # number of new observations
  if (mean.response && m > 1) stop("Only one mean response value allowed.")
  
  # Fit a simple linear regression model and compute necessary components
  z <- lm.fit(cbind(1, x), y)
  b <- unname(z$coefficients)
  n <- length(r <- z$residuals)  # sample size and residuals
  DF <- (DF1 <- n - 2) + (DF2 <- m - 1)  # degrees of freedom
  var1 <- sum(r ^ 2) / z$df.residual  # stage 1 variance estimate
  var2 <- if (m == 1) 0 else var(y0)  # stage 2 variance estimate
  var.pooled <- (DF1 * var1 + DF2 * var2) / DF  # pooled estimate of variance
  sigma.pooled <- sqrt(var.pooled)  # sqrt of pooled variance estimate
  ssx <- sum((x - mean(x))^2)  # sum-of-squares for x, Sxx
  x0.mle <- (eta - b[1L])/b[2L]  # MLE of x0
  
  # Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") return(x0.mle)
  
  # Adjustment for simultaneous intervals
  adjust <- match.arg(adjust)  # FIXME: Does simultaneous work for m > 1?
  crit <- if (m != 1 || adjust == "none") qt((1 + level)/2, n+m-3) else {
    switch(adjust,
          "Bonferroni" = qt((level + 2*k - 1) / (2*k), n+m-3),
          "Scheffe"    = sqrt(k * qf(level, k, n+m-3)))
  }

  # Inversion interval --------------------------------------------------------
  if (interval == "inversion") { 

    c1 <- b[2L]^2 - (sigma.pooled^2 * crit^2)/ssx
    c2 <- if (mean.response) {
      c1/n + (eta - mean(y))^2/ssx
    } else {
      c1*(1/m + 1/n) + (eta - mean(y))^2/ssx
    }
    c3 <- b[2L] * (eta - mean(y))
    c4 <- crit * sigma.pooled
    
    # FIXME: catch errors and throw an appropriate warning
    if (c1 < 0 && c2 <= 0) {
      
      warning("The calibration line is not well determined.", call. = FALSE)
      lwr <- -Inf
      upr <- Inf
      
    } else {
      
      lwr <- mean(x) + (c3 - c4*sqrt(c2))/c1
      upr <- mean(x) + (c3 + c4*sqrt(c2))/c1
      if (c1 < 0 && c2 > 0) {
        stop(paste("The calibration line is not well determined. The resulting \nconfidence region is the union of two semi-infinite intervals:\n(", -Inf, ",", 
                   round(upr, 4), ") U (", round(lwr, 4), ",", Inf, ")"), 
             call. = FALSE)
      }
      
    }
    res <- list("estimate" = x0.mle, 
                "lower"    = lwr, 
                "upper"    = upr, 
                "interval" = interval)
  
  }
    
  # Wald interval  
  if (interval == "Wald") { 
    
    # Compute standard error for Wald interval
    se <- if (mean.response) {
      abs((sigma.pooled/b[2]))*sqrt((1/n + (x0.mle - mean(x))^2/ssx))
    } else {
      abs((sigma.pooled/b[2]))*sqrt((1/m + 1/n + (x0.mle - mean(x))^2/ssx))
    }

    # Store results in a list
    res <- list("estimate" = x0.mle, 
                "lower"    = x0.mle - crit * se,
                "upper"    = x0.mle + crit * se, 
                "se"       = se, 
                "interval" = interval)
    
  } 
  
  # Assign class label and return results
  class(res) <- "invest"
  return(res)
  
} 


#' @rdname calibrate
#' @export
#' @method calibrate formula
calibrate.formula <- function(formula, data = NULL, ..., subset, 
                              na.action = na.fail) {
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, sys.parent()))) {
    m$data <- as.data.frame(data)
  }
  m$... <- NULL
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  Terms <- attr(m, "terms")
  attr(Terms, "intercept") <- 0
  y <- model.extract(m, "response")
  mm <- model.matrix(Terms, m)
  if (ncol(mm) > 1) stop("This function only works for the simple linear regression model (i.e., y ~ x).")
  x <- as.numeric(mm)
  calibrate(cbind(x, y), ...)
} 


#' #' @rdname calibrate
#' #' @export
#' #' @method calibrate lm
#' calibrate.lm <- function(object, ...) {
#'   calibrate(formula(object), data = eval(object$call$data), ...)
#' } 


#' @rdname calibrate
#' @export
#' @method calibrate lm
calibrate.lm <- function(object, y0, interval = c("inversion", "Wald", "none"), 
                         level = 0.95, mean.response = FALSE, 
                         adjust = c("none", "Bonferroni", "Scheffe"), k, ...) {
  
  # Check model formula for correctness
  xname <- all.vars(formula(object)[[3L]])
  yname <- all.vars(formula(object)[[2L]])
  if (length(xname) != 1L) {
    stop("Only one independent variable allowed.")
  }
  if (length(yname) != 1L) {
    stop("Only one dependent variable allowed.")
  }
  
  # Check for intercept using terms object from model fit. Alternatively, this
  # can also be checked by testing if the first column name in model.matrix is 
  # equal to "(Intercept)".
  if (!attr(terms(object), "intercept")) {
    stop(paste(deparse(substitute(object)), "must contain an intercept."))
  }
  
  # Extract x values and y values from model frame
  mf <- model.frame(object)
  if (ncol(mf) != 2) {
    stop("calibrate only works for the simple linear regression model.")
  } 
  x <- model.matrix(object)[, 2]
  y <- model.response(mf)

  # Eta - mean response or mean of observed respone values
  eta <- mean(y0)  # mean of new observations
  m <- length(y0)  # number of new observations
  if (mean.response && m > 1) stop("Only one mean response value allowed.")
  
  # Fit a simple linear regression model and compute necessary components
  b <- unname(object$coefficients)
  n <- length(r <- object$residuals)  # sample size and residuals
  DF <- (DF1 <- n - 2) + (DF2 <- m - 1)  # degrees of freedom
  var1 <- sum(r ^ 2) / object$df.residual  # stage 1 variance estimate
  var2 <- if (m == 1) 0 else var(y0)  # stage 2 variance estimate
  var.pooled <- (DF1 * var1 + DF2 * var2) / DF  # pooled estimate of variance
  sigma.pooled <- sqrt(var.pooled)  # sqrt of pooled variance estimate
  ssx <- sum((x - mean(x))^2)  # sum-of-squares for x, Sxx
  x0.mle <- (eta - b[1L])/b[2L]  # MLE of x0
  
  # Return point estimate only
  interval <- match.arg(interval)
  if (interval == "none") return(x0.mle)
  
  # Adjustment for simultaneous intervals
  adjust <- match.arg(adjust)  # FIXME: Does simultaneous work for m > 1?
  crit <- if (m != 1 || adjust == "none") qt((1 + level)/2, n+m-3) else {
    switch(adjust,
           "Bonferroni" = qt((level + 2*k - 1) / (2*k), n+m-3),
           "Scheffe"    = sqrt(k * qf(level, k, n+m-3)))
  }
  
  # Inversion interval --------------------------------------------------------
  if (interval == "inversion") { 
    
    c1 <- b[2L]^2 - (sigma.pooled^2 * crit^2)/ssx
    c2 <- if (mean.response) {
      c1/n + (eta - mean(y))^2/ssx
    } else {
      c1*(1/m + 1/n) + (eta - mean(y))^2/ssx
    }
    c3 <- b[2L] * (eta - mean(y))
    c4 <- crit * sigma.pooled
    
    # FIXME: catch errors and throw an appropriate warning
    if (c1 < 0 && c2 <= 0) {
      
      warning("The calibration line is not well determined.", call. = FALSE)
      lwr <- -Inf
      upr <- Inf
      
    } else {
      
      lwr <- mean(x) + (c3 - c4*sqrt(c2))/c1
      upr <- mean(x) + (c3 + c4*sqrt(c2))/c1
      if (c1 < 0 && c2 > 0) {
        stop(paste("The calibration line is not well determined. The resulting \nconfidence region is the union of two semi-infinite intervals:\n(", -Inf, ",", 
                   round(upr, 4), ") U (", round(lwr, 4), ",", Inf, ")"), 
             call. = FALSE)
      }
      
    }
    res <- list("estimate" = x0.mle, 
                "lower"    = lwr, 
                "upper"    = upr, 
                "interval" = interval)
    
  }
  
  # Wald interval  
  if (interval == "Wald") { 
    
    # Compute standard error for Wald interval
    se <- if (mean.response) {
      abs((sigma.pooled/b[2]))*sqrt((1/n + (x0.mle - mean(x))^2/ssx))
    } else {
      abs((sigma.pooled/b[2]))*sqrt((1/m + 1/n + (x0.mle - mean(x))^2/ssx))
    }
    
    # Store results in a list
    res <- list("estimate" = x0.mle, 
                "lower"    = x0.mle - crit * se,
                "upper"    = x0.mle + crit * se, 
                "se"       = se, 
                "interval" = interval)
    
  } 
  
  # Assign class label and return results
  class(res) <- "invest"
  return(res)
  
} 

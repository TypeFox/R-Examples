#' @title Maximum Likelihood Estimation for Lambert W\eqn{ \times} F distributions
#' @name MLE_LambertW
#' 
#' @description Maximum Likelihood Estimation (MLE) for Lambert W \eqn{\times F}
#'     distributions computes \eqn{\widehat{\theta}_{MLE}}.
#' 
#' For \code{type = "s"}, the skewness parameter \eqn{\gamma} is estimated and
#' \eqn{\delta = 0} is held fixed; for \code{type = "h"} the one-dimensional
#' \eqn{\delta} is estimated and \eqn{\gamma = 0} is held fixed; and for
#' \code{type = "hh"} the 2-dimensional \eqn{\delta} is estimated and
#' \eqn{\gamma = 0} is held fixed.
#' 
#' By default \eqn{\alpha = 1} is fixed for any \code{type}. If you want to
#' also estimate \eqn{\alpha} (for \code{type = "h"} or \code{"hh"}) 
#' set \code{theta.fixed = list()}.
#' 
#' @inheritParams common-arguments
#' @param y a numeric vector of real values.
#' @param theta.init a list containing the starting values of \eqn{(\alpha,
#'     \boldsymbol \beta, \gamma, \delta)} for the numerical optimization;
#'     default: see \code{\link{get_initial_theta}}.
#' @param hessian indicator for returning the (numerically obtained) Hessian at
#'     the optimum; default: \code{TRUE}. If the \pkg{numDeriv} package is
#'     available it uses \code{numDeriv::hessian()}; otherwise
#'     \code{stats::optim(..., hessian = TRUE)}.
#' @param theta.fixed a list of fixed parameters in the optimization; default
#'     only \code{alpha = 1}.
#' @param return.estimate.only logical; if \code{TRUE}, only a named flattened
#'     vector of \eqn{\widehat{\theta}_{MLE}} will be returned (only the
#'     estimated, non-fixed values). This is useful for simulations where it is
#'     usually not necessary to give a nicely organized output, but only the
#'     estimated parameter. Default: \code{FALSE}.
#' @param optim.fct character; which R optimization function should be
#'     used. Either \code{'optim'} (default), \code{'nlm'}, or \code{'solnp'}
#'     from the \pkg{Rsolnp} package (if available).  Note that if \code{'nlm'}
#'     is used, then \code{not.negative = TRUE} will be set automatically.
#' @param not.negative logical; if \code{TRUE}, it restricts \code{delta} or
#'     \code{gamma} to the non-negative reals. See \code{\link{theta2unbounded}}
#'     for details.
#' 
#' @return 
#' A list of class \code{LambertW_fit}:
#' \item{data}{ data \code{y},}
#' \item{loglik}{scalar; log-likelihood evaluated at the optimum
#'               \eqn{\widehat{\theta}_{MLE}},} 
#' \item{theta.init}{list; starting values for numerical optimization,} 
#' \item{beta}{ estimated \eqn{\boldsymbol \beta} vector of the input distribution via Lambert W MLE (In general this is not exactly identical to \eqn{\widehat{\boldsymbol \beta}_{MLE}} for the input data), } 
#' \item{theta}{list; MLE for \eqn{\theta}, } 
#' \item{type}{see Arguments,} 
#' \item{hessian}{Hessian matrix; used to calculate standard errors (only if \code{hessian = TRUE}, 
#' otherwise \code{NULL}),}
#' \item{call}{function call,} 
#' \item{distname}{see Arguments,} 
#' \item{message}{message from the optimization method. What kind of convergence?,} 
#' \item{method}{estimation method; here \code{"MLE"}.}
#' @keywords optimize
#' @export
#' @examples
#' 
#' # See ?LambertW-package
#' 
MLE_LambertW <- function(y, distname, type = c("h", "s", "hh"), 
                         theta.fixed = list(alpha = 1), 
                         use.mean.variance = TRUE,
                         theta.init = get_initial_theta(y, distname = distname,
                                                        type = type, 
                                                        theta.fixed = theta.fixed,
                                                        use.mean.variance =
                                                            use.mean.variance,
                                                        method = "IGMM"),
                         hessian = TRUE, return.estimate.only = FALSE,
                         optim.fct = c("optim", "nlm", "solnp"),
                         not.negative = FALSE) {
  check_distname(distname)
  type <- match.arg(type)
  optim.fct <- match.arg(optim.fct)
  
  if (optim.fct == "nlm") {
    not.negative <- TRUE
    warning("For optim.fct = 'nlm' not.negative was set to TRUE. ",
            "If you want to allow negative alpha/delta/gamma, ", "
            please use optim.fct = 'optim' instead (default).")
  }
  
  # initialize result ouput
  result <- list(data = y,
                 distname = distname,
                 type = type,
                 theta.fixed = theta.fixed,
                 call = match.call(),
                 method = "MLE")
  if (return.estimate.only) {
    hessian <- FALSE
  }

  if (length(theta.fixed) != 0) {
    # remove fixed values from initial parameters
    for (param.name in names(theta.fixed)) {
      theta.init[[param.name]] <- theta.fixed[[param.name]]
    }
  }
  if (type == "hh") {
    if (length(theta.init$alpha) == 1) {
      theta.init$alpha <- rep(theta.init$alpha, 2)
    }
    if (length(theta.init$delta) == 1) {
      theta.init$delta <- rep(theta.init$delta, 2)
    }
  }
  check_theta(theta.init, distname)
  result$theta.init <- theta.init
  if (not.negative) {
    # change 'optim' to 'nlm'
    # optim.fct <- "nlm"
    if (any(theta.init$delta < 0)) {
      theta.init$delta[theta.init$delta < 0] <- 0
    }
    if (any(theta.init$alpha < 0)) {
      theta.init$alpha[theta.init$alpha < 0] <- 0
    }
    
    theta.init$delta <- theta.init$delta + 0.01 # add eps to avoid log(0)
    theta.init$alpha <- theta.init$alpha + 0.01 # add eps to avoid log(0)
    theta.init <- theta2unbounded(theta.init, distname = distname, type = type)
  }
  if (length(theta.fixed) != 0) {
    # remove fixed values from initial parameters
    for (param.name in names(theta.fixed)) {
      theta.init[[param.name]] <- NULL
    }
  }
  params.init <- flatten_theta(theta.init)
  #####################################################################
  ### Define negative log-likelihood
  #####################################################################
  .neg_loglik_LambertW <- function(param, param.names = names(param)) {
    if (anyNA(param)) {
      neg.loglik <- 10^12
    } else {
      names(param) <- param.names
      theta.tmp <- unflatten_theta(param, distname = distname, type = type)
      # back transform unbounded theta back to original space so that
      # loglik_LambertW receives a valid theta vector
      if (not.negative) {
        theta.tmp <- theta2unbounded(theta.tmp, distname = distname,
                                     inverse = TRUE)
      }
      if (length(theta.fixed) > 0) {
        for (nn in names(theta.fixed)) {
          theta.tmp[[nn]] <- theta.fixed[[nn]]
        }
      }
      neg.loglik <- try(loglik_LambertW(theta.tmp, distname = distname,
                                        y = y, return.negative = TRUE,
                                        type = type, 
                                        flattened.theta.names = param.names,
                                        use.mean.variance = use.mean.variance),
                        silent = TRUE)
      if (inherits(neg.loglik, "try-error")) {
        warning("Parameter search lead to candidates that lead to errors",
                " in the log-likelihood calculation.  Please check the ",
                "estimate theta for feasibility.\n This can happen if",
                " estimates imply non-finite mean or standard deviation.")
        neg.loglik <- NA
      }
    }
    if (is.na(neg.loglik) || any(neg.loglik == c(-Inf, Inf))) { 
      neg.loglik <- 10^12
    }
    return(neg.loglik)
  }
  
  # if not.negative it uses unrestricted optimization; thus no bounds are
  # required
  if (not.negative) {
    lb <- NULL
    ub <- NULL
  } else {
    lower.upper.bounds <- get_theta_bounds(type = type, distname = distname, 
                                           beta = theta.init$beta, 
                                           not.negative = not.negative)
    lb <- lower.upper.bounds$lower + 1e-6
    ub <- lower.upper.bounds$upper - 1e-6
  }
  
  # remove fixed variables from initial parameter and upper/lower bounds
  keep.names <- setdiff(names(params.init), names(flatten_theta(theta.fixed)))
  if (length(keep.names) != length(params.init)) {
    params.init <- params.init[keep.names]
    if (!is.null(lb)) {
      lb <- lb[keep.names]
      ub <- ub[keep.names]
    }
  }
  
  if (any("df" == names(params.init)) && use.mean.variance) {
    lb["df"] <- 2.01
  }
  ########################################################################
  ### Do optimization to obtain theta.hat
  ########################################################################
  
  if (optim.fct == "optim") {
    if (not.negative) {
      fit <- optim(par = params.init,
                   fn = .neg_loglik_LambertW, 
                   param.names = names(params.init),
                   control = list(trace = 1),
                   hessian = FALSE)
    } else {
      fit <- optim(par = params.init,
                   fn = .neg_loglik_LambertW, 
                   param.names = names(params.init), 
                   lower = lb, upper = ub, 
                   control = list(trace = 0), method = "L-BFGS-B",
                   hessian = FALSE)
    }
    params.hat <- fit$par
    conv.message <- fit$message
  } else if (optim.fct == "nlm") {
    fit <- nlm(p = params.init,
               f = .neg_loglik_LambertW, param.names = names(params.init), 
               hessian = FALSE)
    # decode nlm messages
    nlm.code <-
        list("1" = paste0("Relative gradient is close to zero, current ",
                          "iterate is probably solution."),
             "2" = paste0("Successive iterates within tolerance, current",
                          "iterate is probably solution."),
             "3" = paste0("WARNING: Last global step failed to locate a ",
                          "point lower than estimate.\n Either approximate",
                          "local minimum or steptol is too small."),
             "4" = paste0("WARNING: Iteration limit exceeded. Most likely ",
                          "no the correct solution."),
             "5" = paste0("WARNING: Maximum step size stepmax exceeded five ",
                          "consecutive times.  Either the function is ",
                          "unbounded below, becomes asymptotic to a finite ",
                          "value from above in some direction or stepmax is",
                          "too small."))
    conv.message <- nlm.code[[as.character(fit$code)]]
    
    if (grepl("WARNING", conv.message)) {
      warning(conv.message)
    }
    params.hat <- fit$estimate
  } else if (optim.fct == "solnp") {
    requireNamespace("Rsolnp", quietly = TRUE)
    fit <- suppressWarnings(Rsolnp::solnp(params.init,
                                          fun = .neg_loglik_LambertW, 
                                          param.names = names(params.init),
                                          control = list(trace = 0)))
    params.hat <- fit$pars
    conv.message <- fit$elapsed
  } else {
    stop("Something went wrong with the optim.fct.")
  }
  # post-process estimated params hat
  names(params.hat) <- names(params.init)
  
  fit$theta <- unflatten_theta(params.hat, distname = distname, type = type)
  if (not.negative) {
    # transform theta back to original world
    fit$theta <- theta2unbounded(fit$theta, distname = distname, type = type,
                                 inverse = TRUE)
    # overwrite params.hat with original space coordinates (e.g., delta <0 will
    # become exp(delta) > 0)
    params.hat <- flatten_theta(fit$theta)
  }

  params.hat <- round(params.hat, 5)
  if (return.estimate.only) {
    return(params.hat)
  } 
  #############################################################################
  ### Estimate the Hessian
  #############################################################################
  if (hessian) {
    # make a hypercube around theta.hat
    if (requireNamespace("numDeriv", quietly = TRUE)) {
      result$hessian <- -numDeriv::hessian(func = .neg_loglik_LambertW,
                                           x = params.hat,
                                           param.names = names(params.hat))
    } else {
      params.0 <- params.hat
      ub <- params.0 + 1e-5
      lb <- params.0 - 1e-5
      opt.tmp <- optim(par = params.0,
                       fn = .neg_loglik_LambertW, 
                       param.names = names(params.0), 
                       lower = lb, upper = ub,
                       control = list(trace = 0), method = "L-BFGS-B",
                       hessian = TRUE)
      # return negative Hessian since it's the negativ log-likelihood
      result$hessian <- -opt.tmp$hessian  
    }
  } else {
    result$hessian <- NULL
  }

  if (length(theta.fixed) > 0) {
    for (nn in names(theta.fixed)) {
      fit$theta[[nn]] <- theta.fixed[[nn]]
    }
  }
  # make estimated parameters on bounded space
  tau.hat <- theta2tau(fit$theta, distname = distname, 
                       use.mean.variance = use.mean.variance)
  
  result <- c(result,
              list(loglik = loglik_LambertW(fit$theta, y = y, type = type,
                                            distname = distname,
                                            use.mean.variance = 
                                                use.mean.variance)$loglik.LambertW,
                   params.init = params.init,
                   params.hat = params.hat,
                   theta = fit$theta,
                   theta.fixed = theta.fixed,
                   tau = tau.hat,
                   fit = fit,
                   optim.fct = optim.fct,
                   message = conv.message,
                   use.mean.variance = use.mean.variance))
  class(result) <- "LambertW_fit"
  return(result)
} 

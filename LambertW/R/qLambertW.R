#' @rdname LambertW-utils
#' @inheritParams loglik-LambertW-utils
#' @export
qLambertW <- function(p, distname = NULL, theta = NULL, beta = NULL, gamma = 0, 
                      delta = 0, alpha = 1, 
                      input.u = NULL, tau = NULL, 
                      is.non.negative = FALSE,
                      use.mean.variance = TRUE) {
  
  stopifnot(xor(is.null(distname), is.null(input.u)))
  
  if (p < 0 || p > 1) {
    stop("Probability 'p' must lie in [0, 1].")
  }
  
  if (is.null(theta)) {
    warning("Please specify parameters by passing a list",
            "to the 'theta' argument directly.\n",
            "Specifying parameters by alpha, beta, gamma, delta will be",
            "deprecated.")
    theta <- list(beta = beta, alpha = alpha, gamma = gamma, delta = delta)
  } 
  theta <- complete_theta(theta)

  if (is.null(input.u)) {
    is.non.negative <- get_distname_family(distname)$is.non.negative
    check_theta(theta = theta, distname = distname)
    tau <- theta2tau(theta = theta, distname = distname, 
                     use.mean.variance = use.mean.variance)
    q.U <- function(p) qU(p, beta = theta$beta, distname = distname,
                          use.mean.variance = use.mean.variance)
    # support of the LambertW RV; is bounded for skewed Lambert W x F
    # distributions
    rv.support <- get_support(tau, is.non.negative = is.non.negative)
  } else {
    q.U <- input.u
    if (is.null(tau)) {
      stop("You must provide a 'tau' argument if 'input.u' is not NULL.")
    }
  }
    
  type.tmp <- tau2type(tau)
  
  if (type.tmp == "s") {
    # For skewed Lambert W x F distribution it depends whether the
    # family is non-negative or not.  If non-negative, the quantile fct
    # is in closed form; if not, it must be computed numerically.
    dist.family <- get_distname_family(distname)
    
    if (dist.family$is.non.negative) {
      .compute_quantile <- function(prob) {
        u.alpha <- q.U(prob)
        QQ.u <- H_gamma(u.alpha, gamma = theta$gamma)
        QQ <- normalize_by_tau(QQ.u, tau, inverse = TRUE)
        names(QQ) <- NULL
        return(QQ)
      }
    } else {
      #  of a family, the quantile function must be obtained
      # by matching the inverse.  This is slow, but in general no closed form is available
      # (at least not for location-scale family).
      .compute_quantile <- function(prob) {
        if (prob == 0) {
          QQ <- rv.support[1]
        } else if (prob == 1) {
          QQ <- rv.support[2]
        } else if (prob > 0 && prob < 1) {
          # define auxiliary objective function as distance to alpha-level
          # then minimize quantile so it matches this alpha-level
          aux.p <- function(y.a) {
            return(lp_norm(pLambertW(y.a, theta = theta, distname = distname,
                                     use.mean.variance = use.mean.variance) -
                             prob, 2))
          }
          # use default (-10, 10) as values for standard Gaussian
          S.10 <- normalize_by_tau(c(-10, 10), tau, inverse = TRUE)
          if (is.non.negative) {
            # use 10 as maximum (given standardized exp has pexp(10, 1))
            S.10 <- normalize_by_tau(c(0, 10), tau, inverse = TRUE)
          }
          
          intv <- c(max(S.10[1], rv.support[1]), min(S.10[2], rv.support[2]))
          fit <- suppressWarnings(optimize(aux.p, interval = intv))
          QQ <- fit$min
        }  # end of 'if (prob>0 && prob<1)'
        names(QQ) <- NULL
        return(QQ)
      }  # end of 'aux' function
    }
  } else if (type.tmp %in% c("hh", "h")) {
    # For a heavy-tailed Lambert W x  F distribution the quantile function can 
    # be obtained analytically.
    .compute_quantile <- function(prob) {
      u.alpha <- q.U(prob)
      if (type.tmp == "hh") {
        QQ.u <- G_2delta_2alpha(u.alpha, delta = theta$delta, alpha = theta$alpha)
      } else if (type.tmp == "h") {
        QQ.u <- G_delta_alpha(u.alpha, delta = theta$delta, alpha = theta$alpha)
      }
      QQ <- normalize_by_tau(QQ.u, tau, inverse = TRUE)
      names(QQ) <- NULL
      return(QQ)
    }
  }
  # apply aux function to all p values
  return(sapply(p, .compute_quantile))
}

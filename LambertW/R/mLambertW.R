#' @rdname LambertW-utils
#' @description
#' \code{mLambertW} computes the first 4 central/standardized moments of a Lambert W
#' \eqn{\times} F.  Works only for Gaussian distribution.
#' 
#' @export
mLambertW <- function(theta = NULL, distname = c("normal"), beta, gamma = 0, delta = 0, 
                      alpha = 1) {
  if (distname != "normal") {
    stop("Moments calculation for other than Normal distributions are not supported yet!")
  }
  check_distname(distname)
  if (is.null(theta)) {
    warning("Please specify parameters by passing a list",
            "to the 'theta' argument directly.\n",
            "Specifying parameters by alpha, beta, gamma, delta will be",
            "deprecated.")
    theta <- list(beta = beta, alpha = alpha, gamma = gamma, delta = delta)
  } 
  theta <- complete_theta(theta)
  
  if (length(theta$alpha) != 1) {
    stop("'alpha' must be of length 1; 'hh' distributions are not supported yet.")
  } 
  
  if (length(theta$delta) != 1) {
    stop("Parameter 'delta' must be of length 1; 'hh' distributions are not supported yet.")
  } 
  check_theta(theta = theta, distname = distname)
  tau <- theta2tau(theta = theta, distname = distname)
  type.tmp <- tau2type(tau)
  if (tau["gamma"] == 0 && all(theta$delta == 0)) {
    out <- list(mean = tau["mu_x"], sd = tau["sigma_x"], skewness = 0, kurtosis = 3)
  } else if (type.tmp == "s") {
    .mom_LambertW_U_Gauss <- function(gamma) {
      m1 <- gamma * exp((gamma^2)/2)
      m2 <- exp(2 * gamma^2) * (1 + 4 * gamma^2)
      cm2 <- exp(gamma^2) * (exp(gamma^2) * (1 + 4 * gamma^2) - gamma^2)
      m3 <- exp((9/2) * gamma^2) * (3 * gamma) * (9 * gamma^2 + 3)
      cm3 <- m3 - 3 * m1 * m2 + 2 * m1^3
      m4 <- 3 * exp(8 * gamma^2) + 96 * gamma^2 * exp(8 * gamma^2) + 256 * 
        gamma^4 * exp(8 * gamma^2)
      cm4 <- m4 - 4 * m1 * m3 + 6 * m1^2 * m2^2 - 3 * m1^4
      skew <- cm3/cm2^(3/2)
      kurt <- cm4/cm2^2
      return(list(mu.z = m1, sigma2_z = cm2, skew = skew, kurt = kurt))
    }
    mu.x <- tau["mu_x"]
    sigma.x <- tau["sigma_x"]
    mom.z <- .mom_LambertW_U_Gauss(tau["gamma"])
    mu.y <- mu.x + sigma.x * mom.z$mu.z
    sigma2_y <- sigma.x^2 * mom.z$sigma2_z #exp(gamma^2) * ((4 * gamma^2 + 1) * exp(gamma^2) - gamma^2)
    out <- list(mean = mu.y, 
                sd = sqrt(sigma2_y), 
                skewness = mom.z$skew, 
                kurtosis = mom.z$kurt)
  } else if (type.tmp %in% c("h", "hh")) {
    # initialize with NA and Inf (for delta > 1) this is true)
    out <- list(mean = NA, sd = Inf, skewness = NA, kurtosis = Inf)
    
    .moments <- function(n = 1, delta) {
      # odd moments are zero
      if (n %% 2 != 0) {
        return(0)
      } else {
        return(factorial(n) * (1 - n * delta)^(-(n + 1)/2)/(2^(n/2) * factorial(n/2)))
      }
    }
    
    .var_delta <- function(delta) {
      .moments(n = 2, delta = theta$delta)
    }
    
    kurt <- .moments(n = 4, delta = theta$delta)/.var_delta(theta$delta)^2
    if (theta$delta < 1) {
      out$mean <- tau["mu_x"]
    } 
    if (theta$delta < 1/2) {
      out$sd <- tau["sigma_x"] * sqrt(.var_delta(theta$delta))
    }
    if (theta$delta < 1/3) {
      out$skewness <- 0
    }
    if (theta$delta < 1/4) {
      out$kurtosis <- kurt
    }
  }

  names(out$mean) <- NULL
  names(out$sd) <- NULL
  return(out)
} 
#' @rdname tau-utils
#' @description
#' \code{get_initial_tau} provides starting estimates for \eqn{\tau}.
#' @param location.family logical; if \code{FALSE} it sets \code{mu_x} to 0 and only estimates 
#' \code{sigma_x}; if \code{TRUE} (default), it estimates \code{mu_x} as well.
#' @export
#' 

get_initial_tau <- function(y, type = c("h", "hh", "s"), location.family = TRUE) {
  
  type <- match.arg(type)
  
  if (location.family) {
    tau.init <- c(mu_x = median(y), sigma_x = mad(y), alpha = 1, gamma = 0, delta = 0)
  } else {
    tau.init <- c(mu_x = 0, sigma_x = mad(y), alpha = 1, gamma = 0, delta = 0)
  }
  z.init <- (y - tau.init["mu_x"]) / tau.init["sigma_x"]

  if (type == "s") {
    if (location.family) {
      tau.init["gamma"] <- gamma_Taylor(z.init, skewness.x = 0)
    } else {
      # use exponential distribution as baseline for a scale-family
      tau.init["gamma"] <- gamma_Taylor(z.init, skewness.x = 2)
    }
  } else if (type %in% c("h", "hh")) {
      tau.init["delta"] <- delta_Taylor(z.init)
      if (type == "hh") {
        # put a bit more weight to left (right) delta, if data is left (right) skewed
        if (skewness(z.init) > 0) {
          tau.init[c("delta_l", "delta_r")] <- tau.init["delta"] * c(1.1, 0.9)
        } else {
          tau.init[c("delta_l", "delta_r")] <- tau.init["delta"] * c(0.9, 1.1)
        }
        tau.init[c("alpha_l", "alpha_r")] <- c(1, 1)
        # remove 'delta' and 'alpha' from tau
        tau.init <- tau.init[setdiff(names(tau.init), c("delta", "alpha"))]
      }
  }
  x.input <- get_input(y, tau.init)
  # update parameters (and remove 'NA' if exist; happens if 'gamma' is too extreme for the given data)
  if (location.family) {
    tau.init[c("mu_x", "sigma_x")] <- c(mean(x.input, na.rm = TRUE), 
                                        sd(x.input, na.rm = TRUE))
  } else {
    tau.init["sigma_x"] <- sd(x.input, na.rm = TRUE)
  }
  return(tau.init)
}
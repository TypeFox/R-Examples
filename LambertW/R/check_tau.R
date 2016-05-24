#' @rdname tau-utils
#' @description
#' \code{check_tau} checks if \eqn{\tau} is correctly specified (correct names, non-negativity
#' constraints, etc.)
#' @export
check_tau <- function(tau)  {
  
  stopifnot('mu_x' %in% names(tau),
            "sigma_x" %in% names(tau))
  # use the beta check for Normality to test tau mu_x and sigma_x
  check_beta(tau[c("mu_x", "sigma_x")], distname = "normal")
  
  if ("gamma" %in% names(tau)) {
    if (length(tau["gamma"]) != 1) {
      print(tau['gamma'])
      stop("'gamma' must have length 1.")
    }
  } else {
    if (!any(grepl("^[delta]|^[delta_]", names(tau)))) {
      stop("'tau' must either contain 'gamma' or 'delta' (or 'delta_l' and 'delta_r').")
    }
  }
  
  if (any(grepl("^[delta]|^[delta_]", names(tau)))) {
    if (!(any(grepl("^[alpha]|^[alpha_]", names(tau))))) {
      stop("If 'tau' contains 'delta' it must also contain an alpha value.")
    }
  }

  alpha.values <- tau[grepl("alpha", names(tau))]
  delta.values <- tau[grepl("delta", names(tau))]
  stopifnot(length(alpha.values) <= 2,
            length(delta.values) <= 2,
            alpha.values > 0)

  if ("gamma" %in% names(tau)) {
    if (any(delta.values != 0) && tau["gamma"] != 0) {
      print(paste("'gamma':", tau["gamma"]))
      print(paste("'delta'", delta.values))
      stop("Both 'gamma' and 'delta' values can't be non-zero.")
    }
  }
} 

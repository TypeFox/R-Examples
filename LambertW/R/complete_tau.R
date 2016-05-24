#' @rdname tau-utils
#' @description
#' 
#' \code{complete_tau} completes missing values so users don't have to specify
#'     every element of \eqn{\tau} explicitly. \code{'mu_x'} and
#'     \code{'sigma_x'} must be specified, but \code{alpha = 1}, \code{gamma =
#'     0}, and \code{delta = 0} will be set automatically if missing.
#' 
#' @export
complete_tau <- function(tau, type = tau2type(tau)) {
  
  stopifnot("mu_x" %in% names(tau),
            "sigma_x" %in% names(tau))
  
  # check that the type seems right
  if (type == "hh" && ("delta" %in% names(tau))) {
    warning("Changed tau to type 'h', since 'delta' was found in tau.")
    type <- "h"
  } else if (type == 'h' && (any(grepl("delta_", names(tau))))) {
    warning("Changed tau to type 'hh', since 'delta_l' (or 'delta_r') was found in tau.")
    type <- "hh"
  }
  
  if (is.na(tau["gamma"])) {
      tau["gamma"] <- 0
  } 
  if (type == 'h') {
    if (is.na(tau["alpha"])) {
      tau["alpha"] <- 1
    }
    if (is.na(tau["delta"])) {
      tau["delta"] <- 0
    }
  } else if (type == "hh") {
    if (is.na(tau["delta_l"])) {
        tau["delta_l"] <- 0
    }
    if (is.na(tau["delta_r"])) {
      tau["delta_r"] <- 0
    }
    if (is.na(tau["alpha_l"])) {
      tau["alpha_l"] <- 1
    }
    if (is.na(tau["alpha_r"])) {
      tau["alpha_r"] <- 1
    }
    #  remove 'delta' and 'alpha'
    tau <- tau[setdiff(names(tau), "alpha")]
    tau <- tau[setdiff(names(tau), "delta")]
  }
  return(tau)
} 

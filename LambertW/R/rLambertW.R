#' @rdname LambertW-utils
#' @export
rLambertW <- function(n, distname, theta = NULL, beta = NULL, gamma = 0, delta = 0, alpha = 1, 
                      return.x = FALSE, input.u = NULL, tau = NULL,
                      use.mean.variance = TRUE) {
  
  stopifnot(n > 0)
  
  if (is.null(theta)) {
    warning("Please specify parameters by passing a list",
            "to the 'theta' argument directly.\n",
            "Specifying parameters by alpha, beta, gamma, delta will be",
            "deprecated.")
    theta <- list(beta = beta, alpha = alpha, gamma = gamma, delta = delta)
  } 
  theta <- complete_theta(theta)
  
  if (is.null(input.u)) {
    check_distname(distname)
    check_theta(theta = theta, distname = distname)
    tau <- theta2tau(theta = theta, distname = distname, 
                     use.mean.variance = use.mean.variance)
    uu <- rU(n = n, beta = theta$beta, distname = distname,
             use.mean.variance = use.mean.variance)
  } else {
    if (is.numeric(input.u)) {
      uu <- input.u
    } else if (is.function(input.u)) {
      uu <- input.u(n = n, beta = theta$beta)
    } else {
      stop("'input.u' must be either numeric (simulated U values); or a function 
           that draws a random sample from U.")
    }
    if (is.null(tau)) {
      stop("You must provide a 'tau' argument if 'input.u' is not NULL.")
    }
  }
  check_tau(tau)
  
  type.tmp <- tau2type(tau)
  if (all(tau[grepl("delta", names(tau))] == 0) && all(tau[grepl("alpha", names(tau))] == 1) && 
        tau["gamma"] == 0) {
    zz <- uu
  } else if (type.tmp == "h") {
      zz <- G_delta_alpha(uu, delta = tau["delta"], alpha = tau["alpha"])
  } else if (type.tmp == "hh") {
      zz <- G_2delta_2alpha(uu, delta = tau[grepl("delta", names(tau))], 
                            alpha = tau[grepl("alpha", names(tau))])
  } else if (type.tmp == "s") {
    zz <- H_gamma(uu, tau["gamma"])
  } else {
    stop("Something went wrong with the type of the distribution.")
  }
  
  yy <- normalize_by_tau(zz, tau, inverse = TRUE)
  names(yy) <- NULL
  if (return.x) {
    xx <- normalize_by_tau(uu, tau)
    return(list(x = xx, y = yy)) 
  } else {
    return(yy)
  }
} 
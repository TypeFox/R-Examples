#' @rdname LambertW-utils
#' @export
dLambertW <- function(y, distname = NULL, theta = NULL, beta = NULL, gamma = 0, 
                      delta = 0, alpha = 1, input.u = NULL, tau = NULL, 
                      use.mean.variance = TRUE, log = FALSE) {
  
  if (is.null(theta)) {
    warning("Please specify parameters by passing a list",
            "to the 'theta' argument directly.\n",
            "Specifying parameters by alpha, beta, gamma, delta will be",
            "deprecated.")
    theta <- list(beta = beta, alpha = alpha, gamma = gamma, delta = delta)
  } 
  theta <- complete_theta(theta)
  
  if (is.null(input.u)) {
    check_theta(theta = theta, distname = distname)
    tau <- theta2tau(theta = theta, distname = distname,
                     use.mean.variance = use.mean.variance)
    fU <- function(u) dU(u, beta = theta$beta, distname = distname,
                         use.mean.variance = use.mean.variance) 
  } else {
    fU <- input.u
    if (is.null(tau)) {
      stop("You must provide a 'tau' argument if 'input.u' is not NULL.")
    }
  }
  
  fX <- function(x) fU((x - tau["mu_x"])/tau["sigma_x"])/tau["sigma_x"]
  
  type.tmp <- tau2type(tau)
  if (all(tau[grepl("delta", names(tau))] == 0) && 
        all(tau[grepl("alpha", names(tau))] == 1) && 
          tau["gamma"] == 0) {
    gg <- fX(y)
  } else {
    gg <- rep(NA, length(y))
    zz <- normalize_by_tau(y, tau)
    names(zz) <- NULL
    
    ## the heavy-tail version (if delta != 0)
    if (type.tmp == "h") {
      uu <- W_delta_alpha(zz, delta = tau["delta"], alpha = tau["alpha"])
      # g = fU(u) * deriv_W_delta(z, delta=delta) * tau["sigma_x"] g = 1/tau["sigma_x"] * sign(z) *
      # fU(sign(z) * u) * u / (z*(1 + delta*u^2)) #deriv_W_delta(z, delta=delta) *
      # tau["sigma_x"]
      gg <- 1/tau["sigma_x"] * fU(uu) * deriv_W_delta_alpha(zz, delta = tau["delta"], 
                                                            alpha = tau["alpha"])
    } else if (type.tmp == "hh") {
      ind.pos <- (zz > 0)
      theta.l <- list(beta = theta$beta, gamma = 0, 
                      delta = tau["delta_l"], alpha = tau["alpha_l"])
      theta.r <- list(beta = theta$beta, gamma = 0, 
                      delta = tau["delta_r"], alpha = tau["alpha_r"])
      gg[!ind.pos] <- dLambertW(y[!ind.pos], theta = theta.l, distname = distname)
      gg[ind.pos] <- dLambertW(y[ind.pos], theta = theta.r, distname = distname)
    } else if (type.tmp == "s") {
      if (tau["gamma"] < 0) {
        # revert data and signs of gamma and mu_x so that we only have to implement the gamma > 0 case
        y <- -y
        tau["gamma"] <- -tau["gamma"]
        tau["mu_x"] <- -tau["mu_x"]
        zz <- normalize_by_tau(y, tau)
        names(zz) <- NULL
      }
      
      # principal branch for all values (negative and positive)
      r.0 <- W_gamma(zz, gamma = tau["gamma"], branch = 0)
      x.0 <- normalize_by_tau(r.0, tau, inverse = TRUE)
      # for derivative it's good to use 
      # gamma * z = r.0 * gamma (r.0 = W(z * gamma) / gamma)
      g.0 <- fX(x.0) * deriv_W(tau["gamma"] * zz, W.z = r.0 * tau["gamma"])
      gg[zz >= 0] <- g.0[zz >= 0]
      
      zz.neg <- zz[zz < 0]
      if (length(zz.neg) > 0) {
        r.neg.1 <- W_gamma(zz.neg, gamma = tau["gamma"], branch = -1)
        x.neg.1 <- normalize_by_tau(r.neg.1, tau, inverse = TRUE)
        g.neg <- g.0[zz < 0] - 
          fX(x.neg.1) * deriv_W(tau["gamma"] * zz.neg, 
                                branch = -1,
                                W.z = r.neg.1 * tau["gamma"])
        gg[zz < 0] <- g.neg
      }
    }
  } # end of else 
  names(gg) <- NULL
  if (log) {
    return(log(gg))
  } else {
    return(gg)
  } 
}

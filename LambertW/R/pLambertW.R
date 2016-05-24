#' @rdname LambertW-utils
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#' \eqn{P(X \leq x)} otherwise, \eqn{P(X > x)}.
#' @export
pLambertW <- function(q, distname, theta = NULL, beta = NULL, gamma = 0, delta = 0, alpha = 1, 
                      input.u = NULL, tau = NULL, log = FALSE,
                      lower.tail = FALSE, use.mean.variance = TRUE) {
  
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
    FU <- function(u) pU(u, beta = theta$beta, distname = distname,
                         use.mean.variance = use.mean.variance) 
  } else {
    FU <- input.u
    if (is.null(tau)) {
      stop("You must provide a 'tau' argument if 'input.u' is not NULL.")
    }
  }
  y <- q
  
  FX <- function(x) {
    return(FU((x - tau["mu_x"])/tau["sigma_x"]))
  }
  
  type.tmp <- tau2type(tau)
  if (all(tau[grepl("delta", names(tau))] == 0) && all(tau[grepl("alpha", names(tau))] == 1) && 
        tau["gamma"] == 0) {
    GG <- FX(y)
  } else {
    # begin of else for non-trivial calculation of pLambertW
    zz <- normalize_by_tau(y, tau)
    names(zz) <- NULL
    ## the heavy-tail version (if theta$delta != 0)
    if (type.tmp == "h") {
      u <- W_delta_alpha(zz, delta = tau["delta"], alpha = tau["alpha"])
      GG <- FU(u)
    } else if (type.tmp == "hh") {
      ind.pos <- (zz > 0)
      theta.l <- list(beta = theta$beta, gamma = 0, 
                      delta = tau["delta_l"], alpha = tau["alpha_l"])
      theta.r <- list(beta = theta$beta, gamma = 0, 
                      delta = tau["delta_r"], alpha = tau["alpha_r"])
      GG[!ind.pos] <- pLambertW(y[!ind.pos], theta = theta.r, distname = distname,
                                use.mean.variance = use.mean.variance)
      GG[ind.pos] <- pLambertW(y[ind.pos], theta = theta.l, distname = distname,
                               use.mean.variance = use.mean.variance)

    } else if (type.tmp == "s") {
      gamma.negative <- FALSE
      GG <- rep(NA, length(y))
      if (tau["gamma"] < 0) {
        y <- -y
        tau["gamma"] <- -tau["gamma"]
        tau["mu_x"] <- -tau["mu_x"]
        gamma.negative <- TRUE
        zz <- normalize_by_tau(y, tau)
        names(zz) <- NULL
      }
      # principal branch for all values (negative and positive)
      r.0 <- W_gamma(zz, gamma = tau["gamma"], branch = 0)
      x.0 <- normalize_by_tau(r.0, tau, inverse = TRUE)
      # for derivative it's good to use 
      # gamma * z = r.0 * gamma (r.0 = W(z * gamma) / gamma)
      G.0 <- FX(x.0)
      GG[zz >= 0] <- G.0[zz >= 0]
      
      zz.neg <- zz[zz < 0]
      if (length(zz.neg) > 0) {
        r.neg.1 <- W_gamma(zz.neg, gamma = tau["gamma"], branch = -1)
        x.neg.1 <- normalize_by_tau(r.neg.1, tau, inverse = TRUE)
        G.neg <- G.0[zz < 0] - FX(x.neg.1)
        GG[zz < 0] <- G.neg
      }
      if (gamma.negative) {
        GG <- 1 - GG
      }
    }
  }  # end of else 
  names(GG) <- NULL
  
  if (lower.tail) {
    GG <- 1 - GG
  }
  
  if (log) {
    return(log(GG))
  } else {
    return(GG)
  }
} 
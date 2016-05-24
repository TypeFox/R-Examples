#' Prediction Intervals for ARIMA Processes with Exogenous Variables Using Importance Sampling
#'
#' Function \code{arima_pi} computes prediction intervals for ARIMA processes
#' with exogenous variables using importance sampling. For regression coefficients,
#' diffuse (uninformative) prior is used, whereas multiple options for
#' prior distributions for ARMA coefficients are supported.
#'
#' @export
#' @name arima_pi
#
#' @param x vector containing the time series
#' @param xreg matrix or data frame containing the exogenous variables
#' (not including the intercept which is always included for non-differenced series)
#' @param order vector of length 3 with values p,d,q
#' corresponding to the number of AR parameters, degree of differencing and number of MA parameters.
#' @param nsim number of simulations used in importance sampling. Default is 1000.
#' @param n_ahead length of the forecast horizon.
#' @param level desired frequentist coverage probability of the prediction intervals.
#' @param median compute the median of the prediction interval.
#' @param se_limits compute the standard errors of the prediction interval limits.
#' @param prior prior to be used in importance sampling for AR and MA parameters.
#' Defaults to uniform prior. Several Jeffreys' priors are also available (see \code{\link{jeffreys}} for details).
#' If "custom", a user-defined custom prior is used (see next arguments).
#' All priors assume that the ARMA parameters lie in stationarity/invertibility region.
#' @param custom_prior function for computing custom prior.
#' First argument must be a vector containing the AR and MA parameters (in that order).
#' @param custom_prior_args list containing additional arguments to \code{custom_prior}.
#' @param invertibility Logical, should the priors include invertibility constraint? Default is \code{FALSE}.
#' @param last_only compute the prediction intervals only for the last prediction step.
#' @param return_weights Return (scaled) weights used in importance sampling.
#' @param ... additional arguments for \code{\link{arima}}.
#' @seealso \code{\link{tsPI}}, \code{\link{struct_pi}}
#' @return a list containing the prediction intervals.
#'  @references
#' \enumerate{
#' \item{Helske, J. and Nyblom, J. (2015). Improved frequentist prediction
#' intervals for autoregressive models by simulation.
#' In Siem Jan Koopman and Neil Shephard, editors,
#' Unobserved Components and Time Series Econometrics. Oxford University Press.
#' \url{http://urn.fi/URN:NBN:fi:jyu-201603141839}}
#'  \item{Helske, J. and Nyblom, J. (2014). Improved frequentist prediction intervals for
#'  ARMA models by simulation.
#'  In Johan Knif and Bernd Pape, editors,
#' Contributions to Mathematics, Statistics, Econometrics, and Finance:
#' essays in honour of professor Seppo Pynnönen,
#' number 296 in Acta Wasaensia, pages 71--86. University of Vaasa.
#' \url{http://urn.fi/URN:NBN:fi:jyu-201603141836}}
#' }
#' @examples
#'
#' set.seed(123)
#' x <- arima.sim(n = 30, model = list(ar = 0.9))
#'
#' pred_arima <- predict(arima(x, order = c(1,0,0)), n.ahead = 10, se.fit = TRUE)
#' pred_arima <- cbind(pred = pred_arima$pred,
#'   lwr = pred_arima$pred - qnorm(0.975)*pred_arima$se,
#'   upr = pred_arima$pred + qnorm(0.975)*pred_arima$se)
#'
#' pred <- arima_pi(x, order = c(1,0,0), n_ahead = 10)
#'
#' ts.plot(ts.union(x,pred_arima, pred[,1:3]), col = c(1,2,2,2,3,3,3),
#'   lty = c(1,1,2,2,1,2,2))
#'
arima_pi <- function(x, order, xreg = NULL, n_ahead = 1, level = 0.95, median = TRUE, se_limits = TRUE,
  prior = "uniform", custom_prior, custom_prior_args  = NULL, nsim = 1000, invertibility = FALSE, last_only = FALSE,
  return_weights = FALSE, ...){

  distfkt <- function(a, prob, ex, sdx, w){
    sum(w * pnorm(q = a, mean = ex, sd = sdx)) - prob
  }

  if (level >= 1 | level < 0)
    stop("Invalid value of argument 'level'.")
  prior <- match.arg(prior,
    c("uniform", "approx_joint_jeffreys", "approx_marginal_jeffreys",
      "exact_joint_jeffreys", "exact_marginal_jeffreys", "custom"))
  if (prior == "custom" && missing(custom_prior))
    stop("Missing custom prior.")

  n <- length(x)
  fit <- arima(x = x, order = order, xreg = xreg[1:n,, drop = FALSE], ...)
  if (fit$code != 0 || sum(diag(fit$var.coef) < 1e-7) != 0)
  {
    stop("arima function returned non-convergence or coefficient variances smaller than 1e-7.")
  }

  p <- order[1]
  d <- order[2]
  q <- order[3]
  npar <- p + q
  psihat <- as.numeric(fit$coef[1:npar])
  psivar <- matrix(fit$var.coef[1:npar, 1:npar], npar, npar)
  psivarchol <- try(t(chol(psivar)), TRUE)
  if (inherits(psivarchol, "try-error") || any(diag(psivarchol) < 1e-12))
    stop("Covariance matrix fit$var.coef[1:npar, 1:npar] obtained from arima function is not positive definite.")
  psivarinv <- crossprod(solve(psivarchol))

  psisim <- array(rnorm(n = nsim * npar),c(nsim, npar))
  for (i in 1:nsim)
    psisim[i,] <- psihat + psivarchol %*% psisim[i, ]

  sigma2hat <- fit$sigma2
  k <- length(fit$coef) - npar
  sigmasim <- 1/rchisq(nsim, df = n - k)

  #stationarity
  if (p > 0) {
    w <- apply(matrix(abs(apply(cbind(rep(1, nsim),-psisim[, 1:p]), 1, polyroot)) > 1, p, nsim), 2, sum) == p
  } else w <- rep(1, nsim)
  #invertibility?
  if (q > 0 && invertibility) {
    w <- w * apply(matrix(abs(apply(cbind(rep(1, nsim), psisim[, (p + 1):(p + q)]),1, polyroot)) > 1, q, nsim), 2, sum) == q
  }

  x <- window(x, end = end(x) + c(0, n_ahead), extend = TRUE)
  ex <- sdx <- matrix(0,n_ahead, nsim)


  valid <- which(w > 0)[1]
  valid_ar <- if (p > 0) psisim[valid, 1:p] else NULL
  valid_ma <- if ( q > 0) psisim[valid, (p + 1):(p + q)] else NULL
  if (is.null(xreg)) {
    model <- SSModel(x ~ SSMarima(ar = valid_ar, ma = valid_ma, d = order[2], Q = 1), H = 0)
  } else {
    model <- SSModel(x ~ xreg + SSMarima(ar = valid_ar, ma = valid_ma, d = order[2], Q = 1), H = 0)
  }
  kd <- k + d
  m <- max(p, q + 1)
  nd <- (kd + 1):(kd + m)

  for (i in which(w)) {
    if (q > 0)
      model$R[(kd + 2):(kd + 1 + q)] <- psisim[i, (p + 1):npar]
    if (p > 0)
      model$T[(kd + 1):(kd + p),kd + 1,] = psisim[i, 1:p]
    model$P1[nd, nd] <- solve(a = diag(m^2) - as.matrix(kronecker(model$T[nd, nd, ], model$T[nd, nd, ])),
      b = c(model$R[nd] %*% t(model$R[nd])))

    out <- KFS(model, filtering = "mean", smoothing = "none")
    if (out$d > 0) {
      s2 <- sum(out$v[1:out$d][out$Finf == 0]^2/out$F[1:out$d][out$Finf == 0]) +
        sum(out$v[(out$d + 1):n]^2/out$F[(out$d + 1):n])
    } else s2 <- sum(c(out$v[1:n])^2/out$F[1:n])

    sigmasim[i] <-  sqrt(s2 * sigmasim[i])
    ex[1:n_ahead,i] <- out$m[(n + 1):(n + n_ahead)]
    sdx[1:n_ahead,i] <- sqrt(out$P_mu[(n + 1):(n + n_ahead)]) * sigmasim[i]

    #s2 is replaced by (s2/n)/sigma2hat so that w>>0
    detVXinvX <- prod(c(out$Finf[out$Finf > 0], out$F[1:out$d][out$Finf == 0], out$F[(out$d + 1):n]))
    weight <-
      (((s2 / n) / sigma2hat)^(-0.5 * (n - k)) / sqrt(detVXinvX)) /
      exp(-0.5 * t(psisim[i,] - psihat) %*% psivarinv %*% (psisim[i, ] - psihat))

    w[i] <- weight * switch (prior,
      uniform = 1,
      exact_joint_jeffreys =  exact_joint_jeffreys(psisim[i, ], xreg, p, q, n),
      exact_marginal_jeffreys = exact_marginal_jeffreys(psisim[i, ], p, q, n),
      approx_joint_jeffreys = approx_joint_jeffreys(psisim[i, ], xreg, p, q, n),
      approx_marginal_jeffreys = approx_marginal_jeffreys(psisim[i, ], p, q),
      custom = do.call(custom_prior, list(psisim[i, ], custom_prior_args)))
  }
  w <- w/sum(w)
  out <- ts(matrix(NA, n_ahead, 2 + median + 2 * se_limits), end = end(model$y), frequency = frequency(model$y))
  colnames(out) <- c(if (median) "median", "lwr", "upr", if (se_limits) c("se_lwr", "se_upr"))

  if (sum(is.na(w)) == 0) {
    nz_w <- w[w != 0]
    for (j in 1:n_ahead) {
      nz_ex <- ex[j, w != 0]
      nz_sdx <- sdx[j, w != 0]
      interval <- c(mean(nz_ex) + c(-1, 1) * 8 * max(nz_sdx))
      out[j, "lwr"] <- uniroot(distfkt, interval = interval, prob = (1 - level) / 2,
        ex = nz_ex, sdx = nz_sdx,w = nz_w, tol = 1e-12)$root
      out[j, "upr"] <- uniroot(distfkt,  interval = interval, prob = 1 - (1 - level) / 2,
        ex = nz_ex, sdx = nz_sdx,w = nz_w, tol = 1e-12)$root
      if (median) {
        out[j, "median"] <- uniroot(distfkt,  interval = interval, prob = 0.5,
          ex = nz_ex, sdx = nz_sdx,w = nz_w, tol = 1e-12)$root
      }
      if (se_limits) {
        out[j, "se_lwr"] <-
          sqrt(sum((nz_w * ((1 - level) / 2 - pnorm(q = out[j, "lwr"], nz_ex, nz_sdx)))^2) / (nsim - 1)) /
          (sum(nz_w * dnorm(x = out[j, "lwr"], nz_ex, nz_sdx)/sqrt(nsim)))
        out[j, "se_upr"] <-
          sqrt(sum((nz_w * (1 - (1 - level) / 2 - pnorm(q = out[j, "upr"], nz_ex, nz_sdx)))^2) / (nsim - 1)) /
          (sum(nz_w * dnorm(x = out[j, "upr"], nz_ex, nz_sdx)/sqrt(nsim)))
      }
    }
  } else stop("NA values in weights.")


  if (return_weights)
    out <- list(pred = out, weights = w)
  out
}

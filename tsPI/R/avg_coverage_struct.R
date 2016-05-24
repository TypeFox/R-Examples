#' Compute the average coverage of the prediction intervals computed by \code{\link{struct_pi}}
#' and plug-in method
#'
#' Computes expected coverage probabilities of the prediction intervals of structural time series model.
#' Note that for the plug-in method only standard deviations are assumed to be identical to their estimates,
#' but the initial values for the states are still treated as diffuse. Because of this,
#' plug-in method often performs relatively well in case of structural time series models
#' compared to similar type of ARIMA models
#' (local level and local linear trend models are closely related to ARIMA(0,1,1) and ARIMA(0,2,2) models),
#' and in some cases even outperforms the importance sampling approach with uniform prior (see examples).
#' This is not suprising, as local level and local linear trend models are closely related to ARIMA(0,1,1) and ARIMA(0,2,2) models,
#' and the effect of uncertainty in MA components is not as significant as the uncertainty of AR components
#'
#' @export
#' @seealso \code{\link{struct_pi}}.
#
#' @param type Type of model. See \code{\link{struct_pi}}.
#' @param sds vector containing the standard deviations of the model (observation error, level, slope, and seasonal).
#' @param frequency frequency of the series, needed for seasonal component.
#' @param n length of the time series
#' @param n_ahead length of the forecast horizon
#' @param nsim number of simulations used in importance sampling
#' @param nsim2 number of simulations used in computing the expected coverage
#' @param level desired coverage probability of the prediction intervals
#' @param prior prior to be used in importance sampling.
#' @param return_all_coverages return raw results i.e. coverages for each simulations. When \code{FALSE} (default), summary statistics are returned.
#' @param ... additional arguments to \code{\link{struct_pi}}.
#' @return a list containing the coverage probabilities
#' @examples
#' \dontrun{
#' set.seed(123)
#' # takes a while, notice se, increase nsim2 to get more accurate results
#' avg_coverage_struct(type = "level", sds = c(1, 0.1), n = 50, n_ahead = 10, nsim2 = 100)
#' avg_coverage_struct(type = "BSM", sds = c(1, 1, 0.1, 10),
#'  frequency = 4, n = 50, n_ahead = 10, nsim2 = 100)
#' }
avg_coverage_struct <- function(type = c("level", "trend", "BSM"), sds, frequency = 1, n, n_ahead = 1,
  nsim2, nsim = 100, level = 0.95, prior = "uniform", return_all_coverages = FALSE, ...){



  prior <- match.arg(prior, c("uniform", "custom"))

  type <- match.arg(type)

  npar <- switch(type,
    level = 2,
    trend = 3,
    BSM = 4)
  if (npar != length(sds))
    stop("Incorrect number of standard deviations provided.")

  model <- true_model <- switch(type,
    level = SSModel(rep(NA,n) ~ SSMtrend(1, NA), H = NA),
    trend = SSModel(rep(NA,n) ~ SSMtrend(2, list(NA, NA)), H = NA),
    BSM = SSModel(ts(rep(NA,n), frequency = frequency)  ~ SSMtrend(2, list(NA, NA)) + SSMseasonal(frequency, Q = NA), H = NA))
  model$y[1] <- 0



  dx <- 1 + 0:(npar - 2) * npar

  likfn <- function(pars, model){
    # parameters are log(standard deviation)
    model$Q[dx] <- exp(2 * pars[-1])
    model$H[1] <- exp(2 * pars[1])
    - logLik(model)
  }

  true_model$Q[dx] <- exp(2 * log(sds[-1]))
  true_model$H[1] <- exp(2 * log(sds[1]))

  covprobs <- array(0, c(n_ahead, nsim2, 2))
  dimnames(covprobs)[[3]] <- c("plug-in", prior)

  count <- count2 <- 0
  for (i in 1:nsim2) {

    true_model$y[-1] <- NA
    model$y[] <- true_model$y[] <- simulateSSM(true_model, "obs")
    # use true parameters as initial values, in real application multiple initial values would be used
    fit <- optim(fn = likfn, par = log(sds), method = "BFGS", model = model, hessian = TRUE, control = list(reltol= 1e-10))
    if (inherits(try(t(chol(solve(fit$hessian))), TRUE), "try-error")){
      count <- count + 1
      covprobs[, i, 2] <- covprobs[, i, 1] <- NA
      next
    }
    model$Q[dx] <- exp(2 * fit$par[-1])
    model$H[1] <- exp(2 * fit$par[1])

    pred <- predict(model, n.ahead = n_ahead, level = level, interval = "prediction")

    true_pred <- predict(true_model, n.ahead = n_ahead, se.fit = TRUE)


    ipi <- try(struct_pi(model$y, type = type, inits = fit$par, nsim = nsim, level = level, n_ahead = n_ahead,
      prior = prior, median = FALSE, se_limits = FALSE, ...), TRUE)
    if (!inherits(ipi, "try-error")) {
      covprobs[, i, 2] <- pnorm(q = ipi[, "upr"], mean = true_pred[, "fit"], sd = sqrt(true_pred[,"se.fit"]^2 + sds[1]^2)) -
        pnorm(q = ipi[, "lwr"], mean = true_pred[, "fit"], sd = sqrt(true_pred[,"se.fit"]^2 + sds[1]^2))
    } else {
      count2 <- count2 + 1
    }
    covprobs[, i, 1] <- pnorm(q = pred[, 3],mean = true_pred[, "fit"], sd = sqrt(true_pred[,"se.fit"]^2 + sds[1]^2)) -
      pnorm(q = pred[, 2], mean = true_pred[, "fit"], sd = sqrt(true_pred[,"se.fit"]^2 + sds[1]^2))

  }
  if(count > 0)
    warning(paste0("There were ",count, " cases where the simulateSSM generated a series which caused the estimation of the model to fail. These cases were set as NA when computing coverage probabilities. " ))
  if(count2 > 0)
    warning(paste0("There were ",count2, " cases where the the importance sampling method generated error. The coverage probability for these cases were set to zero when computing coverage probabilities. " ))
  if(!return_all_coverages){
    out <- vector("list", 2)
    names(out) <- c("plugin", prior)
    for(i in 1:2){
      out[[i]] <- matrix(0, n_ahead, 8)
      colnames(out[[i]]) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.","Standard error", "failures")
      rownames(out[[i]]) <- paste("n_ahead =", 1:n_ahead)
      for(j in 1:n_ahead)
        out[[i]][j,] <- c(summary(na.exclude(covprobs[j,,i])), sd(covprobs[j,,i], na.rm = TRUE)/sqrt(nsim2-count), sum(is.na(covprobs[j,,i]))+(i==2)*count2)
    }
    return(out)
  } else return(covprobs)
}

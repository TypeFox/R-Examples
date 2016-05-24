#' Create an environmental time series.
#'
#' Generate various types of environmental time series.
#'
#' @param n_t Length of time series.
#' @param type Type of time series to produce.
#' @param sine_params Parameters controlling sine wave time series.
#' @param arma_params Parameters controlling ARMA time series.
#' @param regime_params Parameters controlling regime-shift time series.
#' @param linear_params Parameters controlling warming or cooling time series.
#'   Minimum environmental value, maximum environmental value, environmental
#'   standard deviation, and the year to start the linear trend (useful if
#'   you're going to throw out the early years as burn in).
#' @param linear_arma_params A combination of \code{arma_params} and 
#'   \code{linear_params}.
#' @param constant_params Parameter controlling constant time series.
#' @export
#' @examples
#' types <- c("sine", "arma", "regime", "linear", "linear_arma", "constant")
#' x <- list()
#' for(i in 1:6) x[[i]] <- generate_env_ts(n_t = 100, type = types[i])
#' op <- par(mfrow = c(5, 1), mar = c(3,3,1,0), cex = 0.7)
#' for(i in 1:6) plot(x[[i]], type = "o", main = types[i])
#' par(op)

generate_env_ts <- function(
  n_t,
  type = c("sine", "arma", "regime", "linear", "linear_arma", "constant"),
  sine_params = list(amplitude = 1, ang_frequency = 0.2, phase = 0,
    mean_value = 0, slope = 0, sigma_env = 0.02),
  arma_params = list(mean_value = 0, sigma_env = 0.50, ar = 0.4, ma = 0),
  regime_params = list(break_pts = c(25, 75), break_vals = c(-1, 0, 1)),
  linear_params = list(min_value = -1, max_value = 1, sigma_env = 0.10, start_t = 1),
  linear_arma_params = list(min_value = -1, max_value = 1, sigma_env = 0.10, 
    start_t = 1, ar = 0.4, ma = 0),
  constant_params = list(value = 0)
  ) {
  type <- type[1]
  env_ts <- switch(type,
    arma = as.numeric(with(arma_params, mean_value + arima.sim(model =
          list(ar = ar, ma = ma), n = n_t, sd = sigma_env))),
    sine = {
      x <- seq_len(n_t)
      y <- with(sine_params, amplitude * sin(ang_frequency * x + phase) +
        mean_value + x * slope + rnorm(n_t, mean = 0, sd = sigma_env))
      y
    },
    regime = {
      regime_params$break_pts <- c(1, regime_params$break_pts, n_t)
      y <- vector(length = n_t, mode = "numeric")
      for(i in 2:length(regime_params$break_pts))
      {
        y[regime_params$break_pts[i-1]:regime_params$break_pts[i]] <-
          regime_params$break_vals[i-1]
      }
      y
    },
    linear = {
      trend_years <- n_t - (linear_params$start_t - 1)
      trend <- with(linear_params, seq(min_value, max_value, length.out =
          trend_years) + rnorm(trend_years, mean = 0, sd = sigma_env))
      burnin_and_trend <- with(linear_params, c(rep(min_value, start_t - 1),
          trend))
      burnin_and_trend
    },
    linear_arma = {
      trend_years <- n_t - (linear_arma_params$start_t - 1)
      trend <- with(linear_arma_params,  
        seq(min_value, max_value, length.out = trend_years)) +  
        with(linear_arma_params, 
          as.numeric(arima.sim(model = list(ar = ar, ma = ma),  
            n = trend_years, sd = sigma_env))) 
      burnin_and_trend <- with(linear_arma_params, 
        c(rep(min_value, start_t - 1), trend))
      burnin_and_trend
    },
    constant = {
      rep(constant_params$value, n_t)
    }
    )
  return(env_ts)
}


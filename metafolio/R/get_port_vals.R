#' Get portfolio mean and variance values
#'
#' Takes a list created by \code{\link{meta_sim}} and returns the mean and
#' variance (or risk metric) values. This function is used by other internal
#' functions, but can also be used as its own low-level function.
#'
#' @param x A list object as returned from \code{\link{meta_sim}}
#' @param burn Number of years to throw out as burn in
#' @param risk_fn Type of variance or risk metric. By default takes the variance.
#'   Instead you can supply any function that takes a numeric vector and returns
#'   some single numeric value. E.g. CVaR.
#' @export
#' @seealso \code{\link{get_conserv_plans_mv}}, \code{\link{plot_cons_plans}}
#' @return A data frame with columns for the mean (m) and variance (v).
#' @examples
#' arma_env_params <- list(mean_value = 16, ar = 0.1, sigma_env = 2, ma = 0)
#' base1 <- meta_sim(n_pop = 10, env_params = arma_env_params, env_type =
#'   "arma", assess_freq = 5)
#' get_port_vals(base1)

get_port_vals <- function(x, risk_fn = var, burn = 1:30) {
    port.x <- rowSums(x$A[-burn, ])
    ret.x <- diff(log(port.x))
    var.x <- risk_fn(ret.x)
    mean.x <- mean(ret.x)
    data.frame(m = mean.x, v = var.x)
}

#' Conditional Value at Risk
#'
#' Get the conditional value at risk.
#'
#' @param x A numeric vector
#' @param probs The probability cutoff to pass to the CVaR function.
#' @export

CVaR <- function(x, probs = 0.05) {
  -mean(x[x < VaR(x, probs = probs)], na.rm = TRUE)
}

#' Value at Risk
#'
#' Get the value at risk.
#'
#' @param x A numeric vector
#' @param probs The probability cutoff to pass to the value at risk.
#' @export

VaR <- function(x, probs = 0.05) {
  quantile(x, probs = probs, na.rm = TRUE)
}


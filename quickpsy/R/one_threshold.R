#' Calculates the threshold for one condition
#' \code{one_threshold} calculates the threshold for one condition
#' one_threshold
#' @keywords internal
#' @export
one_threshold <- function(d, prob, log, groups, funname,
                          guess, lapses, curves) {

  if (length(groups) == 0) curves <- curves
  else curves <- semi_join(curves, d, by = as.character(groups))

  if (funname %in%  names(get_functions())) {
    par <- d$par
    if (is.numeric(guess) && is.numeric(lapses))
      q <- (prob - guess) / (1 - guess - lapses)
    if (is.logical(guess) && is.logical(lapses))
      q <- (prob - par[3]) / (1 - par[3] - par[4])
    if (is.logical(guess) && is.numeric(lapses))
      q <- (prob - par[3]) / (1 - par[3] - lapses)
    if (is.numeric(guess) && is.logical(lapses))
      q <- (prob - guess) / (1 - guess - par[3])

    if (q < 0 || q > 1) {
      warning('probabilities not whitin 0 and 1')
      thre <- approx(curves$y,curves$x, xout= prob)$y
    }
    else {
      if (funname == 'cum_normal_fun')
        thre <- inv_cum_normal_fun(q, c(par[1], par[2]))
      if (funname == 'logistic_fun')
        thre <- inv_logistic_fun(q, c(par[1], par[2]))
      if (funname == 'weibull_fun')
        thre <- inv_weibull_fun(q, c(par[1], par[2]))
    }
  }
  else {
    thre <- approx(curves$y,curves$x, xout= prob)$y
  }
  if (log) thre <- exp(thre)

  data.frame(thre, prob)
}






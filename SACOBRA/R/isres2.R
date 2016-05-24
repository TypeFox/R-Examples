isres2<-function (x0, fn, lower, upper, hin = NULL, heq = NULL, maxeval = 10000, 
          pop.size = 20 * (length(x0) + 1), xtol_rel = 1e-06, nl.info = FALSE, 
          ...) 
{
  opts <- list()
  opts$maxeval <- maxeval
  opts$xtol_rel <- xtol_rel
  opts$population <- pop.size
  opts$algorithm <- "NLOPT_GN_ISRES"
  fun <- match.fun(fn)
  fn <- function(x) fun(x, ...)
  if (!is.null(hin)) {
    .hin <- match.fun(hin)
    hin <- function(x) (-1) * .hin(x)
  }
  if (!is.null(heq)) {
    .heq <- match.fun(heq)
    heq <- function(x) .heq(x)
  }
  S0 <- nloptr::nloptr(x0 = x0, eval_f = fn, lb = lower, ub = upper, 
                       eval_g_ineq = hin, eval_g_eq = heq, opts = opts)
  if (nl.info) 
    print(S0)
  S1 <- list(par = S0$solution, value = S0$objective, iter = S0$iterations, 
             convergence = S0$status, message = S0$message)
  return(S1)
}

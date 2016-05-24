AICc.gambin <-
function (object, ...) 
{
  if (length(list(...)) > 0L) 
    warning("additional arguments ignored")
  ll <- logLik(object)
  df <- attr(ll, "df")
  val <- -2 * as.numeric(ll) + 2 * df
  val <- val + 2 * df * (df + 1)/(nobs(object) - df - 1)
  val
}

confint.gambin <-
function(object, parm = "alpha", level = 0.95, ...)
{
  if(!tolower(parm) == "alpha") stop("Only the alpha parameter has confidence intervals")
  .est_confint(object$Alpha, object$logLik, object$Data, level)
}

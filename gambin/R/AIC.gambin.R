AIC.gambin <-
function(object, ...)
{
  if (length(list(...)) > 0L) 
    warning("additional arguments ignored")
  AIC(logLik(object))
}

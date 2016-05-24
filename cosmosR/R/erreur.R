erreur <- function(err, var, ..., ret="ERR")
{
  warning("La variable ", var, " a cause l'erreur :\n", err, ..., call.=F)
  return (ret)
}
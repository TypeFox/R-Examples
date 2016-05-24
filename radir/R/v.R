v <-
function(f, x, beta) 
{
  pars <- get("pars",envir=.GlobalEnv)
  gradfunc <- function(x) {function(x){deriv(f, pars, function.arg="x")}}
  grad <- list()
  for (i in 1:length(beta)) 
  {
    grad[[i]] <- gradfunc(i)
  }
  for (i in 1:length(beta)) 
  {
    assign(eval(parse(text = paste0("pars[", i, "]"))), beta[i])
  }
  return(grad)
}

`compose.ar1` <-
function(arcoef, innov, init = 0, xregcoef = 0, xreg = NULL,
         init.xreg = rep(0, length(xregcoef)))
{
  init.innov <- init
  if(length(xreg) > 0)
    init.innov <- init.innov - as.vector(init.xreg %*% xregcoef)
  
  x <- apply(as.matrix(innov), 2, filter,
             filter = arcoef, method = 'r', init = init.innov)

  if(length(xreg) > 0)
    x <- x + as.vector(as.matrix(xreg) %*% xregcoef)
  
  result <- empty(innov)
  result[] <- x
  return(result)
}

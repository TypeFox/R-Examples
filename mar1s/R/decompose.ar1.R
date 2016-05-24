`decompose.ar1` <-
function(arcoef, data, init = NA, xregcoef = 0, xreg = NULL,
         init.xreg = rep(NA, length(xregcoef)))
{
  x <- as.matrix(data)
  if(length(xreg) > 0)
    x <- x - as.vector(as.matrix(xreg) %*% xregcoef)

  init.x <- init
  if(length(xreg) > 0)
    init.x <- init.x - as.vector(init.xreg %*% xregcoef)

  resid <- apply(rbind(init.x, x), 2, filter,
                 filter = c(1, -arcoef), method = 'c', sides = 1)
  
  result <- empty(data)
  result[] <- tail(resid, -1)
  return(result)
}

Adjust.thin <-
function(y)
{
n <- length(y)
m <- ar.ols(y,aic=F,order.max=1)
b <- m$ar
e <- m$resid
e <- matrix(e[!is.na(e)],nrow=n-1)
r <- 1/(1-b[1]) * e
return(r)
}

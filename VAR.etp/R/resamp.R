resamp <-
function(e)
{
n <- nrow(e)
random <- runif(n)
index <- as.integer(1+n*random)
es <- e[index, , drop=FALSE]
return(es)
}

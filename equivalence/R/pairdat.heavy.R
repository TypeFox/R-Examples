"pairdat.heavy" <-
function(n=100, VAR=4, M=0)
{
   x <- (exp(rnorm(n, 0, sqrt(2))) - exp(1)) * sqrt(VAR/2/(exp(4)-exp(2)))
   y <- (exp(rnorm(n, 0, sqrt(2))) - exp(1)) * sqrt(VAR/2/(exp(4)-exp(2))) + M
   list(x=x, y=y)
}


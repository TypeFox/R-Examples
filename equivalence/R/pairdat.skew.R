"pairdat.skew" <-
function(n=100, VAR=4, M=0)
{
   x <- rexp(n, sqrt(11/VAR)) - sqrt(VAR/11)
   y <- rexp(n, sqrt(11/10/VAR)) - sqrt(10*VAR/11) + M
   list(x=x, y=y)
}


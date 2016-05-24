"pairdat.norm" <-
function(n=100, M=0, SD=1)
{
   x <- rnorm(n, 0, SD/sqrt(2))
   y <- rnorm(n, 0, SD/sqrt(2)) + M
   list(x=x,y=y)
}


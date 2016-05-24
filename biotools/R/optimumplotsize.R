optimumplotsize <-
function(a, b)
{
   pow <- 1 / (2*b + 2)
   x0 <- (a^2 * b^2 * (2*b + 1) / (b + 2))^pow
   return(as.vector(x0))
}

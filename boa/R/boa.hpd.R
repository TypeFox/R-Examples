"boa.hpd" <-
function(x, alpha)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   n <- length(x)
   m <- max(1, ceiling(alpha * n))

   y <- sort(x)
   a <- y[1:m]
   b <- y[(n - m + 1):n]

   i <- order(b - a)[1]

   structure(c(a[i], b[i]), names = c("Lower Bound", "Upper Bound"))
}

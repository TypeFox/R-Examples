"sign.boot" <-
function(x, i, p1=0.2, p2=0.8)
{
   x <- x[i]
   n <- length(x)

   Nplus <- sum(x > 0)
   Nzero <- sum(x == 0)
   Nminus <- n - Nplus - Nzero
   if(sum(Nzero) > 0) cat("WARNING:", sum(Nzero), "ties occurred.\n")
   S1 <- qbinom(0.95, n, p1)
   S2 <- qbinom(0.05, n, p2)

   result <- sum(Nplus >= S1 & Nplus <= S2)
}


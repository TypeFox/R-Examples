TestUNey <- function(x, nrep = 10000, sim = NA, n.min = 30)
{
# This routine tests whether the values in each row of x are unif(0,1). It
# uses the Neyman's smooth test (see e.g., Ledwina 1994, TAS)
# x is a vector
# P-values are computed based on a
# resampling method from unif(0,1).
# All values of $x$ are between 0 and 1

  n <- length(x)
  pi <- LegNorm(x)
  n4 <- (apply(pi$p1, 2, sum) ^ 2 + apply(pi$p2, 2, sum) ^ 2 + 
         apply(pi$p3, 2, sum) ^ 2 + apply(pi$p4, 2, sum) ^ 2) / n 
  if (n < n.min){
      if(is.na(sim[1])) {
        sim <- SimNey(n, nrep)
      }
      pn <- length(which(sim > n4)) / nrep
  } else {
      pn <- pchisq(n4, 4, lower.tail = FALSE)
  }
  list(pn = pn, n4 = n4)
}
SimNey <- function(n, nrep)
{
 x <- matrix(runif(nrep * n), ncol = nrep)
 pi <- LegNorm(x)
 n4sim <- (apply(pi$p1, 2, sum) ^ 2 + apply(pi$p2, 2, sum) ^ 2 + 
          apply(pi$p3, 2, sum) ^ 2 + apply(pi$p4, 2, sum) ^ 2) / n 
 n4sim <- sort(n4sim)
}

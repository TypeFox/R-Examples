library(itertools)
library(foreach)

n <- 10
x <- matrix(rnorm(n * n), n)

# Split matrix x into four submatrices and put them back
# together again
y <-
  foreach(a=iarray(x, c(1,2), chunks=2), .combine='cbind') %:%
    foreach(b=a, .combine='rbind') %do% {
      b
    }
print(identical(x, y))

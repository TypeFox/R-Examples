library(maxLik)
set.seed(0)

## Test standard methods for "lm"
x <- runif(20)
y <- x + rnorm(20)
m <- lm(y ~ x)
print(nObs(m))
print(stdEr(m))

## Test maxControl methods:
set.seed(9)
x <- rnorm(20)
ll <- function(x) dnorm(x, log=TRUE)
for(method in c("NR", "BFGS", "BFGSR")) {
   cat("-- method", method, "--\n")
   m <- maxLik(ll, start=0, method=method, control=list(iterlim=1))
   cat("MaxControl structure:\n")
   show(maxControl(m))
}

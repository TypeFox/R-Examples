## Not really a test:
## test.rosenbrock <- function() {
##   n <- 10
  
##   ## Random optimum in [-50, 50]^n
##   opt <- runif(n, -50, 50)
##   bias <- 0
##   f <- genShiftedRosenbrock(opt, bias)

##   ## Inital parameter values
##   start <- runif(n, -100, 100)  
##   res <- cmaES(start, f, control=list(maxit=500))
##   print(res$par - opt)
## }

## Test what happens when the optimum is in one of the corners of the
## box constraint.
test.corner <- function() {
  f <- function(x, ...)
    drop(crossprod(x))

  n <- 5
  start <- runif(n, 1, 100)
  opt <- rep(1.0, n)

  res <- cma_es(start, f,
                lower=rep(1, n), upper=rep(100, n),
                control=list(maxit=2000))
  
  mse <- drop(sqrt(crossprod(res$par - opt)))
  checkTrue(mse < 0.002)
  checkTrue(res$value < 10.5)
}

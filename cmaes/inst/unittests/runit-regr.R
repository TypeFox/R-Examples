test.regr_stopfitness1 <- function() {
  f <- function(x, ...)
    drop(crossprod(x))
  
  n <- 2
  start <- c(0, 0)
  res <- cma_es(start, f,
                lower=rep(-10, n), upper=rep(10, n),
                control=list(stopfitness=1e-5, maxit=400))

  checkEqualsNumeric(res$convergence, 0)
  checkTrue(res$value < 1e-5, sprintf("res$value = %f >= 1e-5", res$value))
}

test.regr_stopfitness2 <- function() {
  f <- function(x, ...)
    -drop(crossprod(x))

  n <- 2
  start <- c(5, 5)
  res <- cma_es(start, f,
                lower=rep(-10, n), upper=rep(10, n),
                control=list(stopfitness=-1e-5, fnscale=-1, maxit=400))

  checkEqualsNumeric(res$convergence, 0)
  checkTrue(res$value > -1e-5, sprintf("res$value = %f <= -1e-5", res$value))
}

test.regr_bounds <- function() {
  f <- function(x, ...)
    drop(crossprod(x))

  par <- c(2, 2)
  l <- c(0.5, -10)
  u <- c(10, 10)
  ## Optimum lies on bounds. Try several times to reach infeasible
  ## region:
  for (i in 1:10) {
    res <- cma_es(par, f, lower=l, upper=u)
    par <- res$par
    vpar <- pmin(pmax(par, l), u)
    checkEqualsNumeric(drop(crossprod(par - vpar)), 0)
  }
}

resr.regr_names <- function() {
  f <- function(x, ...) {
    if (any(names(x) != c("a", "b", "c")))
      stop("BAM")
    drop(crossprod(x))
  }
  
  par <- c(a=2, b=2, c=3)
  ## Optimum lies on bounds. Try several times to reach infeasible
  ## region:
  res <- cma_es(par, f, lower=-10, upper=10)
  checkEquals(names(par), names(res$par))
}

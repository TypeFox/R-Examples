ifun <- function(fun) {
#	returns inverse of function fun
###########################
  recur <- function(fun, funinv=quote(x)) {
    fun <- as.expression(fun)[[1]]
#   if bracketed drop brackets
    while (length(fun) > 1 && fun[[1]] == as.name('(')) fun <- fun[[2]]
    if ((ne <- length(fun)) > 1) {
#     element of expression containing varname (either 2 or 3, ignoring leading symbol)
      x1 <- which(vapply(fun, function(f) nvars(f) == 1, TRUE))[-1]
#     element not containing varname (only relevant when ne == 3)
      n1 <- 5 - x1
#     expressions matching leading symbol with same length and x arg position
      nf <- which(vapply(fns, function(f)
        f[[1]] == fun[[1]] && length(f) == length(fun) && f[[x1]] == 'x', TRUE))
      if (length(nf) == 0) stop (paste('unrecognised name:', deparse(fun[[1]])))
#     if multiple matches check if numeric args are equal
      if (length(nf) > 1 && ne == 3) {
        nft <- which(vapply(fns[nf], function(f) f[[n1]] == fun[[n1]], TRUE))
        if (length(nft)) nf <- nf[nft]
      }
#     if more than one match use the first
      nf <- nf[[1]]
#     use complement of pair as inverse function
      fn2 <- if (nf %% 2) fns[[nf + 1]] else fns[[nf - 1]]
#     identify position of x arg in inverse function
      x2 <- which(as.list(fn2) == 'x')
#     if length 3 copy n arg
      if (length(fn2) == 3) {
        n2 <- 5 - x2
#       function returns value for n arg
        f <- function(n) {}
        body(f) <- fn2[[n2]]
        fn2[[n2]] <- f(eval(fun[[n1]]))
      }
#     copy x from current inverse function
      fn2[[x2]] <- funinv
#     update function and inverse function and repeat as necessary
      fun <- fun[[x1]]
      if (is.name(fun)) funinv <- fn2 else {
        funinv <- (results <- recur(fun, fn2))$fn
        fun <- results$varname
      }
    }
    return(list(fn=funinv, varname=fun))
  }
###########################
# number of names in function ignoring pi
  nvars <- function(fun) length((av <- all.vars(fun, unique=FALSE))[grep('pi', av, invert=TRUE)])
###########################
  if (nvars(fun) != 1) stop('expression should contain just one instance of one name')
  fns <- quote(c(x+n, x-n, x*n, x/n, x^n, x^(1/n), sqrt(x), x^2, exp(x), log(x),
                 expm1(x), log1p(x), n^x, log(x,n), log10(x), 10^x, log2(x),
                 2^x, n+x, x-n, n-x, n-x, n*x, x/n, n/x, n/x, +x, +x, -x, -x,
                 cos(x), acos(x), sin(x), asin(x), tan(x), atan(x),
                 cosh(x), acosh(x), sinh(x), asinh(x), tanh(x), atanh(x), I(x), I(x)))
  fns[[1]] <- NULL
  fn <- function(x) {}
  results <- with(fns, recur(fun))
  body(fn) <- results$fn
  return(list(fn=fn, varname=deparse(results$varname)))
}

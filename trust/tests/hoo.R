
 library(trust)
 options(digits = 3)

 d <- 5
 mu <- seq(1:d)

 objfun <- function(x) {
     stopifnot(is.numeric(x))
     stopifnot(length(x) == d)
     normxsq <- sum(x^2)
     omnormxsq <- 1 - normxsq
     if (normxsq >= 1) return(list(value = Inf))
     f <- sum(x * mu) - log(omnormxsq)
     g <- mu + 2 * x / omnormxsq
     B <- 4 * outer(x, x) / omnormxsq^2 + 2 * diag(d) / omnormxsq
     list(value = f, gradient = g, hessian = B)
 }

 whoop <- trust(objfun, rep(0, d), 1, 100, blather = TRUE)
 whoop$converged
 ceiling(log10(max(abs(whoop$gradient))))
 length(whoop$r)
 data.frame(type = whoop$steptype, # rho = round(whoop$rho, 2),
     change = whoop$preddiff, accept = whoop$accept, r = whoop$r)

 whoop$argument
 1 - sqrt(sum(whoop$argument^2))
 (whoop$stepnorm / whoop$r)[whoop$accept & whoop$steptype != "Newton"]

 mu <- 10 * mu

 whoop <- trust(objfun, rep(0, d), 1, 100, blather = TRUE)
 whoop$converged
 ceiling(log10(max(abs(whoop$gradient))))
 length(whoop$r)
 # give up.  This is just too ill-determined to be in tests
 if (FALSE) {
 data.frame(type = whoop$steptype, rho = round(whoop$rho, 2),
     change = whoop$preddiff, accept = whoop$accept, r = whoop$r)
 }

 whoop$argument
 1 - sqrt(sum(whoop$argument^2))
 (whoop$stepnorm / whoop$r)[whoop$accept & whoop$steptype != "Newton"]

 try(whoop <- trust(objfun, rep(0.5, d), 1, 100, blather = TRUE))


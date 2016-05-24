
 objfun <- function(x) {
     stopifnot(is.numeric(x))
     stopifnot(length(x) == 2)
     f <- expression(x1^2 - x2^2)
     g1 <- D(f, "x1")
     g2 <- D(f, "x2")
     h11 <- D(g1, "x1")
     h12 <- D(g1, "x2")
     h22 <- D(g2, "x2")
     x1 <- x[1]
     x2 <- x[2]
     f <- eval(f)
     g <- c(eval(g1), eval(g2))
     B <- rbind(c(eval(h11), eval(h12)), c(eval(h12), eval(h22)))
     list(value = f, gradient = g, hessian = B)
 }

 library(trust)

 options(digits = 3)

 goo <- trust(objfun, c(3, 1), 1, 5, blather = TRUE, iterlim = 20)
 goo$converged
 length(goo$r)
 data.frame(type = goo$steptype, rho = goo$rho, change = goo$preddiff,
     accept = goo$accept, r = goo$r)

 ##### note: FAILS to converge because function is unbounded -- optimum
 #####     value does not exist

 (goo$stepnorm / goo$r)[goo$accept & goo$steptype != "Newton"]


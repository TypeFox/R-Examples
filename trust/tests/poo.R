
 objfun <- function(x) {
     ##### Rosenbrock's function #####
     stopifnot(is.numeric(x))
     stopifnot(length(x) == 2)
     f <- expression(100 * (x2 - x1^2)^2 + (1 - x1)^2)
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

 parinit <- c(3, 1)

 out <- trust(objfun, parinit, 1, 1e5, blather = TRUE)
 out$converged
 length(out$r)

 parscale <- c(5, 1)
 shift <- 4
 theta <- parscale * (parinit + shift)

 pobjfun <- function(x) {
     out <- objfun(x / parscale - shift)
     out$gradient <- out$gradient / parscale
     out$hessian <- out$hessian / outer(parscale, parscale)
     return(out)
 }

 pout <- trust(pobjfun, theta, 1, 1e5, blather = TRUE)
 pout$converged
 length(pout$r)
 all.equal(out$argument, pout$argument / parscale - shift)

 qout <- trust(objfun, parinit, 1, 1e5, parscale = parscale, blather = TRUE)
 qout$converged
 length(qout$r)

 all.equal(pout$valpath, qout$valpath)
 transpath <- pout$argpath
 transpath <- sweep(transpath, 2, parscale, "/")
 transpath <- sweep(transpath, 2, shift)
 all.equal(transpath, qout$argpath)


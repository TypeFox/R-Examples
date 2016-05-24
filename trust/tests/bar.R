
 objfun <- function(x) {
     stopifnot(is.numeric(x))
     stopifnot(length(x) == 2)
     f <- expression(x1^2 - x2^2 + (1 / 100) * (x1^4 + x2^4))
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

 tout <- trust(objfun, c(0, 0), 1, 5, blather = TRUE)
 tout

 (tout$stepnorm / tout$r)[tout$accept & tout$steptype != "Newton"]

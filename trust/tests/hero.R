
 library(trust)

 options(error = function() NULL)

 ##### TEST IMMEDIATE ERROR #####

 objfun <- function(x) stop("bogus error in objfun")

 trust(objfun, c(3, 1), 1, 5, blather = TRUE)

 ##### TEST LATER ERROR #####

 objfun <- function(x) {
     kiter <<- kiter + 1
     if (kiter == 5)
         stop("kiter reached 5")
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

 kiter <- 0
 trust(objfun, c(3, 1), 1, 5, blather = TRUE)

 ##### TEST ERROR IN LAST #####

 kiter <- 0
 trust(objfun, c(3, 1), 1, 5, blather = TRUE, iterlim = 3)


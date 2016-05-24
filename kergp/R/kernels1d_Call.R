
k1expFun <- function(x1, x2, par) {
  res <- .Call("k1Exp", x1, x2, par, PACKAGE = "kergp")
}

k1gaussFun <- function(x1, x2, par) {
  res <- .Call("k1Gauss", x1, x2, par, PACKAGE = "kergp")
}

k1powExpFun <- function(x1, x2, par) {
  res <- .Call("k1PowExp", x1, x2, par, PACKAGE = "kergp")
}

k1matern3_2Fun<- function(x1, x2, par) {
  res <- .Call("k1Matern3_2", x1, x2, par, PACKAGE = "kergp")
}

k1matern5_2Fun <- function(x1, x2, par) {
  res <- .Call("k1Matern5_2", x1, x2, par, PACKAGE = "kergp")
}

k1exp <- new("covMan",            
             kernel = k1expFun,
             label = "exponential",
             d = 1L,
             inputNames = NULL,
             parLower = c("range" = 0.0, "var" = 0.0),
             parUpper = c("range" = Inf, "var" = Inf),
             par = c("range" = 1.0, "var" = 1.0),
             parN = 2L,
             kernParNames = c("range", "var"))

k1gauss <- new("covMan",            
               kernel = k1gaussFun,
               label = "gaussian",
               d = 1L,
               inputNames = NULL,
               parLower = c("range" = 0.0, "var" = 0.0),
               parUpper = c("range" = Inf, "var" = Inf),
               par = c("range" = 1.0, "var" = 1.0),
               parN = 2L,
               kernParNames = c("range", "var"))

k1matern3_2 <- new("covMan",            
                   kernel = k1matern3_2Fun,
                   label = "Matern nu = 3/2",
                   d = 1L,
                   inputNames = NULL,
                   parLower = c("range" = 0.0, "var" = 0.0),
                   parUpper = c("range" = Inf, "var" = Inf),
                   par = c("range" = 1.0, "var" = 1.0),
                   parN = 2L,
                   kernParNames = c("range", "var"))

k1matern5_2 <- new("covMan",            
                   kernel = k1matern5_2Fun,
                   label = "Matern nu = 5/2",
                   d = 1L,
                   inputNames = NULL,
                   parLower = c("range" = 0.0, "var" = 0.0),
                   parUpper = c("range" = Inf, "var" = Inf),
                   par = c("range" = 1.0, "var" = 1.0),
                   parN = 2L,
                   kernParNames = c("range", "var"))

k1powExp <- new("covMan",            
                kernel = k1powExpFun,
                label = "power exponential",
                d = 1L,
                inputNames = NULL,
                parLower = c("range" = 0.0, "shape" = 0.0, "var" = 0.0),
                parUpper = c("range" = Inf, "shape" = 2.0, "var" = Inf),
                par = c("range" = 1.0, "shape" = 1.5, "var" = 1.0),
                parN = 3L,
                kernParNames = c("range", "shape", "var"))

setAs("covMan", "function", function(from) from@kernel)

GPkernNames <- c("k1exp", "k1matern3_2", "k1matern5_2", "k1powExp", "k1gauss")

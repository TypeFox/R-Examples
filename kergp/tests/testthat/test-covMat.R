library(testthat)
library(kergp)

set.seet(314159)
myCov <- covTS(inputs = c("Temp", "Humid", "Press"),
               kernel = "k1powExp",
               dep = c(range = "cst", shape = "cst"),
               value = c(shape = 1.8, range = 1.1))
n <- 100; X <- matrix(runif(n*3), nrow = n, ncol = 3)
colnames(X) <- inputNames(myCov)

## check with the same matrix in X and Xnew
CMM <- covMat(myCov, X, X)
CM <- covMat(myCov, X)
err <- max(abs(CM - CMM))

test_that(desc = "covMat at 1 versus 2 identical sets of locations: ",
          code = expect_true(err < 1e-10))

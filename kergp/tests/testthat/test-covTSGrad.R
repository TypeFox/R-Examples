library(testthat)
library(kergp)
library(numDeriv)


## require (numDeriv)  now in 'Depends'
precision <- 1e-6

d <- 3; n <- 10
set.seed(1234)

## utility function
covAsVec <- function(par, object, index, ...) {
  coef(object) <- par
  C <- covMat(object = object, X = X, compGrad = TRUE, index = index)
  grad <- attr(C, "gradient")
  C <- C[row(C) <= col(C)]
  grad <- grad[row(grad) <= col(grad)]
  attr(C, "gradient") <- grad
  C
}

## define a covTS
myCov <- covTS(d = d, kernel = "k1powExp",
               dep = c(range = "input", shape = "input"),
               value = c(range = 0.1, shape = 0.6))

## choose a design and parameter values
X <- array(runif(n * d), dim = c(n, d),
           dimnames = list(NULL, inputNames(myCov)))

Theta <- simulPar(object = myCov, n = 10L)
index <- 1
theta <- Theta[1, ]
res <- covAsVec(theta, myCov, index = index)
grad.check <- jacobian(covAsVec, theta, object = myCov, index = index)

mat <- cbind(check = grad.check[ , index], prog = attr(res, "gradient"))
err <- max(abs(apply(mat, 1, function (x) mean(abs(diff(x))))))

test_that(desc = "gradient of a \"covTS\" structure",
          code = expect_true(err < precision))

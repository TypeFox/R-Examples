context("Fast updating of (I_n - rho * W) for new values of rho")

# faster update of matrix S = (I - rho * W) for new values of rho.
# Problem: matrix subtraction (I - rho * W) is VERY time-consuming 
# because of generic methods and method dispatching.
# Idea: for changing rho we have to update the same matrix elements in S, so lets
# create a template matrix and for each new rho just replace the updated matrix elements directly.
#
# @param S template matrix of (I - rho * W)
# @param ind indizes to replaced
# @param W spatial weights matrix W works only for "dgCMatrix"
# @return (I - rho * W)
update_I_rW <- function(S, ind, rho, W) {
  S@x[ind] <- (-rho*W)@x
  return(S)
}

library(spatialprobit)
n <- 100
I_n <- sparseMatrix(i=1:n,j=1:n,x=1)
rho <- 0.5
W <- kNearestNeighbors(x=rnorm(n), y=rnorm(n), k=6)
# prepare computation of (I_n - rho * W)
if (class(W) == "dgCMatrix") {
  I <- sparseMatrix(i=1:n,j=1:n,x=Inf)  # use Inf as placeholder
  S <- (I - rho * W)
  ind  <- which(is.infinite(S@x))  # Stellen an denen wir 1 einsetzen müssen (I_n)
  ind2 <- which(!is.infinite(S@x))  # Stellen an denen wir -rho*W einsetzen müssen
  S@x[ind] <- 1
} else {
  S <- I_n - rho * W
}

test_that("update_I_rW works correctly", {

  expect_true(class(W) == "dgCMatrix")
  
  # compare update_I_rW vs. (I_n - rho*W)
  rhos <- seq(-1, 1, by=0.01)
  flags <- rep(FALSE, length(rhos))
  for (i in seq(along=rhos)) {
    rho <- rhos[i]
    S1 <- I_n - rho * W
    S2 <- update_I_rW(S, ind=ind2, rho, W)
    flags[i] <- all.equal(S1, S2)
  }
  expect_true(all(flags[i]))
})

test_that("update_I_rW is faster than (I_n - rho*W)", {

  time1 <- system.time(for (rho in seq(-1, 1, by=0.01)) S <- I_n - rho * W)
  time2 <- system.time(for (rho in seq(-1, 1, by=0.01)) S <- update_I_rW(S, ind=ind2, rho, W))

  cond <- time2["elapsed"] < (time1["elapsed"] * 0.7)  # expect at least 30% performance gain
  names(cond) <- NULL                                  # problem with named condition
  expect_true(cond)

})  



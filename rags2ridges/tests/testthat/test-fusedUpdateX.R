context("Unit test of the .armaFusedUpdateX-familiy of functions")

#
# Avoid writing rags2ridges:::
#

.armaFusedUpdateI <- rags2ridges:::.armaFusedUpdateI
.armaFusedUpdateII <- rags2ridges:::.armaFusedUpdateII
.armaFusedUpdateIII <- rags2ridges:::.armaFusedUpdateIII

#
# Define R-versions of the fusedUpdateX functions
#

.fusedUpdateI <- function(g0, Plist, Slist, Tlist, ns, lambda) {
  ##############################################################################
  # - (Internal) "Update" the covariance matrices and use the regular
  #   ridge estimate. The scheme I approach.
  # - g0      > An integer giving the class estimate to be updated.
  # - Plist   > A list of length G of matrices giving the current precision
  #             estimates.
  # - Slist   > A list of length G of sample correlation matrices the same size
  #             as those of Plist.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist
  # - ns      > A vector of length G giving the sample sizes.
  # - lambda  > A G by G symmetric adjacency matrix giving the fused penalty
  #             graph with non-negative entries where lambda[g1, g2]
  #             determine the (rate of) shrinkage between estimates in classes
  #             corresponding to Slist[g1] and Slist[g1].
  #
  #   NOTE: The update function seems to work ok for large lambda.
  #   However, for very large lambda (> 1e154) the exception that the
  #   .armaRidgeP returns the target because of an exception. Which is wrong
  #   in the fused case.
  ##############################################################################

  a <- (sum(lambda[g0, ]))/ns[g0]
  b <- lambda[g0, -g0]/ns[g0]

  OmT <- mapply(`-`, Plist[-g0], Tlist[-g0], SIMPLIFY = FALSE) # Omega - Target
  OmT <- mapply(`*`, b, OmT, SIMPLIFY = FALSE)
  S0 <- Slist[[g0]] - Reduce(`+`, OmT)
  return(rags2ridges:::.armaRidgeP(S0, target = Tlist[[g0]], lambda = a))
}


.fusedUpdateII <- function(g0, Plist, Slist, Tlist, ns, lambda) {
  ##############################################################################
  # - (Internal) "Update" the covariance matrices and use the regular
  #   ridge estimate -- using the alternative II update scheme.
  # - g0      > An integer giving the class estimate to be updated.
  # - Plist   > A list of length G of matrices giving the current precision
  #             estimates.
  # - Slist   > A list of length G of sample correlation matrices the same size
  #             as those of Plist.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist
  # - ns      > A vector of length G giving the sample sizes.
  # - lambda  > A G by G symmetric adjacency matrix giving the fused penalty
  #             graph with non-negative entries where lambda[g1, g2]
  #             determine the (rate of) shrinkage between estimates in classes
  #             corresponding to Slist[g1] and Slist[g1].
  #
  #   NOTE: This update seems to work very poorly for large lambda
  ##############################################################################

  p <- nrow(Plist[[1]])
  lambdaa <- sum(lambda[g0, ])/ns[g0]
  b <- (sum(lambda[g0, ]) - 1)/ns[g0]

  Psum <- Tsum <- matrix(0, p, p)
  for (g in setdiff(seq_along(Plist), g0)) {
    Psum <- Psum + lambda[g0, g]*Plist[[g]]
    Tsum <- Tsum + (lambda[g0, g]/ns[g0])*Tlist[[g]]
  }
  Sbar <- Slist[[g0]] + b*Psum + Tsum
  Tbar <- Tlist[[g0]] + Psum
  return(rags2ridges:::.armaRidgeP(Sbar, target = Tbar, lambda = lambdaa))
}



.fusedUpdateIII <- function(g0, Plist, Slist, Tlist, ns, lambda) {
  ##############################################################################
  # - (Internal) "Update" the covariance matrices and use the regular
  #   ridge estimate -- using the alternative III update scheme.
  # - g0      > An integer giving the class estimate to be updated.
  # - Plist   > A list of length G of matrices giving the current precision
  #             estimates.
  # - Slist   > A list of length G of sample correlation matrices the same size
  #             as those of Plist.
  # - Tlist   > A list of length G of target matrices the same size
  #             as those of Plist
  # - ns      > A vector of length G giving the sample sizes.
  # - lambda  > A G by G symmetric adjacency matrix giving the fused penalty
  #             graph with non-negative entries where lambda[g1, g2]
  #             determine the (rate of) shrinkage between estimates in classes
  #             corresponding to Slist[g1] and Slist[g1].
  #
  #   NOTE: This update function seems to work very well for large lambda.
  #   For very large lambda (> 1e154) the exception triggered in the
  #   .armaRidgeP returns the target because of an exception. However, in this
  #   updating scheme, that is also correct.
  ##############################################################################

  lambdasum <- sum(lambda[g0, ])
  lambdaa <- lambdasum/ns[g0]

  Tbar <- Tlist[[g0]]
  for (g in setdiff(seq_along(Plist), g0)) {
    Tbar <- Tbar + (lambda[g0, g]/lambdasum)*(Plist[[g]] - Tlist[[g]])
  }

  return(rags2ridges:::.armaRidgeP(Slist[[g0]], target = Tbar,
                                   lambda = lambdaa))
}




#
# Actual tests
#

# Create some data
ns. <- c(5, 6, 7)
g0. <- sample(seq_along(ns.) - 1, 1)
Plist. <- createS(n = ns., p = 5, topology = "star", precision = TRUE)
Slist. <- createS(n = ns., Plist = Plist.)
Tlist. <- replicate(length(ns.), diag(5), simplify = FALSE)
lm <- symm(matrix(rchisq(9, df = 3), 3, 3))

# Compute updated
A <- .armaFusedUpdateI(g0 = g0., Plist = Plist., Slist = Slist.,
                       Tlist = Tlist.,ns = ns., lambda = lm)
B <- .armaFusedUpdateII(g0 = g0., Plist = Plist., Slist = Slist.,
                        Tlist = Tlist.,ns = ns., lambda = lm)
C <- .armaFusedUpdateIII(g0 = g0., Plist = Plist., Slist = Slist.,
                         Tlist = Tlist.,ns = ns., lambda = lm)

test_that(".fusedUpdateX returns correctly formatted output", {
  # Test that results are numeric matrices that are symmetric PD
  expect_true(isSymmetricPD(A))
  expect_true(isSymmetricPD(B))
  expect_true(isSymmetricPD(C))

  # Test that they are equal
  expect_equal(A, B)
  expect_equal(A, C)
})


# Note these are index from 1 and forward!
AA <- .fusedUpdateI(g0 = g0. + 1., Plist = Plist., Slist = Slist.,
                    Tlist = Tlist.,ns = ns., lambda = lm)
BB <- .fusedUpdateII(g0 = g0. + 1., Plist = Plist., Slist = Slist.,
                     Tlist = Tlist.,ns = ns., lambda = lm)
CC <- .fusedUpdateIII(g0 = g0. + 1., Plist = Plist., Slist = Slist.,
                      Tlist = Tlist.,ns = ns., lambda = lm)

test_that(".fusedUpdateX return similar results", {
  expect_equal(A, AA)
  expect_equal(B, BB)
  expect_equal(C, CC)
})

test_that(".fusedUpdateX works properly on degenerated data", {
  expect_true(TRUE)  # To be written
})

test_that(".fusedUpdateX works properly on extreme penalties", {
  expect_true(TRUE)  # To be written
})


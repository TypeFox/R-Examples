context("Unit test of the .armaRidgeP.fused function")

armaRidgeP.fused <- rags2ridges:::.armaRidgeP.fused
armaRidgeP       <- rags2ridges:::.armaRidgeP

# Random number of sample, random number of classes
p <- 10
n <- replicate(sample(2:5, 1), sample(3:9, 1))
S <- createS(n = n, p = p)
tgt <- default.target.fused(S, n)
G <- length(n)


test_that("armaRidgeP.fused returns proper format", {

  test.lambdas <- c(1e-200, 1e-100, 1e-50, 1e-10, 1, 1e10, 1e50, 1e100,
                    1e200, 1e300, 1e500, Inf)
  lm <- matrix(1, length(n), length(n))

  for (l in test.lambdas) {

    diag(lm) <- l
    res <- armaRidgeP.fused(S, n, tgt, lambda = lm, Plist = tgt)

    expect_that(res, is_a("list"))
    expect_that(length(res), equals(G))
    for (i in 1:G) {
      expect_that(res[[i]], is_a("matrix"))  # Returns a matrix
      expect_that(typeof(res[[i]]), equals("double"))  # Returns a matrix
      expect_that(dim(res[[i]]), equals(c(p, p)))   # .. of the correct size
    }

  }

})

test_that("armaRidgeP.fused does not ignore the diagonal in lambda", {

  lm <- matrix(1, G, G)
  res.diag <- armaRidgeP.fused(S, n, tgt, lambda = lm, tgt)
  diag(lm) <- 0
  res.nodiag <- armaRidgeP.fused(S, n, tgt, lambda = lm, tgt)

  expect_false(isTRUE(all.equal(res.diag, res.nodiag)))
})


################################################################################
# Test against the old funciton:
################################################################################

# Test against an old, slow implementation
ridgeP.fused.old <- function(Slist, ns, Tlist = default.target.fused(Slist, ns),
                             lambda, lambdaF, lambdaFmat, Plist,
                             maxit = 100L, verbose = TRUE, eps = 1e-4) {

  stopifnot(length(Slist) == length(Tlist))
  G <- length(Slist)  # Number of groups

  # Initialize estimates with the regular ridges from the pooled covariance
  if (missing(Plist)) {
    Spool <- pooledS(Slist, ns, mle = FALSE)
    Plist <- list()
    for (i in seq_len(G)) {
      Plist[[i]] <- .armaRidgeP(Spool, target = Tlist[[i]],
                                lambda = G*lambda/sum(ns))
    }
  }
  stopifnot(length(Slist) == length(Plist))

  if (!missing(lambdaF) && !missing(lambdaFmat)) {
    stop("Supply only either lambdaF or lambdaFmat.")
  } else if (missing(lambdaF) && missing(lambdaFmat)) {
    stop("Either lambdaF or lambdaFmat must be given.")
  } else if (missing(lambdaF) && !missing(lambdaFmat)) {
    lambdaF <- matrix(lambdaF, G, G)
  }

  if (verbose) {
    cat("iteration | difference in Frobenius norm\n")
  }

  tmp.lambda <- lambdaF
  diag(tmp.lambda) <- lambda
  lambda <- tmp.lambda

  lambdasize <- sum(lambdaF)
  tmpPlist <- list()
  diffs <- rep(NA, G)
  i <- 1
  while (i <= maxit) {
    for (g in seq_len(G)) {
      if (lambdasize < 1e50) {
        tmpPlist[[g]] <-
          .armaFusedUpdateI(g0 = g-1, Plist = Plist, Slist = Slist,
                            Tlist = Tlist, ns = ns, lambda = lambda)
      } else {
        tmpPlist[[g]] <-
          .armaFusedUpdateIII(g0 = g-1, Plist = Plist, Slist = Slist,
                              Tlist = Tlist, ns = ns, lambda = lambda)
      }
      diffs[g] <- .FrobeniusLoss(tmpPlist[[g]], Plist[[g]])
      Plist[[g]] <- tmpPlist[[g]]
    }
    mx <- max(diffs)
    if (verbose) {
      cat(sprintf("i = %-3d | max diffs = %0.10f\n", i, mx))
    }
    if (is.nan(mx)) {
      warning("NaNs where introduced likely due to very largs penalties.")
      break
    }
    if (mx < eps) {
      break
    }
    i <- i + 1
  }
  if (i == maxit + 1) {
    warning("Maximum iterations (", maxit, ") hit")
  }
  # Keep dimnames and names
  for (g in seq_along(Slist)) {
    dimnames(Plist[[g]]) <- dimnames(Slist[[g]])
  }
  names(Plist) <- names(Slist)
  return(Plist)
}
environment(ridgeP.fused.old) <- asNamespace('rags2ridges')



# TESTING
test_that("armaRidgeP.fused agrees with the old R implmentation", {

  lambda <- matrix(1, G, G)
  diag(lambda) <- abs(rcauchy(n = 1))
  Spool <- pooledS(S, n, mle = FALSE)
  P <- list()
  for (i in seq_len(G)) {
    P[[i]] <- armaRidgeP(Spool, target = tgt[[i]],
                         lambda = G*lambda[1,1]/sum(n))
  }

  res_old <- ridgeP.fused.old(S, n, tgt, lambda = lambda[1,1], lambdaF = lambda,
                              maxit = 1000, eps = 1e-10, verbose = FALSE)
  res_new <- armaRidgeP.fused(S, n, tgt, lambda, Plist = P,
                              maxit = 1000, eps = 1e-10, verbose = FALSE)

  expect_true(isTRUE(all.equal(res_old, res_new,
                               check.attributes = FALSE,
                               tolerance = 1e-4)))
})






test.varcomp <-
function(fixed, random, test, data = NULL, Sigma = NULL, 
  type = c("RLR", "pseudo", "score"), nsim = 5000L, seed = 130623L,
  keep.matrices = FALSE)
{
  type <- match.arg(type)
  
  mf.fixed <- model.frame(fixed, data = data, drop.unused.levels = TRUE)
  Y <- model.response(mf.fixed)
  mt.fixed <- attr(mf.fixed, "terms")
  X <- model.matrix(mt.fixed, mf.fixed, contrasts.arg = contr.treatment)
  
  mf.random <- model.frame(random, data = data, drop.unused.levels = TRUE)
  mt.random <- attr(mf.random, "terms")
  Z <- lapply(attr(mt.random, "term.labels"), function(term) {
    model.matrix(formula(paste0("~ -1 + ", term)), mf.random,
      contrasts.arg = contr.treatment)
  })
  if (attr(mt.random, "intercept")) {
    Z <- c(list(matrix(1, length(Y), 1L)), Z)
  }
  if (is.null(Sigma)) {
    Sigma <- lapply(sapply(Z, ncol), diag)
  } else {
    if (length(Z) != length(Sigma)) {
      stop(paste0("There are ", length(Z), " random effects design matrices",
        ", but ", length(Sigma), " correlation strutures."))
    }
  }
  Z <- c(Z[-test], Z[test])
  Sigma <- c(Sigma[-test], Sigma[test])
  m0 <- length(Z) - length(test)
  
  if (type == "RLR") {
    result <- rlr.test(Y, X, Z, Sigma, m0, nsim, seed)
  } else if (type == "pseudo") {
    result <- pseudo.rlr.test(Y, X, Z, Sigma, m0, nsim, seed)
    result <- list(pseudo = result)
  } else if (type == "score") {
    result <- score.test(Y, X, Z, Sigma, m0)
    result <- list(score = result)
  }

  if (keep.matrices) {
    result$Y <- Y
    result$X <- X
    result$Z <- Z
    result$Sigma <- Sigma
  }
  
  result
}

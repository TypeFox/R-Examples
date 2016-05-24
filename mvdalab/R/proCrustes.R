proCrustes <- function(X, Y, scaling = TRUE, standardize = FALSE, scale.unit = F, ...) {
  Col.Diff <- (ncol(X) - ncol(Y))
  if (ncol(X) > ncol(Y)) {
    Y <- data.frame(Y, Added = matrix(0, nrow = nrow(X), ncol = Col.Diff))
  } else {
    X <- X
    Y <- Y
  }
  X. <- scale(X, scale = scale.unit)
  Y. <- scale(Y, scale = scale.unit)
  Xmeans <- attr(X., "scaled:center")
  Ymeans <- attr(Y., "scaled:center")
  if (!(standardize)) {
    X. <- X.
    Y. <- Y.
  } else {
    X. <- X.
    Y. <- Y.
    X./sqrt(sum(diag(crossprod(X.))))
    Y./sqrt(sum(diag(crossprod(Y.))))
  }
  SVD <- svd(t(X.) %*% Y.)
  Q <- SVD$v %*% t(SVD$u) #Rotation Matrix
  if (!(scaling)) {
    c. <- 1
  } else {
    c. <- sum(diag(SVD$d)) / sum(diag(Y. %*% t(Y.)))
  }
  M2_min <- sum(diag(X. %*% t(X.))) + (c.^2 * sum(diag(Y. %*% t(Y.)))) - (2 * c. * sum(diag(SVD$d)))
  PRMSE <-   sqrt(M2_min / nrow(X.))
  Yproj <- c. * Y. %*% Q
  Translation <- Xmeans - c. * Ymeans %*% Q
  difference <- X. - t(Q %*% t(Y.))
  residuals. <- sqrt(apply(difference^2, 1, sum))
  MSS <- c.^2 * sum(diag(Y. %*% t(Y.)))
  ESS <- M2_min
  TSS <- MSS + ESS
Results <- list(Rotation.Matrix = Q, Residuals = difference, M2_min = M2_min, Xmeans = Xmeans,
                Ymeans = Ymeans, PRMSE = PRMSE, Yproj = Yproj, 
                scale = c., Translation = Translation, residuals. = residuals., 
                Anova.MSS = MSS, Anova.ESS = ESS, Anova.TSS = TSS)
class(Results) <- "proC"
Results

}
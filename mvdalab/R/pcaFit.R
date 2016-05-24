pcaFit <- function(data, scale = TRUE, ncomp = NULL){
  X <- my.dummy.df(data)
  if(is.null(ncomp)) {
    if(nrow(X) <= ncol(X)) {
      ncomp <- nrow(X) 
    } else ncomp <- ncol(X)
  } else ncomp <- ncomp
  Xmeans <- colMeans(X)
  Scale <- scale
  nobj <- nrow(X)
  npred <- ncol(X)
  X <- X - rep(Xmeans, each = nobj)
    if (Scale == TRUE) {
      scale. <- sqrt(colSums((X - rep(colMeans(X), each = nobj))^2)/(nobj - 1))
      if (any(abs(scale.) < .Machine$double.eps^0.5)) 
        warning("Scaling with (near) zero standard deviation")
      X <- X/rep(scale., each = nobj)
    } else X
SVD <- svd(X)
D.total <- diag((SVD$d^2)) / (nrow(X) - 1)
  if(ncomp == 1) {
    U <- as.matrix(SVD$u[1:ncomp]) %*% diag(as.matrix(SVD$d[1:ncomp]))
    D <- diag(as.matrix(SVD$d[1:ncomp]^2)) / (nrow(X) - 1)
  } else {
    U <- SVD$u[, 1:ncomp] %*% diag((SVD$d[1:ncomp]))
    D <- diag((SVD$d[1:ncomp]^2)) / (nrow(X) - 1)
  }
  loadings.a <- t(sqrt(D) %*% t(SVD$v[, 1:ncomp]))
  loadings <- as.matrix(apply(loadings.a, 2, function(x) x / sqrt(crossprod(x))))
  scores <- as.matrix(U)
  GVC <- sapply(1:ncomp, function(tc) {
    Us <- apply(as.matrix(SVD$u[, 1:tc])^2, 1, sum)
    Vs <- apply(as.matrix(SVD$v[, 1:tc]^2), 1, sum)
    Recon <- (as.matrix(scores[, 1:tc]) %*% t(as.matrix(loadings[, 1:tc])))
      Pre.SS3 <- mean(((nobj * npred - sum(is.na(X))) * (X - Recon)/
                        ((nobj - 1) * npred - sum(is.na(X)) - tc * (nobj + 
                            npred - tc - 1)))^2, na.rm = T)
  }); 
  Zero <- mean((X - rep(colMeans(X, na.rm = TRUE), each = nrow(X)))^2, 
               na.rm = TRUE)
  GVC <- c(Zero, GVC[(is.finite(GVC)) & GVC > 1e-10]) 
  Percents.Explained <- data.frame(Individual.Percent = (diag(D) / sum(diag(D.total))) * 100, 
                                   Cumulative.Percent = cumsum(diag(D) / sum(diag(D.total))) * 100)
  Percents.Explained$Comp <- 1:nrow(Percents.Explained)
  row.names(loadings) <- names(X)
  output <- list(loadings = loadings, scores = scores, D = D,
                 Xdata = X, Percents.Explained = Percents.Explained, 
                 GVC = GVC, ncomp = ncomp, Xmeans = Xmeans,
                 method = "PCA")
  class(output$scores) <- "scores"
  class(output$loadings) <- "loadings"
  class(output) <- "mvdapca"
  output
}





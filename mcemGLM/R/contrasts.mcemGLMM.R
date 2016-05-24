contrasts.mcemGLMM <- function(object, ctr.mat) {
  names0 <- rownames(ctr.mat)
  coef0 <- as.vector(tail(object$mcemEST, n = 1))
  kP <- ncol(object$x)
  kR <- length(coef0) - kP
  for (i in 1:kR) {
    ctr.mat <- cbind(ctr.mat, 0)
  }
  wald0 <- rep(0, nrow(ctr.mat))
  ctr.est0 <- rep(0, nrow(ctr.mat))
  ctr.err0 <- rep(0, nrow(ctr.mat))
  for (i in 1:nrow(ctr.mat)) {
    cm <- t(ctr.mat[i, ])
    ctr.est0[i] <- cm %*% coef0
    ctr.err0[i] <- solve(cm %*% solve(object$iMatrix) %*% t(cm))
    wald0[i] <- ctr.est0[i]^2 * ctr.err0[i]
  }
  std.err <- 1/sqrt(ctr.err0)
  pval0 <- pchisq(wald0, 1, lower.tail = FALSE) * nrow(ctr.mat)
  pval0 <- ifelse(pval0 < 1, pval0, 1)
  tbr <- cbind(ctr.est0, std.err, wald0, pval0)
  colnames(tbr) <- c("Estimate", "Std. Err.", "Wald", "Adj. p-value")
  rownames(tbr) <- rownames(ctr.mat)
  return(tbr)
}
ScoreContrib <- function(object, obs1 = 1, obs2 = NULL) {
  ncomp <- object$ncomp
  S <- as.matrix(object$scores)
  P <- as.matrix(object$loadings)
  X <- as.matrix(object$Xdata)
  if(class(object) == "mvdapca") {
    Ww <- (P %*% solve(t(P) %*% P))
    } else {
    W <- as.matrix(object$weights)
    Ww <- (W %*% solve(t(P) %*% W))
  }
  if(is.null(obs2)) {
    Xrow <- apply(object$Xdata[obs1, ], 2, mean)
  } else {
    Xrow <- apply(object$Xdata[obs1, ], 2, mean) - apply(object$Xdata[obs2, ], 2, mean)
  }
    New.S <- as.numeric(Xrow %*% Ww) 
    Mults <- t(sqrt((New.S / sqrt(diag(cov(S))))^2 * t(((Ww)^2))))
    CPper <- as.matrix(Xrow * Mults)
    CPall <- as.matrix(Xrow * apply(Mults, 1, function(x) sqrt(sum(x^2))))
    Contributions <- data.frame(CPper, Overall = CPall)
    names(Contributions)[1:ncomp] <- c(1:ncomp)
    Contributions$Variable <- rownames(Contributions)
    row.names(Contributions) <- NULL
    Contributions <- list(score.contribution = Contributions)
    class(Contributions) <- "cp"
    Contributions
}

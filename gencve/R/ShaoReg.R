ShaoReg <-
function(n=20, beta=c(3, 1.5, 0, 0, 2, 0, 0, 0), rho=0.5, sig=1) {
    stopifnot(length(beta)>0, abs(rho)<1, sig>0)
    p <- length(beta)
    z <- t(chol(toeplitz(rho^(0:(p-1)))))
    X <- tcrossprod(matrix(rnorm(n*p), nrow=n, ncol=p), z)
    Xy <- cbind(X,(X%*%beta) + rnorm(n, mean=0, sd=sig))
    dimnames(Xy)[[2]] <- c(paste0("x",1:p), "y")
    as.data.frame.matrix(Xy)
  }

theta.HBSYC <-
function(para, y, X1, R){
     p <- ncol(X1)
     m <- nrow(X1)
     beta <- para[1:p]
     Sqsigmav <- para[p+1]
     Sqsigma <- para[(p+2):(p+m+1)]
     lambda <- para[p+m+2]
     diag(R) <- (diag(R) * lambda + 1 - lambda) / Sqsigmav + 1 / Sqsigma
     R[R == -1] <- -lambda / Sqsigmav
     theta_HB <- solve(R, (y - X1 %*% beta) / Sqsigma)
     theta_HB <- theta_HB + X1 %*% beta
}

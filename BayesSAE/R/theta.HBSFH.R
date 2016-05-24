theta.HBSFH <-
function(para, y, X1, R, vardir){
     p <- ncol(X1)
     beta <- para[1:p]
     Sqsigmav <- para[p+1]
     lambda <- para[p+2]
     diag(R) <- (diag(R) * lambda + 1- lambda) / Sqsigmav + 1 / vardir
     R[R == -1] <- -lambda / Sqsigmav
     theta_HB <- solve(R, (y - X1 %*% beta) / vardir)
     theta_HB <- theta_HB + X1 %*% beta
}

#obtain the unbiased estimator of parameters from a MVN dist.
# input
#      X: samples, one row per observation
#note that the MLE of \Sigma is "(n-1)/n * hatSigma"
mvn.ub <- function(X){
    d <- dim(X)
    n <- d[1]
    p <- d[2]
    barX <- colMeans(X)
    S <- rowSums(apply(X, 1, function(x) (x - barX) %*% t(x - barX)))
    S <- matrix(S, nrow=p, ncol=p)
    if(n > 1){
        hatSigma <- S/(n-1)
    }
    else{
        hatSigma <- 0
    }
    list(hatMu=barX, hatSigma=hatSigma)
}

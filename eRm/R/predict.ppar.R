predict.ppar <- function(object, cutpoint = "randomized", ...)
{
# predict method for objects of class ppar
# cutpoint ... either value between 0 and 1, or randomized assignment

Pi <- pmat(object)                            #expected values
X <- object$X.ex
if (max(X, na.rm = TRUE) > 1) stop("Available for dichotomous models only!")

K <- dim(X)[2]
N <- dim(X)[1]

y <- as.vector(t(X))   #observed values
pi.hat <- as.vector(t(Pi))

if (cutpoint == "randomized") {
   pvec <- runif(length(y))
} else {
  pvec <- rep(cutpoint, length(y))
}

classvec <- (pvec < pi.hat)*1                #expected 0/1 vector
classmat <- matrix(classvec, ncol = K, nrow = N, byrow = TRUE)
dimnames(classmat) <- list(rownames(X), colnames(X))
return(classmat)
}

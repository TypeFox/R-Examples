`hyperplane` <-
function(X){
    if (!is.matrix(X)) stop("'X' must be a k x k matrix")
    k <- ncol(X)
    if (k != nrow(X)) stop("'X' must be a k x k matrix")
    d0 <- det(X)
    if (isTRUE( all.equal(d0,0))){
      d <- hyperplaneOrthogonalVector(X)
    }
    else{
      b <- -rep(d0,times=k)
      d <- solve(X,b)
    }
    d <- c(d,d0)
    return(d)
}

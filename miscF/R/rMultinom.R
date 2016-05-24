rMultinom <- function(p){
    if(any(p < 0))
        stop("non-positive probability")
    if(!isTRUE(all.equal(rowSums(p),rep(1, nrow(p)))))
        stop("The sums of probabilities have to be 1.")
    k <- ncol(p)
    n <- nrow(p)
    csm <- ifelse(outer(1:k,1:(k - 1),FUN=`<=`), 1, 0)
    cp <- p %*% csm
    U <- runif(n)
    rowSums(sweep(cp, 1, U, `<=`)) + 1
}

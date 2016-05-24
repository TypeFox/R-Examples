HVK <- function(X, m1=NULL, m2=NULL, ar.order=1) {
    if (!is.numeric(X) | !is.vector(X)) {
        stop("input object should be a vector.")
    }    
    if (any(is.na(X))) {
        stop("input vector should not contain missing values.")
    }
    n <- length(X)
    if(is.null(m1)|is.null(m2)){
        m1 <- round(n^(0.1))
        m2 <- round(n^(0.5))
    }
    m <- c(m1:m2)
    tmp <- sapply(m, function(x) diff(X, lag=x))
    tmp <- sapply(1:length(tmp), function(x) sum(tmp[[x]]^2))
    autocovar0 <- sum(tmp/(2*(n-m)))/(m2-m1+1)
    tmp <- sapply(c(1:ar.order), function(x) diff(X, lag=x))
    if(ar.order==1){
        tmp <- sum(tmp^2)
    }else{
        tmp <- sapply(1:length(tmp), function(x) sum(tmp[[x]]^2))
    }
    autocovarj <- autocovar0 - tmp/(2*(n-c(1:ar.order)))
    autocovar <- c(autocovar0, autocovarj)
    G <- matrix(NA, ar.order, ar.order)
    for (i in 1:(ar.order)){
        for (j in 1:(ar.order)){
            G[i,j] <- autocovar[abs(i-j)+1]
        }
    }
    coeffs <- solve(G)%*%autocovar[2:(ar.order+1)]
    return(as.vector(coeffs))
}
## copied from robustbase
ghq <- function(n = 1, modify = TRUE) {
    ## Adapted from gauss.quad in statmod package
    ## which itself has been adapted from Netlib routine gaussq.f
    ## Gordon Smyth, Walter and Eliza Hall Institute

    n <- as.integer(n)
    if(n<0) stop("need non-negative number of nodes")
    if(n==0) return(list(nodes=numeric(0), weights=numeric(0)))
    ## i <- seq_len(n) # 1 .. n
    i1 <- seq_len(n-1L)

    muzero <- sqrt(pi)
    ## a <- numeric(n)
    b <- sqrt(i1/2)

    A <- numeric(n*n)
    ## A[(n+1)*(i-1)+1] <- a # already 0
    A[(n+1)*(i1-1)+2] <- b
    A[(n+1)*i1] <- b
    dim(A) <- c(n,n)
    vd <- eigen(A,symmetric=TRUE)
    n..1 <- n:1L
    w <- vd$vectors[1, n..1]
    w <- muzero * w^2
    x <- vd$values[n..1] # = rev(..)
    list(nodes=x, weights= if (modify) w*exp(x^2) else w)
}

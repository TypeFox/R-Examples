############
# moran.idx
############
moran.idx <- function(x, prox, addInfo=FALSE){

    ## handle arguments
    if(any(is.na(x))) stop("NA entries in x")
    if(!is.numeric(x)) stop("x is not numeric")

    W <- as.matrix(prox)
    if(!is.matrix(W)) stop("prox is not a matrix")
    if(ncol(W) != nrow(W)) stop("prox is not a square matrix")
    if(any(is.na(W))) stop("NA entries in prox")
    diag(W) <- 0

    n <- nrow(W)


    ## main computations
    x <- x - mean(x)
    sumW <- sum(W)
    num <- n * sum(x * (W %*% x) )
    denom <- sumW * sum(x*x)

    if(denom < 1e-14) stop("denominator equals zero")

    res <- num/denom

    if(addInfo){
        I0 <- -1/(n-1)
        matToDiag <- .5 * (t(W) + W)
        rangeI <- range(eigen(matToDiag)$values)
        attr(res, "I0") <- I0
        attr(res, "Imin") <- rangeI[1]
        attr(res, "Imax") <- rangeI[2]
    }

    return(res)

} # end moran.idx

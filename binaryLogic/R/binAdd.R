#' Binary Addition (+)
#' 
#' @description Adds two binary numbers. (x + y)
#' @details Little-Endian and unsigned is not supported at the moment. No floating point supported.
#' if x or y is signed the return value will also be signed.
#' @usage binAdd(x, y)
#' @param x summand 1 (binary vector)
#' @param y summand 2 (binary vector)
#' @return The sum of x and y. Returns a binary vector.
#' @examples
#' five <- as.binary(5); ten <- as.binary(10);
#' as.numeric(binAdd(ten, five))
#' binAdd(as.binary(c(0,1), logic=TRUE), as.binary(c(1,0), logic=TRUE))
#' @seealso base::as.logical , base::is.logical, base::raw
#' @export
binAdd <- function(x, y) {
    #Should be implemented in C (internal generic). Very slow.
    if (missing(x)) stop("x is missing.")
    if (missing(y)) stop("y is missing.")
    stopifnot(is.binary(x))
    stopifnot(is.binary(y))
    
    x_signed <- attributes(x)$signed
    x_littleEndian <- attributes(x)$littleEndian
    y_signed <- attributes(y)$signed
    y_littleEndian <- attributes(y)$littleEndian
    
    l <- list(class=c("binary","logical"),
              signed=FALSE,
              littleEndian=x_littleEndian)
    
    if (x_signed || y_signed) l$signed <- TRUE
    if (x_littleEndian) y <- switchEndianess(y)
    if (y_littleEndian) x <- switchEndianess(x)    
    
    if (length(x) >= length(y))
    {
        MAX <- length(x)
        a <- rep(FALSE,length(x)-length(y))
        y <- c(a,y)
    } else {
        MAX <- length(y)
        a <- rep(FALSE,length(y)-length(x))
        x <- c(a,x)
    }
    
    ret = binary(MAX)
    temp = binary(MAX+1)
    ret[MAX] <- xor(x[MAX],y[MAX])
    if ((isTRUE(as.logical(x[MAX])) && isTRUE(as.logical(y[MAX])))) temp[MAX+1] <- TRUE
    
    if (MAX > 2)
    {
        for(i in (MAX-1):1)
        {
            ret[i] <- xor(x[i],y[i])
            ret[i] <- xor(ret[i],temp[i+2])
            if(((isTRUE(as.logical(x[i])) && isTRUE(as.logical(y[i]))) || 
                   (isTRUE(as.logical(x[i])) && isTRUE(as.logical(temp[i+2]))) ||
                   (isTRUE(as.logical(y[i])) && isTRUE(as.logical(temp[i+2]))))) temp[i+1] <- TRUE
        }
    }
    if (temp[2] && !l$signed) {
        ret <- c(T,ret)
        #l$signed <- TRUE
        #l$littleEndian <- FALSE
    }
    # reverse big endian.
    return(loadAttributes(ret, l))
}
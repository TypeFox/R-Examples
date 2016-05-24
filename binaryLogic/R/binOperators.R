#' Binary Negation (!)
#' 
#' @description Negates the binary number x. Negation x -> -x or -x -> x
#' @details An »unsigned« number will be returned as »signed« regardless of whether the value is negative. 
#' No floating point supported.
#' @usage negate(x)
#' @param x The number to be negated. A binary vector is expected.
#' @return The negated number of x. Returns a binary vector with signed=TRUE
#' @examples
#' summary(negate(as.binary(5, signed=TRUE)))
#' summary(negate(as.binary(-5, signed=TRUE)))
#' summary(negate(as.binary(5, signed=FALSE)))
#' @seealso \link{switchEndianess} or \link{fillUpToByte}.
#' @export
negate <- function(x) {
    stopifnot(is.binary(x))
    signed <- attributes(x)$signed
    if (all(!x)) return(binary(length(x)))
    #if (!signed) warning("Trying to negate an unsigned digit. treated as signed value. Returns a signed value")
    littleEndian <- attributes(x)$littleEndian

    if (littleEndian) x <- rev(x)
    
    if (length(x) %% byte() != 0) {
        MAX <- (trunc((length(x)/byte())) +1) * byte()
        a <- rep(FALSE, MAX - length(x))
        a <- as.binary(a, littleEndian=littleEndian, logic=TRUE)
        # 'c.binary'(x,y) needs to be implemented
        x <- as.binary(c(a,x), logic=TRUE)
    }
    x <- !x
    x <- binAdd(as.binary(x, logic=TRUE),as.binary(TRUE, logic=TRUE))

    if (littleEndian) x <- rev(x)
    attr(x,"signed") <- TRUE
    attr(x,"littleEndian") <- littleEndian
    class(x) <- c("binary", "logical")
    return(x)
}

#' Binary Left Shift (<<)
#' 
#' @description Logical left shift x << n
#' @usage shiftLeft(x, n)
#' @param x The binary number to shift. (binary or logical vector).
#' @param n The number of bits to shift.
#' @return Pushes 0's(FALSE) to the vector from right(LSB) to left(MSB).
#' Everything on right(MSB) side drops out. Returns a binary/logical vector
#' @examples
#' x <- as.binary(c(1,0,0,1,1,1,0,1), logic=TRUE); x
#' shiftLeft(x,1)
#' shiftLeft(x,2)
#' @seealso \link{shiftRight} and \link{rotate} 
#' @export
shiftLeft <- function(x, n) {
    stopifnot(is.logical(x) || is.binary(x))
    stopifnot(n > 0)
    if (n > length(x)) if (class(x)[1] == "binary") return(binary(length(x))) else return(logical(length(x)))
    l <- saveAttributes(x)
    loadAttributes(c(x[-seq(n)], logical(n)), l)
}

#' Binary Right Shift (>>)
#' 
#' @description Logical right shift 1 >> n
#' @usage shiftRight(x, n)
#' @param x The binary number to shift. (binary or logical vector).
#' @param n The number of bits to shift.
#' @return Pushes 0's(FALSE) to the vector from left(MSB) to right(LSB).
#' Everything on right(LSB) side drops out. Returns a binary/logical vector
#' @examples
#' x <- as.binary(c(1,0,0,1,1,1,0,1), logic=TRUE); x
#' shiftRight(x,1)
#' shiftRight(x,2)
#' @seealso \link{shiftLeft} and \link{rotate} 
#' @export
shiftRight <- function(x, n) {
    stopifnot(is.logical(x) || is.binary(x))
    stopifnot(n > 0)
    if (n > length(x)) if (class(x)[1] == "binary") return(binary(length(x))) else return(logical(length(x)))    
    l <- saveAttributes(x)
    loadAttributes(rev(c(rev(x)[-seq(n)], logical(n))), l)
}

#' Rotate no carry ()
#'
#' @description A circular shift
#' @usage rotate(x, n)
#' @param x The binary number to rotate. (binary or logical vector).
#' @param n The number of bits to rotate.
#' @return rotates the vector from left to right. 
#' The value from MSB is used to fill up the vector at LSB. Returns a binary/logical vector.
#' @examples
#' x <- as.binary(c(1,0,0,1,1,1,0,1), logic=TRUE); x
#' rotate(x,1)
#' rotate(x,2)
#' @seealso \link{shiftLeft} and \link{shiftRight} 
#' @export
rotate <- function(x, n) {
    stopifnot(is.binary(x) || is.logical(x))
    stopifnot(n > 0)
    #save and loadAttributes needs to be done, until I find a way to overload »combine« c(...)
    l <- saveAttributes(x)
    if (n > length(x)) stop("n is larger than length of x")
    loadAttributes(c(x[-seq(n)],x[seq(n)]), l)
}

#' Fill up to Byte (00000000..)
#'
#' @description Fills up the binary number with zeros(0) or ones(1), to the size in Byte.
#' @details No floating point supported.
#' @usage fillUpToByte(x, size=0, value=FALSE)
#' @param x The binary number to fill up with zeros. (Any binary vector).
#' @param value to fill up with FALSE(0) or fill up with TRUE(1).
#' @param size in Byte. 0 = auto (smallest possible Byte).
#' @return binary number. A binary vector with the desired size.
#' @examples
#' fillUpToByte(as.binary(c(1,1), logic=TRUE), size=2)
#' fillUpToByte(as.binary(c(1,0,1), logic=TRUE), size=2, value=FALSE)
#' @seealso \link{fillUpToBit} or \link{bytesNeeded}, \link{negate}, \link{switchEndianess}.
#' @export 
fillUpToByte <- function(x, size=0, value=FALSE) {
    stopifnot(is.binary(x))
    l <- saveAttributes(x)
    if (size == 0 && length(x) %% byte() == 0) return(x)
    if (size > 0 && length(x) >= size * byte()) return(x)

    if (size == 0) {
        append <- binary(((trunc((length(x) / byte())) +1) * byte()) - length(x))
    } else {
        append <- binary(size * byte() - length(x))
    }
    append[1:length(append)] <- value

    if (attributes(x)$littleEndian) {
        x <- c(x,append)
    } else {
        x <- c(append,x)
    }
    x <- loadAttributes(x,l)
    return(x)
}

#' Fill up to bit (000..)
#'
#' @description Fills up the binary number with zeros(0) or ones(1), to the size n in bit.
#' @details No floating point supported.
#' @usage fillUpToBit(x, n, value=FALSE)
#' @param x The binary number to fill up with zeros. (Any binary vector).
#' @param value to fill up with FALSE(0) or fill up with TRUE(1).
#' @param n size in bit.
#' @return binary number. A binary vector with the desired size.
#' @examples
#' fillUpToBit(as.binary(c(1,1), logic=TRUE), n=4)
#' fillUpToBit(as.binary(c(1,0,1), logic=TRUE), n=4, value=FALSE)
#' @seealso \link{fillUpToByte}.
#' @export
fillUpToBit <- function(x, n, value=FALSE) {
    stopifnot(is.binary(x))
    if(missing(n)) stop("n is missing")
    stopifnot(n >= 0)
    
    l <- saveAttributes(x)    
    if (length(x) >= n) return(x)
    
    append <- binary(n - length(x))
    append[seq(length(append))] <- value
    
    if (attributes(x)$littleEndian) {
        x <- c(x,append)
    } else {
        x <- c(append,x)
    }
    x <- loadAttributes(x,l)
    return(x)
}

#' Switch Endianess.
#' 
#' @description Switch little-endian to big-endian and vice versa.
#' @usage switchEndianess(x, stickyBits=FALSE)
#' @param x binary number. Any binary number.
#' @param stickyBits Bits wont change if set TRUE. Only the attribute will be switched. 
#' @return switch little-endian to big-endian and vice versa.
#' @examples
#' x <- as.binary(c(1,1,0,0), logic=TRUE); print(x); summary(x);
#' y <- switchEndianess(x); print(y); summary(y);
#' y <- switchEndianess(x, stickyBits=TRUE); print(y); summary(y);
#' @seealso \link{negate} or \link{fillUpToByte}.
#' @export
switchEndianess <- function(x, stickyBits=FALSE) {
    stopifnot(is.binary(x))
    l <- saveAttributes(x)
    l$littleEndian <- !l$littleEndian
    if (stickyBits == TRUE)
    {
        return(loadAttributes(x,l))
    } else {
        return(loadAttributes(rev(x),l))
    }
}

#' Binary sequence
#' 
#' @description Binary sequence.
#' @usage binSeq(x, ...)
#' @param x a sequence.
#' @param ... used for dec2bin().
#' @return a sequence list of binary digits.
#' @examples
#' binSeq(0:4)
#' @seealso \link{binary}
#' @export
binSeq <- function(x, ...) {
    l <- vector("list", length(x))
    for(i in 1:length(x)) l[[i]] <- dec2bin(x[i], ...)
    if(length(l) == 1){
        return(l[[1]])
    } else {
        return(l)
    }
}

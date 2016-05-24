#' Binary prefix (KiB,MiB,..)
#'
#' @description Num of byte needed to fit in n * KiB, MiB ..etc.
#' @details
#' KiB <- KibiByte
#' MiB <- MebiByte
#' GiB <- GibiByte
#' TiB <- TebiByte
#' PiB <- PebiByte
#' EiB <- ExiByte
#' ZiB <- ZebiByte
#' YiB <- YobiByte 
#' @usage binaryPrefix(n, prefix="KiB")
#' @param n numeric value
#' @param prefix binary prefix * byte. Expeting a »string«
#' @return The number of byte fitting in n * binary prefix * byte
#' @examples
#' #Get the number of byte needed to hold 0.5 and 1:10 KiB
#' binaryPrefix(c(0.5,1:10),"KiB")
#' #Get the number of bit needed to hold 1 KiB
#' binaryPrefix(1,"KiB")*byte()
#' @seealso \link{bytesNeeded} or \link{fillUpToByte} or \link{byte}
#' @export
binaryPrefix <- function(n, prefix="KiB") {
    stopifnot(all(is.numeric(n) || is.na(n)))
    stopifnot(all(n >= 0 || is.na(n)))
    ret <- NULL
    
    prefix <- switch(prefix,
                KiB = 2^10,
                MiB = 2^20,
                GiB = 2^30,
                TiB = 2^40,
                PiB = 2^50,
                EiB = 2^60,
                ZiB = 2^70,
                YiB = 2^80)
    
    ret = n * prefix
    if (any(is.null(ret))) stop("Unknown binary prefix")
    if (isTRUE(all.equal(ret, round(ret)))) {
        return(ret)
    } else {
        message("ceiling called: returns the smallest integer not less than the corresponding element of ret")
        return(ceiling(ret))
    }
}

#' Minimum number of "byte" needed to hold n "bit"
#' 
#' @description A simple helper function
#' that returns the minimum number of byte needed to hold the amount of n bit.
#' @usage bytesNeeded(n)
#' @param n The number of bit.
#' @return The number of minimum byte needed to hold n bit.
#' @examples
#' ten <- as.binary(10)
#' bytesNeeded(length(ten))
#' @seealso \link{fillUpToByte} or \link{binaryPrefix} or \link{byte}
#' @export
bytesNeeded <- function(n) {
    stopifnot(all(is.numeric(n) || is.na(n)))
    stopifnot(all(n >= 0 || is.na(n)))
    
    ifelse(n %% byte() == 0, n %/% byte(), n %/% byte() + 1)
}

#' A simple helper function to return the size of one byte
#' 
#' @description Used to increase readabilaty
#' @usage byte()
#' @return The size of one byte (8)
#' @seealso \link{bytesNeeded} or \link{fillUpToByte} or \link{binaryPrefix} 
#' @export
byte <- function() {
    return(8)
}

#' A gray code converter function
#' 
#' @description This function converts a binary number (base2) to a gray code
#' @usage bin2gray(x)
#' @param x The binary number (base2) or a logical vector.
#' @return The gray code as logical vector.
#' @seealso \link{gray2bin}
#' @export
bin2gray <- function(x) {
    stopifnot(is.logical(x) || is.binary(x))
    max <- length(x)
    g <- logical(max)
    if (max == 1) return(as.logical(x))

    if(is.binary(x) && attributes(x)$littleEndian == TRUE)
    {
        g[max] <- x[max]
        for(i in (max):2) g[i - 1] <- xor(x[i], x[i - 1])
    } else {
        g[1] <- x[1]
        for(i in 1:(max - 1)) g[i + 1] <- xor(x[i], x[i + 1])
    }
    return(g)
}

#' A gray code to binary converter function
#'
#' @description This function converts a gray code to a binary number (base2)
#' @usage gray2bin(x, ...)
#' @param x The gray code as logical vector.
#' @param ... Additional parameter for binary() 
#' @return The binary number (base2).
#' @seealso \link{bin2gray}
#' @export
gray2bin <- function(x, ...) {
    stopifnot(is.logical(x))
    max <- length(x)
    b <- binary(max, ...)
    if (max == 1) return(as.binary(x, logic=TRUE, ...))

    if(is.binary(b) && attributes(b)$littleEndian == TRUE)
    {
        b[max] <- x[max]
        for(i in (max):2) b[i - 1] <- xor(b[i], x[i - 1])
    } else {
        b[1] <- x[1]
        for(i in 1:(max - 1)) b[i + 1] <- xor(b[i], x[i + 1])
    }
    return(b)
}
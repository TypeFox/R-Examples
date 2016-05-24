#' Convert decimal values to binary values
#' 
#' A function to convert decimal values to binary values
#' @param number the decimal number that needs to be converted.
#' @param bits the number of bits of the result
#' 
#' @return a vector contains of binary values 0 and 1.
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' dec2bin(8, 5)

dec2bin <- function(number, bits) {
    if (bits > 1) {
        return(t(sapply(number, function(x) {
            floor(2^c(1 - bits:1) * x)%%2
        })))
    } else {
        return(t(t(sapply(number, function(x) {
            floor(2^c(1 - bits:1) * x)%%2
        }))))
    }
    
} 

#' Generate gray code
#' 
#' A function to generate gray code
#' 
#' @param n an integer to specify the number of bits.
#' @return a vector of the decimal gray code. 
#' @references Dah jyh Guan. (Scientific Note) Generalized Gray Codes with Applications. 1998
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' code<-grays(3)

grays <- function(n) {
    gray_code <- rep(0, 2^n)
    
    # n=1 [0 1]
    gray_code[2] <- 1
    if (n == 1) {
        return(gray_code)
    } else {
        t <- 2
        # for n lager than 2
        for (i in 2:n) {
            t2 <- t + t
            gray_code[(t + 1):t2] <- t + rev(gray_code[1:t])
            t <- t2
        }
        return(gray_code)
    }
} 

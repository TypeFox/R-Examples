#' Graycode Function
#' 
#' A function to generate decimal graycode
#' 
#' @param a the number of targets
#' 
#' @return A list contains the following components:
#' \item{rows}{the number of rows}
#' \item{cols}{the number of columns}
#' \item{dec}{the decimal graycode results}
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @references Dah jyh Guan. (Scientific Note) Generalized Gray Codes with Applications. 1998
#' @examples 
#' code<-graycode2(5)
graycode2 <- function(a) {
    # parameter1 a: the number of targets
    
    rows <- 2^floor(a/2)
    cols <- 2^(a - floor(a/2))
    # get the decimal
    dec <- grays(a)
    # reshape the index to a matrix
    dim(dec) <- c(cols, rows)
    dec <- t(dec)
    # flip every second row of the matrix
    nrow_dec <- nrow(dec)
    if (nrow_dec == 1) {
        return(list(rows, cols, dec))
    } else if (nrow_dec == 2) {
        dec[2, ] <- rev(dec[2, ])
        return(list(rows, cols, dec))
    } else {
        dec[seq(2, nrow_dec, 2), ] <- t(apply(dec[seq(2, nrow_dec, 2), ], 1, rev))  #does not work for nrow<=2
        return(list(rows, cols, dec))
    }
    
} 

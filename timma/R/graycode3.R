#' Gray code function for matrix indexes
#' 
#' A function to generate gray code used for matrix row and column names
#' 
#' @param m an integer to specify the number of bits
#' 
#' @return a list of the following components:
#' \item{gc_row}{binary gray code as row names of the predicted sensitivity matrix}
#' \item{gc_col}{binary gray code as column names of the predicted sensitivity matrix}
#' \item{dec_row}{decimal gray code as row names of the predicted sensitivity matrix}
#' \item{dec_col}{decimal gray code as column names of the predicted sensitivity matrix}
#' 
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' names<-graycode3(3)

graycode3 <- function(m) {
    if (m == 1) {
        gc_row <- c(0, 1)
        gc_col <- c(0, 1)
        dec_row <- gc_row
        dec_col <- gc_col
    } else {
        rows <- floor(m/2)
        cols <- m - rows
        dec_row <- grays(rows)
        dec_col <- grays(cols)
        gc_row <- dec2bin(dec_row, rows)
        gc_col <- dec2bin(dec_col, cols)
    }
    return(list(gc_row = gc_row, gc_col = gc_col, dec_row = dec_row, dec_col = dec_col))
} 

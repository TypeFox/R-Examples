#' Names for the predicted sensitivity matrix
#' 
#' A function to make the target names in the format of gray code for the predected sensitivity matrix
#' 
#' @param m an integer to specify the number of targets
#' @param names a vector of the names of the targets
#' @param gc_row the gray code as row indexes. It can be returned by \code{\link{graycode3}}.
#' @param gc_col the gray code as column indexes. It can be returned by \code{\link{graycode3}}.
#' 
#' @return a list of the following components:
#' \item{nr}{the gray code format target names as row names.}
#' \item{nc}{the gray code format target names as row names.}
#' 
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' k_select<-sffsBinary(tyner_interaction_binary, tyner_sensitivity[, 1])$k_sel
#' gc_timma<-graycode3(length(k_select))
#' select_kinase_names<-dimnames(tyner_interaction_binary)[[2]][k_select]
#' gc_names<-graycodeNames(length(k_select), select_kinase_names, gc_timma$gc_row, gc_timma$gc_col)
#' }
graycodeNames <- function(m, names, gc_row, gc_col) {
    # parameter m: the number of how many selected kinases    
    # parameter names: the names for selected kinases
    dim_row <- dim(gc_row)
    dim_col <- dim(gc_col)
    
    # get the kinase names for rows
    names_row <- names[1:dim_row[2]]
    # get the kinase names for columns
    names_col <- names[dim_row[2] + 1:m]
    
    n_r <- array(NA, dim = dim_row)
    n_c <- array(NA, dim = dim_col)
    for (i in 1:dim_row[1]) {
        for (j in 1:dim_row[2]) {
            if (gc_row[i, j] == 0) {
                n_r[i, j] <- "-"
            } else {
                n_r[i, j] <- names_row[j]
            }
        }
    }
    for (i in 1:dim_col[1]) {
        for (j in 1:dim_col[2]) {
            if (gc_col[i, j] == 0) {
                n_c[i, j] <- "-"
            } else {
                n_c[i, j] <- names_col[j]
            }
        }
    }
    return(list(nr = n_r, nc = n_c))
} 

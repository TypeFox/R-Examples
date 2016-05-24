#' Find the same column from a matrix
#' 
#' A function to seek for the same column from a matrix
#' 
#' @param X a matrix 
#' @param Y a vector with the same length as each column in X. 
#' @return a vector of the column indexes which are the same as vector Y.
#' 
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' x<-data.frame(tyner_interaction_binary)
#' kinase_names<-dimnames(tyner_interaction_binary)
#' float<-sffsBinary(tyner_interaction_binary, tyner_sensitivity[,1])
#' k_select <- float$k_sel
#' select_kinase_names <- findSameSet(x, k_select, kinase_names)
#' }
findSameCol <- function(X, Y) {
    # parameter 1: X is a matrix 
    # parameter 2: Y is a vector with the same length as each column in X
    col_num <- ncol(X)
    same_set <- rep(0, col_num)
    for (i in 1:col_num) {
        if (all(X[, i] == Y)) {
            same_set[i] <- 1
        }
    }
    return(same_set)
} 

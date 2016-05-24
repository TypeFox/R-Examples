#' Number of edits
#' Count the number of edits in a collection of edits.
#' @param E \code{\link{editset}}, \code{\link{editarray}} or \code{\link{editmatrix}}
#' @export
nedits <- function(E){
    if ( is.vector(E) ) return(length(E))
    if (any(class(E) %in% c('editmatrix','editarray'))){ 
        n <- nrow(E)
    } else if ( inherits(E,'editset') ) {
        n <- nrow(E$num)  + nrow(E$mixcat)
    } else {
        stop('Argument must be character, editset, editarray, editmatrix')
    }
    n
}

#' Names of edits
#'
#' Retrieve edit names from editset, -array or -matrix
#' @param E \code{\link{editset}}, \code{\link{editarray}} or \code{\link{editmatrix}}
#' @export
editnames <- function(E){
    
    if (any(class(E) %in% c('editmatrix','editarray'))){ 
        n <- rownames(E)
    } else if ( inherits(E,'editset') ) {
        n <- c(rownames(E$num),rownames(E$mixcat))
    } else {
        stop('Argument must be editset, editarray or editmatrix')
    }
    n

}







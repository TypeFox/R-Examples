#' Returns the derivation history of an edit matrix or array
#'
#' Function \code{\link{eliminate}} tracks the history of edits in a logical array H.
#' H has nrow(E) rows and the number of columns is the number of
#' edits in the \code{\link{editmatrix}} as it was first defined. If 
#' H[i,j1], H[i,j2],...,H[i,jn] are \code{TRUE}, then E[i,] is some 
#' (positive, linear) combination of original edits E[j1,], E[j2,],...,E[jn,]
#'
#' Attributes H and h are used to detect redundant derived edits.
#'
#' @param E \code{\link{editmatrix}}
#' @rdname geth
#' @seealso \code{\link{editmatrix}}, \code{\link{eliminate}}
#'
#'
#' @export
getH <- function(E){
    if ( !class(E) %in% c('editmatrix','editarray') ) 
        stop("E has to be an editmatrix or editarray")
    attr(E,"H")  
}

#' Returns the number of elimination steps performed on an edit matrix or array
#'
#' h records the number of variables eliminated from E by \code{\link{eliminate}}
#'
#' @rdname geth
#' @seealso \code{\link{editmatrix}}, \code{\link{eliminate}}
#' @export
geth <- function(E){
    if ( !class(E) %in% c('editmatrix','editarray') ) 
        stop("E has to be an editmatrix or editarray")
    attr(E,"h")  
}


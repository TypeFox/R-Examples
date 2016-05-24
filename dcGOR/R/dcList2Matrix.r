#' Function to convert a list into a matrix containing three columns
#'
#' \code{dcList2Matrix} is supposed to convert a list into a matrix containing three columns
#'
#' @param x a list, its each component must be a named vector
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' a matrix containing three columns: 1st for the input list names (if exist, otherises an increasing integer), 2nd for the vector names of each list component, and 3rd for the vector values of each list component
#' @note
#' none
#' @export
#' @seealso \code{\link{dcAlgoPropagate}}
#' @include dcList2Matrix.r
#' @examples
#' \dontrun{
#' # load an object 'HIS'
#' Feature2GOMF.sf <- dcRDataLoader(RData='Feature2GOMF.sf')
#' # get a list
#' x <- Feature2GOMF.sf$hscore
#' # convert the list into a matrix
#' res <- dcList2Matrix(x)
#' dim(res)
#' res[1:10,]
#' }

dcList2Matrix <- function(x, verbose=T)
{
    
    if(!is.list(x)){
        stop("The input 'x' must be a list!\n")
    }else{
        tmp <- x[[1]]
        if(!is.null(tmp)){
            if (!is.vector(tmp) | is.null(names(tmp))){
                stop("The components of input list 'x' must be a named vector!\n")
            }
        }
    }
    
    if(is.null(names(x))){
        names(x) <- 1:length(x)
    }
    
    x_names <- names(x)
    output_list <- lapply(1:length(x), function(i){
        tmp <- x[[i]]
        cbind(C1=rep(x_names[i],length(tmp)), C2=names(tmp), C3=tmp)
    })
    x_mat <- base::do.call(base::rbind, output_list)
    rownames(x_mat) <- NULL
    colnames(x_mat) <- NULL

    if(verbose){
        message(sprintf("The input list has been converted into a matrix of %d X %d.", dim(x_mat)[1], dim(x_mat)[2]), appendLF=T)
    }
    
    invisible(x_mat)
}

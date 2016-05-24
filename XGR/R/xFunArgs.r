#' Function to assign (and evaluate) arguments with default values for a given function
#'
#' \code{xFunArgs} is supposed to assign (and evaluate) arguments with default values for a given function.
#'
#' @param fun character specifying the name of the function
#' @param action logical to indicate whether the function will act as it should be (with assigned values in the current environment). By default, it sets to FALSE
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to TRUE for display
#' @return
#' a list containing arguments and their default values
#' @note
#' This function is potentially useful when debugging as it frees developers from specifying default values for all arguments except those arguments of interest
#' @export
#' @seealso \code{\link{xFunArgs}}
#' @include xFunArgs.r
#' @examples
#' fun <- "xRDataLoader"
#' xFunArgs(fun)

xFunArgs <- function(fun, action=F, verbose=TRUE)
{
    
    args_list <- base::formals(fun)
    args_names <- names(args_list)
    for(i in 1:length(args_list)){
        lft <- args_names[[i]]
        rgt <- paste(base::deparse(args_list[[i]]),collapse='')
        if(rgt!=''){
            tmp <- paste(lft, '<-', rgt, sep=' ')
            if(action==T){
                base::eval(base::parse(text=tmp), envir=parent.frame())
            }else{
                base::eval(base::parse(text=tmp))
            }
        }
    }
    
    if(verbose){
        if(action==T){
            message(sprintf("For the function '%s', %d arguments have been assigned with default values in the current environment.", fun, length(args_names)), appendLF=T)
        }else{
            message(sprintf("For the function '%s', there are %d arguments.", fun, length(args_names)), appendLF=T)
        }
    }
    
    invisible(args_list)
}

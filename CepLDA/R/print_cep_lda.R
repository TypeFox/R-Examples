#' @name print.ceplda
#' @aliases print.ceplda
#' @title Print Method for Class 'ceplda'
#' 
#' @description
#' Print the results from \code{cep.lda}
#' @param x an object of class "ceplda".
#' @param ... additional arguments.
#' @method print ceplda
#' @export
#' 
#' @seealso
#'  \code{\link{cep.lda}} 


print.ceplda <- function(x,...){
        if(!inherits(x,"ceplda"))
                stop("Object must be of class 'ceplda'")
        cat("Optimal L selected: \n")
        print(x$Lopt)
        cat("\n")
        cat("Linear Discriminat Analysis results: \n")
        print(x$C.lda,...)
}
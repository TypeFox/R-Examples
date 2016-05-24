#================================================================
#' @title Print summaries of \code{meteRelat} objects
#'
#' @description S3 method for class \code{meteRelat}
#'
# @details
#' 
#' @param x an object of class meteRelat
#' @param ... arguments to be passed to methods
#' 
#' @export
#' 
# @examples
#'
#' @return x silently
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.
#' @importFrom graphics plot points


## print method for objects of class `meteRelat'

print.meteRelat <- function(x,...) {
    cat(switch(attr(x$pred, 'type'), 
               'sar' = 'Species area relationship',
               'ear' = 'Endemcis area relationship',
               'damuth' = 'Abundance metabolic rate relationship'),
        sprintf('predicted using %s', ifelse(is.null(x$obs),
                                             'state variables only',
                                             'raw data')),
        '\n')

    print(x$pred)
    if(!is.null(x$obs)) print(x$obs)
    
    invisible(x)
}
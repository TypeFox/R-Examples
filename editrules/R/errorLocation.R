#' The errorLocation object
#'
#' Object storing information on error locations in a dataset.
#' 
#' The \code{errorlocation} objects consists of the following slots wich can be 
#' accessed with the dollar operator, just like with lists. Right now the only
#' functions creating such objects are \code{\link{localizeErrors}} and \code{\link{checkDatamodel}}.
#'
#' \itemize{
#'      \item{\code{adapt} a \code{logical} array where each row/column shows which record/variable should be adapted.}
#'      \item{\code{status} A \code{data.frame} with the same number of rows as \code{adapt}. It contains the following
#'          columns
#'          \itemize{
#'              \item{\code{weight} weight of the found solution}
#'              \item{\code{degeneracy} number of equivalent solutions found}
#'              \item{\code{user} user time used to generate solution (as in \code{sys.time})}
#'              \item{\code{system} system time used to generate solution (as in \code{sys.time})}
#'              \item{\code{elapsed} elapsed time used to generate solution (as in \code{sys.time})}
#'              \item{\code{maxDurationExceeded} Was the maximum search time reached?}
#'              \item{\code{memfail} Indicates whether a branch was broken off due to memory allocation failure (branch and bound only)}
#'          }
#'      }
#'      \item{\code{method} The error localization method used, can be "mip", "localizer" or "checkDatamodel".}
#'      \item{\code{call} The R calls to the function generating the object.}
#'      \item{\code{user} \code{character} user who generated the object.}
#'      \item{\code{timestamp} \code{character} timestamp.}
#' }
#'
#' It is possible to \code{plot} objects of class \code{errorLocation}. An overview containing
#' three or four graphs will be plotted in a new window. Axes in scatterplots are set to logarithmic
#' if their scales maxima exceed 50.
#'
#' @rdname errorLocation
#' @name errorLocation
#' @seealso \code{\link{localizeErrors}}, \code{\link{checkDatamodel}}
{}


#' Generate new errorlocation object
#'
#' @param adapt Logical array of same dimensions of dataset, indicating which variable have to be changed per record
#' @param status data.frame 
#' @param call A \code{call} object
#' @param user \code{character} username
#' @param timestamp \code{character} returned by \code{date()}.
#'
#' @keywords internal
newerrorlocation <- function(adapt, status, call=sys.call(-1), method, user=Sys.info()["user"], timestamp=date()){
    structure(
        list(
            adapt  = adapt,
            status = status,
            call = call,
            method = method,
            user = user,
            timestamp = timestamp
        ),
        class=c('errorLocation')
    ) 
}


# empty errorLocation object for a data.frame adapt array
emptyerrorlocation <- function(dat, ...){
    m <- nrow(dat)
    n <- ncol(dat)
    newerrorlocation(
        array(FALSE,dim=c(m,n),dimnames=dimnames(dat)),
        status = emptyStatus(m),
        ...
    )
}


#' Print object of class errorLocation 
#'
#' @param x object of class errorLocation
#' @param ... arguments to be passed to other methods
#' @method print errorLocation
#' @keywords internal
#' @export
print.errorLocation <- function(x,...){
    cat("Object of class 'errorLocation' generated at",x$timestamp,'\n')
    cat("call :", as.character(as.expression(x$call)),'\n')
    cat("method :", x$method,'\n')
    cat("slots:",paste("$",names(x),sep=''),'\n\n')
    
    if ( nrow(x$adapt) <= 10 ){ 
        cat('Values to adapt:\n')
        print(x$adapt)
        cat('\n','Status:\n')
        print(x$status)
    } else {
        cat('Values to adapt:')
        print(x$adapt[1:10,,drop=FALSE])
        cat("...print truncated",'\n\n')
        cat('\n','Status:\n')
        print(x$status[1:10,,drop=FALSE])
        cat("...print truncated",'\n')
    }

}




# Combine two objects
#
# @param x object 1
# @param y object 2
# @param ... extra parameters
# @keywords internal
# @name plus
# @rdname plus
`%+%` <- function(x,y,...){
    if ( is.null(y) ) return(x)
    UseMethod("%+%")
}

# Combine two objects of class errorLocation
# @method `%+%` errorLocation
# @keywords internal
# 
# @rdname plus
`%+%.errorLocation` <- function(x,y,...){
    stopifnot( class(y) == 'errorLocation' )
    call <- unique(c(x$call,y$call))
    if (length(call) == 1) call <- call[[1]]
    newerrorlocation(    
        adapt=x$adapt | y$adapt,
        status=data.frame(
            weight     = x$status$weight + y$status$weight,
            degeneracy = x$status$degeneracy * y$status$degeneracy,
            user       = x$status$user + y$status$user,
            system     = x$status$system + y$status$system,
            elapsed    = x$status$elapsed + y$status$elapsed,
            maxDurationExceeded = x$status$maxDurationExceeded | y$status$maxDurationExceeded,
            memfail    = x$status$memfail | y$status$memfail
            ),
        call = call,
        method=unique(c(x$method, y$method)),
        user = unique(c(x$user,y$user)),
        timestamp=date()
    )

}


# NULL method
# @method `%+%` NULL
# @keywords internal
# @rdname plus
`%+%.NULL` <- function(x,y,...){
    y
}


# create an empty Status data.frame with n rows
emptyStatus <- function(n, weight=0, degeneracy=1, user=0, system=0, elapsed=0, maxDurationExceeded=FALSE,memfail=FALSE){
    data.frame(
        weight=rep(weight,n),
        degeneracy=rep(degeneracy,n),
        user=rep(user,n),
        system=rep(system,n),
        elapsed=rep(elapsed,n),
        maxDurationExceeded=rep(maxDurationExceeded,n),
        memfail=rep(memfail,n)
    )
}







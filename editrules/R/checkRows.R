# Check rows of data.frame against edits.
#
# This is an S3 generic function for checking rows of a \code{data.frame} against
# a number of edit restrictions. The edits can be entered either in \code{character}
# \code{data.frame} or \code{editmatrix} format.
#
# If edits are represented as a \code{character} vector, the entries of \code{E} are parsed
# and evaluated in the environment of \code{dat}
#
# If the edits are  represented in a \code{data.frame}, the \code{data.frame} must have the format
# described in \code{\link{editmatrix}}. The edits are coerced to a character vector, and passed
# to \code{checkRows.character}.
#
# If the edits are represented by an \code{\link{editmatrix}} (representing linear (in)equalities)
# verbose edits are extracted and passed on to \code{checkRows.character}
#
#
# @aliases checkRows.character checkRows.data.frame checkRows.editmatrix
#
# @param E Edits, in \code{character}, \code{data.frame} or \code{\link{editmatrix}} representation
# @param dat The data to check.
# @return a logical vector with \code{length} equal to \code{nrow(dat)}. If a row is violates 
#      no edit restrictions, \code{TRUE} otherwise \code{FALSE}
#
# @seealso violatedEdits
# @example ../examples/checkRows.R
# 
checkRows <- function(E, dat){
    UseMethod("checkRows")
}



# @rdname checkRows
# @method checkRows editmatrix
# 
checkRows.editmatrix <- function(E, dat){
    stopifnot(is.data.frame(dat))
    vars <- getVars(E) %in% names(dat)
    if (!all(vars)){
       stop("Edits contain variable(s):", paste(colnames(E)[!vars], collapse=","), 
            ", that are not available in the data.frame")
    }
    
    check <- !logical(nrow(dat))
    ed <- as.expression(E)
    for ( i in 1:length(ed)){
        check <- check & tryCatch(eval(ed[[i]], envir=dat), error=function(e){
            stop(paste("Edit",ed[[i]],"can not be checked. Evaluation returned",e$message,sep="\n" ))
        })
    }
    return(check)    
} 

# @rdname checkRows
# @method checkRows character
# 
checkRows.character <- function(E, dat){
   ed <- parseEdits(E)
    check <- !logical(nrow(dat))
    for ( i in 1:length(E)){
        check <- check & tryCatch(eval(ed[[i]], envir=dat), error=function(e){
            stop(paste("Edit",ed[[i]],"can not be checked. Evaluation returned",e$message,sep="\n" ))
        })
    }
    return(check)
}

# @rdname checkRows
# @method checkRows data.frame
# 
checkRows.data.frame <- function(E, dat){
    checkRows(as.character(E$edit),dat)
}




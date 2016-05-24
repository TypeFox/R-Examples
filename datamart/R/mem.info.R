#' Information on objects in R environment
#'
#' This function creates a list of objects that are currently in a given R environment (default the global workspace). 
#' Hence it is an extended version of \code{ls}.
#'
#' @param envir            the environment to inspect, default is .GlobalEnv
#' @param sortBy           the result will be decreasingly sorted by this column. Possible values "Type", "Size" (default), "Rows", "Columns"
#'
#' @return data.frame with object information
#' @references 
#' \href{http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session}{Stackoverflow}
#' @export
#' @author Petr Pikal, David Hinds and Dirk Eddelbuettel
mem.info <- function(envir=.GlobalEnv, sortBy="Size") {
    napply <- function(names, fn) sapply(names, function(x) fn(get(x, envir=envir)))
    names <- ls(envir=envir)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x) as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    out <- out[order(out[[sortBy]], decreasing=TRUE), ]
    return(out)
}


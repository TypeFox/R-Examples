
#' get index list from editmatrix
#' 
#' The 'ind' attribute is a named list of named integer vectors. The list names are the 
#' variable names. The vectors in the list index the columns in the editarray associated with the
#' variables. The names of the vectors are the names of the columns of the editarray.
#' 
#' @param E \code{\link{editarray}}
#' @return named list, indexing category levels in the editarray (columns)
#' @keywords internal
getInd <- function(E) attr(E,"ind")


#' get seprator used to seperate variables from levels in editarray
#' @param E \code{\link{editarray}}
#' @return character
#' @keywords internal
getSep <- function(E) attr(E,"sep")

#' Get named logical array from editarray
#' @param E \code{\link{editarray}}
#' @return logical array
#' @keywords internal
getArr <- function(E) unclass(E)[,,drop=FALSE]

#' retrieve level names from editarray
#' @param editarray \code{\link{editarray}}
#' @return character vector
#' @keywords internal
getlevels <- function(E) colnames(E)

#' retrieve edit names from editarray
#' @param E \code{\link{editarray}}
#' @return character vector
#' @keywords internal
getnames <- function(E) rownames(E)

#' Summarize data model of an editarray in a data.frame
#'
#' @param E \code{\link{editarray}}
#' @return \code{data.frame} describing the categorical variables and their levels.
#' @example ../examples/datamodel.R
#' @seealso \code{\link{checkDatamodel}}
#' @export
datamodel <- function(E){
    if (ncol(E) == 0 ) return(data.frame(variable=character(0),value=character(0)))
    st <- stack(getInd(E))
    vals <- do.call(c,lapply(getInd(E),names))
    data.frame(variable=as.character(st[,2]),value=vals,row.names=1:nrow(st))
}


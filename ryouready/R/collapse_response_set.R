# allow rec to be a character or numeric vector that matches 
# the number of variables
# rec     the recoding. May be a string as used in recode in car package
#         or a vector whose length matches the number of variables to recode
# x       the vector after collapsing the variables
#
update_rec <- function(rec, x) 
{
  if (is.null(rec))
    return(rec)
  if (length(rec) > 1) {             # is rec a vector with several entries?
    if (length(x) != length(rec))    # length must match the number of vars
      stop("the length of 'rec' must match the number of variables to recode", call.=FALSE)
    if (is.character(rec))           # add ticks around character vector elements
      rec <- paste0("'", rec, "'")
    rec <- paste0(1:length(x), "=", rec, collapse="; ")   # create recodes as required by car::recode
  }
  rec
}


# collapse variables into one and assign values according to recordings
# x       dataframe containing the variables of the multiple response set
# rec     the recoding. May be a string as used in recode in car package
#         or a vector whose length matches the number of variables to recode
#
collapse_vars <- function(x, rec=NULL) 
{
  x[is.na(x)] <- 0
  if (sum(x) != 1)                  # if 0 or more than 1 entry
    r <- NA
  else {
    r <- sum(x * 1:length(x))       # multiply each column by different number
  }
  if (!is.null(rec)) {
    r <- car::recode(r, update_rec(rec, x))
  } 
  r
}


#' @export
#' @method collapse_responseset data.frame
#' @rdname collapse_responseset
#'  
collapse_responseset.data.frame <- function(x, vars=NULL, rec=NULL, ...) 
{
  if (is.null(vars))
    vars <- 1:ncol(x)
  apply(x[vars], 1, collapse_vars, rec=rec)  
}


#' @export
#' @method collapse_responseset default
#' @rdname collapse_responseset
#' 
collapse_responseset.default <- function(..., rec=NULL)
{
  d <- as.data.frame(list(...))
  collapse_responseset(d, rec=rec)
}


#' Collapse multiple response sets to single variable
#' 
#' This functions allows to collapse several multiple response set varables into
#' one variable. It can be applied either to a \code{dataframe} or 
#' within the \code{transform} function.
#'  
#' @param  x A dataframe.
#' @param ... Several vector of the same length (for default method).
#' @param vars The names or indexes of the dataframe columns that contain the 
#'   multi response set. By default all variables from dataframe are used.
#' @param rec A vector of the same length as the number of variables specifying 
#'   the new values for each column.
#' @return A vector with the with the new values.
#' @rdname collapse_responseset
#' @export
#' @author Mark Heckmann
#' @examples 
#'   
#'  d <- data.frame(t1=c(1,0,NA,0,0),
#'                  t2=c(0,1,0,NA,0),
#'                  t3=c(0,0,1,0,0) )
#'
#'  # collapse all variables of a dataframe
#'  collapse_responseset(d)
#'  
#'  # collapse columns 1 to 3 (which is all in this case as well)
#'  collapse_responseset(d, vars=1:3)
#'  collapse_responseset(d, vars=c("t1", "t2", "t3"))
#'  
#'  # use letters instead fo numbers for recoding
#'  collapse_responseset(d, vars=1:3, rec=letters[1:3])
#'  
#'  # use with several vectors
#'  collapse_responseset(d$t1, d$t2, d$t3)
#'  
#'  # use inside of transform
#'  transform(d, new=collapse_responseset(t1, t2, t3))
#'
#  # use inside of transform with letters as recodes
#'  transform(d, new=collapse_responseset(t1, t2, t3, rec=letters[1:3]))
#'  
collapse_responseset <- function(x, ...) 
{
  UseMethod("collapse_responseset")
}





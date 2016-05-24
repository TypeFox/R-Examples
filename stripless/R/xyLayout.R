#' @title xyLayout
#' 
#' @description Unexported generic and methods used by \code{strucplot} to check and
#' fix or create a \code{xyLayout} list.
#' 
#' @details The various methods provide convenience and flexibility for specifying
#' the \code{xyLayout} list argument that controls the format of strucplot
#' displays. Essentially any sensible way of specifying the xyLayout should
#' work. See the Help page for \code{\link{strucplot}} for details.
#' 
#' @param xylay An appropriate list, matrix, or vector for determining the
#' plot structure. Can also be missing, length 0, etc.
#' 
#' @param n The number of conditioning factors. The integers in the combined list
#' that is created will be a permutation of 1, 2, \ldots ,n .
#' 
#' @param \ldots x and/or y vector for xylay list. The remaining component will be
#' constructed as needed.
#' 
#' @return A list of class "xyLayout" suitable for the xyLayout argument of 
#' \code{strucplot}.
#' 
#' @seealso \code{\link{strucplot}}
#' 
xyLayout <- function(xylay, ...)UseMethod("xyLayout")
#
#'@describeIn xyLayout List method.
xyLayout.list <- function(xylay = list(), n = 
                     stop("Number of conditioning factors is missing"))
  ## xylay is the layout list. 
  ## n is the number of conditioning variables
{
  ### check length, names arguments for agreement, duplicates,  missing x or y
  ### etc. and "fix" if possible
  if(n < 1 || floor(n) != n)
    stop("n must be a nonnegative integer = the number of conditioning factors")
  xy <- c("x","y")
  ix <- unlist(xylay)
  nm <- tolower(names(xylay))
  sq.n <- seq_len(n)
  make_even_layout <- function(n)
  {
    if(n== 1)list(x=1, y = integer(0)) else {
      k <- ceiling(n/2)
      list(x=seq_len(k),y=seq.int(k+1,n)) ## return 
    }
  }
  if(!length(ix)) xylay <- make_even_layout(n)
  else if(!((length(xylay) == 2) && all(xy %in% nm) &&
          !anyDuplicated(ix) && all(sort(ix) == sq.n))) {
    if(length(xylay) > 2)stop("Layout list cannot have length > 2")
    if(anyNA(ix))stop("NA's not permitted")
    else if(!is.numeric(ix))stop("Layout lists must be numeric")
    else if(any(!ix %in% sq.n) ||anyDuplicated(ix) ) 
      stop(sprintf("Layout lists must be non-duplicated integers between 1 and %d",
                   n))
    else {
      if(identical(character(0),nm))names(xylay)<- nm <- xy[seq_along(xylay)]
      else {
        nm[is.na(nm)] <- ""
        if(any(!nm %in% c(xy,""))) stop("Bad names")
        else if(anyDuplicated(nm[nm != ""]))stop("Duplicated names not permitted")
      }
      if(length(xylay) == 2){
        if(any(nm == ""))names(xylay)[nm == ""]<- setdiff(xy,nm)
        hasvals <- sapply(xylay,length)
          if(all(hasvals) && any(!sq.n %in% ix))
            stop("Layout lists do not match number of conditioning factors")
          else if(any(!hasvals))xylay[[which(!hasvals)]]<- setdiff(sq.n, ix)
        } else {
          if(!length(xylay[[1]]))xylay <- make_even_layout(n)
          else {
            xylay[[2]] <- setdiff(sq.n, ix)
            if(nm == "")names(xylay) <- xy
            else names(xylay)[2] <- setdiff(xy,nm)
          }
        }
      }
    }
  structure(xylay[xy], class = c("xyLayout",class(xylay)))
}
  
#' @describeIn xyLayout Matrix method.
xyLayout.matrix <- function(xylay,n)
{
    nc <- ncol(xylay)
    if(!(nc %in% 1:2)) stop("Matrix must have 1 or 2 columns")
    if(!is.numeric(xylay))stop("Matrix must be numeric")
    if(anyNA(xylay))stop("NA's not allowed")
    xy <- c("x","y","",NA)
    nm <- tolower(dimnames(xylay)[[2]])
    if(identical(nm,character(0)) || all(nm %in% xy[3:4]))
      dimnames(xylay) <- list(NULL,xy[seq_len(nc)])
    else if(any(!nm %in% xy)) stop(
      "Column names must be absent,NA,'x', or 'y'")
    else dimnames(xylay)[[2]][!nm %in% xy[1:2]] <-setdiff(xy[1:2], nm) 
    xyLayout(as.list(data.frame(xylay)),n)
}
#' @describeIn xyLayout Data frame method.
xyLayout.data.frame <- function(xylay,n)xyLayout(as.list(xylay),n)

#' @describeIn xyLayout Default methods in which the individual components
#' are given as 1 or 2 vectors of integers.
xyLayout.default <- function(xylay,n) xyLayout(list(x = xylay), n)

    
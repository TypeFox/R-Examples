#' Categorize (semi-)continuous variables
#' 
#' Categorize continuous or semi-continuous variables.  This is a utility
#' function that is useful for writing custom wrapper functions such as
#' \code{\link{simEUSILC}}.
#' 
#' If \code{zeros} is \code{TRUE}, 0 is added to the break points and treated
#' as its own factor level.  Consequently, intervals for negative values are
#' left-closed and right-open, whereas intervals for positive values are
#' left-open and right-closed.
#' 
#' @name getCat
#' @param x a numeric vector to be categorized.
#' @param breaks a numeric vector of two or more break points.
#' @param zeros a logical indicating whether \code{x} is semi-continuous, i.e.,
#' contains a considerable amount of zeros.  See \dQuote{Details} on how this
#' affects the behavior of the function.
#' @param right logical; if \code{zeros} is not \code{TRUE}, this indicates
#' whether the intervals should be closed on the right (and open on the left)
#' or vice versa.
#' @return A \code{\link{factor}} containing the categories.
#' @author Andreas Alfons
#' @export
#' @seealso \code{\link{getBreaks}}, \code{\link{cut}}
#' @keywords manip
#' @examples
#' 
#' data(eusilcS)
#' 
#' ## semi-continuous variable
#' breaks <- getBreaks(eusilcS$netIncome, 
#'     weights=eusilcS$rb050, equidist = FALSE)
#' netIncomeCat <- getCat(eusilcS$netIncome, breaks)
#' summary(netIncomeCat)
#' 
getCat <- function(x, breaks, zeros = TRUE, right = FALSE) {
  # initializations
  if(!is.numeric(x)) stop("'x' must be a numeric vector")
  checkBreaks(breaks)
  if(isTRUE(zeros)) {
    pos <- which(x > 0)
    neg <- which(x < 0)
    # positive values (also works if none exist)
    bpos <- c(0, breaks[breaks > 0])
    if(length(bpos) == 1) lpos <- NULL
    else {
      cpos <- cut(x[pos], bpos)
      lpos <- levels(cpos)
    }
    # negative values (also works if none exist)
    bneg <- c(breaks[breaks < 0], 0)
    if(length(bneg) == 1) lneg <- NULL
    else {
      cneg <- cut(x[neg], bneg, right=FALSE)
      lneg <- levels(cneg)
    }
    # put it all together
    categories <- factor(ifelse(is.na(x), NA, 0), levels=c(lneg, 0, lpos))
    if(length(pos) > 0 && length(bpos) > 1) categories[pos] <- cpos
    if(length(neg) > 0 && length(bneg) > 1) categories[neg] <- cneg
    # return vector
    categories
  } else cut(x, breaks, include.lowest=TRUE, right=right)
}
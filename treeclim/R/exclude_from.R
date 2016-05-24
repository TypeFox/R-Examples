##' Exclude months from analysis
##' 
##' Single months or ranges of months can be excluded from
##' analysis. This is helpful for e.g. excluding winter months without
##' cambial activity.
##' @param month range of numeric month ids
##' @param exclude range or set of months to exclude
##' @return a reduced set of numeric month ids
##' @keywords manip
##' @examples
##' exfr(-5:10, -10:3)
##' @seealso \code{link{.range}}, \code{link{.mean}},
##' \code{link{.sum}}
##' @export
exclude_from <- function(month, exclude = NULL) {
  ## check if exclude is given
  if (is.null(exclude))
    return(month)

  ## check if exclude is supplied as a range
  if (is_continuous(exclude)) {
    exclude <- correct_continuous(exclude)
  } else {
    if (any(exclude == 0)) {
      stop("It is not possible to mix ranges through zero with other specifications.")
    }
  }
  
  ## month must be supplied as a range!
  if (!is_continuous(month)) {
    stop("`month` has to be a continuous range.")
  } else {
    month <- correct_continuous(month)
    ## exclude all excluded months (if available)
    matching <- na.omit(match(exclude, month))
    if (length(matching) > 0) {
      month <- month[-na.omit(matching)]
    } 
  }
  month
}

## (aliased for convenience)

##' @rdname exclude_from
##' @export
exfr <- exclude_from

#' Group a numerical variable into an n-level factor
#'
#' \code{GroupNumeric} converts a numerical variable x
#' into a factor variable with n groups for categorical
#' association analysis.
#'
#' This function uses the classIntervals function from
#' the classInt package to compute the breakpoints
#' that define the groups.  The style parameter is
#' passed to the classIntervals function to specify
#' the grouping method.  Note that some methods may
#' return a different number of groups than that
#' requested via the n parameter.  If groupNames is
#' specified consistently with n, this different
#' number of returned groups will cause an error.
#' The recommended approach in this case is to either
#' change the style parameter or to re-run without
#' groupNames specified.
#'
#' @param x Numeric vector to be grouped.
#' @param n Integer number of groups; if NULL (the
#' default), the number of groups will be inferred
#' from the groupNames parameter.
#' @param groupNames Character vector of names for
#' the levels of the factor variable created; if
#' NULL (the default), the default names from R's
#' cut function will be used.
#' @param orderedFactor Logical, specifying whether
#' the factor returned is ordered or not; default
#' is FALSE.
#' @param style Character string, passed to the
#' classIntervals function from the classInt
#' package as its style parameter (see help
#' file for the classInterval function for
#' details); default is "quantile".
#' @param \dots Optional parameters passed to the
#' classIntervals function from the classInt package.
#' @return Factor variable with n distinct levels,
#' named according to groupNames (if specified).
#' @author Ron Pearson
#' @export
#'
GroupNumeric <- function(x, n = NULL, groupNames = NULL,
                         orderedFactor = FALSE,
                         style = "quantile", ...){
  #
  #  Check consistency of n and groupNames
  #
  nNames <- length(groupNames)
  if (nNames > 0){
    if (is.null(n)){
      n <- nNames
    } else {
      if (n != nNames){
        stop("Length of groupNames is not equal to n")
      }
    }
  } else {
    if (is.null(n)){
      stop("Either n or groupNames must be specified")
    }
  }
  #
  #  Call classIntervals for group cutoffs
  #
  y <- classInt::classIntervals(x, n,
                                style = style, ...)
  #
  #  For some classInt styles, number of groups may differ
  #  from n - if this happens and groupNames is specified,
  #  the call to cut will fail - test for this and give
  #  the user an informative error message
  #
  cutoffs <- y$brks
  nGrps <- length(cutoffs) - 1
  if ((nNames > 0) & (nGrps != n)){
    errMsg <- paste("Internal classInt call with style = ",
                    style, "returned", nGrps,
                    "groups instead of", n,
    "\n To see the resulting groups, rerun without groupNames")
    stop(errMsg)
  }
  #
  #  Call cut with group cutoffs to create factor
  #
  z <- cut(x, breaks = cutoffs, labels = groupNames,
           include.lowest = TRUE,
           ordered_result = orderedFactor)
  return(z)
}

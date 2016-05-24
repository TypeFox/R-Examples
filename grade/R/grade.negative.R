###########################################################################
#
# checks a single number to see if it's negative.  correctans may be null in
# this context, as it has no meaning.
#
###########################################################################
grade.negative <- function(correctans=NULL, studentans, tolerance=0.01,
                           useeval=TRUE, usena=FALSE, useinf=FALSE,
                           quiet=TRUE) {
  stopifnot(grade.allnotnull(studentans, tolerance, useeval, usena, useinf))

  tol.num <- grade.gettolerance(tolerance)
  if (is.null(tol.num)) return(FALSE)

  if (usena) {
    warning("usena is not used in grade.negative, setting to FALSE")
    usena <- FALSE
  }
  #################################################################
  # simply convert studentans to a number, then pass it off to
  # grade.orderedset
  sans <- grade.parse(studentans, useeval, usena, useinf, quiet)
  if (is.null(sans)) return(FALSE)

  if (length(sans) != 1) return(FALSE)

  return( sans < -tol.num )
}

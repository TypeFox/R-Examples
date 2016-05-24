###########################################################################
#
# matches a correct interval to a student interval.  First it converts
# the student
#
###########################################################################
grade.interval <- function(correctans,studentans,tolerance=0.01,
                           useeval=TRUE, usena=FALSE, useinf=FALSE,
                           quiet=TRUE) {
  #
  # No null input allowed
  #
  stopifnot(grade.allnotnull(correctans, studentans, tolerance,
                             useeval, usena, useinf, quiet))
  stopifnot(all(is.logical(useeval), is.logical(usena), is.logical(useinf),
                is.logical(quiet)))

  if (usena) {
    warning("grade.interval does not use usena, setting to FALSE")
    usena <- FALSE
  }

  #################################################################
  # correctans should be a numeric vector of length 2
  correctans <- grade.parse(correctans, useeval, usena, useinf, quiet)
  if (is.null(correctans)) {
    if (quiet == FALSE)
      warning("correctans failed to parse in grade.interval")
    return(FALSE)
  }
  if (length(correctans) != 2) {
    warning("correctans is not a vector of length 2")
    return(FALSE)
  }
  if (correctans[1] > correctans[2]) {
    # I guess this is okay, but we'll whine about it
    if (quiet == FALSE)
      warning("correctans interval in REVERSE order")
  }

  tol.num <- grade.gettolerance(tolerance)
  if (is.null(tol.num)) return(FALSE)

  studentans <- grade.parse(studentans, useeval, usena, useinf, quiet)
  if (is.null(studentans)) return (FALSE)
  #################################################################
  return(isTRUE(all.equal(correctans, studentans, tolerance=tol.num)))
}

###########################################################################
#
# matches an unordered set to a students
#
###########################################################################
grade.set <- function(correctans,studentans,tolerance=0.01,
                      useeval=TRUE, usena=FALSE, useinf=FALSE,
                      quiet=TRUE) {
  #
  # No null input allowed
  #
  stopifnot(grade.allnotnull(correctans, studentans, tolerance,
                             useeval, usena, useinf))
  stopifnot(all(is.logical(useeval), is.logical(usena), is.logical(useinf),
                is.logical(quiet)))

  #################################################################
  correctans <- grade.parse(correctans, useeval, usena, useinf, quiet)
  stopifnot(is.null(correctans) == FALSE)
  
  studentans <- grade.parseset(studentans, useeval, usena, useinf)
  if (is.null(studentans)) return(FALSE)
  
  tol.num <- grade.gettolerance(tolerance)
  if(is.null(tol.num)) return(FALSE)

  # sort and compare the lists
  correctans <- sort(correctans, na.last=FALSE)
  studentans <- sort(studentans, na.last=FALSE)
  return(isTRUE(all.equal(correctans, studentans, tolerance=tol.num)))
}

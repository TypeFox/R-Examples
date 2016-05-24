###########################################################################
#
# checks an ordered set of elements against a correctans.  has an
# extra option, 'converted' which defaults to FALSE, if true, it
# expects input to already have been verified by others
#
###########################################################################
grade.orderedset <- function(correctans, studentans, tolerance=0.01,
                             useeval=TRUE, usena=FALSE, useinf=FALSE,
                             quiet=TRUE) {
  #################################################################
  #
  # NO NULL INPUT
  #
  #################################################################
  stopifnot(grade.allnotnull(correctans, studentans, tolerance,
                             useeval, usena, useinf, quiet))

  stopifnot(all(is.logical(useeval), is.logical(usena), is.logical(useinf),
                is.logical(quiet)))
  #################################################################
  tol.num <- grade.gettolerance(tolerance)
  if(is.null(tol.num)) return(FALSE)

  correctans <- grade.parse(correctans, useeval, usena, useinf, quiet)
  if (is.null(correctans)) return(FALSE)
  
  studentans <- grade.parse(studentans, useeval, usena, useinf, quiet)
  if (is.null(studentans)) return(FALSE)
  #################################################################
  # the input is verified, so all we have left is a simple check
  return(isTRUE(all.equal(studentans, correctans, tolerance=tol.num)))
}

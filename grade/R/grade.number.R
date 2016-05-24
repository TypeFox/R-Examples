###########################################################################
#
# checks a single number for equality
#
###########################################################################
grade.number <- function(correctans, studentans, tolerance=0.01,
                         useeval=TRUE, usena=FALSE, useinf=FALSE, quiet=TRUE) {
  stopifnot(grade.allnotnull(correctans, studentans, tolerance))

  #################################################################
  # convert studentans to a number and make sure correctans
  # is a scalar
  correctans <- grade.parsechunk(correctans, useeval, usena, useinf, quiet)
  stopifnot(is.null(correctans) == FALSE)
  if (length(correctans) != 1)
    stop("Vector of length > 1 passed for correctans in grade.number")

  tol.num <- grade.gettolerance(tolerance)
  if (is.null(tol.num)) return(FALSE)
  
  sans <- grade.parsechunk(studentans, useeval, usena, useinf)
  if (is.null(sans)) return(FALSE)
  
  return(isTRUE(all.equal(correctans, sans, tolerance=tol.num)))
}

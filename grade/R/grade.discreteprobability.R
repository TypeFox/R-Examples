###########################################################################
#
# checks a set to see if it can be a discrete probability.  It has
# an extra option 'ordered', TRUE/FALSE to indicate if order matters
# checkcorrect -- TRUE/FALSE indicates whether or not to check studentans
#                 against correctans
#
###########################################################################
grade.discreteprobability <- function(correctans, studentans, tolerance=.01,
                                      useeval=TRUE, usena=FALSE, useinf=FALSE,
                                      quiet=TRUE, ordered=FALSE,
                                      checkcorrect=TRUE) {
  #################################################################
  # get input, stop on null
  stopifnot(grade.allnotnull(studentans, tolerance, ordered, checkcorrect,
                             useeval, usena, useinf, quiet))
  stopifnot(all(is.logical(useeval), is.logical(usena), is.logical(useinf),
                is.logical(quiet), is.logical(ordered),
                is.logical(checkcorrect)))

  if (usena == TRUE) {
    warning("usena is ignored in grade.discreteprobability, setting FALSE")
    usena <- FALSE
  }
  if (useinf == TRUE) {
    warning("useinf is ignored in grade.discreteprobability, setting FALSE")
    useinf <- FALSE
  }
  #################################################################
  tol.num <- grade.gettolerance(tolerance)
  if (is.null(tol.num)) return(FALSE)
  
  if (checkcorrect == TRUE) {
    correctans <- grade.parse(correctans, useeval, usena, useinf, quiet)
    stopifnot(is.null(correctans) == FALSE)
    stopifnot(isTRUE(all.equal(sum(correctans), 1, tolerance=tol.num)))
    stopifnot(all(correctans > -tol.num))
  }
  
  #################################################################
  sans <- grade.parseset(studentans, useeval, usena, useinf, quiet)
  if (is.null(sans)) return(FALSE)
  
  # 1.  Probability should sum to onen
  if (abs(sum(sans) - 1) > tol.num) return(FALSE)
  # 2.  Elements >= 0
  if (all(sans > -tol.num) == FALSE) return(FALSE)
  
  # if we aren't checking against correctans, we are done
  if (checkcorrect == FALSE) return(TRUE)

  # if order doesn't matter, sort both lists before passing
  if (ordered == FALSE) {
    sans <- sort(sans)
    correctans <- sort(correctans)
  }

  # 3.  Match the correct distribution
  return (isTRUE(all.equal(correctans, sans, tolerance=tol.num)))
}


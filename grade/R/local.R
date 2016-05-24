###########################################################################
#
# simple function, checks to see if all of it's arguments are null or not
# returns true if no arguments, or they are all not null
#
###########################################################################
grade.allnotnull <- function(...) {
  for(x in list(...))
    if (is.null(x)) return(FALSE)
  return(TRUE)
}
###########################################################################
#
# tolerance handling argument, converts tolerance for functions.  Possibly
# from a string.
#
###########################################################################
grade.gettolerance <- function(tolerance) {
  ## this could be updated to use eval etc like everything else
  tol.num <- as.numeric(tolerance)
  if (grade.isscalar(tol.num) == FALSE) return(NULL)
  return(tol.num)
}

###########################################################################
#
# checks to see if the argument is a finite scalar, returns true or false as
# appropriate
#
###########################################################################
grade.isscalar <- function(x, usena=FALSE, useinf=FALSE, quiet=TRUE) {
  # these go first to force it into a scalar
  if (is.vector(x) == FALSE) {
    if (quiet == FALSE)
      warning("Object passed to grade.isscalar is not a vector")
    return(FALSE)
  }
  if (length(x) != 1) {
    if (quiet == FALSE)
      warning("Object passed to grade.isscalar has length != 1")
    return(FALSE)
  }
  ## scalar NA's aren't numeric.
  if (is.logical(x)) {
    if (usena == FALSE) {
      if (quiet == FALSE)
        warning("Object passed to grade.isscalar is logical but usena=FALSE")
      return(FALSE)
    }
    if (is.na(x) == FALSE) {
      if (quiet == FALSE)
        warning("Object passed to grade.isscalar is logical but not na") 
      return(FALSE)
    }
    return(TRUE)
  }
  if (is.numeric(x) == FALSE) {
    if (quiet == FALSE)
      warning("Inappropriate object passed to grade.isscalar")
    return(FALSE)
  }
  if ((usena==FALSE) && is.na(x) == TRUE) {
    if (quiet == FALSE)
        warning("NA passed to grade.isscalar, but usena=FALSE")
    return(FALSE)
  }
  if ((useinf==FALSE) && is.infinite(x) == TRUE) {
    if (quiet == FALSE)
      warning("Inf or -Inf passed to grade.isscalar but useinf=FALSE")
    return(FALSE)
  }
  return(TRUE)
}

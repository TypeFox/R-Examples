###########################################################################
#
# A function to help parse true false elements
#
###########################################################################
grade.parsetruefalse <- function(chunk, useeval=TRUE, quiet=TRUE) {
  if (is.null(chunk)) return(NULL)
  if (is.character(chunk) == FALSE) {
    if (is.vector(chunk) == TRUE && length(chunk) > 1) {
      grade.warn("non vector non character passed to grade.parsetruefalse", quiet)
      return(NULL)
    }
    if (is.logical(chunk) == FALSE) {
      grade.warn("non logical non character passed to grade.parsetruefalse", quiet)
      return(NULL)
    }
    if (is.na(chunk) == TRUE) {
      grade.warn("NA passed to grade.parsetruefalse", quiet)
      return(NULL)
    }
    return(chunk)
  }
  if (grade.safechunk(chunk, quiet) == FALSE) {
    return(NULL)
  }
  ans <- NULL
  if (useeval) {
    ans <- eval(parse(text=toupper(chunk)))
  } else {
    ans <- as.logical(chunk)
  }
  if (is.logical(ans)) {
    if (is.na(ans)) {
      ans <- NULL
      grade.warn("grade.parsetruefalse received text that parsed to NA", quiet)
    }
  } else {
    grade.warn("grade.parsetruefalse received text that parsed to a non-logical", quiet)
    ans <- NULL
  }
  return(ans)
}

###########################################################################
#
# matches a correct interval to a student interval.  First it converts
# the student
#
###########################################################################
grade.truefalse <- function(correctans,studentans,tolerance=0.01,
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
    warning("grade.truefalse does not use usena, setting to FALSE")
    usena <- FALSE
  }
  if (useinf) {
    warning("grade.truefalse does not use useinf, setting to FALSE")
    useinf <- FALSE
  }

  #################################################################
  # correctans should be a numeric vector of length 2
  correctans <- grade.parsetruefalse(correctans, useeval, quiet)
  studentans <- grade.parsetruefalse(studentans, useeval, quiet)
  if (is.null(correctans)) {
    grade.warn("correctans failed to parse in grade.truefalse", quiet)
    return(FALSE)
  }
  if (is.null(studentans)) return (FALSE)
  #################################################################
  return(correctans==studentans)
}

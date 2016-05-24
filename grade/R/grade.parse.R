###########################################################################
#
# grade.warn -- issue a 'warning' of the supplied argument if quiet is false
#
###########################################################################
grade.warn <- function(ans, quiet=TRUE) {
  if (quiet == FALSE)
    warning(ans)
}
###########################################################################
#
# grade.safechunk -- takes a chunk of text and says if it's safe or not
#                    to eval.  The goal is to not eval things that would
#                    be function calls or assignments.  I.e. try to contain
#                    the students answers 
#
###########################################################################
grade.safechunk <- function(chunk, quiet=TRUE) {
  if (is.character(chunk) == FALSE) {
    warning("non character passed to grade.safechunk")
    return(FALSE)
  }
  check.expr <- "[],[()><=]"
  check.res <- regexpr(check.expr, chunk)
  if (check.res[1] != -1) {
    grade.warn("An element passed to grade.parsechunk contained a forbidden character", quiet)
    return(FALSE)
  }
  return(TRUE)
}
###########################################################################
#
# grade.parse -- takes a student's answer string, and returns either
#                NULL -- input not valid
#                a scalar -- input valid, no brackets
#                a vector -- input valid, bracketed by [[(] and [])]
#
###########################################################################
grade.parse <- function(ans, useeval=TRUE, usena=FALSE, useinf=FALSE,
                        quiet=TRUE) {
  ##################################################################
  # check to see if the input doesn't need any processing
  if (is.null(ans)) return(NULL)
  if (is.character(ans) == FALSE) {
    if (is.vector(ans) == FALSE) {
      grade.warn("non vector non character passed to grade.parse", quiet)
      return(NULL)
    }
    for(x in ans) {
      if (grade.isscalar(x, usena, useinf, quiet) == FALSE)
        return(NULL)
    }
    return(ans)
  }
  #################################################################
  # branch to either parsechunk or parseset
  check.expr <- "[][(),]"
  check.res <- regexpr(check.expr, ans)
  if (check.res[1] == -1) # no brackets or commas present
    return(grade.parsechunk(ans, useeval, usena, useinf, quiet))
  return(grade.parseset(ans, useeval, usena, useinf, quiet))
}
###########################################################################
#
# string in, number in it out (or null if it is not valid)
#
###########################################################################
grade.parsechunk <- function(ans, useeval=TRUE, usena=FALSE, useinf=FALSE,
                             quiet=TRUE) {
  #################################################################
  # make sure our input is ready to go
  if (is.null(ans)) return(NULL)

  cur <- 0
  if (is.character(ans)) {
    # refuse to eval any strings with these characters
    # this prevents assignments and function calls
    if (grade.safechunk(ans, quiet) == FALSE) {
      return(NULL)
    }
    #################################################################
    # we are ready to evaluate this string,
    if (useeval == TRUE) {
      try.res <- try(cur <- eval(parse(text=ans)), silent=TRUE)
      if (class(try.res) == "try-error") {
        grade.warn("An element passed to grade.parsechunk failed to eval", quiet)
        return(NULL)
      }
    } else {
      cur <- as.numeric(ans)
    }
  } else {
    cur <- ans
  }
  if (grade.isscalar(cur, usena, useinf, quiet)) return(cur)
  return(NULL)
}
###########################################################################
#
# takes a string in, returns a vector of the set it represents if that makes
# sense, NULL otherwise
#
###########################################################################
grade.parseset <- function(ans, useeval=TRUE, usena=FALSE, useinf=FALSE,
                           quiet=TRUE)  {
  stopifnot(grade.allnotnull(useeval, usena, useinf, quiet))
  stopifnot(all(is.logical(useeval), is.logical(usena), is.logical(useinf),
                is.logical(quiet)))
  # bail right away if input passed input is already done
  if (is.null(ans)) return(NULL)
  
  #################################################################
  if (is.character(ans)) {
    #################################################################
    # get rid of leading and trailing whitespace
    ans <- sub("^\\s*", "", ans)
    ans <- sub("\\s*$", "", ans)

    # make sure string starts/ends with '[' or '(' and ']' or '('
    # and none of those in the middle
    brack.expr <- "^[[(][^][()]*[])]$"
    brack.res <- regexpr(brack.expr, ans)
    if (brack.res[1] == -1) {
      return(NULL);
    }

    #################################################################
    # iterate through the string, converting each piece to numeric as we go
    # return NULL if we run into problems in the middle, otherwise
    # return a numeric vector of what we found
    ret.interval <- c()
    comma.expr <- ","
    subans <- substr(ans, 2, nchar(ans)-1) # cut out the brackets
    repeat {
      comma.res <- regexpr(comma.expr, subans)
      chunk <- substr(subans, 1, ifelse(comma.res[1] == -1, nchar(subans),
                                        comma.res[1]-1))
      cur <- grade.parsechunk(chunk, useeval, usena, useinf, quiet)
      if (is.null(cur)) return(NULL)
      ret.interval <- c(ret.interval, cur)
      if (comma.res[1] == -1) {
        break; # hit the end
      } else {
        subans <- substr(subans, comma.res[1]+1, nchar(subans)) # next piece
      }
    }
    return(ret.interval)
  }
  #################################################################
  # Getting here means that it was not a character passed, so we just need to
  # check to see if the elements are okay or not
  if (is.vector(ans) == FALSE) {
    grade.warn("non vector non character passed to grade.parseset", quiet)
    return(NULL)
  }
  for (x in ans) {
    if (grade.isscalar(x, usena, useinf, quiet) == FALSE)
      return(NULL)
  }
  return(ans)
}

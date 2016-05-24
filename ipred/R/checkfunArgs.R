# $Id: checkfunArgs.R,v 1.1 2003/02/17 09:49:31 hothorn Exp $

checkfunArgs <- function(fun, type=c("model", "predict")) {

  # check for appropriate arguments of user-supplied function "fun" 
  # this will not work for generics in R < 1.7.0 and therefore not used by
  # now

  type <- match.arg(type)

  if (!is.function(fun)) {
    warning("fun is not a function")
    return(FALSE)
  }

  funargs <- formals(fun)

  switch(type, "model"={
    if (!all(names(funargs)[1:2] %in% c("formula", "data"))) {
      warning("fun is not a function with at least 'formula' and 'data' arguments")
      return(FALSE)
    } else {
      return(TRUE)
    }
  }, "predict"={
    if (length(funargs) < 2) {
      warnings("fun is not a function with at least 'object' and 'newdata' arguments")
      return(FALSE)
    } else {
      return(TRUE)
    }
  })
}

# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.




if (FALSE)
family.vglm <- function(object, ...) 
    object$vfamily


if (FALSE)
print.vfamily <- function(x, ...) {
  f <- x$vfamily
  if (is.null(f))
    stop("not a VGAM family function")

  nn <- x$blurb
  if (is.null(nn))
    invisible(return(x))

  cat("Family: ", f[1], "\n") 
  if (length(f)>1)
    cat("Classes:", paste(f, collapse=", "), "\n")
  cat("\n")

  for (ii in 1:length(nn))
    cat(nn[ii])
  cat("\n")
  invisible(return(x))
}




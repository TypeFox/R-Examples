########################################################################
#
#   Extensions to mgcv package allowing a projection Slepian basis
#   to be used. 
#
########################################################################

########################################################################
#
#   summary.gam
#
#   When computing a summary(object) for a gam object, if the sp basis
#   is used we require the frequentist version of the confidence interval
#   calculation. To ensure this happens, this function overloads the
#   summary.gam() from mgcv, checks the gam object to see if it contains
#   a sp.smooth-class smooth object, and then calls mgcv::summary.gam
#   with the right parameters.
#
#   If there are no sp.smooth-class objects, everything is passed 
#   through unchanged to mgcv::summary.gam.
#
#   mgcv::summary.gam is called as ...
#     summary.gam <- function (object, dispersion = NULL, freq = FALSE, p.type=0, ...) {}
#
########################################################################
summary.gam <- function(object, ...) {
  sp <- unlist(lapply(X = object[["smooth"]], FUN = function(X) {
    test <- "sp.smooth" %in% attr(X, "class")
    test
  }))

  if(length(which(sp==TRUE)) > 0) { # at least one 'sp.smooth' object
    summ.obj <- mgcv::summary.gam(object, freq=TRUE, p.type=5, ...)
  } else { # any other case
    summ.obj <- mgcv::summary.gam(object, ...)
  }
  summ.obj
}


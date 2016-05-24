# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.






nobs.vlm <- function(object, type = c("lm", "vlm"), ...) {


  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("lm", "vlm"))[1]


  if (type == "lm") {
    object@misc$n
  } else {
    object@misc$nrow.X.vlm
  }
}



if (!isGeneric("nobs"))
  setGeneric("nobs", function(object, ...)
             standardGeneric("nobs"),
             package = "VGAM")


setMethod("nobs", "vlm",
         function(object, ...)
         nobs.vlm(object, ...))











nvar.vlm <- function(object, type = c("vlm", "lm"), ...) {

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("vlm", "lm"))[1]


  if (type == "lm") {
    object@misc$p
  } else {
    object@misc$ncol.X.vlm
  }
}



nvar.vgam <- function(object, type = c("vgam", "zz"), ...) {

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("vgam", "zz"))[1]

  stop("function nvar.vgam() has not been written yet")

  if (type == "vgam") {
    object@misc$p
  } else {
    object@misc$ncol.X.vlm
  }
}


nvar.rrvglm <- function(object, type = c("rrvglm", "zz"), ...) {

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("rrvglm", "zz"))[1]

  stop("function nvar.rrvglm() has not been written yet")

  if (type == "vgam") {
    object@misc$p
  } else {
    object@misc$ncol.X.vlm
  }
}



nvar.qrrvglm <- function(object, type = c("qrrvglm", "zz"), ...) {

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("qrrvglm", "zz"))[1]

  stop("function nvar.qrrvglm() has not been written yet")

  if (type == "qrrvglm") {
    object@misc$p
  } else {
    object@misc$ncol.X.vlm
  }
}



nvar.rrvgam <- function(object, type = c("cao", "zz"), ...) {

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("rrvglm", "zz"))[1]

  stop("function nvar.rrvgam() has not been written yet")

  if (type == "cao") {
    object@misc$p
  } else {
    object@misc$ncol.X.vlm
  }
}



nvar.rcim <- function(object, type = c("rcim", "zz"), ...) {

  if (mode(type) != "character" && mode(type) != "name")
    type <- as.character(substitute(type))
  type <- match.arg(type,
                    c("rcim", "zz"))[1]

  stop("function nvar.rcim() has not been written yet")

  if (type == "rcim") {
    object@misc$p
  } else {
    object@misc$ncol.X.vlm
  }
}





if (!isGeneric("nvar"))
  setGeneric("nvar", function(object, ...)
             standardGeneric("nvar"),
             package = "VGAM")


setMethod("nvar", "vlm",
         function(object, ...)
         nvar.vlm(object, ...))



setMethod("nvar", "vgam",
         function(object, ...)
         nvar.vgam(object, ...))


setMethod("nvar", "rrvglm",
         function(object, ...)
         nvar.rrvglm(object, ...))



setMethod("nvar", "qrrvglm",
         function(object, ...)
         nvar.qrrvglm(object, ...))



setMethod("nvar", "rrvgam",
         function(object, ...)
         nvar.rrvgam(object, ...))



setMethod("nvar", "rcim",
         function(object, ...)
         nvar.rcim(object, ...))






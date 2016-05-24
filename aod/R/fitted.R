if(!isGeneric("fitted"))
  setGeneric(name = "fitted", def = function(object, ...) standardGeneric("fitted"))

## fitted values for models of class glimML (functions betabin and betapois)
setMethod(f = "fitted", signature = "glimML", definition = function(object, ...) {
  mf <- object@CALL
  mb <- match(c("formula", "data", "na.action"), names(mf), 0)
  mfb <- mf[c(1, mb)]
  mfb$drop.unused.levels <- TRUE
  mfb[[1]] <- as.name("model.frame")
  names(mfb)[2] <- "formula"
  mfb <- eval(mfb, parent.frame())
  mt <- attr(mfb, "terms")
  Y <- model.response(mfb, "numeric")
  X <- if(!is.empty.model(mt)) model.matrix(mt, mfb, contrasts) else matrix(, NROW(Y), 0)
  offset <- model.offset(mfb)
  b <- coef(object)
  eta <- as.vector(X %*% b)
  eta <- if(is.null(offset)) eta else eta + offset
  invlink(eta, type = object@link)
  })
  
## fitted values for models of class glimQL (functions quasibin and quasipois)
setMethod(f = "fitted", signature(object = "glimQL"), definition = function(object, ...) fitted(object@fm))

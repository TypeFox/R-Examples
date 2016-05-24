if(!isGeneric("residuals"))
  setGeneric(name = "residuals", def = function(object, ...) standardGeneric("residuals"))

## Residuals for objects of formal class "glimML" (functions betabin and negbin)
setMethod(f = "residuals", signature = "glimML", definition = function(object, type = c("pearson", "response"), ...){
  type <- match.arg(type)
  method <- object@method
  mf <- object@CALL
# get phi.fit
  param <- object@param
  nb.b <- length(coef(object))
  nb.phi <- length(param) - nb.b
  mr <- match(c("random", "data", "na.action"), names(mf), 0)
  mr <- mf[c(1, mr)]
  mr$drop.unused.levels <- TRUE
  mr[[1]] <- as.name("model.frame")
  names(mr)[2] <- "formula"
  mr <- eval(mr, parent.frame())
  modmatrix.phi <- model.matrix(object = object@random, data = mr)
  phi.fit <- modmatrix.phi %*% param[(nb.b + 1):(nb.b + nb.phi)]
# get observed and fitted response
  mb <- match(c("formula", "data", "na.action"), names(mf), 0)
  mfb <- mf[c(1, mb)]
  mfb$drop.unused.levels <- TRUE
  mfb[[1]] <- as.name("model.frame")
  names(mfb)[2] <- "formula"
  mfb <- eval(mfb, parent.frame())
  Y <- model.response(mfb, "numeric")
  y <- if(NCOL(Y) == 2) Y[,1] else Y
  y.fit <- fitted(object)
# compute residuals
  if(object@method == "BB"){
    n <- rowSums(Y)
    res <- switch(type,
                  pearson = ifelse(n == 0, 0, (y - n * y.fit) / sqrt(n * y.fit * (1 - y.fit) * (1 + (n - 1) * phi.fit))),
                  response = ifelse(n == 0, 0, y/n - y.fit))
    }
  if(object@method == "NB"){
    nam <- names(y)
    v <- y.fit + phi.fit * (y.fit^2)
    r <- as.vector((y - y.fit) / sqrt(v))
    names(r) <- nam
    res <- switch(type,
                  pearson = r,
                  response = y - y.fit)
    }
  if(!is.null(object@na.action)) 
    res <- naresid(object@na.action, res)
  res
  })

## residuals for objects of formal class "glimQL" (functions quasibin and quasipois)
setMethod(f = "residuals", signature = "glimQL", definition = function(object, type = c("pearson", "response"), ...){
  type <- match.arg(type)
  residuals(object@fm, type = type, ...)
  })

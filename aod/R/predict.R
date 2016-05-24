if(!isGeneric("predict"))
  setGeneric(name = "predict", def = function(object, ...) standardGeneric("predict"))

## predicted values for models of class glimML (functions betabin and betapois)
setMethod(f = "predict", signature = "glimML",
          definition = function(object, newdata = NULL, type = c("response", "link"), se.fit = FALSE, ...){
  type <- match.arg(type)
  mf <- object@CALL
  b <- coef(object)
  f <- object@formula[-2]
  data <- object@data
  offset <- NULL
  if(is.null(newdata)){
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
    }
  else{
    mfb <- model.frame(f, newdata)
    offset <- model.offset(mfb)
    X <- model.matrix(object = f, data = newdata)
    }
  eta <- as.vector(X %*% b)
  eta <- if(is.null(offset)) eta else eta + offset
  varparam <- object@varparam
  varb <- as.matrix(varparam[seq(length(b)), seq(length(b))])
  vareta <- X %*% varb %*% t(X)
  if(type == "response"){
    p <- invlink(eta, type = object@link)
    J <- switch(object@link,
               logit = diag(p * (1 - p), nrow = length(p)),
               cloglog = diag(-(1 - p) * log(1 - p), nrow = length(p)),
               log = diag(p, nrow = length(p)))
    varp <- J %*% vareta %*% J
    se.p <- sqrt(diag(varp))
    }
  se.eta <- sqrt(diag(vareta))
  if(!se.fit)
    res <- switch(type, response = p, link = eta)
  else
    res <- switch(type, response = list(fit = p, se.fit = se.p), link = list(fit = eta, se.fit = se.eta))
  res
  })

## predicted values for models of class glimQL (functions quasibin and quasipois)
setMethod(f = "predict", signature(object = "glimQL"),
          function(object, newdata = NULL, type = c("response", "link"), se.fit = FALSE, ...) 
            predict(object@fm, newdata = newdata, type = type, se.fit = se.fit, ...))

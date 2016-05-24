setClass("FLXMRmgcv",
         representation(G = "list",
                        control = "list"),
         contains="FLXMRglm")

FLXMRmgcv <- function(formula = .~., family = c("gaussian", "binomial", "poisson"),
                      offset = NULL, control = NULL, optimizer = c("outer", "newton"),
                      in.out = NULL, eps = .Machine$double.eps, ...)
{
  if (is.null(control)) control <- mgcv::gam.control()
  family <- match.arg(family)

  am <- if (family == "gaussian" && get(family)()$link == "identity") TRUE else FALSE
  z <- new("FLXMRmgcv", FLXMRglm(formula = formula, family = family, offset = offset), 
           name=paste("FLXMRmgcv", family, sep=":"), control = control)

  scale <- if (family %in% c("binomial", "poisson")) 1 else -1

  gam_fit <- function(G, w) {
    G$family <- get(family)()
    G$am <- am
    G$w <- w
    G$conv.tol <- control$mgcv.tol
    G$max.half <- control$mgcv.half
    zero_weights <- any(w < eps)
    if (zero_weights) {
      ok <- w >= eps
      w <- w[ok]
      G$X <- G$X[ok,,drop=FALSE]
      if (is.matrix(G$y)) G$y <- G$y[ok,,drop=FALSE] else G$y <- G$y[ok]
      G$mf <- G$mf[ok,,drop=FALSE]
      G$w <- G$w[ok]
      G$offset <- G$offset[ok]

      if (G$n.paraPen > 0) {
        OMIT <- which(colSums(abs(G$X)) == 0)
        if (length(OMIT) > 0) {
          Ncol <- ncol(G$X)
          Assign <- unique(G$assign[OMIT])
          G$assign <- G$assign[-OMIT]
          G$nsdf <- G$nsdf - length(OMIT)
          G$X <- G$X[,-OMIT,drop=FALSE]
          G$mf$Grouping <- G$mf$Grouping[,-which(colSums(abs(G$mf$Grouping))==0),drop=FALSE]
          if (length(G$off) > 1) G$off[2] <- G$off[2] - length(OMIT)
          for (i in seq_along(G$smooth)) {
            G$smooth[[i]]$first.para <- G$smooth[[i]]$first.para - length(OMIT)
            G$smooth[[i]]$last.para <- G$smooth[[i]]$last.para - length(OMIT)
          }
          G$S[[1]] <- G$S[[1]][-c(OMIT-sum(G$assign != Assign)),
                               -c(OMIT-sum(G$assign != Assign))]
        }
      }
    }
    z <- mgcv::gam(G = G, method = "ML", optimizer = optimizer, control = control, scale = scale,
                   in.out = in.out, ...)
    if (zero_weights) {
      residuals <- z$residuals
      z$residuals <- rep(0, length(ok))
      z$residuals[ok] <- residuals
      if (G$n.paraPen > 0 && length(OMIT) > 0) {
        coefficients <- z$coefficients
        z$coefficients <- rep(0, Ncol)
        z$coefficients[-OMIT] <- coefficients
      }
    }
    z
  }

  if (family=="gaussian"){
    z@fit <- function(x, y, w, G){
      gam.fit <- gam_fit(G, w)
      with(list(coef = gam.fit$coefficients, df = sum(gam.fit$edf)+1,
                sigma = sqrt(sum(w * gam.fit$residuals^2 /
                  mean(w))/ (nrow(x)-sum(gam.fit$edf)))),
           eval(z@defineComponent))
    }
  }
  else if(family %in% c("binomial", "poisson")){
    z@fit <- function(x, y, w, G){
      gam.fit <- gam_fit(G, w)
      with(list(coef = gam.fit$coefficients, df = sum(gam.fit$edf)),
                eval(z@defineComponent))
    }
  }
  else stop(paste("Unknown family", family))
  z
}

setMethod("FLXmstep", signature(model = "FLXMRmgcv"), function(model, weights, ...)
{
  apply(weights, 2, function(w) model@fit(model@x, model@y, w, model@G))
})

setMethod("FLXgetModelmatrix", signature(model="FLXMRmgcv"), function(model, data, formula, lhs=TRUE,
                                             paraPen = list(), ...)
{
  formula <- RemoveGrouping(formula)
  
  if (length(grep("\\|", deparse(model@formula)))) stop("no grouping variable allowed in the model")
  if(is.null(model@formula))
    model@formula <- formula

  model@fullformula <- update(terms(formula, data=data), model@formula)
  gp <- mgcv::interpret.gam(model@fullformula)
  if (lhs) {
    model@terms <- terms(gp$fake.formula, data = data)   
    mf <- model.frame(model@terms, data=data, na.action = NULL, drop.unused.levels = TRUE)
    response <- as.matrix(model.response(mf, "numeric"))
    model@y <- model@preproc.y(response)
  }
  else {
    model@terms <- terms(gp$fake.formula, data = data)   
    mf <- model.frame(delete.response(model@terms), data=data, na.action = NULL, drop.unused.levels = TRUE)
  }
  model@G <- mgcv::gam(model@fullformula, data = data, fit = FALSE)
  model@x <- model@G$X
  model@contrasts <- attr(model@x, "contrasts")
  model@x <- model@preproc.x(model@x)
  model@xlevels <- .getXlevels(delete.response(model@terms), mf)
  model
})

setMethod("predict", signature(object="FLXMRmgcv"), function(object, newdata, components, ...)
{
  predict_gam <- function (object, newdata, ...) {
    nn <- names(newdata)
    mn <- colnames(object$model)
    for (i in 1:length(newdata)) if (nn[i] %in% mn && is.factor(object$model[, nn[i]])) {
      newdata[[i]] <- factor(newdata[[i]], levels = levels(object$model[, nn[i]]))
    }
    if (length(newdata) == 1) 
      newdata[[2]] <- newdata[[1]]
    n.smooth <- length(object$smooth)
    Terms <- delete.response(object$pterms)
    X <- matrix(0, nrow(newdata), length(object$coefficients))
    Xoff <- matrix(0, nrow(newdata), n.smooth)
    mf <- model.frame(Terms, newdata, xlev = object$xlevels)
    if (!is.null(cl <- attr(object$pterms, "dataClasses"))) 
      .checkMFClasses(cl, mf)
    Xp <- model.matrix(Terms, mf, contrasts = object$contrasts)
    if (object$nsdf) 
      X[, 1:object$nsdf] <- Xp
    if (n.smooth) 
      for (k in 1:n.smooth) {
        Xfrag <- mgcv::PredictMat(object$smooth[[k]], newdata)
        X[, object$smooth[[k]]$first.para:object$smooth[[k]]$last.para] <- Xfrag
        Xfrag.off <- attr(Xfrag, "offset")
        if (!is.null(Xfrag.off)) {
          Xoff[, k] <- Xfrag.off
        }
      }
    X
  }
  object@G$model <- object@G$mf
  z <- list()
  for(k in seq_along(components)) {
    object@G$coefficients <- components[[k]]@parameters$coef
    X <- predict_gam(object@G, newdata)
    z[[k]] <- components[[k]]@predict(X, ...)
  }
  z
})
                           

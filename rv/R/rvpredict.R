
rvpredict <- function (object, ...) {
  UseMethod("rvpredict")
}

.X.and.offset <- function (object, newdata) {
  offset <- object$offset
  tt <- terms(object)
  if (missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
    offset <- object$offset
  } else {
    Terms <- delete.response(tt)
    m <- rvmodel.frame(Terms, newdata, na.action = na.fail, xlev = object$xlevels)
    ##
    if (!is.null(cl <- attr(Terms, "dataClasses"))) {
      stats::.checkMFClasses(cl, m)
    }
    X <- rvmodel.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- rep(0, nrow(X))
    if (!is.null(off.num <- attr(tt, "offset"))) {
      for (i in off.num) {
        offset <- offset + eval(attr(tt, "variables")[[i + 1]], newdata)
      }
      if (!is.null(object$call$offset)) {
        offset <- offset + eval(object$call$offset, newdata)
      }
    }
  }
  list(X=X, offset=offset)
}

rvpredict.lm <- function (object, newdata, ...) {
  L <- .X.and.offset(object, newdata)
  X <- L$X
  ## end code from predict.lm
  post <- posterior(object)
  X.names <- colnames(X)
  all.coeff <- post$beta
  if (! all(ok <- (X.names %in% names(all.coeff)))) {
    stop("Some names are missing: ", paste(X.names[! ok], collapse=", "))
  }
  beta <- all.coeff[X.names]
  mu.pred <- (X %**% beta) ## rv-compatible multiplication
  mu.pred <- (mu.pred + L$offset)
  sigma <- post$sigma
  y.pred <- rvnorm(mean=mu.pred, sd=sigma)
  return(y.pred)
}

rvmodel.frame <- function (formula, data=NULL, ...) {
  if (is.null(data)) {
    return(model.frame(formula, data=data, ...))
  }
  L <- as.list(data)
  if (! any(sapply(L, is.rv))) {
    m <- match.call()
    m[[1L]] <- as.name("model.frame")
    return(eval(m, envir=parent.frame()))
  }
  terms <- NULL
  .mf <- function (data.matrix, formula, ...) {
    df <- as.data.frame(data.matrix)
    mf <- model.frame(formula, data=df, ...)
    terms <<- attr(mf, "terms")
    return(mf)
  }
  W <- lapply(L, as.rv)
  W <- do.call(cbind.rv, W)
  colnames(W) <- names(L)
  x <- simapply(W, .mf, formula=formula, ...)
  attr(x, "terms") <- terms
  return(x)
}

rvmodel.matrix <- function(object, data, ...) {
  L <- as.list(data)
  if (! any(sapply(L, is.rv))) {
    m <- match.call()
    m[[1L]] <- as.name("model.matrix")
    return(eval(m, envir=parent.frame()))
  }
  terms <- attr(data, "terms")
  mm.names <- NULL
  .mm <- function (data.matrix, object, ...) {
    df <- structure(as.data.frame(data.matrix), terms=terms)
    mm <- model.matrix(object, data=df, ...)
    mm.names <<- colnames(mm)
    return(mm)
  }
  W <- lapply(L, as.rv)
  W <- do.call(cbind.rv, W)
  colnames(W) <- names(L)
  x <- simapply(W, .mm, object=object, ...)
  colnames(x) <- mm.names
  return(x)
}

##

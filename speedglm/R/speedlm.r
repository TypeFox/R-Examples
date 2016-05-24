speedlm.fit <- function(y, X, intercept = FALSE, offset = NULL, row.chunk = NULL, 
                        sparselim = 0.9, camp = 0.01, eigendec = TRUE, 
                        tol.solve = .Machine$double.eps, sparse = NULL, tol.values = 1e-07, 
                        tol.vectors = 1e-07, method=c('eigen','Cholesky','qr'), ...) 
{
  method <- match.arg(method) 
  nvar <- ncol(X)
  nobs <- nrow(X)
  if (is.null(offset)) 
  offset <- rep(0, nobs)
  colnam <- colnames(X)
  if (is.null(sparse)) sparse <- is.sparse(X = X, sparselim, camp)
  if (sparse) X <- as(X, "dgCMatrix")
  A <- if (sparse | is.null(row.chunk)) crossprod(X) else cp(X, , row.chunk)
  y <- y - offset
  Xy <- if (sparse) crossprod(X, y) else t(crossprod(y, X))
  X1X <- colSums(X)
  names(X1X) <- colnam
  yy <- crossprod(y)
  
  switch(method,
         eigen={
           ris <- if (eigendec) control(A, , tol.values, tol.vectors, , method) else 
                                list(XTX = A, rank = nvar, pivot = 1:nvar)
           ris$XTX <- if (sparse) as(ris$XTX, "dgCMatrix") else 
                                  as(ris$XTX, "matrix")	
           ok <- ris$pivot[1:ris$rank]
           coef<-as(solve(ris$XTX, Xy[ok]), "numeric")
           coefficients <- rep(NA, nvar)
           coefficients[ok] <- coef				
           RSS <- yy - 2 * crossprod(coef, Xy[ok]) + t(coef) %*% ris$XTX %*% coef
         },
         Cholesky={
           ris <- if (eigendec) 
             control(A, , tol.values, tol.vectors, , method)
           else list(XTX = A, rank = nvar, pivot = 1:nvar)
           ris$XTX <- if (sparse) 
             as(ris$XTX, "dgCMatrix")
           else as(ris$XTX, "matrix")	
           ok <- ris$pivot[1:ris$rank]
           coef<-as(solve(ris$XTX, Xy[ok]), "numeric")
           coefficients <- rep(NA, nvar)
           coefficients[ok] <- coef				
           RSS <- yy - 2 * crossprod(coef, Xy[ok]) + t(coef) %*% ris$XTX %*% coef
         },
         qr={
           if(class(A)=='dsCMatrix') {
             A<-as(A,'matrix')
             Xy<-as(Xy,'matrix')
           }
           C_Cdqrls<-getNativeSymbolInfo('Cdqrls', PACKAGE=getLoadedDLLs()$stats)
           ris<-c(list(XTX=A), .Call(C_Cdqrls,A, Xy, tol.values, FALSE))
           coefficients<-ris$coefficients
           coef<-coefficients[ris$pivot[1:ris$rank]]
           ord <- ris$pivot
           RSS <- yy - 2 * crossprod(coefficients, Xy[ris$pivot]) + t(coefficients[ord]) %*% ris$XTX %*% coefficients[ord]
           ok<-ris$pivot[1:ris$rank]
           if (ris$rank < nvar) 
             coefficients[(ris$rank + 1L):nvar] <- NA
           coefficients<-coefficients[ord]						 	
         },
         stop('speedlm.fit: Unknown method value')
  )
  names(coefficients) <- colnames(X)
  dfr <- nrow(X) - ris$rank
  rval <- list(coefficients = coefficients, coef = coef, df.residual = dfr, 
               XTX = as(ris$XTX, "matrix"), Xy = Xy, nobs = nrow(X), 
               nvar = nvar, ok = ok, A = as(A, "matrix"), RSS = as.numeric(RSS), 
               rank = ris$rank, pivot = ris$pivot, sparse = sparse, 
               yy = yy, X1X = X1X, intercept = intercept,method=method)
  class(rval) <- "speedlm"
  rval
}

speedlm.wfit <- function(y, X, w, intercept = FALSE, offset = NULL, row.chunk = NULL, 
                         sparselim = 0.9, camp = 0.01, eigendec = TRUE, 
                         tol.solve = .Machine$double.eps, sparse = NULL, tol.values = 1e-07, 
                         tol.vectors = 1e-07, method=c('eigen','Cholesky','qr'), ...) 
{
  nvar <- ncol(X)
  nobs <- nrow(X)
  if (is.null(offset)) 
    offset <- rep(0, length(y))
  if (any(is.na(w))) 
    stop("some weights are NA")
  if (any(w < 0)) 
    stop("weights must not be negative")
  colnam <- colnames(X)
  if (is.null(sparse)) 
    sparse <- is.sparse(X = X, sparselim, camp)
  if (sparse) 
    X <- as(X, "dgCMatrix")
  pw <- prod(w[w != 0])
  sqw <- sqrt(w)
  sqwX <- sqw * X
  SW <- sum(w)
  yy <- crossprod(sqw * y)
  X1X <- colSums(sqwX)
  names(X1X) <- colnam
  XW1 <- crossprod(w, X)
  A <- if (sparse) 
    crossprod(sqwX)
  else cp(sqwX, , row.chunk)
  y <- y - offset
  Xy <- if (sparse) 
    crossprod(X, (w * y))
  else t(crossprod((w * y), X))
  
  switch(method,
         eigen={
           ris <- if (eigendec)
             control(A, , tol.values, tol.vectors, , method)
           else 
             list(XTX = A, rank = nvar, pivot = 1:nvar)
           ris$XTX <- if (sparse)
             as(ris$XTX, "dgCMatrix")
           else as(ris$XTX, "matrix")
           ok <- ris$pivot[1:ris$rank]
           coef <- as(solve(ris$XTX, Xy[ok]), "numeric")
           coefficients <- rep(NA, nvar)
           coefficients[ok] <- coef
           RSS <- yy - 2 * crossprod(coef, Xy[ok]) + t(coef) %*% ris$XTX %*% coef
         },
         Cholesky={
           ris <- if (eigendec)
             control(A, , tol.values, tol.vectors, , method)
           else 
             list(XTX = A, rank = nvar, pivot = 1:nvar)
           ris$XTX <- if (sparse)
             as(ris$XTX, "dgCMatrix")
           else as(ris$XTX, "matrix")
           ok <- ris$pivot[1:ris$rank]
           coef <- as(solve(ris$XTX, Xy[ok]), "numeric")
           coefficients <- rep(NA, nvar)
           coefficients[ok] <- coef
           RSS <- yy - 2 * crossprod(coef, Xy[ok]) + t(coef) %*% ris$XTX %*% coef
         },
         qr={
           if(class(A)=='dsCMatrix') {
             A<-as(A,'matrix')
             Xy<-as(Xy,'matrix')
           }
           C_Cdqrls<-getNativeSymbolInfo('Cdqrls', PACKAGE=getLoadedDLLs()$stats)
           ris<-c(list(XTX=A), .Call(C_Cdqrls,A, Xy, tol.values, FALSE))
           coefficients<-ris$coefficients
           coef<-coefficients[ris$pivot[1:ris$rank]]
           ord<-order(ris$pivot)
           RSS <- yy - 2 * crossprod(coefficients, Xy[ris$pivot]) + t(coefficients[ord]) %*% ris$XTX %*% coefficients[ord]
           ok<-ris$pivot[1:ris$rank]
           if (ris$rank < nvar) 
             coefficients[(ris$rank + 1L):nvar] <- NA
           coefficients<-coefficients[ord]						 	
         },
         stop('speedlm.fit: Unknown method value')
  )
  
  names(coefficients) <- colnames(X)
  zero.w <- sum(w == 0)
  dfr <- nrow(X) - ris$rank - zero.w
  rval <- list(coefficients = coefficients, coef = coef, weights = w, 
               df.residual = dfr, XTX = as(ris$XTX, "matrix"), Xy = Xy, 
               nobs = nrow(X), nvar = nvar, ok = ok, A = as(A, "matrix"), 
               RSS = as.numeric(RSS), rank = ris$rank, pivot = ris$pivot, 
               sparse = sparse, yy = yy, X1X = X1X, SW = SW, XW1 = XW1, 
               zero.w = zero.w, pw = pw, intercept = intercept,method=method)
  class(rval) <- "speedlm"
  rval
}	


speedlm <- function(formula, data, weights = NULL, offset = NULL, sparse = NULL, 
                  set.default = list(), method = c('eigen','Cholesky','qr'), model = FALSE,  
                  y = FALSE, fitted = FALSE, subset = NULL,...) 
{
  target <- y
  call <- match.call()
  M <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(M), 0L)
  M <- M[c(1L, m)]
  M$drop.unused.levels <- TRUE
  M[[1L]] <- quote(stats::model.frame)
  M <- eval(M, parent.frame())
  y <- M[[1]]
  tf <- terms(formula, data=data)
  X <- model.matrix(tf, M)
  offset <- model.offset(M)
  if (is.null(offset)) 
    offset <- rep(0, length(y))
  set <- list(sparselim = 0.9, camp = 0.01, eigendec = TRUE, 
              row.chunk = NULL, tol.solve = .Machine$double.eps, tol.values = 1e-07, 
              tol.vectors = 1e-07, method = method)
  nmsC <- names(set)
  set[(namc <- names(set.default))] <- set.default
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in set.default: ", paste(noNms, 
                                                    collapse = ", "))
  if (!is.null(weights)) 
    rval <- speedlm.wfit(y, X, offset = offset, w = weights, 
                         camp = set$camp, sparselim = set$sparselim, eigendec = set$eigendec, 
                         row.chunk = set$row.chunk, tol.solve = set$tol.solve, 
                         sparse = sparse, tol.values = set$tol.values, tol.vectors = set$tol.vectors, 
                         method = set$method, intercept = attr(tf, "intercept"))
  else rval <- speedlm.fit(y, X, offset = offset, row.chunk = set$row.chunk, 
                           camp = set$camp, sparselim = set$sparselim, eigendec = set$eigendec, 
                           tol.solve = set$tol.solve, sparse = sparse, tol.values = set$tol.values, 
                           tol.vectors = set$tol.vectors, method = set$method, intercept = attr(tf, 
                                                                                                "intercept"))
  rval$terms <- tf
  rval$call <- call
  if (ncol(M)>1) {
    for (i in 2:ncol(M)) {
      if (is.factor(M[, i])) 
        eval(parse(text = paste("rval$levels$'", names(M)[i], 
                                "'", "<-levels(M[,i])", sep = "")))
    }
  }
  if(model) rval$model <- M
  if(target) rval$y <- y
  class(rval) <- "speedlm"
  if(fitted) rval$fitted.values <- predict.speedlm(rval,M) 
  rval
}

formula.speedlm <- function (x, ...) 
{
    form <- formula(x$terms)
    environment(form) <- environment(x$formula)
    form
}

update.speedlm <- function(object, formula, data, add=TRUE, evaluate=FALSE, 
                           subset=NULL, ...) {
  if (!inherits(object, "speedlm")) 
    stop("object must be of class speedlm")
  if ((!missing(formula))&(!missing(data))&(add)) 
    stop('cannot specify a formula while adding new data')
  if ((!missing(data))&(add)) updateWithMoreData(object, data, subset=subset,...) else{
    if (missing(data)) update.default(object, formula, evaluate=evaluate) else 
      update.default(object, formula, data=data, evaluate=evaluate,subset=subset, ...)
  }  
}

updateWithMoreData <- function(object, data, weights = NULL, offset = NULL,
                               sparse = NULL, all.levels = FALSE, 
                               set.default = list(), subset=NULL,...){
  if (!inherits(object, "speedlm")) 
    stop("object must be of class speedlm")
  if (is.null(call <- getCall(object))) 
    stop("need an object with call component")
  M <- match.call(expand.dots = F)
  formula <- as.formula(object$call[[2]])
  M$formula <- formula 
  m <- match(c("formula", "data", "subset"), names(M), 0L)
  M <- M[c(1L, m)]
  M$drop.unused.levels <- TRUE
  M[[1L]] <- quote(stats::model.frame)
  M <- eval(M, parent.frame())
  y <- M[[1]]    
  set <- list(sparselim = 0.9, camp = 0.01, eigendec = TRUE, 
              row.chunk = NULL, tol.solve = .Machine$double.eps, tol.values = 1e-07, 
              tol.vectors = 1e-07, method = object$method)
  nmsC <- names(set)
  set[(namc <- names(set.default))] <- set.default
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in set.default: ", paste(noNms, 
                                                    collapse = ", "))
  fa <- which(attributes(attributes(M)$terms)$dataClasses=="factor")
  if ((!(all.levels))&(length(fa)>0)) {
    flevels <- list()
    j <- 0
    for (i in 1:length(fa)) {
      j <- j + 1
      eval(parse(text = paste("flevels$'", names(M)[fa[i]],
                              "'", "<-levels(M[,fa[i]])", sep = "")))
      a <- c(object$levels[[j]][!(object$levels[[j]] %in% flevels[[j]])], flevels[[j]])
      flevels[[j]] <- a
    }
    if (!is.null(subset)){
      M <- model.frame(formula, data, subset=subset, drop.unused.levels=T,xlev = flevels)
      X <- model.matrix(formula, M)
      object$levels <- flevels
   } else {
     M <- model.frame(formula, data, xlev = flevels)
     X <- model.matrix(formula, M, xlev = flevels)
     object$levels <- flevels
   }
  } else {
    X <- model.matrix(object$terms, M)
    flevels <- object$levels
  }
  pw <- if (is.null(weights)) weights else 
    prod(weights[weights != 0]) + object$pw
  w <- weights
  zero.w <- sum(w == 0)
  offset <- model.offset(M)
  colnam <- colnames(X)
  if (is.null(sparse)) 
    sparse <- is.sparse(X = X, set$sparselim, set$camp)
  if (sparse) 
    X <- as(X, "dgCMatrix")
  if (!is.null(w)) {
    if (object$rank == ncol(X)) {
      XTX <- if (sparse) crossprod(sqrt(w) * X) else cp(X, w, set$rowchunk)
      XTX <- as(XTX, "matrix") + object$XTX
      A <- XTX
      Xy <- as(object$Xy, "matrix") + as(crossprod(X, (w * y)), "matrix")
      rank <- object$rank
      ok <- 1:ncol(X)
      nvar <- object$nvar
    } else {
      A <- if (sparse) crossprod(sqrt(w) * X) else cp(X, w, set$rowchunk)
      A <- as(A, "matrix")
      A[rownames(object$A), colnames(object$A)] <- A[rownames(object$A), 
                                                     colnames(object$A)] + object$A
      ris <- control(A, , set$tol.values, set$tol.vectors,, set$method)
      XTX <- ris$XTX
      ok <- ris$pivot[1:ris$rank]
      Xy <- as(crossprod(X, (w * y)), "matrix")
      Xy[rownames(object$A), ] <- Xy[rownames(object$A), ] + as(object$Xy, "matrix")
      rank <- ris$rank
      nvar <- length(ris$pivot)
    }
    sqw <- sqrt(w)
    coef <- solve(XTX, Xy[ok], tol = set$tol.solve)
    yy <- object$yy + crossprod(sqw * y)
    X1X <- crossprod(sqw, X)
    names(X1X) <- colnam
    X1X[names(object$X1X)] <- X1X[names(object$X1X)] + object$X1X
    XW1 <- crossprod(w, X)
    names(XW1) <- colnam
    XW1[names(object$X1X)] <- XW1[names(object$X1X)] + object$XW1
    SW <- sum(w) + object$SW
    dfr <- object$df.residual + nrow(X) - zero.w
  } else {
    if (object$rank == ncol(X)) {
      XTX <- if (sparse) 
        crossprod(X) else cp(X, w = NULL, set$rowchunk)
      XTX <- as(XTX, "matrix") + object$XTX
      A <- XTX
      Xy <- as(object$Xy, "matrix") + as(crossprod(X, y), "matrix")
      rank <- object$rank
      ok <- 1:ncol(X)
      nvar <- object$nvar
    } else {
      A <- if (sparse) crossprod(X) else cp(X, w = NULL, set$rowchunk)
      A <- as(A, "matrix")

      A[rownames(object$A), colnames(object$A)] <- A[rownames(object$A), 
                                                     colnames(object$A)] + object$A
      ris <- control(A, , set$tol.values, set$tol.vectors,,set$method)
      XTX <- ris$XTX
      ok <- ris$pivot[1:ris$rank]
      Xy <- as(crossprod(X, y), "matrix")
      Xy[rownames(object$A), ] <- Xy[rownames(object$A), 
                                     ] + as(object$Xy, "matrix")
      rank <- ris$rank
      nvar <- length(ris$pivot)
    }
    coef <- solve(XTX, Xy[ok], tol = set$tol.solve)
    yy <- object$yy + crossprod(y)
    X1X <- colSums(X)
    names(X1X) <- colnam
    X1X[names(object$X1X)] <- X1X[names(object$X1X)] + object$X1X
    dfr <- object$nobs + nrow(X) - rank
    XW1 <- SW <- NULL
  }
  RSS <- yy - 2 * crossprod(coef, Xy[ok]) + crossprod(coef, 
                                                      XTX) %*% coef
  coefficients <- rep(NA, nvar)
  coefficients[ok] <- coef
  names(coefficients) <- colnames(X)
  rval <- list(coefficients = coefficients, coef = coef, df.residual = dfr, 
               X1X = X1X, Xy = Xy, XW1 = XW1, SW = SW, yy = yy, 
               A = as(A,"matrix"), nobs = object$nobs + nrow(X), RSS = as.numeric(RSS), 
               rank = rank, ok = ok, nvar = nvar, weights = weights, 
               zero.w = zero.w + object$zero.w, pw = pw, XTX = as(XTX,"matrix"), 
               sparse = sparse, "intercept"=object$intercept,method=set$method)
  rval$terms <- object$terms
  rval$call <- call
  rval$levels <- flevels
  class(rval) <- "speedlm"
  rval
}

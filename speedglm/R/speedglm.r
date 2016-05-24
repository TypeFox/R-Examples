speedglm <- function(formula,data,family=gaussian(),weights=NULL,start=NULL,
                     etastart=NULL,mustart=NULL,offset=NULL,maxit=25, k=2, 
                     sparse=NULL,set.default=list(), trace=FALSE,
                     method=c('eigen','Cholesky','qr'), model=FALSE, y=FALSE, 
                     fitted=FALSE,...){
  call <- match.call()
  target <- y
  M <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset"), names(M), 0L)
  M <- M[c(1L, m)]
  M$drop.unused.levels <- TRUE
  M[[1L]] <- quote(stats::model.frame)
  M <- eval(M, parent.frame())
  y  <- M[[1]] 
  tf <- terms(formula,data=data)
  X  <- model.matrix(tf,M)
  offset <- model.offset(M)
  intercept <- attributes(tf)$intercept    
  set <- list(sparselim=.9,camp=.01,eigendec=TRUE,row.chunk=NULL,
              tol.solve=.Machine$double.eps,acc=1e-8,tol.values=1e-7,
              tol.vectors=1e-7, method = match.arg(method))
  nmsC <- names(set)
  set[(namc <- names(set.default))] <- set.default
  if (length(noNms <- namc[!namc %in% nmsC]) > 0) 
    warning("unknown names in set.default: ", paste(noNms, collapse = ", "))  
  rval <- speedglm.wfit(y=y,X=X,family=family,weights=weights,start=start,
                        etastart=etastart,mustart=mustart,offset=offset, 
                        intercept=intercept,row.chunk=set$row.chunk,maxit=maxit,
                        k=k,acc=set$acc,sparselim=set$sparselim,camp=set$camp,
                        eigendec=set$eigendec,tol.solve=set$tol.solve,
                        sparse=sparse,tol.values=set$tol.values,trace=trace,
                        tol.vectors=set$tol.vectors, method=set$method) 
  rval$terms <- tf 
  rval$call <- call
  class(rval)<- c("speedglm","speedlm")
  if (model) rval$model <- M
  if (fitted) rval$linear.predictors <- predict.speedlm(rval, newdata=M)
  if (target) rval$y <- y 
  if ((rval$iter==maxit)&(!rval$convergence)) 
    warning("Maximum number of iterations reached without convergence")     
  rval  
}

speedglm.wfit <- function(y, X, intercept=TRUE, weights=NULL,row.chunk=NULL,
                          family=gaussian(), start=NULL, etastart=NULL,
                          mustart=NULL, offset=NULL, acc=1e-08, maxit=25, k=2,
                          sparselim=.9,camp=.01, eigendec=TRUE, tol.values=1e-7,
                          tol.vectors=1e-7, tol.solve=.Machine$double.eps,
                          sparse=NULL,method = c('eigen','Cholesky','qr'), 
                          trace=FALSE,...){
  nobs <- NROW(y)
  nvar <- ncol(X) 
  if (missing(y)) stop("Argument y is missing")
  if (missing(X)) stop("Argument X is missing") 
  if (is.null(offset)) offset <- rep.int(0, nobs)  
  if (is.null(weights)) weights <- rep(1, nobs)  
  col.names <- dimnames(X)[[2]]
  method <- match.arg(method)
  fam  <- family$family
  link <- family$link  
  variance <- family$variance
  dev.resids <- family$dev.resids
  aic <- family$aic
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta 
  if (is.null(sparse)) sparse <- is.sparse(X=X,sparselim,camp) 
  if (is.null(start)) {
    if (is.null(mustart)) eval(family$initialize)
    eta <- if (is.null(etastart)) family$linkfun(mustart) else etastart
    mu <- mustart
    start <- rep(0,nvar)
  } else { 
    eta <- offset + as.vector(if (nvar == 1) X*start else {
      if (sparse) X%*%start else tcrossprod(X,t(start))})
    mu <- linkinv(eta)
  }  
  iter <- 0
  dev <- sum(dev.resids(y, mu, weights)) 
  tol <- 1
  if ((fam=="gaussian")&(link=="identity")) maxit <- 1
  C_Cdqrls <- getNativeSymbolInfo('Cdqrls', PACKAGE=getLoadedDLLs()$stats)
  while((tol>acc)&(iter<maxit)){
    iter <- iter + 1
    beta <- start
    dev0 <- dev
    varmu <- variance(mu)
    mu.eta.val <- mu.eta(eta)    
    z <- (eta - offset) + (y - mu)/mu.eta.val    
    W <- (weights*mu.eta.val*mu.eta.val)/varmu 
    XTX <- cp(X,W,row.chunk,sparse)
    XTz <- if (sparse) t(X)%*%(W*z) else t(crossprod((W*z),X))    
    if (iter==1 & method != 'qr') {
      variable <- colnames(X) 
      ris <- if (eigendec) control(XTX,,tol.values,tol.vectors,,method) else 
        list("rank"= nvar,"pivot"=1:nvar)               
      ok <- ris$pivot[1:ris$rank]
      if (eigendec) {
        XTX <- ris$XTX
        X <- X[,ok]
        XTz<-XTz[ok]
        start <- start[ok]
      }     
      beta <- start           
    }  
    if(method=='qr') {
      ris<-.Call(C_Cdqrls,XTX, XTz, tol.values, FALSE)
      start <- if (ris$rank < nvar) ris$coefficients[ris$pivot]
      else ris$coefficients
    } else {
      start <- solve(XTX,XTz,tol=tol.solve)
    }
    eta <- if (sparse) drop(X%*%start) else drop(tcrossprod(X,t(start)))
    mu <- linkinv(eta <- eta + offset)
    dev <- sum(dev.resids(y, mu, weights)) 
    tol <- max(abs(dev0-dev)/(abs(dev)+0.1))
    if (trace) cat("iter",iter,"tol",tol,"\n")
  } 
  wt <- sum(weights)
  wtdmu <- if (intercept) sum(weights * y)/wt  else linkinv(offset)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  n.ok <- nobs - sum(weights == 0)
  nulldf <- n.ok - as.integer(intercept) 
  rank <- ris$rank
  dfr <- nobs - rank - sum(weights == 0)
  aic.model <- aic(y, nobs, mu, weights, dev) + k * rank
  ll.nuovo <- ll.speedglm(fam,aic.model,rank)
  res <- (y-mu)/mu.eta(eta)
  resdf <- n.ok - rank
  RSS <- sum(W*res*res)
  var_res <- RSS/dfr  
  dispersion <- if (fam %in% c("poisson","binomial")) 1 else var_res   
  if(method=='qr'){
    coefficients<-start
    coefficients[coefficients==0]=NA
    ok <- ris$pivot[1:rank]
  } 
  else {
    coefficients <- rep(NA,nvar)
    start<- as(start,"numeric")
    coefficients[ok] <- start
  }
  names(coefficients) <- col.names
  rval <- list("coefficients"=coefficients,"logLik"=ll.nuovo,"iter"=iter,
               "tol"= tol,"family"=family,"link"=link,"df"=dfr,"XTX"=XTX,
               "dispersion"=dispersion,"ok"=ok,"rank"=rank,"RSS"=RSS, method=method,
               "aic"=aic.model, "sparse"=sparse,"deviance"=dev,"nulldf"=nulldf, 
               "nulldev"=nulldev,"ngoodobs"=n.ok,"n"=nobs,"intercept"=intercept,
               "convergence"=(!(tol>acc)))   
  class(rval) <- "speedglm"
  rval
}





shglm <- function(formula, datafun, family = gaussian(), weights.fo = NULL,
                  start = NULL, etastart = NULL, mustart = NULL, offset = NULL,
                  maxit = 25, k = 2, chunksize = 5000, sparse = NULL,trace=FALSE, 
                  all.levels = FALSE, set.default = list(), ...){
  if (!is.null(start)) stop("Sorry, code for argument start is not implemented yet")
  if (!is.null(mustart)) stop("Sorry, code for argument mustart is not implemented yet")
  if (!is.null(etastart)) stop("Sorry, code for argument etastart is not implemented yet")
  dati <- datafun(reset = TRUE)
  dati <- datafun(reset = FALSE)
  call <- match.call()
  tf <- terms(formula, data=dati)
  M <- model.frame(tf, dati)
  y <- M[[1]]
  X <- model.matrix(tf, M)
  offset <- model.offset(M)
  set <- list(sparselim = 0.9, camp = 0.01, eigendec = TRUE,
              row.chunk = NULL, tol.solve = .Machine$double.eps, acc = 1e-08,
              tol.values = 1e-07, tol.vectors = 1e-07, method = "eigen")
  nmsC <- names(set)
  set[(namc <- names(set.default))] <- set.default
  if (length(noNms <- namc[!namc %in% nmsC]) > 0)
    warning("unknown names in set.default: ", paste(noNms,
                                                    collapse = ", "))
  obj <- list()
  obj$terms <- tf
  fa <- NULL
  if (!all.levels) {
    fa <- which(attributes(attributes(M)$terms)$dataClasses=="factor")
    if (length(fa)>0){
      for (i in 1:length(fa)) {
        eval(parse(text = paste("obj$levels$'", names(M)[fa[i]],
                                "'", "<-levels(M[,fa[i]])", sep = "")))
      }
    }
  }
  nomicol <- colnames(dati)
  variable <- colnames(X)
  nobs <- length(y)
  weights <- if (is.null(weights.fo))
    rep(1, nobs) else model.frame(weights.fo, dati)[[1]]
  intercept <- attributes(tf)$intercept
  nvar <- ncol(X)
  if (is.null(offset))
    offset <- rep.int(0, nobs)
  if (is.null(sparse)) {
    sp <- NULL
    sparse <- is.sparse(X = X, set$sparselim, set$camp)
  } else sp <- sparse
  ok <- 1:ncol(X)
  rank <- ncol(X)
  dfr <- nobs - rank - sum(weights == 0)
  fam <- family$family
  link <- family$link
  linkfun <- family$linkfun
  variance <- family$variance
  dev.resids <- family$dev.resids
  linkinv <- family$linkinv
  mu.eta <- family$mu.eta
  iter <- 0
  tol <- 1
  dev <- nobs
  weights.cum <- wt <- wy <- nulldev <- tot.obs <- zero.weights <- 0
  block.cum <- nrow(X)
  end.iter <- FALSE
  
  while ((tol > set$acc) & (iter < maxit)) {
    dev.prec <- dev
    aic.model <- dev <- Wy <- RSS <- 0
    iter <- iter + 1
    XTX <- matrix(0, ncol(X), ncol(X))
    XTz <- matrix(0, ncol(X), 1)
    eof <- FALSE
    iter2 <- 0
    while (!eof) {
      iter2 <- iter2 + 1
      if (iter > 1) {
        if (is.null(sp))
          sparse <- is.sparse(X = X, set$sparselim, set$camp)
        if (nvar == 1) {
          eta <- (offset + as.vector(X * start))[1:length(y)]
        } else {
          eta <- (offset + if (sparse)
            drop(X[, ok] %*% start) else drop(tcrossprod(X[, ok], t(start))))[1:length(y)]
        }
        mu <- linkinv(eta <- eta)
      } else {            
        if (is.null(mustart)) eval(family$initialize)
        if (is.null(start)) {
          eta <- (if (is.null(etastart)) linkfun(mustart) else etastart)[1:length(y)]
          mu <- mustart[1:length(y)]
        } else {
          if (nvar == 1) {
            eta <- (offset + as.vector(X * start))[1:length(y)]
          } else {
            eta <- (offset + if (sparse)
              drop(X %*% start) else as.vector(tcrossprod(X, t(start))))[1:length(y)]
          }
          mu <- linkinv(eta)
        }
      }
      if (trace) cat("iter",iter,"elaborating chunk",iter2,"dim",dim(X),"\n")
      dev <- dev + sum(dev.resids(y, mu, weights))
      if (iter>1) aic.model <- aic.model + aic.shglm(fam,y,wt,mu,weights,dev.prec)
      varmu <- variance(mu)
      mu.eta.val <- mu.eta(eta)
      z <- (eta - offset) + (y - mu)/mu.eta.val
      W <- (weights * mu.eta.val * mu.eta.val)/varmu
      if (iter2 == 1) {
        XTX <- XTX + cp(X, W, set$row.chunk, sparse)
        XTz <- XTz + if (sparse) t(X)%*%(W*z) else t(crossprod((W*z),X)) 
      } else {
        Ax <- XTX
        XTX <- cp(X, W, set$row.chunk, sparse)
        XTX[rownames(Ax),colnames(Ax)]<-XTX[rownames(Ax),colnames(Ax)]+Ax
        Az <- XTz
        XTz <- if (sparse) t(X)%*%(W*z) else t(crossprod((W*z), X))
        XTz[rownames(Az), ] <- XTz[rownames(Az), ] + Az
      }
      res <- (y - mu)/mu.eta(eta)
      RSS <- RSS + sum(W * res * res)
      if (iter == 1) weights.cum <- weights.cum + sum(weights == 0)           
      if (is.null(dati <- datafun(reset = FALSE))) {
        if (iter == 1) {
          ris <- if (set$eigendec)
            control(XTX,,set$tol.values,set$tol.vectors,out.B=FALSE,set$method) else
              list(rank = nvar, pivot = 1:nvar)
          ok <- ris$pivot[1:ris$rank]
        }
        if (is.null(start)) start <- rep(0, rank)
        beta <- if (iter == 1) start[ok] else start
        start <- solve(XTX[ok, ok], XTz[ok], tol = set$tol.solve)
        if (trace) cat("coef",start,"tol",tol,"\n")
        tol <- max(abs(dev.prec - dev)/(abs(dev) + 0.1))        
        if ((tol > set$acc) & (iter < maxit)) {
          dati <- datafun(reset = TRUE)
          dati <- datafun(reset = FALSE)
          eof <- TRUE
        } else break
      }
      colnames(dati) <- nomicol
      M <- model.frame(tf, dati)
      y <- M[[1]]
      if ((!(all.levels))&(length(fa)>0)) {
        flevels <- list()
        j <- 0
        for (i in 1:length(fa)) {
          j <- j + 1
          eval(parse(text = paste("flevels$'", names(M)[fa[i]],
                                  "'", "<-levels(M[,fa[i]])", sep = "")))
          a <- c(obj$levels[[j]][!(obj$levels[[j]] %in%
                                     flevels[[j]])], flevels[[j]])
          flevels[[j]] <- sort(a)
        }
        M <- model.frame(obj$terms, dati, xlev = flevels)
        X <- model.matrix(obj$terms, M, xlev = flevels)
        obj$levels <- flevels
      } else {
        X <- model.matrix(obj$terms, M)
        flevels <- obj$levels
      }
      offset <- model.offset(M)
      nobs <- length(y)
      if (is.null(offset))
        offset <- rep.int(0, nobs)
      weights <- if (is.null(weights.fo))
        rep(1, nobs) else model.frame(weights.fo, dati)[[1]]
      if (iter == 1) {
        tot.obs <- tot.obs + nobs
        wt <- wt + sum(weights)
        wy <- wy + crossprod(weights, y)
        zero.weights <- zero.weights + sum(weights == 0)
        block.cum <- block.cum + nrow(X)
      }
      if (iter == 2) {
        wtdmu <- if (intercept) wy/wt  else linkinv(offset)
        nulldev <- nulldev + sum(dev.resids(y, wtdmu, weights))
      }
    }
  }
  datafun(reset = TRUE)
  rank <- ris$rank
  n.ok <- tot.obs - zero.weights
  nulldf <- n.ok - as.integer(intercept)
  aic.rest <- ifelse((fam %in% c("Gamma", "inverse.gaussian",
                                 "gaussian")), 2, 0)
  aic.model <- aic.model + k * rank + aic.rest
  ll.nuovo <- ll.speedglm(fam, aic.model, rank)
  resdf <- n.ok - rank
  var_res <- RSS/resdf
  dispersion <- if (fam %in% c("poisson", "binomial")) 1  else var_res
  coefficients <- rep(NA, ncol(X))
  start <- as(start, "numeric")
  coefficients[ok] <- start
  names(coefficients) <- colnames(X)
  rval <- list(coefficients = coefficients, logLik = ll.nuovo,
               iter = iter, tol = tol, family = family, link = link, df = resdf,
               XTX = XTX[ok, ok], dispersion = dispersion, nok = ris$nok,
               ok = ok, RSS = RSS, ncoll = ris$ncoll, sparse = sparse,
               nulldev = nulldev, rank = rank, deviance = dev,
               nulldf = nulldf, ngoodobs = n.ok, n = tot.obs, intercept = intercept,
               aic = aic.model, convergence = (!(tol> set$acc)),method="eigen")
  rval$tf <- tf
  rval$call <- call
  if ((rval$iter==maxit)&(!rval$convergence))
    warning("Maximum number of iterations reached without convergence")
  class(rval) <- "speedglm"
  rval
}


#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

momentEstim <- function(object, ...)
  {
  UseMethod("momentEstim")
  }

momentEstim.baseGmm.twoStep <- function(object, ...)
  {
  P <- object
  x <- P$x
  if (P$optfct == "optimize")
    {
    n = nrow(P$g(P$t0[1], x))
    q = ncol(P$g(P$t0[1], x))
    k = 1
    }
  else
    {
    n = nrow(P$g(P$t0, x))
    q = ncol(P$g(P$t0, x))
    k = length(P$t0)
    }
  k2 <- k
  df <- q - k
  w=diag(q)
  if (P$optfct == "optim")
    {
    if (P$gradvf)
      {
      gradvOptim <- P$gradv
      gr2 <- function(thet, x,  w, gf, INV)
		{
		gt <- gf(thet, x)
		Gbar <- gradvOptim(thet, x) 
		gbar <- as.vector(colMeans(gt))
		if (INV)		
		  	obj <- crossprod(Gbar, solve(w, gbar))
		else
			obj <- crossprod(Gbar,w)%*%gbar
		return(obj*2)
		}
      argDots <- list(...)
      allArgOptim <- list(par = P$t0, fn = .obj1, gr = gr2, x = P$x, w = w, gf = P$g, INV = TRUE)
      allArgOptim <- c(allArgOptim,argDots)
      res <- do.call(optim,allArgOptim)
      } else {
      res <- optim(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
	}
    }
  if (P$optfct == "constrOptim")
	{
	if (!any(c("ui","ci") %in% names(list(...))))
		stop("You must specify ui and ci when optfct is set to constrOptim")
	argDots <- list(...)
	ui <- argDots$ui
	ci <- argDots$ci
	argDots$ui <- NULL
	argDots$ci <- NULL
      	allArgOptim <- list(theta = P$t0, f = .obj1, grad = NULL, ui = ui, ci = ci, x = P$x, w = w, gf = P$g, INV = TRUE)
        allArgOptim <- c(allArgOptim,argDots)
	res <- do.call(constrOptim,allArgOptim)
	}

  if (P$optfct == "nlminb")
    {
    res <- nlminb(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    res$value <- res$objective
    }
  if (P$optfct == "optimize")
    {
    res <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
    res$par <- res$minimum
    res$value <- res$objective
    }	
  if (q == k2 | P$wmatrix == "ident")
    {
    z = list(coefficients = res$par, objective = res$value, k=k, k2=k2, n=n, q=q, df=df)	
    if (P$optfct == "optim")
	z$algoInfo <- list(convergence = res$convergence, counts = res$counts, message = res$message)
    else if(P$optfct == "nlminb")
	z$algoInfo <- list(convergence = res$convergence, counts = res$evaluations, message = res$message)
    }
  else
    {
    if (P$vcov == "iid")
      w <- P$iid(res$par, P$x, P$g, P$centeredVcov)

    if (P$vcov == "HAC")
     {
     gmat <- P$g(res$par, P$x)
     w <- .myKernHAC(gmat, P)
     }

    if (P$optfct == "optim")
      {
      argDots <- list(...)
      if (P$gradvf)
        {
        gradvOptim <- P$gradv
        gr2 <- function(thet, x,  w, gf, INV)
		{
		gt <- gf(thet, x)
		Gbar <- gradvOptim(thet, x) 
		gbar <- as.vector(colMeans(gt))
		if (INV)		
		  	obj <- crossprod(Gbar, solve(w, gbar))
		else
			obj <- crossprod(Gbar,w)%*%gbar
		return(obj*2)
		}
	argDots <- list(...)
        allArgOptim <- list(par = P$t0, fn = .obj1, gr = gr2, x = P$x, w = w, gf = P$g, INV = TRUE)
        allArgOptim <- c(allArgOptim,argDots)
        res2 <- do.call(optim,allArgOptim)
        }
      else
        res2 <- optim(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
      }
  if (P$optfct == "constrOptim")
	{
	if (!any(c("ui","ci") %in% names(list(...))))
		stop("You must specify ui and ci when optfct is set to constrOptim")
	argDots <- list(...)
	ui <- argDots$ui
	ci <- argDots$ci
	argDots$ui <- NULL
	argDots$ci <- NULL
      	allArgOptim <- list(theta = P$t0, f = .obj1, grad = NULL, ui = ui, ci = ci, x = P$x, w = w, gf = P$g, INV = TRUE)
        allArgOptim <- c(allArgOptim,argDots)
	res2 <- do.call(constrOptim,allArgOptim)
	}
    if (P$optfct == "nlminb")
      {
      res2 <- nlminb(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
      res2$value <- res2$objective
      }

    if (P$optfct == "optimize")
      {
      res2 <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
      res2$par <- res2$minimum
      res2$value <- res2$objective
      }	

    z = list(coefficients = res2$par, objective = res2$value, k=k, k2=k2, n=n, q=q, df=df, initTheta = res$par)	
    if (P$optfct == "optim")
	z$algoInfo <- list(convergence = res2$convergence, counts = res2$counts, message = res2$message)
    else if(P$optfct == "nlminb")
	z$algoInfo <- list(convergence = res2$convergence, counts = res2$evaluations, message = res2$message)
    }

  z$dat <- P$x
  z$gt <- P$g(z$coefficients, P$x)
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  z$WSpec <- P$WSpec

  names(z$coefficients) <- P$namesCoef
  if (is.null(colnames(z$gt)))
	colnames(z$gt) <- paste("gt",1:ncol(z$gt),sep="")
  class(z) <- paste(P$TypeGmm,".res",sep="")	
  z$specMod <- P$specMod 
  z$w0 <- w
  return(z)
  }

momentEstim.baseGmm.twoStep.formula <- function(object, ...)
  {
  P <- object
  g <- P$g
  dat <- P$x
  
  x <- dat$x
  k <- dat$k
  k2 <- k*dat$ny
  n <- nrow(x)
  q <- dat$ny*dat$nh
  df <- q-k*dat$ny
  w <- diag(q)

  if (q == k2 | P$wmatrix == "ident")
    {
    res2 <- .tetlin(dat, w, P$gradv, P$g)
    z = list(coefficients = res2$par, objective = res2$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df, df.residual = (n-k))
    }
  else
    {  
    if (P$vcov == "iid")
   	{            
      res1 <- .tetlin(dat, w, P$gradv, P$g, type="2sls")
      initTheta <- res1$par
      gmat <- g(res1$par, dat)
      w <- crossprod(gmat)/n
      res2 <- .tetlin(dat, w, P$gradv, g)	
      res2$firstStageReg <- res1$firstStageReg
      res2$fsRes <- res1$fsRes
      }
    if (P$vcov == "HAC")
      {          
      res1 <- .tetlin(dat, w, P$gradv, P$g, type="2sls")
      initTheta <- res1$par
      gmat <- g(res1$par, dat)
      w <- .myKernHAC(gmat, P)
      res2 <- .tetlin(dat, w, P$gradv, g)
      }
    z = list(coefficients = res2$par, objective = res2$value, dat=dat, k=k, k2=k2, n=n, q=q, df=df, initTheta = initTheta, df.residual = (n-k))	
    }
  z$gt <- g(z$coefficients, dat) 
  b <- z$coefficients
  y <- as.matrix(dat$x[,1:dat$ny])
  b <- t(matrix(b, nrow = dat$ny))
  x <- as.matrix(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
  yhat <- x %*% b
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  z$WSpec <- P$WSpec
  z$w0 <- w
  names(z$coefficients) <- P$namesCoef
  colnames(z$gt) <- P$namesgt
 
  if (P$vcov == "iid" & P$wmatrix != "ident")
	z$fsRes <- res2$fsRes

  class(z) <- paste(P$TypeGmm,".res",sep="")
  z$specMod <- P$specMod
  return(z)	
  }

momentEstim.baseGmm.iterative.formula <- function(object, ...)
  {
  P <- object
  g <- P$g
  dat <- P$x
  x <- dat$x
  k <- dat$k
  k2 <- k*dat$ny
  n <- nrow(x)
  q <- dat$ny*dat$nh
  df <- q-k*dat$ny
  w <- diag(q)
  
  if (q == k2 | P$wmatrix == "ident")
    {
    res <- .tetlin(dat, w, P$gradv, g)
    z = list(coefficients = res$par, objective = res$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df, df.residual = (n-k))
    }
  else
    {
    res <- .tetlin(dat, w, P$gradv, g, type="2sls")
    fsRes <- res$fsRes
    initTheta <- res$par
    ch <- 100000
    j <- 1
    while(ch > P$crit)
      {
      tet <- res$par
      if (P$vcov == "HAC")
	{
        gmat <- g(tet, dat)
	if (!is.null(attr(w,"Spec")))
		P$WSpec$sandwich$bw <- attr(w,"Spec")$bw
	w <- .myKernHAC(gmat, P)
        }
      if(P$vcov == "iid")
	{
        gmat <- g(tet, dat)
	w <- crossprod(gmat)/n
	}
      res <- .tetlin(dat, w, P$gradv, g)
      ch <- crossprod(abs(tet- res$par)/tet)^.5
      if (j>P$itermax)
        {
        cat("No convergence after ", P$itermax, " iterations")
        ch <- P$crit
        }
	if(P$traceIter)
		cat("Iter :",j,": value=",res$value,", Coef = ", res$par,"\n") 
        j <- j+1	
      }
    z = list(coefficients = res$par, objective = res$value, dat=dat, k=k, k2=k2, n=n, q=q, df=df, initTheta=initTheta, df.residual = (n-k))	
   }
  z$gt <- g(z$coefficients, dat) 
  b <- z$coefficients
  y <- as.matrix(dat$x[,1:dat$ny])
  b <- t(matrix(b, nrow = dat$ny))
  x <- as.matrix(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
  yhat <- x %*% b
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  z$WSpec <- P$WSpec
  z$w0 <- w

  names(z$coefficients) <- P$namesCoef
  colnames(z$gt) <- P$namesgt
  if (P$vcov == "iid" & P$wmatrix != "ident")
	z$fsRes <- fsRes
  class(z) <- paste(P$TypeGmm,".res",sep="")
  z$specMod <- P$specMod
  return(z)	
  }

momentEstim.baseGmm.iterative <- function(object, ...)
  {
  P <- object
  x <- P$x
  if (P$optfct == "optimize")
    {
    n = nrow(P$g(P$t0[1], x))
    q = ncol(P$g(P$t0[1], x))
    k = 1
    }
  else
    {
    n = nrow(P$g(P$t0, x))
    q = ncol(P$g(P$t0, x))
    k = length(P$t0)
    }
  k2 <- k
  df <- q - k
  w=diag(q)
  if (P$optfct == "optim")
    {
    if (P$gradvf)
        {
        gradvOptim <- P$gradv
        gr2 <- function(thet, x,  w, gf, INV)
		{
		gt <- gf(thet, x)
		Gbar <- gradvOptim(thet, x) 
		gbar <- as.vector(colMeans(gt))
		if (INV)		
		  	obj <- crossprod(Gbar, solve(w, gbar))
		else
			obj <- crossprod(Gbar,w)%*%gbar
		return(obj*2)
		}
      argDots <- list(...)
      allArgOptim <- list(par = P$t0, fn = .obj1, gr = gr2, x = P$x, w = w, gf = P$g, INV = TRUE)
      allArgOptim <- c(allArgOptim,argDots)
      res <- do.call(optim,allArgOptim)
      }
    else
      res <- optim(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    }
 if (P$optfct == "constrOptim")
    {
	if (!any(c("ui","ci") %in% names(list(...))))
		stop("You must specify ui and ci when optfct is set to constrOptim")
	argDots <- list(...)
	ui <- argDots$ui
	ci <- argDots$ci
	argDots$ui <- NULL
	argDots$ci <- NULL
      	allArgOptim <- list(theta = P$t0, f = .obj1, grad = NULL, ui = ui, ci = ci, x = P$x, w = w, gf = P$g, INV = TRUE)
        allArgOptim <- c(allArgOptim,argDots)
	res <- do.call(constrOptim,allArgOptim)
    }
 
  if (P$optfct == "nlminb")
    {
    res <- nlminb(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    res$value <- res$objective
    }
  if (P$optfct == "optimize")
    {
    res <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
    res$par <- res$minimum
    res$value <- res$objective
    }	

  if (q == k2 | P$wmatrix == "ident")
    {
    z <- list(coefficients = res$par, objective = res$value, k=k, k2=k2, n=n, q=q, df=df)
    if (P$optfct == "optim")
	z$algoInfo <- list(convergence = res$convergence, counts = res$counts, message = res$message)
    else if(P$optfct == "nlminb")
	z$algoInfo <- list(convergence = res$convergence, counts = res$evaluations, message = res$message)
    }	
  else
    {
    initTheta = res$par
    ch <- 100000
    j <- 1
    while(ch > P$crit)
      {
      tet <- res$par
      if (P$vcov == "iid")
        w <- P$iid(tet, P$x, P$g, P$centeredVcov)
      if (P$vcov == "HAC")
	{
	gmat <- P$g(tet, P$x)
	if (!is.null(attr(w,"Spec")))
		P$WSpec$sandwich$bw <- attr(w,"Spec")$bw
	w <- .myKernHAC(gmat, P)
        }

      if (P$optfct == "optim")
        {
        if (P$gradvf)
          {
          gradvOptim <- P$gradv
          gr2 <- function(thet, x,  w, gf, INV)
		{
		gt <- gf(thet, x)
		Gbar <- gradvOptim(thet, x) 
		gbar <- as.vector(colMeans(gt))
		if (INV)		
		  	obj <- crossprod(Gbar, solve(w, gbar))
		else
			obj <- crossprod(Gbar,w)%*%gbar
		return(obj*2)
		}
   	  argDots <- list(...)
          allArgOptim <- list(par = tet, fn = .obj1, gr = gr2, x = P$x, w = w, gf = P$g, INV = TRUE)
          argDots$gr <- NULL
          allArgOptim <- c(allArgOptim,argDots)
          res <- do.call(optim,allArgOptim)
          }
        else
          res <- optim(tet, .obj1, x = P$x, w = w, gf = P$g, ...)
        }
      if (P$optfct == "constrOptim")
        {
	if (!any(c("ui","ci") %in% names(list(...))))
		stop("You must specify ui and ci when optfct is set to constrOptim")
	argDots <- list(...)
	ui <- argDots$ui
	ci <- argDots$ci
	argDots$ui <- NULL
	argDots$ci <- NULL
      	allArgOptim <- list(theta = tet, f = .obj1, grad = NULL, ui = ui, ci = ci, x = P$x, w = w, gf = P$g, INV = TRUE)
        allArgOptim <- c(allArgOptim,argDots)
	res <- do.call(constrOptim,allArgOptim)
        }
      if (P$optfct == "nlminb")
        {
        res <- nlminb(tet, .obj1, x = P$x, w = w, gf = P$g, ...)
        res$value <- res$objective
        }
      if (P$optfct == "optimize")
        {
        res <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
        res$par <- res$minimum
        res$value <- res$objective
        }	
        ch <- crossprod(tet-res$par)^.5/(1+crossprod(tet)^.5)
        if (j>P$itermax)
          {
          cat("No convergence after ", P$itermax, " iterations")
          ch <- P$crit
          }
	if(P$traceIter)
		cat("Iter :",j,": value=",res$value,", Coef = ", res$par,"\n") 
        j <- j+1	
      }
    z = list(coefficients = res$par, objective = res$value,k=k, k2=k2, n=n, q=q, df=df, initTheta=initTheta)	
    if (P$optfct == "optim")
	z$algoInfo <- list(convergence = res$convergence, counts = res$counts, message = res$message)
    else if(P$optfct == "nlminb")
	z$algoInfo <- list(convergence = res$convergence, counts = res$evaluations, message = res$message)

    }

  z$dat <- P$x
  z$gt <- P$g(z$coefficients, P$x)
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  z$WSpec <- P$WSpec
  z$w0 <- w

  names(z$coefficients) <- P$namesCoef
  if (is.null(colnames(z$gt)))
	colnames(z$gt) <- paste("gt",1:ncol(z$gt),sep="")
  z$specMod <- P$specMod
  class(z) <- paste(P$TypeGmm,".res",sep="")	
  return(z)
  }

momentEstim.baseGmm.cue.formula <- function(object, ...)
  {

  P <- object
  g <- P$g
  
  dat <- P$x
  x <- dat$x
  k <- dat$k
  k2 <- k*dat$ny
  n <- nrow(x)
  q <- dat$ny*dat$nh
  df <- q-k*dat$ny
  w <- diag(q)
  P$weightMessage <- "Weights for kernel estimate of the covariance are fixed and based on the first step estimate of Theta"	

  if (q == k2 | P$wmatrix == "ident")
    {
    res <- .tetlin(dat, w, P$gradv, g)
    z = list(coefficients = res$par, objective = res$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df, df.residual = (n-k))
    P$weightMessage <- "No CUE needed because the model is just identified"
    }
  else
    {
    if (is.null(P$t0))
	{
	P$t0 <- .tetlin(dat, w, P$gradv, g, type="2sls")$par
	initTheta <- P$t0
	}
    else
	{
	initTheta <- P$t0
	}

    gt0 <- g(P$t0,dat)
    w <- .myKernHAC(gt0, P)
    P$WSpec$sandwich$bw <- attr(w,"Spec")$bw

    if (P$optfct == "optim")
      res2 <- optim(P$t0,.objCue, x = dat, P = P, ...)
    if (P$optfct == "constrOptim")
      {
	if (!any(c("ui","ci") %in% names(list(...))))
		stop("You must specify ui and ci when optfct is set to constrOptim")
	argDots <- list(...)
	ui <- argDots$ui
	ci <- argDots$ci
	argDots$ui <- NULL
	argDots$ci <- NULL
      	allArgOptim <- list(theta = P$t0, f = .objCue, grad = NULL, ui = ui, ci = ci, x = dat, P = P)
        allArgOptim <- c(allArgOptim,argDots)
	res <- do.call(constrOptim,allArgOptim)
      }
    if (P$optfct == "nlminb")
      {
      res2 <- nlminb(P$t0,.objCue, x = dat, P = P, ...)
      res2$value <- res2$objective
      }
    if (P$optfct == "optimize")
      {
      res2 <- optimize(.objCue,P$t0, x = dat, P = P, ...)
      res2$par <- res2$minimum
      res2$value <- res2$objective
      }
    z = list(coefficients = res2$par, objective = res2$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df, initTheta=initTheta, df.residual = (n-k))
    if (P$optfct == "optim")
	z$algoInfo <- list(convergence = res2$convergence, counts = res2$counts, message = res2$message)
    else if(P$optfct == "nlminb")
	z$algoInfo <- list(convergence = res2$convergence, counts = res2$evaluations, message = res2$message)
    }

  z$gt <- g(z$coefficients, dat) 
  b <- z$coefficients
  y <- as.matrix(dat$x[,1:dat$ny])
  b <- t(matrix(b, nrow = dat$ny))
  x <- as.matrix(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
  yhat <- x %*% b
  z$dat <- dat 
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  z$specMod <- P$specMod
  z$cue <- list(weights=P$fixedKernW,message=P$weightMessage)
  z$WSpec <- P$WSpec
  z$w0 <- w
  names(z$coefficients) <- P$namesCoef
  colnames(z$gt) <- P$namesgt

  class(z) <- paste(P$TypeGmm,".res",sep="")
  return(z)	
  }

momentEstim.baseGmm.cue <- function(object, ...)
  {
  P <- object
  x <- P$x

  res <- try(gmm(P$g,P$x,P$t0,wmatrix="ident",optfct=P$optfct, ...))
  if(class(res)=="try-error")
	stop("Cannot get a first step estimate to compute the weights for the Kernel estimate of the covariance matrix; try different starting values")

  initTheta <- res$coef
  n <- nrow(res$gt)
  q <- ncol(res$gt)
  w <- diag(q)

  if (P$optfct == "optimize")
    k = 1
  else
    k = length(P$t0)

  k2 <- k
  df <- q - k

  if (q == k2 | P$wmatrix == "ident")
    {
    z <- list(coefficients = res$coef, objective = res$objective, algoInfo = res$algoInfo, k=k, k2=k2, n=n, q=q, df=df)
    P$weightMessage <- "No CUE needed because the model if just identified or you set wmatrix=identity"
    }	
  else
    {
    gt0 <- P$g(res$coef,x) #Should we compute the weigths with the initial value provided?
    w <- .myKernHAC(gt0, P)
    P$WSpec$sandwich$bw <- attr(w,"Spec")$bw
    P$weightMessage <- "Weights for kernel estimate of the covariance are fixed and based on the first step estimate of Theta"

    if (P$optfct == "optim")
      {
      res2 <- optim(P$t0, .objCue, x = x, P = P, ...)
      }
   if (P$optfct == "constrOptim")
      {
	if (!any(c("ui","ci") %in% names(list(...))))
		stop("You must specify ui and ci when optfct is set to constrOptim")
	argDots <- list(...)
	ui <- argDots$ui
	ci <- argDots$ci
	argDots$ui <- NULL
	argDots$ci <- NULL
      	allArgOptim <- list(theta = P$t0, f = .objCue, grad = NULL, ui = ui, ci = ci, x = x, P = P)
        allArgOptim <- c(allArgOptim,argDots)
	res2 <- do.call(constrOptim,allArgOptim)
      }
    if (P$optfct == "nlminb")
      {
      res2 <- nlminb(P$t0, .objCue, x = x, P = P, ...)
      res2$value <- res2$objective
      }
    if (P$optfct == "optimize")
      {
      res2 <- optimize(.objCue,P$t0, x = x, P = P, ...)
      res2$par <- res2$minimum
      res2$value <- res2$objective
      }
    z = list(coefficients=res2$par,objective=res2$value, k=k, k2=k2, n=n, q=q, df=df, initTheta=initTheta)	
    if (P$optfct == "optim")
	z$algoInfo <- list(convergence = res2$convergence, counts = res2$counts, message = res2$message)
    else if(P$optfct == "nlminb")
	z$algoInfo <- list(convergence = res2$convergence, counts = res2$evaluations, message = res2$message)

    }

  z$dat <- P$x
  z$gradv <- P$gradv
  z$gt <- P$g(z$coefficients, P$x)
  z$iid <- P$iid
  z$g <- P$g
  z$cue <- list(weights=P$fixedKernW,message=P$weightMessage)
  names(z$coefficients) <- P$namesCoef
  if (is.null(colnames(z$gt)))
	colnames(z$gt) <- paste("gt",1:ncol(z$gt),sep="")
  z$WSpec <- P$WSpec

  z$specMod <- P$specMod
  class(z) <- paste(P$TypeGmm, ".res", sep = "")	
  return(z)
  }

momentEstim.baseGel.modFormula <- function(object, ...)
  {
  P <- object
  g <- P$g
  l0Env <- new.env()
  assign("l0",rep(0,P$q),envir=l0Env)

  n <- nrow(P$dat$x)
  if (!P$constraint)
    {
    if (P$optfct == "optim")
      res <- optim(P$tet0, .thetf, P = P, l0Env = l0Env, ...)
    if (P$optfct == "nlminb")
      res <- nlminb(P$tet0, .thetf, P = P, l0Env = l0Env, ...)
	
    if (P$optfct == "optimize")
      { 
      res <- optimize(.thetf, P$tet0, P = P, l0Env = l0Env, ...)
      res$par <- res$minimum
      res$convergence <- "There is no convergence code for optimize"
      }
    }

  if(P$constraint)
    res <- constrOptim(P$tet0, .thetf, grad = NULL, P = P, l0Env = l0Env, ...)

  All <- .thetf(res$par, P, "all",l0Env = l0Env)
  gt <- All$gt
  rlamb <- All$lambda

  z <- list(coefficients = res$par, lambda = rlamb$lambda, conv_lambda = rlamb$conv, conv_par = res$convergence, dat=P$dat)
  rho1 <- .rho(gt, z$lambda, derive = 1, type = P$typel, k = P$k1/P$k2)
  z$foc_lambda <- crossprod(colMeans(rho1*gt))

  z$type <- P$type
  z$gt <- gt
  rhom <- .rho(z$gt, z$lambda, type = P$typet, k = P$k1/P$k2)
  z$pt <- -.rho(z$gt, z$lambda, type = P$typet, derive = 1, k = P$k1/P$k2)/n
  # Making sure pt>0
  if (P$type=="CUE")
	{
	eps <- -length(z$pt)*min(min(z$pt),0)
	z$pt <- (z$pt+eps/length(z$pt))/(1+eps)
	}
  ###################

  z$conv_moment <- colSums(z$pt*z$gt)
  z$conv_pt <- sum(as.numeric(z$pt))
  z$objective <- sum(rhom - .rho(1, 0, type = P$typet))/n

  namex <- colnames(P$dat$x[, (P$dat$ny+1):(P$dat$ny+P$dat$k)])
  nameh <- colnames(P$dat$x[, (P$dat$ny+P$dat$k+1):(P$dat$ny+P$dat$k+P$dat$nh)])


  if (P$dat$ny > 1)
    {
    namey <- colnames(P$dat$x[, 1:P$dat$ny])
    names(z$coefficients) <- paste(rep(namey, P$dat$k), "_", rep(namex, rep(P$dat$ny, P$dat$k)), sep = "")
    colnames(z$gt) <- paste(rep(namey, P$dat$nh), "_", rep(nameh, rep(P$dat$ny,P$dat$nh)), sep = "")
    names(z$lambda) <- paste("Lam(",rep(namey,P$dat$nh), "_", rep(nameh, rep(P$dat$ny,P$dat$nh)), ")", sep = "")
    }
  if (P$dat$ny == 1)
    {
    names(z$coefficients) <- namex
    colnames(z$gt) <- nameh
    names(z$lambda) <- nameh
    }

 if(P$gradvf)
    G <- P$gradv(z$coefficients, P$dat)
  else
    G <- P$gradv(z$coefficients, P$dat, g = P$g, z$pt)

  khat <- crossprod(c(z$pt)*z$gt,z$gt)/(P$k2)*P$bwVal
  G <- G/P$k1 

  kg <- solve(khat, G)
  z$vcov_par <- solve(crossprod(G, kg))/n
  p_temp <- solve(khat,G)
  z$vcov_lambda <- solve(khat, ( diag(ncol(khat)) - G %*% (z$vcov_par*n) %*% t(p_temp) ))/n*P$bwVal^2

  z$weights <- P$w
  z$bwVal <- P$bwVal
  names(z$bwVal) <- "Bandwidth"
  dimnames(z$vcov_par) <- list(names(z$coefficients), names(z$coefficients))
  dimnames(z$vcov_lambda) <- list(names(z$lambda), names(z$lambda))
  b <- z$coefficients
  y <- as.matrix(model.response(P$dat$mf, "numeric"))
  ny <- P$dat$ny
  b <- t(matrix(b, nrow = ny))
  x <- as.matrix(model.matrix(P$dat$mt, P$dat$mf, NULL))
  yhat <- x%*%b
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- P$dat$mt
  if(P$model) z$model <- P$dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$call <- P$call
  z$k1 <- P$k1
  z$k2 <- P$k2
  z$khat <- khat
  z$CGEL <- P$CGEL
  class(z) <- paste(P$TypeGel, ".res", sep = "")
  return(z)
  }

momentEstim.baseGel.mod <- function(object, ...)
  {
  P <- object
  x <- P$dat
  n <- ifelse(is.null(dim(x)),length(x),nrow(x))
  l0Env <- new.env()
  assign("l0",rep(0,P$q),envir=l0Env)
  if (!P$constraint)
    {
    if (P$optfct == "optim")
      res <- optim(P$tet0, .thetf, P = P, l0Env = l0Env, ...)
    if (P$optfct == "nlminb")
      res <- nlminb(P$tet0, .thetf, P = P, l0Env = l0Env, ...)
	
    if (P$optfct == "optimize")
      { 
      res <- optimize(.thetf, P$tet0, P = P, l0Env = l0Env, ...)
      res$par <- res$minimum
      res$convergence <- "There is no convergence code for optimize"
      }
    }

  if(P$constraint)
    res <- constrOptim(P$tet0, .thetf, grad = NULL, P = P,l0Env = l0Env, ...)

  All <- .thetf(res$par, P, "all",l0Env = l0Env)
  gt <- All$gt
  rlamb <- All$lambda

  z <- list(coefficients = res$par, lambda = rlamb$lambda, conv_lambda = rlamb$conv, conv_par = res$convergence, dat=P$dat)
  rho1 <- .rho(gt, z$lambda, derive = 1, type = P$typel, k = P$k1/P$k2)
  z$foc_lambda <- crossprod(colMeans(rho1*gt))

  z$type <- P$type
  z$gt <- gt
  rhom <- .rho(z$gt, z$lambda, type = P$typet, k = P$k1/P$k2)
  z$pt <- -.rho(z$gt, z$lambda, type = P$typel, derive = 1, k = P$k1/P$k2)/n

# making sure pt>0
  if (P$type=="CUE")
	{
	eps <- -length(z$pt)*min(min(z$pt),0)
	z$pt <- (z$pt+eps/length(pt))/(1+eps)
	}
##################
  z$conv_moment <- colSums(as.numeric(z$pt)*z$gt)
  z$conv_pt <- sum(as.numeric(z$pt))
  z$objective <- sum(as.numeric(rhom) - .rho(1, 0, type = P$typet, k = P$k1/P$k2))/n

  if(P$gradvf)
    G <- P$gradv(z$coefficients, x)
  else
    G <- P$gradv(z$coefficients, x, g = P$g, z$pt)
 
  khat <- crossprod(c(z$pt)*z$gt, z$gt)/(P$k2)*P$bwVal
  G <- G/P$k1 

  kg <- solve(khat, G)
  z$vcov_par <- solve(crossprod(G, kg))/n
  p_temp <- solve(khat,G)
  z$vcov_lambda <- solve(khat, ( diag(ncol(khat)) - G %*% (z$vcov_par*n) %*% t(p_temp) ))/n*P$bwVal^2
	
  z$weights <- P$w
  z$bwVal <- P$bwVal
  names(z$bwVal) <- "Bandwidth"
 
  if(is.null(names(P$tet0)))
    names(z$coefficients) <- paste("Theta[" ,1:P$k, "]", sep = "")
  else
    names(z$coefficients) <- names(P$tet0)

  colnames(z$gt) <- paste("gt[",1:ncol(z$gt),"]", sep = "")
  names(z$lambda) <- paste("Lambda[",1:ncol(z$gt),"]", sep = "")
  dimnames(z$vcov_par) <- list(names(z$coefficients), names(z$coefficients))
  dimnames(z$vcov_lambda) <- list(names(z$lambda), names(z$lambda))
  if(P$X) z$x <- x
  z$call <- P$call
  z$k1 <- P$k1
  z$k2 <- P$k2
  z$khat <- khat
  z$CGEL <- P$CGEL

  class(z) <- paste(P$TypeGel, ".res", sep = "")
  return(z)
  }



momentEstim.fixedW.formula <- function(object, ...)
  {
  P <- object
  g <- P$g

  dat <- P$x
  x <- dat$x
  k <- dat$k
  k2 <- k*dat$ny
  n <- nrow(x)
  q <- dat$ny*dat$nh
  df <- q-k*dat$ny
  w <- P$weightsMatrix
  if(!all(dim(w) == c(q,q)))
    stop("The matrix of weights must be qxq")
  eigenW <- svd(w)$d
  if(!is.double(eigenW))
    warning("The matrix of weights is not strictly positive definite")
  if(is.double(eigenW))
    {
    if(any(eigenW<=0))
      warning("The matrix of weights is not strictly positive definite")
    }
  
  res2 <- .tetlin(dat, w, P$gradv, g, inv=FALSE)
  z = list(coefficients = res2$par, objective = res2$value, dat=dat, k=k, k2=k2, n=n, q=q, df=df, df.residual = (n-k))	

  z$gt <- g(z$coefficients, dat) 
  b <- z$coefficients
  y <- as.matrix(dat$x[,1:dat$ny])
  b <- t(matrix(b, nrow = dat$ny))
  x <- as.matrix(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
  yhat <- x %*% b
  z$dat <- dat 
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  z$WSpec <- P$WSpec

  names(z$coefficients) <- P$namesCoef
  colnames(z$gt) <- P$namesgt

  z$specMod <- P$specMod
  class(z) <- paste(P$TypeGmm,".res",sep="")
  return(z)	
  }

momentEstim.fixedW <- function(object, ...)
  {
  P <- object
  x <- P$x
  if (P$optfct == "optimize")
    {
    n = nrow(P$g(P$t0[1], x))
    q = ncol(P$g(P$t0[1], x))
    k = 1
    }
  else
    {
    n = nrow(P$g(P$t0, x))
    q = ncol(P$g(P$t0, x))
    k = length(P$t0)
    }
  k2 <- k
  df <- q - k
  w <- P$weightsMatrix
  if(!all(dim(w) == c(q,q)))
    stop("The matrix of weights must be qxq")
  eigenW <- svd(w)$d
  if(!is.double(eigenW))
    warning("The matrix of weights is not strictly positive definite")
  if(is.double(eigenW))
    {
    if(any(eigenW<=0))
      warning("The matrix of weights is not strictly positive definite")
    }

  if (P$optfct == "optim")
    {
    if (P$gradvf)
        {
        gradvOptim <- P$gradv
        gr2 <- function(thet, x,  w, gf, INV)
		{
		gt <- gf(thet, x)
		Gbar <- gradvOptim(thet, x) 
		gbar <- as.vector(colMeans(gt))
		if (INV)		
		  	obj <- crossprod(Gbar, solve(w, gbar))
		else
			obj <- crossprod(Gbar,w)%*%gbar
		return(obj*2)
		}
      argDots <- list(...)
      allArgOptim <- list(par = P$t0, fn = .obj1, gr = gr2, x = P$x, w = w, gf = P$g, INV = FALSE)
      argDots$gr <- NULL
      allArgOptim <- c(allArgOptim,argDots)
      res2 <- do.call(optim,allArgOptim)
      }
    else
      res2 <- optim(P$t0, .obj1, x = P$x, w = w, gf = P$g, INV = FALSE,  ...)
    }

  if (P$optfct == "constrOptim")
	{
	if (!any(c("ui","ci") %in% names(list(...))))
		stop("You must specify ui and ci when optfct is set to constrOptim")
	argDots <- list(...)
	ui <- argDots$ui
	ci <- argDots$ci
	argDots$ui <- NULL
	argDots$ci <- NULL
      	allArgOptim <- list(theta = P$t0, f = .obj1, grad = NULL, ui = ui, ci = ci, x = P$x, w = w, gf = P$g, INV = FALSE)
        allArgOptim <- c(allArgOptim,argDots)
	res2 <- do.call(constrOptim,allArgOptim)
	}
  
  if (P$optfct == "nlminb")
    {
    res2 <- nlminb(P$t0, .obj1, x = P$x, w = w, gf = P$g, INV = FALSE, ...)
    res2$value <- res2$objective
    }
  if (P$optfct == "optimize")
    {
    res2 <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, INV = FALSE, ...)
    res2$par <- res2$minimum
    res2$value <- res2$objective
    }	
  z = list(coefficients = res2$par, objective = res2$value, k=k, k2=k2, n=n, q=q, df=df)	
  if (P$optfct == "optim")
     z$algoInfo <- list(convergence = res2$convergence, counts = res2$counts, message = res2$message)
  else if(P$optfct == "nlminb")
     z$algoInfo <- list(convergence = res2$convergence, counts = res2$evaluations, message = res2$message)

  z$dat <- P$x
  z$gt <- P$g(z$coefficients, P$x)
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  z$WSpec <- P$WSpec

  names(z$coefficients) <- P$namesCoef
  if (is.null(colnames(z$gt)))
	colnames(z$gt) <- paste("gt",1:ncol(z$gt),sep="") 
  z$specMod <- P$specMod
  class(z) <- paste(P$TypeGmm,".res",sep="")	
  return(z)
  }


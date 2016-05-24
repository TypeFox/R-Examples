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

getModel <- function(object, ...)
  {
  UseMethod("getModel")
  }

getModel.constGmm <- function(object, ...)
  {
  class(object) <- "baseGmm"
  obj <- getModel(object)
  if (!is.null(object$t0))
	{
	if (!is.null(dim(object$eqConst)))
		stop("When t0 is provided, eqConst must be a vector")
	if (length(object$eqConst)>=length(object$t0))
		stop("Too many constraints")
        if (is.character(object$eqConst))
		{
		if (is.null(names(object$t0)))
		   stop("t0 must be a named vector if you want eqConst to be names")
		if (any(!(object$eqConst %in% names(object$t0))))
		   stop("Wrong coefficient names in eqConst")
		object$eqConst <- sort(match(object$eqConst,names(object$t0)))
		}
	restTet <- object$t0[object$eqConst]
	obj$t0 <- object$t0[-object$eqConst]
	object$eqConst <- cbind(object$eqConst,restTet)	
	}  else {
	if (is.null(dim(object$eqConst)))
		stop("When t0 is not provided, eqConst must be a 2xq matrix")
	}
  rownames(object$eqConst) <- obj$namesCoef[object$eqConst[,1]]
  if(is(object$g, "formula"))
	{
	if (obj$x$ny>1)
	  stop("Constrained GMM not implemented yet for system of equations")
	obj$g2 <- function(tet, dat)
	      {
	      x <- attr(dat,"eqConst")$Xunc
	      y <- attr(dat,"eqConst")$Yunc
	      ny <- 1
	      nh <- dat$nh
	      tet <- matrix(tet, ncol = ncol(x))
	      e <- y - x %*% t(tet)
	      gt <- e * dat$x[, ny+dat$k+1]
	      if(nh > 1)
			for (i in 2:nh)	  gt <- cbind(gt, e*x[, (ny+dat$k+i)])
	      return(gt)
	      }
    	obj$gradv2 <- function(dat)
	      {
	      x <- attr(dat,"eqConst")$Xunc
	      y <- attr(dat,"eqConst")$Yunc
	      nh <- dat$nh
	      k <- ncol(x)
	      ny <- 1
	      dat$x <- cbind(y,x,dat$x[, dat$ny+dat$k+(1:nh)])
	      dgb <- -(t(dat$x[,(ny+k+1):(ny+k+nh)]) %*% dat$x[,(ny+1):(ny+k)]) %x% diag(rep(1,ny))/nrow(x)
	      return(dgb)
	      }
	attr(obj$x,"eqConst") <- list(unConstg = obj$g2, unConstgradv = obj$gradv2, eqConst = object$eqConst, 
				Yunc = obj$x$x[,1], Xunc = as.matrix(obj$x$x[,1+(1:obj$x$k)]))
	x <- as.matrix(obj$x$x[,1+(object$eqConst[,1])])%*%object$eqConst[,2]
	obj$x$x <- obj$x$x[,-(1+(object$eqConst[,1]))]
	obj$x$x[,1] <- obj$x$x[,1]-x
	obj$x$k <- obj$x$k-nrow(object$eqConst)
	if (obj$x$k<=0)
		stop("Nothing to estimate")
        } else {
	attr(obj$x,"eqConst") <- list(unConstg = obj$g, unConstgradv = obj$gradv, eqConst = object$eqConst)   	   	
	obj$g <- function(tet, dat)
		{
		resTet <- attr(dat,"eqConst")$eqConst
		tet2 <- vector(length=length(tet)+nrow(resTet))
		tet2[resTet[,1]] <- resTet[,2]
		tet2[-resTet[,1]] <- tet
		attr(dat,"eqConst")$unConstg(tet2, dat)
		}
	obj$gradv <- function(tet, dat)
		{
		resTet <- attr(dat,"eqConst")$eqConst
		tet2 <- vector(length=length(tet)+nrow(resTet))
		tet2[resTet[,1]] <- resTet[,2]
		tet2[-resTet[,1]] <- tet
		if (!is.null(as.list(args(attr(dat,"eqConst")$unConstgradv))$g))
			attr(dat,"eqConst")$unConstgradv(tet2, dat, g=attr(dat,"eqConst")$unConstg)[,-resTet[,1]]
		else
			attr(dat,"eqConst")$unConstgradv(tet2, dat)[,-resTet[,1]]
		}
	}
  
  obj$eqConst <- object$eqConst
  obj$namesCoef <- obj$namesCoef[-object$eqConst[,1]]
  obj$type <- paste(obj$type,"(with equality constraints)",sep=" ")	
  mess <- paste(rownames(object$eqConst), " = " , object$eqConst[,2], "\n",collapse="")
  mess <- paste("#### Equality constraints ####\n",mess,"##############################\n\n",sep="")
  obj$specMod <- mess                                   
  return(obj)
  }

getModel.baseGmm <- function(object, ...)
  {
  object$allArg <- c(object, list(...))
  if(is(object$g, "formula"))
    {
    object$gradvf <- FALSE
    if (is.null(object$data))
    	dat <- getDat(object$g, object$x)
    else
    	dat <- getDat(object$g, object$x, object$data)
    
    if(is.null(object$weightsMatrix))
      {
      clname <- paste(class(object), ".", object$type, ".formula", sep = "")
      }
    else
      {
      clname <- "fixedW.formula"
      object$type <- "One step GMM with fixed W"
      }
    object$gform<-object$g
    g <- function(tet, dat)
      {
      x <- dat$x
      ny <- dat$ny
      nh <- dat$nh
      k <- dat$k
      tet <- matrix(tet, ncol = k)
      e <- x[,1:ny] - x[,(ny+1):(ny+k)] %*% t(tet)
      gt <- e * x[, ny+k+1]
      if(nh > 1)
	for (i in 2:nh)	  gt <- cbind(gt, e*x[, (ny+k+i)])
      return(gt)
      }
    gradv <- function(dat)
      {
      x <- dat$x
      ny <- dat$ny
      nh <- dat$nh
      k <- dat$k
      dgb <- -(t(x[,(ny+k+1):(ny+k+nh)]) %*% x[,(ny+1):(ny+k)]) %x% diag(rep(1,ny))/nrow(x)
      return(dgb)
      }

    namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
    nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
 
    if (dat$ny > 1)
      {
      namey <- colnames(dat$x[,1:dat$ny])
      object$namesCoef <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
      object$namesgt <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny, dat$nh)), sep = "")
      } else {
      object$namesCoef <- namex
      object$namesgt <- nameh
      }
    object$g <- g
    object$x <- dat
    attr(object$x,"ModelType") <- "linear"
    }
  else
    {
    attr(object$x,"ModelType") <- "nonlinear"
    k <- length(object$t0)
    if(is.null(names(object$t0)))
       object$namesCoef <- paste("Theta[" ,1:k, "]", sep = "")
    else
       object$namesCoef <- names(object$t0)

    if(is.null(object$weightsMatrix))
      {
      clname <- paste(class(object), "." ,object$type, sep = "")
      }
    else
	{
        clname <- "fixedW"
	object$type <- "One step GMM with fixed W"
	}
    if (!is.function(object$gradv))
      { 
      gradv <- .Gf
      object$gradvf <- FALSE
      }
    else
      {
      gradv <- object$gradv
      object$gradvf <- TRUE
      }
    }
	
  iid <- function(thet, x, g, centeredVcov)
    {
    gt <- g(thet,x)
    if(centeredVcov) gt <- residuals(lm(gt~1))
    n <- ifelse(is.null(nrow(gt)), length(gt), nrow(gt))
    v <- crossprod(gt,gt)/n
    return(v)
    }
 
  
  object$iid<-iid
  object$TypeGmm <- class(object)
  object$gradv <- gradv	
  object$WSpec <- list(vcov = object$vcov, sandwich = list(kernel = object$kernel, bw = object$bw, prewhite = object$prewhite, 
			ar.method = object$ar.method, approx = object$approx, tol = object$tol))
  class(object)  <- clname
  return(object)
  }

getModel.baseGel <- function(object, ...)
  {

  P <- object
  if (P$type == "ETEL")
    {
    P$typel <- "ET"
    P$typet <- "EL"	
    }
  else
    {
    P$typel <- P$type
    P$typet <- P$type
    }
  if(P$optfct == "optim" | P$optfct == "nlminb")
    P$k <- length(P$tet0)
  else
    P$k <- 1
  
  if (is(P$g, "formula"))
    {
    clname <- paste(class(P), ".modFormula", sep = "")
    if (is.null(P$data))
    	dat <- getDat(P$g, P$x)
    else
    	dat <- getDat(P$g, P$x, P$data)
    if (P$k != dat$k)
	stop("The length of tet0 does not match the number of regressors")

    g <- function(tet, dat)
      {
      x <- dat$x
      ny <- dat$ny
      nh  <- dat$nh
      k <- dat$k
      tet <- matrix(tet, ncol = k)
      e <- x[,1:ny] -  x[, (ny+1):(ny+k)]%*%t(tet)
      gt <- e*x[, ny+k+1]
      if (nh > 1)
        {	
        for (i in 2:nh)
          {
          gt <- cbind(gt, e*x[,(ny+k+i)])
          }
        }
      return(gt)
      }
    gradv <- function(tet, dat, pt = NULL)
      {
      x <- dat$x
      ny <- dat$ny
      nh  <- dat$nh
      k <- dat$k
      tet <- matrix(tet, ncol = k)
      if (is.null(pt))
	      dgb <- -(t(x[,(ny+k+1):(ny+k+nh)])%*%x[,(ny+1):(ny+k)])%x%diag(rep(1, ny))/nrow(x)
      else
	      dgb <- -(t(c(pt)*x[,(ny+k+1):(ny+k+nh)])%*%x[,(ny+1):(ny+k)])%x%diag(rep(1, ny))

      return(dgb)
      }
    P$dat <- dat
    P$gform <- P$g
    P$g <- g
    if (P$smooth)
	{
	P$gradv <- .Gf
	P$gradvf <- FALSE
    } else {
	P$gradv <- gradv
	P$gradvf <- TRUE
	}
    }	
  else
    {
    P$dat <- P$x
    clname <- paste(class(P), ".mod", sep = "")
    P$gform <- NULL
    if (!is.function(object$gradv))
      { 
      P$gradv <- .Gf
      P$gradvf <- FALSE
      }
    else
      {
      P$gradvf <- TRUE
      }

    }
  if (P$smooth)
    {
    if (is.function(object$gradv))
	warning("Since you provided gradv, smooth does not apply to it. Do not provide it if you want it to be based on the smoothed g()",
		call. = FALSE)		
    if(P$kernel == "Truncated")
        {
        P$wkernel <- "Bartlett"
        P$k1 <- 2
        P$k2 <- 2
        }
    if(P$kernel == "Bartlett")
        {
        P$wkernel <- "Parzen"
        P$k1 <- 1
        P$k2 <- 2/3
        }
    P$g1 <- P$g

    gt <- gmm(P$g, P$dat, P$tet0, wmatrix = "ident")$gt
    gt <- sweep(gt,2,colMeans(gt),FUN="-")
    class(gt) <- "gmmFct"	
    if (is.function(P$bw))
	    P$bwVal <- P$bw(gt, kernel = P$wkernel, prewhite = P$prewhite, 
	               ar.method = P$ar.method, approx = P$approx)
    else
	    P$bwVal <- P$bw

    P$w <- smoothG(gt, bw = P$bwVal)$kern_weights
    attr(P$dat,"smooth") <- list(bw=P$bwVal, w=P$w, g = P$g1,kernel=P$kernel)

    P$g <- function(thet, x)
      {
      gf <- attr(x,"smooth")$g(thet, x)
      bw <- attr(x,"smooth")$bw
      gt <- smoothG(gf, bw = bw)$smoothx
      return(gt)
      }
    P$q <- ncol(gt)
    }
  else
   {
   P$q <- ncol(P$g(P$tet0, P$dat))
   P$k1 <- 1
   P$k2 <- 1
   P$w <- kernel(1)
   P$bwVal <- 1
   }	
  P$CGEL <- P$alpha
  class(P) <- clname
  return(P)
  }


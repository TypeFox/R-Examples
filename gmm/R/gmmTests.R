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


# This function compute what is needed for the K statistics of KleiBergen (2005)
#####################################################################################

.BigCov <- function(obj,theta0)
	{
	insertRC <- function(A,w,v)
		{
		NewA <- matrix(ncol=ncol(A)+length(w),nrow=nrow(A)+length(w))
		NewA[-w,-w] <- A
		NewA[w,] <- v
		NewA[,w] <- v
		NewA
		}
	dg <- function(obj)
		{
		dat <- obj$dat
		if (!is.null(attr(dat,"eqConst")))
			x <- attr(dat,"eqConst")$Xunc
		else
			x <- dat$x[,(dat$ny+1):(dat$ny+dat$k)]
		k <- ncol(x)
		h <- dat$x[,(dat$ny+dat$k+1):ncol(dat$x)]
		qt <- array(dim=c(dim(obj$gt),k))
		for (i in 1:k)
			qt[,,i] <- -x[,i]*h
		qt}

	if (attr(obj$dat,"ModelType") == "nonlinear")
		{
		Myenv <- new.env()
		assign("obj", obj, envir=Myenv)
		assign("theta", theta0, envir=Myenv)
		gFunct <- if (!is.null(attr(obj$dat,"eqConst")))
				attr(obj$dat,"eqConst")$unConstg
			  else
				obj$g
		assign("g",gFunct,envir=Myenv)
		res <- numericDeriv(quote(g(theta,obj$dat)),"theta",Myenv)
		qT <- attr(res,"gradient")
		} else {
		qT <- dg(obj)}

	qTmat <- apply(qT,3,colSums)
	qT <- matrix(qT,nrow=dim(qT)[1])
	gt <- obj$g(theta0,obj$dat)
	fT <- colSums(gt)
	n <- nrow(gt)
	q <- ncol(gt)
	All <- cbind(gt,qT)
	All <- sweep(All,2,colMeans(All),FUN="-")
	f <- function(x)
		all(abs(x)<1e-7)
	w <- which(apply(All,2,f))
	if (length(w) != 0)
		All <- All[,-w]
	if (dim(All)[2] >= dim(All)[1])
		stop("Too many moment conditions. Cannot estimate V")
	if (obj$WSpec$vcov == "iid") 
		{
		V <- crossprod(All)/nrow(All) 
	} else {
		class(All) <- "gmmFct"
		argSand <- obj$WSpec$sandwich
		argSand$x <- All
		argSand$sandwich <- FALSE
		V <- do.call(kernHAC,argSand)
		}
	if (length(w) != 0)
		V <- insertRC(V,w,0)
	Vff <- V[1:q,1:q]
	Vthetaf <- V[(q+1):nrow(V),1:q]
	list(Vff=Vff,Vthetaf=Vthetaf,qT=qTmat,fT=fT,n=n,q=q)
	}		

KTest <- function(obj, theta0=NULL, alphaK = 0.04, alphaJ = 0.01)
	{
	if (class(obj) != "gmm")
		stop("KTest is only for gmm type objects")

	if (!is.null(attr(obj$dat,"eqConst")))
		{
		if (!is.null(theta0))
			warning("setting a value for theta0 has no effect when the gmm is already constrained")
		resTet <- attr(obj$dat,"eqConst")$eqConst
		tet <- obj$coefficients
		theta0 <- vector(length=length(tet)+nrow(resTet))
		theta0[resTet[,1]] <- resTet[,2]
		theta0[-resTet[,1]] <- tet
		testName <- paste(rownames(resTet), " = ", resTet[,2], collapse="\n")
		if (is.list(obj$dat))		
			{
			x <- model.matrix(obj$dat$mt,obj$dat$mf)
			y <- model.response(obj$dat$mf)
			obj$dat$x <- cbind(y,x,obj$dat$x[,(obj$dat$ny+obj$dat$k+1):ncol(obj$dat$x)])
			obj$dat$k <- ncol(x)
			} else {
			obj$g <- attr(obj$dat,"eqConst")$unConstg
			}	
		dfK <- nrow(resTet)
		which <- resTet[,1]
	} else {
		if (is.null(theta0))
			stop("You must either estimate a restricted model first or set theta0 under H0")		
		if (length(theta0) != length(obj$coef))
			stop("theta0 is only for tests on the whole vector theta when obj is an unrestricted GMM")		
		dfK <- length(theta0)
		testName <- paste(names(obj$coef), " = ", theta0, collapse="\n")
		which <- 1:length(theta0)
		}
	V <- .BigCov(obj, theta0)
	Vff <- V$Vff
	Vtf <- V$Vthetaf
	qT <- V$qT
	fT <- V$fT
	dfJ <- V$q-length(theta0)
	# the following is vec(D)
	D <- c(qT)-Vtf%*%solve(Vff,fT)
	D <- matrix(D,ncol=length(theta0))
	meat <- t(D)%*%solve(Vff,D)
	bread <- t(fT)%*%solve(Vff,D)
	K <- bread%*%solve(meat,t(bread))/V$n
	pv <- 1-pchisq(K,dfK)
	S <- t(fT)%*%solve(Vff,fT)/V$n
	J <- S-K
        dfS <- dfK +  dfJ
	pvJ <- 1-pchisq(J,dfJ)
	type <- c("K statistics","J statistics", "S statistics")
	test <- c(K,J,S)
        pvS <- 1-pchisq(S,dfS)
	test <- cbind(test,c(pv,pvJ,pvS),c(dfK,dfJ,dfS))
	dist <- paste("Chi_sq with ", c(dfK,dfJ,dfS), " degrees of freedom", sep="")
	if(dfJ>0)
		ans <- list(test=test,dist=dist,type=type,testName=testName)
	else
		ans <- list(test=matrix(test[1,],nrow=1),dist=dist[1],type=type[1],testName=testName)
	if (pvJ<alphaJ)	{
		message <- "reject"
	} else {
		if (pv < alphaK)
			message <- "reject"
		else
			message <- "do not reject"
	}
	ans$KJ <- (message == "do not reject")
	ans$Matrix <- list(D=D,bread=bread,meat=meat,qT=qT,fT=fT)	
	ans$message <- paste("KJ-test result: We ", message, " H0 (alphaJ = ", alphaJ, ", alphaK = ", alphaK, ")", sep="")	
	class(ans) <- "gmmTests"
	ans
	}	

print.gmmTests <- function(x, digits = 5, ...)
	{
	lab <- paste("%0.",digits,"f",sep="")
	cat("\nTest robust to weak identification\n")
	cat("**********************************\n\n")
	cat("The Null Hypothesis\n")
	cat(x$testName,"\n\n")
	for (i in 1:length(x$type))
		{
		cat(x$type[i],"\n")
		cat("Test: ", sprintf(lab,x$test[i,1]), "(",x$dist[i],")\n")
		cat("P-value: ", sprintf(lab,x$test[i,2]),"\n\n")
		}
	cat(x$message,"\n")
	}

gmmWithConst <- function(obj, which, value)
	{
	argCall <- obj$allArg
	argCall$call = NULL
	if (!is.null(attr(obj$w0,"Spec")))
		if (is.function(argCall$bw))
			argCall$bw <- attr(obj$w0,"Spec")$bw

	if (length(which)>=length(obj$coefficients))
		stop("Too many constraints")
        if (is.character(which))
		{
		if (any(!(which %in% names(obj$coefficients))))
		   stop("Wrong coefficient names in eqConst")
		if (attr(obj$dat,"ModelType") == "linear")
			which <- match(which,names(obj$coefficients))
		}  
	if (!is.null(argCall$t0))
		{
		argCall$t0 <- obj$coefficients
		argCall$t0[which] <- value
		}
		
	if (attr(obj$dat,"ModelType") == "nonlinear")
		{
		eqConst <- which
	} else {
		eqConst <- cbind(which,value)		
		}
	argCall$eqConst <- eqConst
	res <- do.call(gmm,argCall)
	res$call <- match.call()
	return(res)
	}



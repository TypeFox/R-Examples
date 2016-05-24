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


summary.gmm <- function(object, ...)
	{
	z <- object
	se <- sqrt(diag(z$vcov))
	par <- z$coefficients
	tval <- par/se
	ans <- list(met=z$met,kernel=z$kernel,algo=z$algo,call=z$call)
	names(ans$met) <- "GMM method"
	names(ans$kernel) <- "kernel for cov matrix"
		
	ans$coefficients <- cbind(par,se, tval, 2 * pnorm(abs(tval), lower.tail = FALSE))
    	dimnames(ans$coefficients) <- list(names(z$coefficients), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
	ans$stest <- specTest(z)
        ans$algoInfo <- z$algoInfo
	if(z$met=="cue")
		ans$cue <- object$cue
	if (!is.null(object$initTheta))
		{
		ans$initTheta <- object$initTheta
		names(ans$initTheta) <- names(z$coefficients)
		}
        ans$specMod <- object$specMod
	ans$bw <- attr(object$w0,"Spec")$bw
	ans$weights <- attr(object$w0,"Spec")$weights
	if(object$infVcov == "iid")
		ans$kernel <- NULL
	class(ans) <- "summary.gmm"
	ans
	}

summary.tsls <- function(object, vcov = NULL, ...)
	{
	if (!is.null(vcov))
		object$vcov=vcov
	else
		object$vcov=vcov(object)
	ans <- summary.gmm(object)
	ans$met <- paste(ans$met, "(Meat type = ", attr(object$vcov, "vcovType"), ")",sep="")
	k <- object$dat$k
	if (!is.null(object$fsRes))
		{
		fstat <- vector()
                fsRes <- object$fsRes
                if (class(fsRes) == "listof")
                    {
                    nendo <- length(fsRes)
                } else {
                    nendo <- 1
                }
                if (nendo == 1)
                    { 
                    fstat[1] <- fsRes$fstatistic[1]
                    df1 <- fsRes$fstatistic[2]
                    df2 <- fsRes$fstatistic[3]
                } else {
                    fstat[1] <- fsRes[[1]]$fstatistic[1]
                    df1 <- fsRes[[1]]$fstatistic[2]
                    df2 <- fsRes[[1]]$fstatistic[3]
                }
                if (nendo > 1){
                    for (i in 2:nendo)	
			fstat[i] <- fsRes[[i]]$fstatistic[1]
                }
		pvfstat <- 1-pf(fstat,df1, df2)
		names(fstat) <- attr(fsRes,"Endo")
		ans$fstatistic <- list(fstat = fstat, pvfstat = pvfstat, df1 = df1, df2 = df2)
		}
        ans$specMod <- object$specMod
	class(ans) <- "summary.tsls"
	return(ans)
	}


print.summary.tsls <- function(x, digits = 5, ...)
	{
	print.summary.gmm(x,digits)
	if (!is.null(x$fstatistic))
		{
		cat("\n First stage F-statistics: \n")
		if(names(x$fstatistic$fstat)[1]=="(Intercept)")
			start=2
		else
			start=1
		for (i in start:length(x$fstatistic$fstat))
			cat(names(x$fstatistic$fstat)[i], 
			": F(",x$fstatistic$df1,", ",x$fstatistic$df2,") = ",x$fstatistic$fstat[i], 
			" (P-Vavue = ",x$fstatistic$pvfstat[i],")\n")
	} else {
		cat("\n No first stage F-statistics (just identified model)\n")
		}
	}

print.summary.gmm <- function(x, digits = 5, ...)
	{
	cat("\nCall:\n")
	cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
	cat("\nMethod: ", x$met,"\n")
	if (x$met=="cue")
		cat("         (",x$cue$message,")\n\n")
	else
		cat("\n")
	if( !is.null(x$kernel))
		{
		cat("Kernel: ", x$kernel)
		if (!is.null(x$bw))
			cat("(with bw = ", round(x$bw,5),")\n\n")
		else
			cat("\n\n")	
		}
	cat("Coefficients:\n")
	print.default(format(x$coefficients, digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	cat(x$stest$ntest,"\n")
	print.default(format(x$stest$test, digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	if(!is.null(x$initTheta))
		{
		cat("Initial values of the coefficients\n")
		print(x$initTheta)
		cat("\n")
		}
        cat(x$specMod)
	if(!is.null(x$algoInfo))
		{	
		cat("#############\n")
 		cat("Information related to the numerical optimization\n")
		}
	if(!is.null(x$algoInfo$convergence))
		cat("Convergence code = ", x$algoInfo$convergence,"\n")
	if(!is.null(x$algoInfo$counts))
		{	
		cat("Function eval. = ",x$algoInfo$counts[1],"\n")
		cat("Gradian eval. = ",x$algoInfo$counts[2],"\n")
		}	
	if(!is.null(x$algoInfo$message))
		cat("Message: ",x$algoInfo$message,"\n")
	invisible(x)
	}


formula.gmm <- function(x, ...)
{
    if(is.null(x$terms))
	stop("The gmm object was not created by a formula")
    else
	formula(x$terms)
}

confint.gmm <- function(object, parm, level=0.95, ...)
		{
		z <- object
		se <- sqrt(diag(z$vcov))
		par <- z$coefficients
			
		zs <- qnorm((1-level)/2,lower.tail=FALSE)
		ch <- zs*se
		ans <- cbind(par-ch,par+ch)
		dimnames(ans) <- list(names(par),c((1-level)/2,0.5+level/2))
		if(!missing(parm))
			ans <- ans[parm,]
		ans
		}
		
residuals.gmm <- function(object,...) 
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$residuals
	}

fitted.gmm <- function(object,...)
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$fitted.value
	}

print.gmm <- function(x, digits=5, ...)
	{
	cat("Method\n", x$met,"\n\n")
	cat("Objective function value: ",x$objective,"\n\n")
	print.default(format(coef(x), digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	if(!is.null(x$algoInfo$convergence))
		cat("Convergence code = ", x$algoInfo$convergence,"\n")
	cat(x$specMod)
	invisible(x)
	}

coef.gmm <- function(object,...) object$coefficients

vcov.gmm <- function(object,...) object$vcov

estfun.gmmFct <- function(x, y = NULL, theta = NULL, ...)
	{
	if (is(x, "function"))
		{
		gmat <- x(theta, y)
		return(gmat)
		}
	else
		return(x)
	}
estfun.tsls <- function(x, ...)
	{
	model.matrix(x)*c(residuals(x))
	}
model.matrix.tsls <- function(object, ...)
{
dat <- object$dat
ny <- dat$ny
nh <- dat$nh
k <- dat$k
x <- dat$x
n <- nrow(x)
hm <- as.matrix(x[,(ny+k+1):(ny+k+nh)])
xm <- as.matrix(x[,(ny+1):(ny+k)])
xhat <- lm(xm~hm-1)$fitted
assign <- 1:ncol(xhat)
if (attr(object$terms,"intercept")==1)
	assign <- assign-1
attr(xhat,"assign") <- assign
xhat
}
vcov.tsls <- function(object, type=c("Classical","HC0","HC1","HAC"), hacProp = list(), ...)
	{
	type <- match.arg(type)
	if (type == "Classical")
		{
		sig  <- sum(c(residuals(object))^2)/(nrow(object$dat$x)-object$dat$k)
  		ny <- object$dat$ny
		nh <- object$dat$nh
		k <- object$dat$k
		n <- nrow(object$dat$x)
		hm <- as.matrix(object$dat$x[,(ny+k+1):(ny+k+nh)])
		Omega <- crossprod(hm)*sig/nrow(object$dat$x)
		vcovType <- "Classical"
		V <- solve(crossprod(object$G,solve(Omega,object$G)))/nrow(object$dat$x)
		}
	else if (strtrim(type,2) == "HC")
		{
		meat <- meatHC(object, type)
		bread <- bread(object)
		vcovType <- paste("HCCM: ", type, sep="")
		V <- crossprod(bread, meat%*%bread)/nrow(object$dat$x)
		}
	else
		{
		object$centeredVcov <- TRUE
		gt <- model.matrix(object)*c(residuals(object))
		gt <- lm(gt~1)
		arg <- c(list(x=gt,sandwich=FALSE),hacProp)
		meat <- do.call(kernHAC, arg)
		KType <- ifelse(is.null(hacProp$kernel),  formals(kernHAC)$kernel[[2]], hacProp$kernel)
		vcovType <- paste("HAC: ", KType, sep="")
		bread <- bread(object)
		V <- crossprod(bread, meat%*%bread)/nrow(object$dat$x)
		}
	attr(V, "vcovType") <- vcovType
	return(V)
	}

estfun.gmm <- function(x, ...)
  {
  foc <- x$gt %*% x$w %*% x$G
  return(foc)
  }

bread.gmm <- function(x, ...)
  {
  GWG <- crossprod(x$G, x$w %*% x$G)
  b <- try(solve(GWG), silent = TRUE)
  if (class(b) == "try-error")
    stop("The bread matrix is singular")
  return(b)
  }
bread.tsls <- function(x, ...)
	{
	dat <- x$dat
  	ny <- dat$ny
	nh <- dat$nh
	k <- dat$k
	x <- dat$x
	n <- nrow(x)
	hm <- as.matrix(x[,(ny+k+1):(ny+k+nh)])
	xm <- as.matrix(x[,(ny+1):(ny+k)])
	xhat <- lm(xm~hm-1)$fitted
	solve(crossprod(xhat)/n)
	}




		



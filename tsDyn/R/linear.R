#' @export

##Copyright (C) 2005/2006  Antonio, Fabio Di Narzo
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

#linear model fitter (via OLS)
#str: call to nlar.struct
linear <- function(x, m, d=1, steps=d, series,include = c("const", "trend","none", "both"), type=c("level", "diff", "ADF")) {
	str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
	type<-match.arg(type)
	
	###Build regressor matrix
	if(type=="level"){
	  xx <- getXX(str)
	  yy <- getYY(str)
	}
	else{ 
	  if(type=="diff"){
	    xx <- getdXX(str)
	    yy <- getdYY(str)
	  }
	  else if(type=="ADF"){
	    xx <- cbind(getdX1(str),getdXX(str))
	    yy <- getdYY(str)
	  }
	  str$xx<-xx
	  str$yy<-yy
	}
	
	constMatrix<-buildConstants(include=include, n=nrow(xx)) #stored in miscSETAR.R
	incNames<-constMatrix$incNames #vector of names
	const<-constMatrix$const #matrix of none, const, trend, both
	ninc<-constMatrix$ninc #number of terms (0,1, or 2)
	xx <- cbind(const,xx)
	###name the regressor matrix
	phi<-ifelse(type=="level", "phi", "Dphi")
	dX1<- if(type=="ADF") "phi.1" else NULL
	nlags<-if(type=="ADF") ncol(xx)-ninc-1 else ncol(xx)-ninc
	colnames(xx) <- c(incNames, dX1, paste(phi,1:nlags, sep="."))
	
	#evaluate the model
	res <- lm.fit(xx, yy)
	res$incNames<-incNames
	#check if unit root lie outside unit circle
	if(type=="level")
	  is<-isRoot(coef(res), regime=".", lags=seq_len(m))
	
	return(extend(nlar(str,
	  coefficients=res$coefficients,
	  fitted.values=res$fitted.values,
	  residuals=res$residuals,
	  k=res$rank,
	  model=data.frame(yy,xx),
	  model.specific=res),
		"linear"))
}

#' @S3method print linear
print.linear <- function(x, ...) {
	NextMethod(...)
	cat("\nAR model\n")
	cat("Coefficients:\n")
	print(x$coef, ...)
	invisible(x)
}

#' @S3method summary linear
summary.linear <- function(object, ...) {
	ans <- list()
	obj <- c(object, object$model.specific)
	Qr <- obj$qr
	n <- nrow(Qr$qr)
	p <- obj$rank
	resvar <- mse(object)*n/(n-p)
	p1 <- 1:p
	R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
	se <- sqrt(diag(R) * resvar)
	est <- obj$coefficients[Qr$pivot[p1]]
	tval <- est/se
	coef <- cbind(est, se, tval, 2*pt(abs(tval), n-p, lower.tail = FALSE))
	dimnames(coef) <- list(names(est), c(" Estimate"," Std. Error"," t value","Pr(>|t|)"))
	ans$coef <- coef
	return( extend(summary.nlar(object, ...), "summary.linear", listV=ans) )
}

#' @S3method print summary.linear
print.summary.linear <- function(x, digits=max(3, getOption("digits") - 2),
	signif.stars = getOption("show.signif.stars"), ...) {
	NextMethod(...)
	cat("\nCoefficient(s):\n")
	printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)
	invisible(x)
}

oneStep.linear <- function(object, newdata, ...) {
	if((names(coef(object))[1]!="const" |names(coef(object))[2]=="trend"))
	  stop("Currently, only arg const is implemented for predict method")
	cbind(1,newdata) %*% object$coefficients
}


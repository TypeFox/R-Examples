## Copyright (C) 2006  Antonio, Fabio Di Narzo
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
#' @export
#AAR fitter
aar <- function(x, m, d=1, steps=d, series){
	str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
	xx <- str$xx
	yy <- str$yy
  dat <- data.frame(cbind(xx,y=yy))
  predNames <- paste("s(",names(dat)[1:NCOL(xx)],", bs=\"cr\")")
  predNames <- paste(predNames, collapse="+")
  form <- as.formula(paste("y ~", predNames))
  model <- gam(form, data=dat)
	return( extend(nlar(str,
		coefficients = coef(model),
		fitted.values = model$fitted.values,
		residuals = model$residuals,
		k = model$rank,
		model= model$model,
		model.specific=model
		), "aar") )
}

#' @S3method print aar
print.aar <- function(x, ...) {
	NextMethod(...)
	cat("\nAAR model\n")
  print.gam(x$model.specific, ...)
  invisible(x)
}

#' @S3method summary aar
summary.aar <- function(object, ...) {
	extend(summary.nlar(object), "summary.aar", internals=summary(object$model.specific, ...))
}

#' @S3method print summary.aar
print.summary.aar <- function(x, digits=max(3, getOption("digits") - 2),
	signif.stars = getOption("show.signif.stars"), ...) {
	NextMethod(digits=digits, signif.stars=signif.stars, ...)
  print.summary.gam(x$internals, digits=digits, signif.stars=signif.stars, ...)
  invisible(x)
}

#' @S3method plot aar
plot.aar <- function(x, ask=interactive(), ...) {
	op <- par(no.readonly=TRUE)
	par(ask=ask)
	NextMethod(ask=ask, ...)
  m <- x$str$m
  x <- x$model.specific
  pd <- list()
  for(i in 1:m) {
    raw<-x$model[x$smooth[[i]]$term]
    xx<-seq(min(raw),max(raw),length=100)   # generate x sequence for prediction
    dat<-data.frame(x=xx)
    names(dat)<-x$smooth[[i]]$term
    X <- PredictMat(x$smooth[[i]],dat)   # prediction matrix from this term
    first<-x$smooth[[i]]$first.para;
    last<-x$smooth[[i]]$last.para
    p<-x$coefficients[first:last]      # relevent coefficients 
    fit<-X%*%p                         # fitted values
    edf<-sum(x$edf[first:last])
    xterm <- x$smooth[[i]]$term
    xlabel <- xterm
    ylabel<-paste("s(",xterm,",",as.character(round(edf,2)),")",sep="")
    pd.item<-list(fit=fit,dim=1,x=xx,ylab=ylabel,xlab=xlabel,raw=raw[[1]])
    pd[[i]]<-pd.item
  }

  for (i in 1:m)
    plot(pd[[i]]$x,pd[[i]]$fit,type="l",xlab=pd[[i]]$xlab,ylab=pd[[i]]$ylab,...)

  par(op)
  invisible(x)
}

oneStep.aar <- function(object, newdata, ...){
  newdata <- data.frame(newdata)
  names(newdata) <- attr(object$model.specific$terms, "term.labels")
  predict.gam(object$model.specific, data.frame(newdata))
}


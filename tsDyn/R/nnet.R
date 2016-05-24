#' @export

## Copyright (C) 2005/2006  Antonio, Fabio Di Narzo
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

#Neural Network fitter
nnetTs <- function(x, m, d=1, steps=d, series, size, control=list(trace=FALSE)) {
	str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
	args <- c(list(getXX(str), getYY(str), size=size, linout=TRUE), control)
	res <- do.call(nnet::nnet, args)
	res$k <- length(res$wts)
	res$fitted <- res$fitted.values	
	if(length(unique(res$fitted))/length(res$fitted)<0.1) warning("Fitted model seems to lead poor result. Please check result. ")
	return(extend(nlar(str,
		coefficients = res$wts,
		residuals=res$residuals,
		fitted.values=res$fitted,
		k=res$k,
		model=NULL,
		model.specific=res), "nnetTs"))
}

#' @S3method print nnetTs
print.nnetTs <- function(x, ...) {
	NextMethod(...)
	cat("\nNNET time series model\n")
	print(x$model.specific, ...)
	invisible(x)
}

oneStep.nnetTs <- function(object, newdata, ...)
	predict(object$model.specific, newdata)

#' @export
selectNNET <- function(x, m, d=1, steps=d, size=1:(m+1), maxit=1e3, trace=FALSE) {
	IDS <- as.matrix( size )
	colnames(IDS) <- c("size")

	lis_models <- list()
	for(i in 1:length(size)){
	  lis_models[[i]] <- nnetTs(x=x, m=m, d=d, steps=steps, size=i, control=list(maxit=maxit, trace=trace))
	}
	

	computedAIC <- sapply(lis_models, AIC)
	computedBIC <- sapply(lis_models, BIC)

	res <- cbind(IDS, AIC = computedAIC, BIC = computedBIC)
	idSel <- sort(computedAIC, index=TRUE)$ix
	idSel <- idSel[1:min(10, length(idSel))]
	res <- data.frame(res[idSel,], row.names=NULL)
	return(res)
}



plot.clogitL1 = function (x, logX = T, add.legend=F, add.labels=T, lty=1:ncol(x$beta), col=1:ncol(x$beta), ...) 
{
    if (logX) 
        horiz = log(x$lambda)
    else horiz = x$lambda
    matplot(x = horiz, y = x$beta, type = "l", xlab = "Regularisation parameter", 
        ylab = "Parameter estimate", lty=lty, col=col, ...)

	# add the legend, if required	
	if(add.legend){
		# first find the variable names
		if(is.null(dimnames(x$X_c))){
			var.names = paste("Variable", 1:ncol(x$beta))
		}else{
			var.names = dimnames(x$X_c)[[2]]
		}

		# now add the legend
		legend ("topright", legend=var.names, lty=lty, col=col, ...)
	}

	# label the curves
	if (add.labels){
		# first find the variable names
		if(is.null(dimnames(x$X_c))){
			plot.names = 1:ncol(x$beta)
		}else{
			plot.names = dimnames(x$X_c)[[2]]
		}

		text(x=min(horiz),y=x$beta[nrow(x$beta),], labels=plot.names, lty=lty, col=col, ...)
	}
}

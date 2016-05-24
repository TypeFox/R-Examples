
print.circsizer<-function(x, digits=NULL, ...) {
	cat("\nCall:\n\t", deparse(x$call), "\n\nngrid: ", x$ngrid, 
	"\tlog.scale: ", x$log.scale, "\n\nBandwidths = ", "\n\t", sep = "")
	print(x$bw, digits=digits, ...)
}


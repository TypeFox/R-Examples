print.clogitL1 = function (x, digits = 6, ...){
	# prints a summary of the regularisation path of a clogitL1 object
	#	output:
	#		- 3 column matrix containing: non-zero beta, %deviance explained, lambda

	out = data.frame(Df=x$nz_beta, DevPerc=round(x$dev_perc, digits), Lambda=round(x$lambda, digits))
	print(out, ...)
}
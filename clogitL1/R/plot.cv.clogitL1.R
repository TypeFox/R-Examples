plot.cv.clogitL1 = function (x, ...) 
{
	matplot(x=x$lambda, y=cbind(x$mean_cv + x$se_cv, x$mean_cv - x$se_cv), type="p", pch=c("-","-"), col=c("gray", "gray"), xlab="log(lambda)", ylab="Conditional likelihood deviance", ...)
	matlines (x=rbind(x$lambda, x$lambda), y=rbind(x$mean_cv + x$se_cv, x$mean_cv - x$se_cv), col="gray", lty=2, ...)
	points(x=x$lambda, y=x$mean, pch="o", col="red", ...)
    	axis(side = 3, at = x$lambda, labels = x$nz_beta, tick = FALSE, ...)
    	abline(v = x$minCV_lambda, lty = 2, ...)
   	abline(v = x$minCV1se_lambda, lty = 2, ...)
}

residualsStudy <- function(model)
{
	residuals <- model$model$residuals		

	if (model$type=="PolyMARS"){
		Yfit	<- model$model$fitted
	} else Yfit	<- model$model$fitted.values

	op <- par(mfrow = c(1, 3),pty="s")
	plot(residuals ,ylab = "residuals")

	plot(Yfit,residuals,xlab = "fitted values",ylab = "residuals")

	hist(residuals, freq=FALSE,xlab = "residuals", ylab="density",main= "",
		ylim=c(0,max(density(residuals)$y, hist(residuals,plot=FALSE)$density)))
  	lines(density(residuals) ,col="red")
  
  	# qqnorm(residuals,xlab = "theorical quantiles",ylab = "sample quantiles",main="")
  	# qqline(residuals)
	par(op)
	#mtext("Residuals study", side=3, line=0, font=1, cex=1.3)
        title("Residuals study", cex.main=1.3)
}

plot.CADFtest <- function(x, plots=(1:4), ...)
{
#	x    : an x of class `CADFtest'
#	plots:  specify the plots to be produced
	switch(length(plots),
			par(mfrow=c(1,1)),
			par(mfrow=c(2,1)),
			layout(matrix(c(1,1,2,3),2,2, byrow=TRUE)),
			par(mfrow=c(2,2))
	)
	r <- residuals(x)

    if (1 %in% plots)
	{
		sr <- rstandard(x$est.model)
		plot(sr, type="h", main="standardized residuals", 
			ylab="", xlab="Time", ...)
        abline(h=0)
	}
	if (2 %in% plots)
	{
		jb <- jarque.bera.test(r)
		plot(density(r), main="residuals density", 
			xlab=paste("p-value of the Jarque-Bera test = ", round(jb$p.value,4), sep=""), ...)
	}
	if (3 %in% plots)
	{
		acf(r, main="residuals ACF", ...)
	}
    if (4 %in% plots)
	{
		pacf(r, main="residuals PACF", ...)
	}
}

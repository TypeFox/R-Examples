cdf.plot <-
function(z)
{
	plot(z$x.info[,'z'],z$x.info[,'Emp'],xlab='Sorted observations',las=2,cex.axis=.8,ylab='Cumulative distribution function',
		type='n')
	lines(z$x.info[,'z'],z$x.info[,'Fz'],col='red')
	points(z$x.info[,'z'],z$x.info[,'Emp'],cex=.5,pch=19)
}


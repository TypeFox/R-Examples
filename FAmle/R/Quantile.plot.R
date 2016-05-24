Quantile.plot <-
function(z,ci=FALSE,alpha=.05)
{
	plot(z$x.info[,'zF'],z$x.info[,'z'],xlab='Theoretical quantile',ylab='Observed quantile',las=2,
		cex.axis=.8,type='n')
	abline(0,1,col='red')
	points(z$x.info[,'zF'],z$x.info[,'z'],pch=19,cex=.5)
	if(ci)
	{
		Q.int <- delta.QQ(z,alpha)
		lines(z$x.info[,'zF'],Q.int[1,],col='red',lty='dashed')
		lines(z$x.info[,'zF'],Q.int[2,],col='red',lty='dashed')
		legend('topleft',paste('Approx. ',100*(1-alpha),'%-CI',sep=''),col='red',lty='dashed',bty='n')
	}
}


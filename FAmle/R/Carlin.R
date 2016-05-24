Carlin <-
function(x)
{
	sims <- x$sims
	names. <- colnames(sims)
	d <- ncol(sims)
	lag.1 <- sapply(as.list(1:d),function(g) acf(sims[,g],plot=FALSE)$acf[2])
	p <- c(.025,.5,.975)
	Qp <- sapply(as.list(1:d),function(g) quantile(sims[,g],p))
	layout(t(matrix(1:(2*d),ncol=d)))
	for(i in 1:d)
	{
		plot.ts(sims[,i],xlab='Iteration',ylab=names.[i],main='',
		cex.axis=1)
		name.i <- paste(names.[i],': lag 1 acf = ',round(lag.1[i],3))
		title(main=name.i)
		hist(sims[,i],freq=FALSE,col='gray',border='white',xlab='',ylab='',
		cex.axis=1,main='')
		name.i.2 <- paste(paste(p,collapse=', '),' quantiles are ',
			paste(round(Qp[,i],3),collapse=', '),sep='')
		title(main=name.i.2,cex.main=.85)
	}
	layout(matrix(1,1,1))
}
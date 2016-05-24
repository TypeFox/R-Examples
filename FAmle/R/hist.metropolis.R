hist.metropolis <-
function(x,density=TRUE,...)
{
	k <- x$input$k
	layout(matrix(1:k,ncol=k))
	for(i in 1:k)
	{
		hist(x$sims[,i],main=colnames(x$sims)[i],freq=FALSE,xlab='',...)
		if(density) lines(density(x$sims[,i]),col='red')
	}
	layout(matrix(1))
}


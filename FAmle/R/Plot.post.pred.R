Plot.post.pred <-
function(x,...)
{
	x.new <- distr(x=x$iter,dist=x$input$dist,param=x$sims,type='r')
	hist(x$input$x.info[,'x'],freq=FALSE,main='Posterior predictive distribution',xlab='',...)
	lines(density(x.new),col='red')
}


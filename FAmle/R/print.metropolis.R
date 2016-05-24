print.metropolis <-
function(x,stats.fun=NULL,...)
{
	if(is.null(stats.fun)) stats.fun <- function(x) c(mean=mean(x),sd=sd(x),quantile(x,c(.01,.5,.99)))


	cat('-----------------------------------------\n')
	cat('       Bayesian Posterior Summary        \n')
	cat('-----------------------------------------\n')
	cat('Data object: ',x$input$data.name,'\n')
	cat('Likelihood: ',x$input$dist,'\n')
	if(x$prior=='no') cat('Prior distribution specified: No\n')
	else cat('Prior distribution specified: Yes\n')
	cat('-----------------------------------------\n')
	cat('Marginal posterior distributions\n')
	cat('---\n')
	print(apply(x$sims,2,stats.fun),digits=3)
	cat('-----------------------------------------\n')
	cat('Posterior correlations\n')
	cat('---\n')
	print(cor(x$sims),digits=2)
	cat('----------------------------------------\n')
	cat('Number of iterations:',x$iter,'\n')
	cat('Burnin: first',x$burn,'iterations\n')
	cat('Acceptance rate:',x$rate,'\n')
	cat('Required time:',round(x$total.time,3),'minutes\n')
	cat('-----------------------------------------\n')
}



vorob_optim_parallel2 <- function(x,other.points, integration.points,integration.weights=NULL,
		intpoints.oldmean,intpoints.oldsd,precalc.data,
		model, T, new.noise.var=NULL,batchsize,alpha,current.vorob){
	
	x.complete <- c(x,other.points)
	return(vorob_optim_parallel(
		x = x.complete, integration.points = integration.points, integration.weights = integration.weights,
		intpoints.oldmean = intpoints.oldmean,intpoints.oldsd = intpoints.oldsd,precalc.data = precalc.data,
		model = model,T = T,new.noise.var = new.noise.var,batchsize=batchsize,alpha=alpha,current.vorob=current.vorob
		)
	)
	
}

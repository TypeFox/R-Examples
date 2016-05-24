## -----------------------------------------------------------------------------
## Fonction generateWithMH
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

generateWithlrmM = function(seeds,seeds_eval=NA,N,lambda=1,limit_f,burnin=30,thinning=4,p=Normal,q=Uniform,modified=TRUE,VA_function=NULL,VA_esp=NULL,VA_var=NULL) {

	#q is the proposed PDF and y=q(x) should generate a random y|x while q(x,y) should return P(y|x)
	#p is the target PDF
	#VA_function is an optional function to add if samples are to be used for a MC estimator as sum(VA_function(x)) to get correlation between samples
	#VA_esp & VA_var are the esperance and the variance of the random variable VA_function(x)

	tic = proc.time()[3]

	#set variables
	seeds = as.matrix(seeds)
	dim = dim(seeds)
	n_seeds = dim[2];
	dimension = dim[1]
	chain_length = ceiling(N/n_seeds)
	U = matrix(NA,dimension,chain_length*n_seeds)
	G = 1:chain_length*n_seeds
	Ncall = 0;

	if(missing(limit_f)) {
		limit_fun = function(x) {-1}
	}
	else {
		limit_fun = function(x) {tryCatch(limit_f(x)$mean, error = function(cond) {return(limit_f(x))})}
	}

	#core loop
	for (j in 1:n_seeds){
		start = (j-1)*chain_length+1
		end = j*chain_length
		seed = seeds[,j]
		eval_seed = seeds_eval[j]
		MH = MetropolisHastings(x0=seed,eval_x0=eval_seed,chain_length=(chain_length-1),modified=modified,modif_parameter=lambda,limit_fun=limit_fun,burnin=burnin,thinning=thinning,p=p,q=q)
		U[,start:end] = MH$points
		G[start:end] = MH$eval
		Ncall = Ncall + MH$Ncall
	}
	
	if(dim(U)[2]>N) {U = U[,1:N]; G = G[1:N]}
	colnames(U) <- NULL
	toc = proc.time()[3]-tic
	cat("",(burnin + (thinning+1)*chain_length)*n_seeds,"points generated in",toc,"sec. with ",n_seeds,"seeds,",N,"points kept : burnin =",burnin,"thinning =",thinning,"\n")

	if(!is.null(VA_function)) {
		cat("====== Beginning of Monte-Carlo estimation ======\n")
		cat("#Calculate Monte-Carlo estimator\n")
		VA_values = apply(U,2,VA_function)
		MC_est = mean(VA_values)
		cat(" MC_est =",MC_est,"\n")

		cat("#Calculate the covariance between samples \n")
		stat = MCMCcovariance(n_seeds=n_seeds,chain_length=chain_length,VA_values=VA_values,VA_esp=MC_est)
		MC_gamma = stat$gamma
		MC_var = stat$var
		MC_delta = stat$cov

		res = list(points=U,
               eval=G,
               chain_length=chain_length,
               Ncall=Ncall,
               VA_values=VA_values,
               est=MC_est,
               var=MC_var,
               delta=MC_delta,
               gamma=MC_gamma)
	}
	else{	res = list(points=U,
                   eval=G,
                   chain_length=chain_length,
                   Ncall=Ncall)
	}

	return(res)
}
plot_all_marginals <-
function(MCMC.output,percent.burnin=0,thinning=1,population.names=NULL){
	MCMC.output.list <- load_MCMC_output(MCMC.output)
	with(MCMC.output.list, {
            for (k in 1:nrow(aE)) {
                plot_marginal(aE[k,]/aD,percent.burnin,thinning,param.name=paste(sprintf("aE_%s",k),"aD",sep="/")); 
                devAskNewPage(ask=TRUE) 
            }
			plot_marginal(a0,percent.burnin,thinning)
			plot_marginal(aD,percent.burnin,thinning)
            for (k in 1:nrow(aE)) { plot_marginal(aE[k,],percent.burnin,thinning,param.name=paste("aE_",k)) }
			plot_marginal(a2,percent.burnin,thinning)
			plot_marginal(beta,percent.burnin,thinning)
				if(exists("phi_mat")){
					plot_all_phi_marginals(phi_mat,percent.burnin,thinning,population.names)
				}
			plot_marginal(LnL_thetas,percent.burnin,thinning)
			plot_marginal(LnL_counts,percent.burnin,thinning)
			plot_marginal(Prob,percent.burnin,thinning)
		})
	}

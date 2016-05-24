plot_all_trace <-
function(MCMC.output,percent.burnin=0,thinning=1,population.names=NULL){
	MCMC.output.list <- load_MCMC_output(MCMC.output)
	with(MCMC.output.list, {
            for (k in 1:nrow(aE)) {
                plot_trace(aE[k,]/aD,percent.burnin,thinning,param.name=paste(sprintf("aE_%s",k),"aD",sep="/")); 
                devAskNewPage(ask=TRUE) 
            }
			plot_trace(a0,percent.burnin,thinning)
			plot_trace(aD,percent.burnin,thinning)
            for (k in 1:nrow(aE)) { plot_trace(aE[k,],percent.burnin,thinning,param.name=paste("aE_",k)) } 
			plot_trace(a2,percent.burnin,thinning)
			plot_trace(beta,percent.burnin,thinning)
				if(exists("phi_mat")){
					plot_all_phi_trace(phi_mat,percent.burnin,thinning,population.names)
				}		
			plot_trace(LnL_thetas,percent.burnin,thinning)
			plot_trace(LnL_counts,percent.burnin,thinning)
			plot_trace(Prob,percent.burnin,thinning)
		})
	}

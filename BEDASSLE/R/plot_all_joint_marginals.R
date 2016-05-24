plot_all_joint_marginals <-
function(MCMC.output,percent.burnin=0,thinning=1){
	MCMC.output.list <- load_MCMC_output(MCMC.output)
	with(MCMC.output.list, {
			plot_joint_marginal(a0,aD,percent.burnin,thinning)
				devAskNewPage(ask=TRUE)
            for (k in 1:nrow(aE)) { plot_joint_marginal(a0,aE[k,],percent.burnin,thinning,param.name2=paste("aE_",k)) }
			plot_joint_marginal(a0,a2,percent.burnin,thinning)
			plot_joint_marginal(a0,beta,percent.burnin,thinning)
            for (k in 1:nrow(aE)) { plot_joint_marginal(aD,aE[k,],percent.burnin,thinning,param.name2=paste("aE_",k))  }
			plot_joint_marginal(aD,a2,percent.burnin,thinning)
			plot_joint_marginal(aD,beta,percent.burnin,thinning)
            for (k in 1:nrow(aE)) { 
                plot_joint_marginal(aE[k,],a2,percent.burnin,thinning,param.name1=paste("aE_",k)) 
                plot_joint_marginal(aE[k,],beta,percent.burnin,thinning,param.name1=paste("aE_",k))
                plot_joint_marginal(aE[k,]/aD,a0,percent.burnin,thinning,param.name1=paste("aE_",k)) 
                plot_joint_marginal(aE[k,]/aD,a2,percent.burnin,thinning,param.name1=paste("aE_",k))
                plot_joint_marginal(aE[k,]/aD,beta,percent.burnin,thinning,param.name1=paste("aE_",k))
            }
			plot_joint_marginal(a2,beta,percent.burnin,thinning)
		})
	}

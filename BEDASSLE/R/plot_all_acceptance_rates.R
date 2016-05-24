plot_all_acceptance_rates <-
function(MCMC.output){
	MCMC.output.list <- load_MCMC_output(MCMC.output)
	with(MCMC.output.list, {
			plot_acceptance_rate(aD_accept,aD_moves)
				devAskNewPage(ask=TRUE)
			plot_acceptance_rate(aE_accept,aE_moves)
			plot_acceptance_rate(a2_accept,a2_moves)
			plot_acceptance_rate(mu_accept,mu_moves)
			plot_acceptance_rate(thetas_accept,thetas_moves)
				if(exists("phi_accept")){
					plot_acceptance_rate(phi_accept,phi_moves)
				}
		})
	}

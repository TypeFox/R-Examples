plot_all_phi_marginals <-
function(phi_mat,percent.burnin=0,thinning=1,population.names=NULL,pop.index=NULL,histogram=TRUE,density=TRUE){
		k <- nrow(phi_mat)
		for(i in 1:k){
			plot_phi_marginal(phi_mat[i,],percent.burnin,thinning,population.names=population.names[i],pop.index=i,histogram,density)
			devAskNewPage(ask=TRUE)
		}
	}

plot_all_phi_trace <-
function(phi_mat,percent.burnin=0,thinning=1,population.names=NULL){
		k <- nrow(phi_mat)
		for(i in 1:k){
			plot_phi_trace(phi_mat[i,],percent.burnin,thinning,population.names=population.names[i],pop.index=i)
			devAskNewPage(ask=TRUE)
		}
	}

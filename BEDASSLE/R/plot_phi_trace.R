plot_phi_trace <-
function(phi,percent.burnin=0,thinning=1,population.names=NULL,pop.index=NULL){
		burnin <- (percent.burnin/100)*length(which(phi!=0))
		x <- seq(from = burnin,to = length(which(phi!=0)),by = thinning)
			Fk <- 1/(1+phi[x])
			if(is.null(population.names)){
				plot(Fk,
					pch=19,
					ylim=c(0,1),
					col=adjustcolor(1,alpha.f=0.4),
					main=sprintf("Trace Plot of F parameter, population %s",pop.index),
					xlab="MCMC sampled generations",
					ylab=sprintf("F parameter, population %s",pop.index)
                )
                abline(h=median(Fk),col="red")	
			}
				if(!is.null(population.names)){
				plot(Fk,
					pch=19,
					ylim=c(0,1),
					col=adjustcolor(1,alpha.f=0.4),
					main=paste("Trace Plot of F parameter,",population.names,sep=" "),
					xlab="MCMC sampled generations",
					ylab=paste("F parameter,",population.names,sep=" ")
                )
                abline(h=median(Fk),col="red")	
			}
	}

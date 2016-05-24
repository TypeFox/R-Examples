plot_phi_marginal <-
function(phi,percent.burnin=0,thinning=1,population.names=NULL,pop.index=NULL,histogram=TRUE,density=TRUE){
		burnin <- (percent.burnin/100)*length(which(phi!=0))
		x <- seq(from = burnin,to = length(which(phi!=0)),by = thinning)
			Fk <- 1/(1+phi[x])
			marginal <- density(Fk,adj=1)
			if(is.null(population.names)){
				plot(1,type="n",
					ylim=c(0,max(marginal$y)+max(marginal$y)/5),
					xlim=c(min(marginal$x),max(marginal$x)),
					main=sprintf("Marginal density of F parameter, population %s",pop.index),
					ylab="density",
					xlab=sprintf("F parameter, population %s",pop.index)
                )
			}
			if(!is.null(population.names)){
				plot(1,type="n",
					ylim=c(0,max(marginal$y)+max(marginal$y)/5),
					xlim=c(min(marginal$x),max(marginal$x)),
					main=paste("Marginal density of F parameter,",population.names,sep=" "),
					ylab="density",
					xlab=paste("F parameter,",population.names,sep=" ")
                )
			}			
			if(histogram){
				hist(Fk,
					freq=FALSE,
					col="gray",
					add=TRUE)
			}
			if(density){
				lines(marginal,adj=1)
				polygon(x=c(0,marginal$x,0),y=c(0,marginal$y,0),col=adjustcolor("blue",0.6))					
			}
				segments(x0=median(Fk),
					y0=0,
					x1=median(Fk),
					y1=max(marginal$y+marginal$y/20),
					col="red",
					lwd=3)
	}

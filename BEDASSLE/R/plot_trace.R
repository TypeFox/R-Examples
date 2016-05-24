plot_trace <-
function(parameter,percent.burnin=0,thinning=1,param.name=deparse(substitute(parameter))){
		burnin <- (percent.burnin/100)*length(which(parameter!=0))
		x <- seq(from = burnin,to = length(which(parameter!=0)),by = thinning)
				plot(parameter[x],
					pch=19,
					col=adjustcolor(1,alpha.f=0.4),
					main=paste("Trace Plot of",param.name,sep=" "),
					xlab="MCMC sampled generations",
					ylab=param.name
                )
                abline(h=median(parameter[x]),col="red")	
	}

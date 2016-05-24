plot_joint_marginal <-
function(parameter1,parameter2,percent.burnin=0,thinning=1,param.name1=deparse(substitute(parameter1)),param.name2=deparse(substitute(parameter2))){
		burnin <- (percent.burnin/100)*length(which(parameter1!=0))
		x <- seq(from = burnin,to = length(which(parameter1!=0)),by = thinning)
			plot(parameter1[x],parameter2[x],
				col=rainbow(start=0.5,end=1.0,length(parameter1[x]),alpha=0.4),
				pch=19,
				main=paste("Joint Marginal Density of",paste(param.name1,param.name2,sep=" and "),sep=" "),
				xlab=param.name1,
				ylab=param.name2,
				ylim=c(min(parameter2[x])-abs(min(parameter2[x])-max(parameter2[x]))/5,max(parameter2[x])))	
					legend(x="bottomright",
						cex=0.6,
						pt.cex=1,
						pch=21,
						col=1,
						pt.bg=c(rainbow(start=0.5,end=1.0,length(parameter1[x]),alpha=0.8)[c(1,floor(length(parameter1[x])/2),length(parameter1[x]))]),
						legend=c("MCMC sampled generation 1",
								paste("MCMC sampled generation",floor(length(parameter1[x])/2)),
								paste("MCMC sampled generation",length(parameter1[x]))))
	}

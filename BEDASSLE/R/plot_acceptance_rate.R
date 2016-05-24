plot_acceptance_rate <-
function(accepted.moves,proposed.moves,param.name=deparse(substitute(accepted.moves))){
        param.name <- strsplit(param.name,split="_")[[1]][1]
		x <- seq(1,length(which(proposed.moves!=0)))
		acceptance.rate <- accepted.moves[x]/proposed.moves[x]
			plot(acceptance.rate,
				pch=20,
				col=adjustcolor(1,alpha.f=0.7),
				xlab="MCMC sampled generations",
				ylab="acceptance rate",
				main=paste(param.name,"acceptance rate",sep=" "),
				ylim=c(0,1))
			abline(h=c(0.2,0.7),col="green",lty="dashed",lwd=2)
			
			if(median(acceptance.rate) > 0.7 || median(acceptance.rate) < 0.2){
				polygon(x=c(0-0.04*length(acceptance.rate),
							0-0.04*length(acceptance.rate),
							length(acceptance.rate)+0.04*length(acceptance.rate),
							length(acceptance.rate)+0.04*length(acceptance.rate)),
						y=c(0-0.04*length(acceptance.rate),
						1+0.04*length(acceptance.rate),
						1+0.04*length(acceptance.rate),
						0-0.04*length(acceptance.rate)),
						col=adjustcolor("red",alpha.f=0.2))
			}
	}

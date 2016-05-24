densityDIFplot <-
function(x,xlim,ylim,xlab,ylab,main,...){

		


		if(missing(main)){main="Observed Score Distribution\n"}
		
		if(missing(xlim)){xlim=c(min(x$subjscore),max(x$subjscore))}
		
		if(missing(ylab)){ylab="Density of Score"}
		xlab <- "Observed Scores"
		
		

		ymax<-max(density(x$subjscore,from=0,to=max(x$subjscore))$y)

		if(missing(ylim)){ylim=c(0,ymax+.02)}


		plot(xlim,ylim,type="n",xlim=xlim, ylim=ylim ,xlab=xlab,ylab=ylab,main = main,...)
		ngrps<-length(x$groups)

		plot_colors <- c("blue","red","forestgreen","black","yellow","orange")

		if(ngrps>6){plot_colors<-c(plot_colors,sample(colors(),ngrps-6))}


		line_type<-rep(1:6,ngrps)


		for(i in 1:ngrps){

			cgrp<-x$DIF[[i]]
			lines(density(cgrp$subjscore,from=0,to=max(cgrp$subjscore)),col=plot_colors[i],lty=line_type[i],...)
		
		}

		legend(min(x$subjscore), ymax+.02,x$groups, cex=0.8, col=plot_colors[1:ngrps],lty=line_type[1:ngrps], lwd=2, bty="n");

		axis(3,at=x$subjscoresummary, labels=labels(x$subjscoresummary),tck=0)
		abline(v=x$subjscoresummary,col="blue",lty=2)
		box()

}


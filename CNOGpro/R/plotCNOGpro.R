plotCNOGpro <- function(experiment){
  par(mfrow=c(1,1))
  if (experiment$is_GC_normalized){obs <- experiment$CorrReadsprWindow}
  else{obs <- experiment$ReadsprWindow}
  plot(x=seq(1,experiment$chrlength, experiment$windowlength), y=obs, xlab="Chromosome coordinate", ylab="Coverage", main=paste("Coverage along chromosome", experiment$accession),cex=0.1)
  par(ask=T)
  if(!is.null(experiment$HMMtable)){
    membership <- getMembership(obs=obs,windowlength=experiment$windowlength,HMMtable=experiment$HMMtable)
    membership <- membership[as.integer(membership) != 0]
    colors <- rainbow(n=length(unique(membership)),alpha=0.5)
    d <- density(obs[membership==1],from=0,to=max(obs))
    plot(d, xlab="Coverage", main="Read count distribution",xlim=c(0,max(obs)))
    polygon(d, col=colors[1])
    if(!max(membership)==1){
      for (state in 2:max(membership)){
        if (length(obs[membership==state]) < 10){cat("Not enough data points for copy number state ", state, ", skipping...\n");next}
        d <- density(obs[membership==state],from=0,to=max(obs))
        lines(d)
        polygon(d,col=colors[state])
      }
      legend("topright",legend=unique(membership),fill=colors)
    }
  }
  par(mfrow=c(1,2))
  boxplot(experiment$ReadsprWindow ~ experiment$GCperwindow, varwidth=T, outline=F,main="Before GC normalization",ylab="Coverage",xlab="GC percentage")
  if(!is.null(experiment$CorrReadsprWindow)){
    boxplot(experiment$CorrReadsprWindow ~ experiment$GCperwindow,varwidth=T,outline=F, main="After GC normalization",ylab="Coverage",xlab="GC percentage")
  }
  par(ask=F)
}
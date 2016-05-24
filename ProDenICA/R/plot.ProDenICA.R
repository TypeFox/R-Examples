plot.ProDenICA=function(x,...){
  null.check=is.null(x$density[[1]])
  if(null.check){
    warning("Nothing to plot; only implemented with 'Gpois'")
  }
  else{
    dens=x$density
    for(i in seq(dens)){
      plot(dens[[i]],xlab=paste("S",i,sep=""),ylab=paste("Density(S",i,")",sep=""),main=paste("Source Component",i),type="l")
    }
  }
  invisible()
    
}
      
      

## plot.nsCosinor.R
## Plots results from nsCosinor using ggplot2
plot.nsCosinor<-function(x, ...){

  ## basic variables
  cycles<-x$cycles;
  k<-length(cycles);
  month<-(12*(x$time-floor(x$time)))+1
  
  # season
  smat<-as.matrix(x$season)
  # loop through seasons
  for (index in 1:k){
    mean<-smat[,(index*3)-2]
    lower<-smat[,(index*3)-1]
    upper<-smat[,(index*3)]
    type<-paste("Season, cycle=",cycles[index],sep="")
    this.frame=data.frame(time=x$time,mean=mean,lower=lower,upper=upper,type=type)
    if(index==1){season.frame=this.frame}else{season.frame=rbind(season.frame,this.frame)}
  }

  # trend
  trend.frame=data.frame(time=x$time,mean=x$trend$mean,lower=x$trend$lower,upper=x$trend$upper,type='Trend')
  plot.frame=rbind(trend.frame,season.frame)
  
  # plot with ribbon
  gplot = ggplot(plot.frame, aes(time, mean)) +
    geom_ribbon(aes(ymin=lower, ymax=upper, alpha=5),show_guide = FALSE)+
    geom_line()+
    theme_bw()+
    xlab('Time') +        
    ylab(' ') +        
    facet_grid(type~.,scales='free_y')
  print(gplot)
  
  # return
  return(gplot)
}

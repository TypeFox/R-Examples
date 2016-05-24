# plot.Monthmean.R
plot.Monthmean<-function(x,...){
## Check
  if (class(x) != "Monthmean"){
    stop("Object must be of class 'Monthmean'")
  } 
## Plot
  par(lwd=2)
  plot(x$mean,type='o',bty='n',xaxt='n',xlab='Month',ylab='Mean',...)
  box(lwd=1)
  x.ticks=seq(1,12,1)
  x.labels=c('J','F','M','A','M','J','J','A','S','O','N','D')
  axis(side=1,at=x.ticks,labels=x.labels)
}


Pchart <-
function(behavior,groupX,bandX,ABxlab,ABylab, ABmain){
  
  sumg<-tapply(behavior,groupX,sum,na.rm=T)
  tn<-table(groupX)
  len<-length(tn)
  A<-seq(1:len)
  p=sumg/tn
  end<-which(is.na(behavior))
  
  e<-end[1]
  e=e-1
  e<-groupX[e]
  np<-length(end)
  A<-p[1:e]
  meanA<-mean(A,na.rm=T)
  sdA<-sd(A,na.rm=T)
  Uband<-meanA+bandX*sdA
  Lband<-meanA-bandX*sdA
  x=(1:length(p))
  
  j=1
  while (j <= np){
    e<-end[j]
    e<-e-1
    e<-groupX[e]
    e<-e+j
    
    p<-insert(p,NA,e)
    x<-insert(x,NA,e)
    j=j+1
  }
  min=Lband-.3
  max=Uband+.3
  xmax<-length(p)+1
  
  graphics.off()
  plot.new()
  layout(rbind(1,2), heights=c(6,1))
 plot(x,p,ylim=c(min,max),xlim=c(1,xmax),type="o",col="red",bty='L',main=ABmain,xlab=ABxlab,ylab=ABylab)
  
  abline(h=meanA,col="green")
  abline(h=Uband,col="blue")
  abline(h=Lband,col="orange")
  
  #par(mar=c(1, 1, 1, 1))
  #plot.new()
  #legend("center", c("behavior","Uband","mean","Lband"), col = c("red","blue", "green","orange"), lwd = 1,ncol=4,bty ="n")
  
  puband<-c("Uband=",round(Uband,3))
  pmean<-c("mean= ",round(meanA,3))
  plband<-c("Lband=",round(Lband,3))
  print(puband)
  print(pmean)
  print(plband)
  ab<-NULL
  
  ab<<-recordPlot()
}

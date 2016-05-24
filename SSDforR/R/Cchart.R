Cchart <-
function(behavior,phaseX,v1,bandX,ABxlab,ABylab, ABmain){
  
  maxy=which.max(behavior)
  max<-behavior[maxy]
  miny=which.min(behavior)
  min<-behavior[miny]
  
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  bmean=mean(A)
  
  btsd<-sqrt(bmean)
  
  Uband=round(btsd*bandX+bmean)
  Lband=round(bmean-btsd*bandX)
  #min=Lband-3
  #max=Uband+3
  y<-na.omit(behavior)
  total=length(y)
  x=(1:total)
  
  end<-which(is.na(phaseX))
  np<-length(end)
  j=1
  while (j <= np){
    e<-end[j]
    
    
    
    y<-insert(y,NA,e)
    x<-insert(x,NA,e)
    j=j+1
  }
  f1=Uband >maxy
  f2=Lband <miny
  
  if (f1==TRUE)
  {maxy=Uband+2}
  
  
  if (f2==TRUE)
  {miny=Lband-2}
  
  graphics.off()
  layout(rbind(1,2), heights=c(6,1))
  plot(x,y,ylim=c(miny,maxy),type="o",col="red",ylab=ABylab,xlab=ABxlab,main=ABmain,bty="l")
  
  abline(h=bmean,col="green")
  abline(h=Uband,col="blue")
  abline(h=Lband,col="orange")
  
  # par(mar=c(1, 1, 1, 1))
  # plot.new()
  #legend("center", c("behavior","Uband","mean","Lband"), col = c("red","blue", "green","orange"), lwd = 1,ncol=4,bty ="n")
  
  puband<-c("Uband=",Uband)
  pmean<-c("mean=",round(bmean,2))
  plband<-c("Lband=",Lband)
  print(puband)
  print(pmean)
  print(plband)
  ab<-NULL
  
  ab<<-recordPlot()
}

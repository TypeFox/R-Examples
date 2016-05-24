ABstat <-
function(behavior,phaseX,v,statX){
 
  writeLines("-------------------------------------------------------------------------------------")
  writeLines("Click the mouse in the beginning of the phase you want the line in")
  writeLines("-------------------------------------------------------------------------------------")
  abmedian<-tapply(behavior, phaseX,statX)
  tphase<-table(phaseX)
  mlin<- tphase[names(tphase)==v]
  omedian<- abmedian[names(abmedian)==v]
  stats<-c(statX,round(omedian,3))
 print(stats)
  l<-locator(1)
  mlin<-l$x+mlin
  if (statX=="mean") {cl="blue"}   
  if (statX=="median") {cl="darkseagreen"}  
  if (statX=="sd") {cl="darkorange"} 
  if (statX=="mean") {tl=1}   
  if (statX=="median") {tl=2}  
  if (statX=="sd") {tl=4}
  segments(x0=l$x,x1=mlin,y0=omedian,col=cl,lty=tl,lwd=3)
  u<-readline("accept line? (y/n) ")
  if (u=="n")
  {replayPlot(ab)}
  
    ab<-NULL
    
  ab<<-recordPlot()
}
ABlineD <-
function(behavior){
  writeLines("-------------------------------------------------------------------------------------")
  writeLines("Click the mouse in the gap between the phases you want the line in.")
  writeLines("-------------------------------------------------------------------------------------")
 
  l<-locator(1)
  
 
  ymin=min(behavior,na.rm = TRUE)
  ymax=max(behavior,na.rm = TRUE)
  ymax=ymax+1
  segments(x0=l$x,y0=ymin,y1=rep(ymin:ymax),lty=2)
  
 
  
  
  u<-readline("accept line? (y/n) ")
 if (u=="n")
  {replayPlot(ab)}
 
    ab<-NULL
    ab<<-recordPlot()
  {assign('ab', recordPlot())}
  
}

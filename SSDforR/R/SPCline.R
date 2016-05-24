SPCline<-function(){
  instruct<-"Click the mouse on the top of the graph where the vertical line should begin. Repeat at the bottom of the graph for where the line should end.  "
  
  l<-locator(1)
  
  l1<-locator(1)
  
  # check<-l$y > l1$y
  
  # swap=l1$y
  #  if(check==TRUE) {l1$y=l$y;l$y=swap} 
  
  
  
  segments(x0=l$x,y0=l$y,y1=rep(l$y:l1$y))
  
 
  u<-readline("accept line? (y/n) ")
  if (u=="n")
  {replayPlot(ab)}
  
  ab<-NULL
  ab<<-recordPlot()

}

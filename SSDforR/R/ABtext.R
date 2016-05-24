
ABtext <-
function(textx){
  text<-NULL
  writeLines("------------------------------------------------------------------------")
  writeLines("Click the mouse where you want the text to begin.")
  writeLines("------------------------------------------------------------------------")
  
  text(locator(1),c(textx),cex=.9)
  u<-readline("accept text? (y/n) ")
  if (u=="n")
  {replayPlot(ab)}
  
    
   
    ab<-NULL
  ab<<-recordPlot()
  
}

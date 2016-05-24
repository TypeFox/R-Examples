ABplot <-
function(behavior,phaseX,ABxlab,ABylab, ABmain){
 
  maxy=which.max(behavior)
  max<-behavior[maxy]+1
  numx<-sum(!is.na(behavior))+3
  miny=which.min(behavior)
  min<-behavior[miny]
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
  
  
  
  graphics.off()
  plot(x,y, ylim=c(min,max),xlim=c(0,numx),lwd=2,type="o",col="red", xlab=ABxlab, ylab=ABylab, main=ABmain,bty='L' )
 
  
  ab<-NULL
  ab<<-recordPlot()
  
  
  writeLines("-------------------------------------------------------------------------------------")
  writeLines("1-You can add lines between phases by using the ABlines() function.")
  writeLines("2-You can add text by using the ABtext() function.")
  writeLines("3-You can add a mean, median, sd line by using the ABstat() function.")
  writeLines("4-You can also add a goal line using the Gline() function.")
  writeLines("-------------------------------------------------------------------------------------")
}

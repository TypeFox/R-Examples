plotINF <- function(x, ...) # NEW NAME
# plot Inf curves
{
  
  THE  <- x$Catinf$thetas
  GROUPnames <- names(x$Catinf$TestInfGROUPS)
  
  par("ask"=TRUE) 


  for(gru in 1:length(x$Catinf$catinfG)) # loops groups 
  {  # loops the groups
    
    
    for(ITEM in 1:length(x$Catinf$catinfG[[1]])) # loops items
    {
      
      MAIN <- paste("Group",GROUPnames[gru],"- Item",ITEM)
    #x$Catinf$catinfG$A$Item1
    
    NAMES <- c("overall",gsub("dt_","cat",x$reshOBJ$aDD[[ITEM]]$categ))
    
    plot(THE,rowSums(x$Catinf$catinfG[[gru]][[ITEM]]),type="l", xlab=expression(theta),ylab="Information", lty=2, lwd=2, main = MAIN, ...)
    
    for(i in 1:ncol(x$Catinf$catinfG[[gru]][[ITEM]]))
    {
      lines(THE,x$Catinf$catinfG[[gru]][[ITEM]][,i],type="l",col=i+1) 
    }
    
        legend("topright",legend=NAMES,lty=1,col=1:(ncol(x$Catinf$catinfG[[gru]][[ITEM]])+1),x.intersp=0.1,y.intersp=0.4,seg.len=0.5,bty="n")
    
    }
    
    
  }  
  
################################################################
  
  for(gru in 1:length(x$Catinf$catinfG)) # loops groups 
  {  # loops the groups
  plot(THE, x$Catinf$TestInfGROUPS[[gru]], main=paste("Test Information Function \n group:", GROUPnames[gru]), type="l", lwd=1.3, xlab=expression(theta),ylab="Information")
  }

  
}






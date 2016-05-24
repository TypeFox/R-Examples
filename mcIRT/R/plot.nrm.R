plot.nrm <-
function(x, numbpoints=100, fromto=c(-4,4),...)
{

  
par("ask"=TRUE)


  gnames <- names(x$ZLpar)
  for(gru in 1:length(x$ZLpar)) # loops groups 
    {
    
    for(ITEM in 1:length(x$ZLpar[[gru]])) # loops items
    {
      
      Km  <- matrix(c(rep(1,numbpoints),seq(fromto[1],fromto[2],length=numbpoints)),ncol=2)
      mat <- matrix(x$ZLpar[[gru]][[ITEM]] ,nrow=2,byrow=TRUE)
      
      Z <- Km %*% mat
      ez <- exp(Z)
      ezrs <- rowSums(ez)        
      ZQstern <- ez / ezrs
      
      MAIN <- paste("Group",gnames[gru],"- Item",ITEM)

      for(i in 1:ncol(mat)) # loops categories
      {
        if(i == 1)
        {
          plot(seq(fromto[1],fromto[2],length=numbpoints),ZQstern[,i],type="l",col=i,main=MAIN,xlab=expression(theta),ylab="Prob",...)
          
        } else {
          MAIN <- paste("Group",gnames[gru],"- Item",i)
          lines(seq(fromto[1],fromto[2],length=numbpoints),ZQstern[,i],col=i)
                }

        inames <- x$reshOBJ$aDD[[ITEM]]$categ
        legend("right",legend=inames,lty=1,col=1:ncol(mat),x.intersp=0.1,y.intersp=0.7,seg.len=0.4,bty="n")
      }
      
    }
    
  }



}

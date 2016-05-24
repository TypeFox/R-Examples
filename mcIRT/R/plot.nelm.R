plot.nelm <-
function(x,numbpoints=100,fromto=c(-4,4),...)
{
  
  par("ask"=TRUE)
  

  for(gru in 1:nlevels(x$reshOBJ$gr)) # loops groups 
  {
    
    for(ITEM in 1:length(x$reshOBJ$aDD)) # loops items
    {
      
      
    
      Km  <- matrix(c(rep(1,numbpoints),seq(fromto[1],fromto[2],length=numbpoints)),ncol=2)
      abpar <- c(x$ZLpar$beta[[gru]][ITEM],x$ZLpar$alpha[[gru]][ITEM])
      
      solit <- twoplpart(Km=Km, abpar=abpar)
      dosolit <- 1-solit
      tplpart <- cbind(solit,dosolit)
      
      pitemNRM <-  unlist(x$ZLpar$nrmpar[[gru]][[ITEM]])

      ZQstern <- coP_nrm(pitemNRM,Km)
      
      ZQstern_nlm <- ZQstern * as.vector(dosolit)      
      ZQstern_all <- cbind(solit,ZQstern_nlm)

      MAIN <- paste("Group",levels(x$reshOBJ$gr)[gru],"- Item",ITEM)

      anzcateg <- length(x$reshOBJ$aDD[[ITEM]]$tabcat)
      for(i in 1:anzcateg) # loops categories
      {
        if(i == 1)
        {
          plot(seq(fromto[1],fromto[2],length=numbpoints),ZQstern_all[,i],type="l",col=i,main=MAIN,xlab=expression(theta),ylab="Prob",...)
          
        } else {
          MAIN <- paste("Group",levels(x$reshOBJ$gr)[gru],"- Item",ITEM)
          lines(seq(fromto[1],fromto[2],length=numbpoints),ZQstern_all[,i],col=i)
        }
        inames <- x$reshOBJ$aDD[[ITEM]]$categ
        legend("right",legend=inames,lty=1,col=1:anzcateg,x.intersp=0.05,y.intersp=0.7,seg.len=0.4,bty="n",cex=0.8)
      }
      
    }
    
  }
  
  
  
}

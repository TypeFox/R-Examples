CentMnrm <-
function(reshOBJ)
{
  
  
  catanz <- rep(rep(sapply(reshOBJ$aDD,function(x)x$anz_cat),each=2),each=nlevels(reshOBJ$gr))
  
  matList <- lapply(catanz,function(scr)
  {
    toD <- matrix( - 1/scr,ncol=scr,nrow=scr)
    diag(toD) <- diag(toD) + 1
    toD
  })  
  
  
  to2 <- cumsum(catanz)
  from1   <- c(1,to2[-length(to2)]+1)
  
  cem1 <- matrix(0,ncol=sum(catanz),nrow=sum(catanz))
  
  
  for(i in 1:length(from1))
  {
    fr <- from1[i]
    tt <- to2[i]
    
    cem1[fr:tt,fr:tt] <- matList[[i]]
  }
  
  cem1  
}

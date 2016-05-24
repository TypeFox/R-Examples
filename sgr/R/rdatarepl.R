#rm(list=ls())
#load("~/lavori/proveR/SGR/L.rda")

#Dx <- matrix(L[,2])
#RM <- replacement.matrix(2,p=c(0,.3))
#j <- 1; i <- 1

rdatarepl <- function(Dx,RM,printfp=TRUE) {
  
  Fx <- Dx
  for (j in 1:ncol(Dx)) {
    for (i in 1:nrow(RM)) {
      (quali <- which(Dx[,j]==i) )
      Fx[quali,j] <- sample(1:ncol(RM),length(quali),TRUE,prob=RM[i,])
    }
  }
  
  Delta <- Dx-Fx
  Delta[Delta!=0] <- 1
  Fperc <- sum(Delta)/(prod(dim(Delta)))*100
  
  if (printfp) cat(paste(round(Fperc,2),"% of data replaced.",sep=""),"\n")
  return(list(Fx=Fx,Fperc=Fperc))
  
}

#rdatarepl(Dx+1,RM)
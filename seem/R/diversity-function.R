diversity <- function(output,param,fileout,pdfout){

  npar <- length(param$plab)

  if(npar==1) {
    labrun <- param$plab; valrun <- param$pval   
    nruns <- length(param$pval); param$pval <- matrix(param$pval)  
  } else {
    nruns <- dim(param$pval)[1]; valrun <- array()
    labrun <- param$plab[1] 
    for(i in 2:npar) labrun <- paste(labrun,param$plab[i],sep=";")
    for(j in 1:nruns){
     valrun[j] <- param$pval[j,1]
     for(i in 2:npar) valrun[j] <- paste(valrun[j],param$pval[j,i],sep=";")  
     }
  }
  legrun <- paste(labrun,"=",valrun,sep="")
  t <- output[,1]; nt <- length(t)
  nv <- (dim(output)[2]-1)/nruns

  # variables by run
  if(nv>1){
   nvars <- matrix(nrow=nruns, ncol=nv)
   outX <- output[,-1]
   for(i in 1:nruns){for(j in 1:nv) nvars[i,j] <- i+nruns*(j-1)}
   Y <- list(); for(i in 1:nruns) Y[[i]] <- outX[,nvars[i,]]
  }
  # allocate list to store proportions and matrix to store Shannon's
  prop <- Y
  H <- matrix(nrow=nt,ncol=nruns)
  # loop thru runs, time, and variables
  for(i in 1:nruns){
   for(j in 1:nt){
     for(k in 1:nv){
       # calculate proportions
       prop[[i]][j,k] <- Y[[i]][j,k]/sum(Y[[i]][j,])
     }
     # calculate Shannon's index
     H[j,i] <- sum(-prop[[i]][j,]*log(prop[[i]][j,]))
   }
  }
  if(pdfout==T) pdf(paste(fileout,".pdf",sep=""))
   mat<- matrix(1:2,2,1,byrow=T)
   layout(mat,c(7,7),c(3.5,3.5),respect=TRUE)
   par(mar=c(4,4,1,.5), xaxs="i", yaxs="r")

   xext <- 1.3*max(t) 
   matplot(t,H,type="l",col=1,xlim=c(0,xext), ylab="Shannon Diversity") 
   legend("topright",legend=legrun,lty=1:3,col=1)
  if(pdfout==T) dev.off()
   
   output <- data.frame(t,H)
   return(output)
}

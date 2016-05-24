 vars.out <- function(prefix, lab.out="", input, param, output, tXlab, pdfout=F){
 # produces output files for multiple runs and variables

# form labels
  if(is.null(param)) {legrun<-"Nominal";nruns<-1;lab.out<-legrun} else {
   npar <- length(param$plab)
   lab.out <- paste(lab.out,param$plab[1],sep=""); if(length(param$plab)>1)
   for(i in 2:length(param$plab)) lab.out <- paste(lab.out,param$plab[i],sep="")

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
  }
  legvar <- tXlab[-1]
  nv <- length(tXlab)-1
  X <- list()
  t <- output[,1]
  # runs by variable
  for(i in 1:nv) X[[i]] <- output[,((i-1)*nruns+2):(i*nruns+1)]

  # variables by run
  if(nv>1){
   nvars <- matrix(nrow=nruns, ncol=nv)
   outX <- output[,-1]
   for(i in 1:nruns){for(j in 1:nv) nvars[i,j] <- i+nruns*(j-1)}
   Y <- list(); for(i in 1:nruns) Y[[i]] <- outX[,nvars[i,]]
  }
  ly <- list(m=matrix(1:2,2,1,byrow=T), w=c(7,7), h=c(3.5,3.5), r=TRUE)

# produce state graphs
  if(pdfout==T) {
   pdf(file=paste(prefix,"-",lab.out,"-out.pdf",sep=""))
   layout(ly$m,ly$w,ly$h,ly$r); par(mar=c(4,4,1,.5), xaxs="i", yaxs="i")

   for(i in 1:nv){
    xext <- 1.3*max(t) 
    matplot(t, X[[i]], type="l", xlim=c(0,xext), xlab=tXlab[1], ylab=tXlab[i+1], col=1,lwd=1.3)
    legend ("topright", legend=legrun, lty=1:nruns,col=1, lwd=1.3,cex=0.7)
    }
  
# produce phase plane for each pair
   layout(ly$m,ly$w,ly$h,ly$r); par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")
   for(i in 1:(nv-1)){
    for(j in (i+1):nv){
    xext <- 1.3*max(X[[i]])  
    matplot(X[[i]], X[[j]], type="l", xlim=c(0,xext), 
          xlab= tXlab[i+1], ylab=tXlab[j+1], col=1,lwd=1.3)
    legend ("topright", legend=legrun, lty=1:nruns,col=1, lwd=1.3,cex=0.7)
   }
  } 

# produce multivariate plots for each run
  if(nv>1){
   layout(ly$m,ly$w,ly$h,ly$r); par(mar=c(4,4,1,.5), xaxs="i", yaxs="i")
   for(i in 1:nruns){
    xext <- 1.3*max(t) 
    matplot(t, Y[[i]], type="l", xlim=c(0,xext), xlab=tXlab[1], ylab="Variables", col=1,lwd=1.3)
    legend ("topright", legend=legvar, lty=1:nv,col=1, lwd=1.3,cex=0.7)
    text(1.2*max(t), min(Y[[i]])+0.1*max(Y[[i]]), legrun[i])
    }
  }
  dev.off()

 } # end of if pdfout=T 

# produce state graphs
  if(pdfout==F){
   win.graph();layout(ly$m,ly$w,ly$h,ly$r); par(mar=c(4,4,1,.5), xaxs="i", yaxs="i")
   for(i in 1:nv){
    xext <- 1.3*max(t) 
    matplot(t, X[[i]], type="l", xlim=c(0,xext), xlab=tXlab[1], ylab=tXlab[i+1], col=1,lwd=1.3)
    legend ("topright", legend=legrun, lty=1:nruns,col=1, lwd=1.3,cex=0.7)
    if(i%%2==0 && i<nv) {win.graph(); layout(ly$m,ly$w,ly$h,ly$r);par(mar=c(4,4,1,.5), xaxs="i", yaxs="i")}
    }
  
# produce phase plane for each pair
   k <- 0
   win.graph();layout(ly$m,ly$w,ly$h,ly$r);par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")
   for(i in 1:(nv-1)){
    for(j in (i+1):nv){
    xext <- 1.3*max(X[[i]]) 
    matplot(X[[i]], X[[j]], type="l", xlim=c(0,xext),
          xlab= tXlab[i+1], ylab=tXlab[j+1], col=1,lwd=1.3)
    legend ("topright", legend=legrun, lty=1:nruns,col=1, lwd=1.3,cex=0.7)
    k <- k+1
    if(k%%2==0 && i <(nv-1)) {win.graph(); layout(ly$m,ly$w,ly$h,ly$r);par(mar=c(4,4,1,.5), xaxs="i", yaxs="i")}
   }
  } 

# produce multivariate plots for each run
   if(nv>1){
   win.graph();layout(ly$m,ly$w,ly$h,ly$r); par(mar=c(4,4,1,.5), xaxs="i", yaxs="i")
   for(i in 1:nruns){
    xext <- 1.3*max(t) 
    matplot(t, Y[[i]], type="l", xlim=c(0,xext), xlab=tXlab[1], ylab="Variables", col=1,lwd=1.3)
    legend ("topright", legend=legvar, lty=1:nv,col=1, lwd=1.3,cex=0.7)
    text(1.2*max(t), min(Y[[i]])+0.1*max(Y[[i]]), legrun[i])
    if(i%%2==0 && i<nruns) {win.graph(); layout(ly$m,ly$w,ly$h,ly$r);par(mar=c(4,4,1,.5), xaxs="i", yaxs="i")}
    }
   }

 } # end of if pdfout=F  

 # writes output file as csv
  fileout.csv=paste(prefix,"-",lab.out,"-out.csv",sep="")
  write.table(input,fileout.csv, sep=",", row.names=F, quote=F)
  name.param <-t(c("Runs",c(legrun))) 
  write.table(name.param, fileout.csv, sep=",",
              col.names=F, row.names=F, append=T, quote=F)
  name.out <- t(names(output)) 
  write.table(name.out,fileout.csv, sep=",", col.names=F, 
              row.names=F, append=T, quote=F)
  write.table(output,fileout.csv, sep=",", col.names=F, row.names=F, 
              append=T, quote=F)

} # end of function

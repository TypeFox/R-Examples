# plot output file from semi

read.plot.semi.out <- function(fileout,spp,label,pdfout=F,plotmarkov=F,ctd=F){

 nspp <- length(spp);nvar <- nspp + 2
 init <- scan(fileout, skip=4,nlines=1)
 P <- matrix(scan(fileout, skip=6,nlines=nspp),ncol=nspp,byrow=T)
 # label for units
 timeunit <- "Years"; yunit <- "Fraction of Cover"
 X <- matrix(ncol=nspp,nrow=20)
 X[1,] <-init; yrs <- array(); yrs[1] <- 0
 for (i in 2:20){
  yrs[i] <- i
  X[i,] <- P%*%X[i-1,]
}

if(pdfout==T) pdf(paste(fileout,".pdf",sep=""))
if(ctd==F){
 mat<- matrix(1:2,2,1,byrow=T)
 layout(mat,c(7,7),c(3.5,3.5),respect=TRUE)
 par(mar=c(4,4,1,.5), xaxs="r", yaxs="r")
}

# plot cover by spp
if(plotmarkov==T){
 matplot(yrs, X, type="l", xlab=timeunit, ylab=yunit, ylim=c(0,1),
         xlim=c(0,1.5*max(yrs)),lty=c(1:nspp), col=1)
 legend(1.05*max(yrs),1, spp[1:nspp], lty=1:nspp, col=1)
 mtext(side=3,line=-1,"Markov Chain")
}

out <- matrix(scan(fileout,skip=(7+nspp)), byrow=T, ncol=nvar)
#rename output
time <- out[,1]
cover<- out[,2:(nvar-1)]
# plot cover by spp
matplot(time*10, cover, type="l", xlab=timeunit, ylab=yunit, ylim=c(0,1),
         xlim=c(0,1.5*max(time*10)),lty=1:nspp, col=1)
legend(1.05*max(time*10),1, spp[1:nspp], lty=1:nspp, col=1)
mtext(side=3,line=-1,label)

if(pdfout==T) dev.off()

return(list(X=X,out=out))
}

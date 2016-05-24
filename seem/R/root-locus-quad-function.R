# ------------------------------
# quadratic characteristic eqn root locus

root.locus.quad <- function(det,Tr){

# dim and arrays
ndet <- length(det);nTr <- length(Tr)
D <- matrix(nrow=nTr,ncol=ndet) 
ReL <- structure(1:(2*nTr*ndet),dim=c(nTr,ndet,2))
ImL <- structure(1:(2*nTr*ndet),dim=c(nTr,ndet,2))

# calculations
for(j in 1:ndet){
 for(i in 1:nTr){
  D[i,j] <- Tr[i]^2-4*det[j]
  if(D[i,j] <=0) {
   ReL[i,j,1] <- Tr[i]/2
   ReL[i,j,2] <- Tr[i]/2
   ImL[i,j,1] <- sqrt(abs(D[i,j]))/2
   ImL[i,j,2] <- -sqrt(abs(D[i,j]))/2
  } else {
   ReL[i,j,1] <- Tr[i]/2 + sqrt(D[i,j])/2
   ReL[i,j,2] <- Tr[i]/2 - sqrt(D[i,j])/2
   ImL[i,j,1] <- 0
   ImL[i,j,2] <- 0
  }
 }
}

return(list(det=det, Tr=Tr, D=D, ReL=ReL, ImL=ImL))
} # end function

# plots
D.plot <- function(x,i){
 matplot(x$Tr,x$D[,],type="l",col=1,xlab="Trace",ylab="Discriminant")
 abline(h=0)
 if(i==3 || i==4) pos <- "topleft"
 if(i==1 || i==2) pos <- "topright"
 legend(pos,legend=paste("Det=",x$det),lty=1:length(x$det))
}

ReIm.plot <- function(x,j,k){

 ReL= x$ReL; ImL=x$ImL; det = x$det; Tr=x$Tr
 ndet=length(det); nTr <- length(Tr)
 x1 <- min(ReL); x2 <- max(ReL)
 y1 <- min(ImL); y2 <- max(ImL)

 plot(ReL[,1,k],ImL[,1,k],type="n",xlab="Re",ylab="Im",
      xlim=c(x1,x2),ylim=c(y1,y2))
 abline(h=0,col="grey");abline(v=0,col="grey")
 for(i in 1:ndet){
  lines(ReL[,i,k],ImL[,i,k],lty=i,lwd=1.7)
  tail <- c(ReL[1,i,k],ImL[1,i,k]) 
  head <- c(ReL[30,i,k],ImL[30,i,k]) 
  if ((head[1]-tail[1])==0 && (head[2]-tail[2])==0) a=1 
  else  arrows(tail[1],tail[2],head[1],head[2],length=.1,lwd=1.7)
 }

 if(j==1 || j==2) pos <- "topleft"
 if(j==3 || j==4) pos <- "topright"
 legend(pos,legend=paste("Det=",det),lty=1:ndet,cex=0.7)
 mtext(line=-1,paste("Tr from",Tr[1],"to",Tr[nTr]),cex=0.7)
 title(paste("Eigenvalue", k),cex.main=0.7)

 } # end of function






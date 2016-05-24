"PCdiagplot" <-
function(x,PCobj,crit=c(0.975,0.99,0.999),ksel=NULL,plot=TRUE,
                    plotbw=TRUE,raw=FALSE,colgrid="black", ...){
#
#  raw ... if raw==TRUE the raw SDist is computed, without using the median correction
#
# ksel ... select range for values of k
# colgrid ... color of the grid

if (is.null(PCobj$scores)) stop("No PC scores have been computed! Provide scores!")

n <- nrow(PCobj$sco)
k <- ncol(PCobj$loadings)
if (is.null(ksel)) {ksel <- seq(1,k)} 
kl <- length(ksel)
if (ksel[kl]>k) {
	ksel <- seq(1,k)
	warning("Not so many PCs available as specified -> set to number of available PCs!")
}

# initialize output matrices
SDist <- matrix(NA,nrow=n,ncol=kl)
ODist <- matrix(NA,nrow=n,ncol=kl)
critOD <- matrix(NA,kl,length(crit))
critSD <- matrix(NA,kl,length(crit))

for (k in 1:kl){
  # compute score distances:
  if (ksel[k]==1){
    SDist[,k]=abs(PCobj$scores[,1])/PCobj$sdev[1]
  }
  else {
    SDist[,k]=sqrt(apply(t(t(PCobj$scores[,1:ksel[k]]^2)/PCobj$sdev[1:ksel[k]]^2),1,sum))
  }
  if (!raw){
    SDist[,k]=SDist[,k]*sqrt(qchisq(0.5,ksel[k]))/median(SDist[,k])
  }

  # compute orthogonal distances:
  xc <- scale(x,center=PCobj$center,scale=PCobj$scale)
  ODist[,k]=sqrt(apply((xc-PCobj$scores[,1:ksel[k]]%*%t(PCobj$loadings[,1:ksel[k]]))^2,1,sum))

  # compute critical values
  critSD[k,]=sqrt(qchisq(crit,ksel[k]))
  critOD[k,]=(median(ODist[,k]^(2/3))+mad(ODist[,k]^(2/3))*qnorm(crit))^(3/2)
}

# plot results:
if (plot){

require(graphics)
if (plotbw) pg=rev(gray(seq(0,0.7,length=ncol(critSD))))  # gray level to plot
else pg=rev(rainbow(ncol(critSD),start=0.9,end=0.1))  # color level to plot
par(mfrow=c(1,2))

plot(0,0,xlim=c(1,n),ylim=c(min(ksel),max(ksel)),type="n",xlab="Index of the observation",
       ylab="Number of PCs", ...)
title("Orthogonal distance")
for (i in 1:n){
  for (j in 1:length(ksel)){
    howbig=sum(ODist[i,j]>critOD[j,])
    if (howbig==0) { rect(i-0.5,ksel[j]-0.5,i+0.5,ksel[j]+0.5,border=colgrid) }
    else { rect(i-0.5,ksel[j]-0.5,i+0.5,ksel[j]+0.5,col=pg[howbig]) }
  }
}

plot(0,0,xlim=c(1,n),ylim=c(min(ksel),max(ksel)),type="n",xlab="Index of the observation",
       ylab="Number of PCs", ...)
title("Score distance")
for (i in 1:n){
  for (j in 1:length(ksel)){
    howbig=sum(SDist[i,j]>critSD[j,])
    if (howbig==0) { rect(i-0.5,ksel[j]-0.5,i+0.5,ksel[j]+0.5,border=colgrid) }
    else { rect(i-0.5,ksel[j]-0.5,i+0.5,ksel[j]+0.5,col=pg[howbig]) }
  }
}

}
list(ODist=ODist,SDist=SDist,critOD=critOD,critSD=critSD)
}


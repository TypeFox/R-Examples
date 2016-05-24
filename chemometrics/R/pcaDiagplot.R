"pcaDiagplot" <-
function(X,X.pca,a=2,quantile=0.975,scale=TRUE,plot=TRUE, ...){
#
# INPUT:
# X ... mean centered data matrix
# X.pca ... PCA object
# a ... dimension of the PCA space (no. of PCs)
# quantile ... critical quantile
# scale ... if TRUE then X will be scaled - and X.pca should be from scaled data too
# plot ... TRUE if plot should be made, otherwise FALSE
# ... additional plotting arguments
#
# OUTPUT:
# SDist ... Score distances
# ODist ... Orthogonal distances
# critSD ... critical cut-off value for the score distances
# critOD ... critical cut-off value for the orthogonal distances


if (is.null(a)) {a=ncol(X.pca$loa)} 

# compute score distances:
if (is.null(X.pca$sdev)){
  sdev <- apply(X.pca$scores,2,sd)
  SDist=sqrt(apply(t(t(X.pca$sco[,1:a]^2)/sdev[1:a]^2),1,sum))
}
else {
  SDist=sqrt(apply(t(t(X.pca$sco[,1:a]^2)/X.pca$sdev[1:a]^2),1,sum))
}

# compute orthogonal distances:
if (scale){
  if (is.null(X.pca$scale)){
    Xs=scale(X,TRUE,TRUE)
  }
  else {
    Xs=scale(X,center=X.pca$center,scale=X.pca$scale)
  }
  ODist=sqrt(apply((Xs-X.pca$sco[,1:a]%*%t(X.pca$loa[,1:a]))^2,1,sum))
}
else {
  ODist=sqrt(apply((scale(X,TRUE,FALSE)-X.pca$sco[,1:a]%*%t(X.pca$loa[,1:a]))^2,1,sum))
}

# compute critical values
critSD=sqrt(qchisq(quantile,a))
critOD=(median(ODist^(2/3))+mad(ODist^(2/3))*qnorm(quantile))^(3/2)

# plot results:
if (plot){
par(mfrow=c(1,2))
plot(SDist,ylim=c(0,max(SDist)),
     ylab="Score distance SD",xlab="Object number",cex.lab=1.2,...)
abline(h=critSD,lty=2)
plot(ODist,ylim=c(0,max(ODist)),
     ylab="Orthogonal distance OD",xlab="Object number",cex.lab=1.2,...)
abline(h=critOD,lty=2)
}
list(SDist=SDist,ODist=ODist,critSD=critSD,critOD=critOD)
}


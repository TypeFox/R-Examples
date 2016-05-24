locoutNeighbor <-
function(dat,X,Y,propneighb=0.1,variant=c("dist","knn"),usemax=1/3,npoints=50,chisqqu=0.975, 
         indices=NULL, xlab=NULL,ylab=NULL,colall=gray(0.7),colsel=1, ...)
{

if (is.null(ylab)){
  ylab <- paste("Degree of isolation from ",(1-propneighb)*100,"% of the neighbors",sep="")
}
n <- nrow(dat)
p <- ncol(dat)
covr <- covMcd(dat)
cinv <- solve(covr$cov) # inverse of the robust covariance matrix
MDglobal <- sqrt(mahalanobis(dat, covr$center, cinv, inverted=TRUE))

# TRUE/FALSE for non-outlying/outlying:
qchi <- sqrt(qchisq(chisqqu,p))
MDglobalTF <- (MDglobal<qchi)

if (!is.null(indices)){
  indices.reg=indices[MDglobalTF[indices]]
  indices.out=indices[!MDglobalTF[indices]]
}

# pairwise robust Mahalanobis distances:
idx <- matrix(1:n,n,n)
sel <- as.vector(idx[lower.tri(idx)])
hlp <- as.matrix(dat[rep(1:(n-1),seq((n-1),1)),]-dat[sel,])
MDij <- sqrt(rowSums((hlp%*%cinv)*hlp))
MDpair <- matrix(0,n,n)
MDpair[lower.tri(MDpair)] <- MDij
MDpair <- t(MDpair)
MDpair[lower.tri(MDpair)] <- MDij

# Distance matrix:
Xmat <- matrix(rep(t(X),length(X)),ncol=dim(t(X))[2],byrow=TRUE)
Ymat <- matrix(rep(t(Y),length(Y)),ncol=dim(t(Y))[2],byrow=TRUE)
EuclD <- sqrt((c(rep(X,length(X)))-Xmat)^2+(c(rep(Y,length(Y)))-Ymat)^2)
sortD <- apply(EuclD,1,sort,index.return=TRUE)
neighbmatrix <- matrix(unlist(unlist(sortD,recursive=FALSE)[seq(from=2,to=2*n,by=2)]),
  ncol=n,nrow=n,byrow=TRUE)
# neighbmatrix tells for each observation (rows) the indices of the next neighbs (cols)

if (variant=="dist"){
# search the alpha-neighbors (in coordinate space) and compute
# mean pairwise MD (variable space) among them
maxD <- max(EuclD)*usemax
neigbound <- matrix(NA,nrow=n,ncol=npoints)
vec <- seq(from=0,to=maxD,length=npoints)
for (i in 1:n){
  for (j in 1:npoints){
    MDneig <- sort(MDpair[i,EuclD[i,]<=vec[j]])
    if (length(MDneig)>1) {MDneig <- MDneig[-1]}
    MDbound <- MDneig[ceiling(length(MDneig)*propneighb)]
    neigbound[i,j] <- pchisq(MDbound^2,p,MDglobal[i]^2)
  }
}
}
else {
# search the alpha-neighbors (in coordinate space) and compute
# mean pairwise MD (variable space) among them
npoints <- min(n,npoints)
maxn <- trunc(n*usemax)
neigbound <- matrix(NA,nrow=n,ncol=npoints)
vec <- ceiling(seq(from=1,to=maxn,length=npoints))
for (i in 1:n){
  for (j in 1:npoints){
    MDneig <- sort(MDpair[i,neighbmatrix[i,1:vec[j]]])
    if (length(MDneig)>1) {MDneig <- MDneig[-1]}
    MDbound <- MDneig[ceiling(length(MDneig)*propneighb)]
    neigbound[i,j] <- pchisq(MDbound^2,p,MDglobal[i]^2)
  }
}
}

###############################################################################################
#Plots:
par(mfrow=c(1,2))
# plot all observations:
if (variant=="dist"){
  xrange <- c(0,maxD)  
  if (is.null(xlab)) xlab=c("Spatial distance")
}
else {
  xrange <- c(0,maxn)  
  if (is.null(xlab)) xlab=c("Number of neighbors")
}


# Plot only regular observations:
idxreg <- (1:n)[MDglobalTF]
plot(0,0,ylim=c(min(neigbound[idxreg,]),max(neigbound[idxreg,])),xlim=xrange,type="n",
   xlab=xlab, ylab=ylab,...)
title("Regular observations")
for (i in 1:length(idxreg)){
  lines(vec,neigbound[idxreg[i],],col=colall,...)
}

if (!is.null(indices)){ # indices were provided to plot
  if (length(indices.reg)>0){ # some of the indices are from regular observations, so plot them:
    for (i in 1:length(indices.reg)){
      lines(vec,neigbound[indices.reg[i],],col=colsel, ...)
    }
  }
}

# plot only outlying observations:
idxout <- (1:n)[!MDglobalTF]
plot(0,0,ylim=c(min(neigbound[idxout,]),max(neigbound[idxout,])),xlim=xrange,type="n",
   xlab=xlab, ylab=ylab,...)
title("Outlying observations")
for (i in 1:length(idxout)){
  lines(vec,neigbound[idxout[i],],col=colall, ...)
}

if (!is.null(indices)){ # indices were provided to plot
  if (length(indices.out)>0){ # some of the indices are from regular observations, so plot them:
    for (i in 1:length(indices.out)){
      lines(vec,neigbound[indices.out[i],],col=colsel, ...)
    }
  }
}

if (!is.null(indices)){ # indices were provided to plot
  indicesreg <- indices.reg
  indicesout <- indices.out
}
else {
  indicesreg <- idxreg
  indicesout <- idxout
}
list(indices.reg=indicesreg,indices.out=indicesout)
}

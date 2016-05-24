locoutSort <-
function(dat,X,Y,distc=NULL,k=10,propneighb=0.1,chisqqu=0.975,sel=NULL,...)
{

#################################################################################################

neighbdist <-
function(dat,X,Y,dist=NULL,k,propneighb=0.1,chisqqu=0.975,
   xlab="Sorted index",ylab="Distance to neighbor",cex=0.1,pch=3,col=gray(0.6), ...)
{
# select neighbors (by distance of knn) and compute pairwise distances
# plot them according to the fraction of neighbors within a tolerance ellipse
# Separately plot values for small and large global MDs
#
# dat .... data set (without coordinates)
# X ... x-coordinate
# Y ... y-coordinate
# dist ... maximum distance to search for neighbors
# k ... number of nearest neighbors to seach - not taken if !is.null(dist)
# propneighb ... proportion of all neighbors that should be within the
#                tolerance ellipse (at least 1/(#neighbors); 0.1 by default)
# chisqqu ... quantile of the chisquare distribution for splitting the plot
# ... additional parameters for plotting

n <- nrow(dat)
p <- ncol(dat)
covr <- covMcd(dat)
cinv <- solve(covr$cov) # inverse of the robust covariance matrix
MDglobal <- sqrt(mahalanobis(dat, covr$center, cinv, inverted=TRUE))

# TRUE/FALSE for non-outlying/outlying:
qchi <- sqrt(qchisq(chisqqu,p))
MDglobalTF <- (MDglobal<qchi)

# pairwise robust Mahalanobis distances:
idx <- matrix(1:n,n,n)
sel <- as.vector(idx[lower.tri(idx)])
hlp <- as.matrix(dat[rep(1:(n-1),seq((n-1),1)),]-dat[sel,])
MDij <- sqrt(rowSums((hlp%*%cinv)*hlp))
MDpair <- matrix(0,n,n)
MDpair[lower.tri(MDpair)] <- MDij
MDpair <- t(MDpair)
MDpair[lower.tri(MDpair)] <- MDij

# indices of neighbors (in list):
if (!is.null(dist)){
  Tp1 <- matrix(rep(t(X),length(X)),ncol=dim(t(X))[2],byrow=TRUE)
  Tp2 <- matrix(rep(t(Y),length(Y)),ncol=dim(t(Y))[2],byrow=TRUE)
  H <- sqrt((c(rep(X,length(X)))-Tp1)^2+(c(rep(Y,length(Y)))-Tp2)^2)
  diag(H) <- NA
  resdist <- apply(H <= dist,1,which)
}
else {
  n <- length(X)
  nnlist <- vector("list",n)
  distanc <- matrix(0, nrow=n, ncol=n)
  uns <- rep(1,n)

  for (i in 1:n)
  {
    distanc[i,] <- sqrt((X - X[i]*uns)*(X - X[i]*uns) + (Y - Y[i]*uns)*(Y - Y[i]*uns))
    d1 <- distanc[i,]
    names(d1) <- 1:n
    d2 <- sort(d1,index.return=TRUE)
    ind <- which(d2$x[2:(length(X)-1)]==d2$x[3:length(d2$x)])
    pl<-length(ind)
    if(pl!=0)
     {for (j in 1:pl)
       d2$ix[which(d2$x==d2$x[ind[j]+1])]=d2$ix[sample(which(d2$x==d2$x[ind[j]+1]))]
     }
    nnlist[[i]] <- d2$ix[2:(k+1)]
  }
  resdist <- nnlist
}

# guarantee that there is at least 1 neighbor to each data point:
if (any(lapply(resdist,length)==0)) 
  {stop("No neighbor to at least 1 observation -> choose larger distance!")}

# create list with pairwise MDs only for neighbors
MDpairN <- vector("list", n)

# boundary that should include required proportion of neighbors:
chibound <- rep(NA,n)
for (i in 1:n){
  MDpairN[[i]] <- MDpair[i,resdist[[i]]]
  nn <- max(1,round(length(MDpairN[[i]])*propneighb)) # number of neighbors in tolerance ellipse
  nval <- MDpairN[[i]][order(MDpairN[[i]])][nn]  # value of largest neighbor to be included
  chibound[i] <- pchisq(nval^2,p,MDglobal[i]^2)
}

# sort according to values of "chibound" - separately for outliers and non-outliers
idx1 <- order(chibound[MDglobalTF])
idx0 <- order(chibound[!MDglobalTF])
# plot pairwise distances
plot(0,0,xlim=c(1,n),ylim=c(min(unlist(MDpairN)),max(unlist(MDpairN))),
     xlab=xlab,ylab=ylab,type="n",cex.lab=1.2)

MDpairN1 <- MDpairN[(1:n)[MDglobalTF]]
MDpairN0 <- MDpairN[(1:n)[!MDglobalTF]]
for (i in 1:sum(MDglobalTF)){
  points(rep(i,length(MDpairN1[[idx1[i]]])),MDpairN1[[idx1[i]]],cex=cex,pch=pch,col=col,...)
}
for (i in 1:sum(!MDglobalTF)){
  points(rep(i+sum(MDglobalTF),length(MDpairN0[[idx0[i]]])),MDpairN0[[idx0[i]]],cex=cex,pch=pch,col=col,...)
}
abline(v=sum(MDglobalTF)+0.5)         

# provide indices referring to original data with quantile ordering
ind1 <- (1:n)[MDglobalTF]
indices.regular <- ind1[idx1]
ind0 <- (1:n)[!MDglobalTF]
indices.outliers <- ind0[idx0]

#indices.regular=order(chibound)[MDglobalTF]
#indices.outliers=order(chibound)[!MDglobalTF]

list(resdist=resdist,idx1=idx1,idx0=idx0,mduse=MDpairN,mdusual=MDglobal,mdusual01=MDglobalTF,
     chibound=chibound,indices.regular=indices.regular,indices.outliers=indices.outliers)
}
#################################################################################################

xy <- cbind(X,Y)

par(mfrow=c(1,2))


if (!is.null(dist)){
  res <- neighbdist(dat,X,Y,k=k,dist=distc,propneighb=0.1,chisqqu=0.975)
}
else {
  res <- neighbdist(dat,X,Y,k=k,propneighb=0.1,chisqqu=0.975)
}

mduse <- res$mduse
idx1 <- res$idx1
idx0 <- res$idx0
resdist <- res$resdist
mdusual <- res$mdusual
mdusual01 <- res$mdusual01

if (is.null(sel)){
        options(locatorBell=FALSE)
        cat("Enter the points defining the polygon of the subsamples:\n")
        sel <- locator(type="l")
}
polygon(sel$x,sel$y,border=1)

n <- nrow(dat)
inpol <- vector("list",n)

mduse1 <- mduse[(1:n)[mdusual01]]
mduse0 <- mduse[(1:n)[!mdusual01]]
for (i in 1:sum(mdusual01)){
  l <- length(unlist(mduse1[[idx1[i]]]))
  inpol[[idx1[i]]] <- in.polygon(rep(i,l),mduse1[[idx1[i]]],sel$x,sel$y)
  if (sum(inpol[[idx1[i]]])>0){
    points(rep(i,sum(inpol[[idx1[i]]])),mduse1[[idx1[i]]][inpol[[idx1[i]]]],cex=0.4,pch=3,col=1)
  }
}

for (i in 1:sum(!mdusual01)){
  l <- length(unlist(mduse0[[idx0[i]]]))
  inpol[[idx0[i]+sum(mdusual01)]] <- in.polygon(rep(i+sum(mdusual01),l),mduse0[[idx0[i]]],sel$x,sel$y)
  if (sum(inpol[[idx0[i]+sum(mdusual01)]])>0){
    points(rep(i+sum(mdusual01),sum(inpol[[idx0[i]+sum(mdusual01)]])),mduse0[[idx0[i]]][inpol[[idx0[i]+sum(mdusual01)]]],cex=0.4,pch=3,col=1)
  }
}


plot(X,Y,col=gray(0.6),pch=1,cex=0.3,xlab="X coordinate",ylab="Y coordinate",cex.lab=1.2)

ind <- (1:n)[mdusual01]
for (i in 1:sum(mdusual01)){
  if (sum(inpol[[idx1[i]]])>0){
    points(X[ind[idx1[i]]],Y[ind[idx1[i]]],col=1,pch=3,cex=0.4)
    hlp0 <- resdist[[ind[idx1[i]]]][inpol[[idx1[i]]]]
    if (length(hlp0)==1){hlp1 <- as.matrix(t(xy[hlp0,]))} else {hlp1 <- xy[hlp0,]}
    hlp2 <- matrix(unlist(rep(xy[ind[idx1[i]],],length(hlp0))),nrow=length(hlp0),byrow=TRUE)
    segments(hlp1[,1],hlp1[,2],hlp2[,1],hlp2[,2],col=gray(0.2),lwd=0.3)
  }
}

ind <- (1:n)[!mdusual01]
for (i in 1:sum(!mdusual01)){
  if (sum(inpol[[idx0[i]+sum(mdusual01)]])>0){
    points(X[ind[idx0[i]]],Y[ind[idx0[i]]],col=1,pch=3,cex=0.4)
    hlp0 <- resdist[[ind[idx0[i]]]][inpol[[idx0[i]+sum(mdusual01)]]]
    if (length(hlp0)==1){hlp1 <- as.matrix(t(xy[hlp0,]))} else {hlp1 <- xy[hlp0,]}
    hlp2 <- matrix(unlist(rep(xy[ind[idx0[i]],],length(hlp0))),nrow=length(hlp0),byrow=TRUE)
    segments(hlp1[,1],hlp1[,2],hlp2[,1],hlp2[,2],col=gray(0.2),lwd=0.3)
  }
}


for (i in 1:sum(!mdusual01)){
  if (sum(inpol[[idx0[i]+sum(mdusual01)]])>0){
    points(X[ind[idx0[i]]],Y[ind[idx0[i]]],col=1,pch=3,cex=0.4)
    hlp0 <- resdist[[ind[idx0[i]]]][inpol[[idx0[i]+sum(mdusual01)]]]
    if (length(hlp0)==1){hlp1 <- as.matrix(t(xy[hlp0,]))} else {hlp1 <- xy[hlp0,]}
    hlp2 <- matrix(unlist(rep(xy[ind[idx0[i]],],length(hlp0))),nrow=length(hlp0),byrow=TRUE)
    segments(hlp1[,1],hlp1[,2],hlp2[,1],hlp2[,2],col=gray(0.2),lwd=0.3)
  }
}

invisible()
list(sel=sel,index.regular=res$indices.regular,index.outliers=res$indices.outliers)
}

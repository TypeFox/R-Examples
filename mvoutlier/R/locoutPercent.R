locoutPercent <-
function(dat,X,Y,dist=NULL,k=10,chisqqu=0.975,sortup=10,sortlow=90,nlinesup=10,nlineslow=10,
   indices=NULL,xlab="(Sorted) Index",ylab="Distance to neighbor",col=gray(0.7), ...)
{

n=nrow(dat)
p=ncol(dat)
covr <- covMcd(dat)
cinv <- solve(covr$cov) # inverse of the robust covariance matrix
MDglobal <- sqrt(mahalanobis(dat, covr$center, cinv, inverted=TRUE))

# TRUE/FALSE for non-outlying/outlying:
qchi <- sqrt(qchisq(chisqqu,p))
MDglobalTF <- (MDglobal<qchi)

if (!is.null(indices)){
  indices.reg <- indices[MDglobalTF[indices]]
  indices.out <- indices[!MDglobalTF[indices]]
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

# indices of neighbors (in list):
if (!is.null(dist)){
  Tp1 <- matrix(rep(t(X),length(X)),ncol=dim(t(X))[2],byrow=TRUE)
  Tp2 <- matrix(rep(t(Y),length(Y)),ncol=dim(t(Y))[2],byrow=TRUE)
  H <- sqrt((c(rep(X,length(X)))-Tp1)^2+(c(rep(Y,length(Y)))-Tp2)^2)
  diag(H) <- NA
  resdist <- apply(H <= dist,1,which)
  usexlab <- c("Percentage of next neighbors")
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
  usexlab <- paste("Percentage of next",k,"neighbors")
}

# guarantee that there is at least 1 neighbor to each data point:
if (any(lapply(resdist,length)==0)) 
  {stop("No neighbor to at least 1 observation -> choose larger distance!")}

# create list with pairwise MDs only for neighbors
MDpairN <- vector("list", n)

# boundary that should include required proportion of neighbors:
propneighb <- seq(from=0.01,to=1,by=0.01)
chibound <- matrix(NA,nrow=n,ncol=length(propneighb))
dimnames(chibound) <- list(c(1:n),paste("p",propneighb,sep=""))
for (j in 1:length(propneighb)){
  for (i in 1:n){
    MDpairN[[i]] <- MDpair[i,resdist[[i]]]
    nn <- max(1,floor(length(MDpairN[[i]])*propneighb[j])) # number of neighbors in tolerance ellipse
    nval <- MDpairN[[i]][order(MDpairN[[i]])][nn]  # value of largest neighbor to be included
    chibound[i,j] <- pchisq(nval^2,p,MDglobal[i]^2)
  }
}

###############################################################################################
#Plots:

par(mfrow=c(1,2))
# plot only regular observations:
plot(0,0,xlim=c(0,100),ylim=c(min(chibound),max(chibound)),type="n",
   xlab=usexlab, ylab="Degree of isolation", ...)
title("Regular observations")

# plot only curves for regular observations:
idxreg <- (1:n)[MDglobalTF]
for (i in 1:length(idxreg)){
  lines(propneighb*100,chibound[idxreg[i],],col=gray(0.7), ...)
}

if (!is.null(indices)){ # indices were provided to plot
  if (length(indices.reg)>0){ # some of the indices are from regular observations, so plot them:
    for (i in 1:length(indices.reg)){
      lines(propneighb*100,chibound[indices.reg[i],],col=1,lwd=1.2, ...)
    }
  }
}
else { # plot lines according to desired number 
  oUpT <- order(chibound[MDglobalTF,sortup])
  selreg.up <- oUpT[c((length(oUpT)-nlinesup+1):length(oUpT))]
  idxreg.up <- idxreg[selreg.up]
  for (i in 1:length(selreg.up)){
    lines(propneighb*100,chibound[idxreg.up[i],],col=1,lwd=1.2, ...)
  }

  oLowT <- order(chibound[MDglobalTF,sortlow])
  selreg.low <- oLowT[1:nlineslow]
  idxreg.low <- idxreg[selreg.low]
  for (i in 1:length(selreg.low)){
    lines(propneighb*100,chibound[idxreg.low[i],],col=1,lwd=1.2,lty=2, ...)
  }
}

####################################################################################
# plot only outlying observations:
plot(0,0,xlim=c(0,100),ylim=c(min(chibound),max(chibound)),type="n",
   xlab=usexlab, ylab="Degree of isolation", ... )
title("Outlying observations")

# plot only curves for outlying observations:
idxout <- (1:n)[!MDglobalTF]
for (i in 1:length(idxout)){
  lines(propneighb*100,chibound[idxout[i],],col=gray(0.7), ...)
}

if (!is.null(indices)){ # indices were provided to plot
  if (length(indices.out)>0){ # some of the indices are from outlying observations, so plot them:
    for (i in 1:length(indices.out)){
      lines(propneighb*100,chibound[indices.out[i],],col=1,lwd=1.2, ...)
    }
  }
}
else { # plot lines according to desired number
  oUpF <- order(chibound[!MDglobalTF,sortup])
  selout.up <- oUpF[c((length(oUpF)-nlinesup+1):length(oUpF))]
  idxout.up <- idxout[selout.up]
  for (i in 1:length(selout.up)){
    lines(propneighb*100,chibound[idxout.up[i],],col=1,lwd=1.2, ...)
  }

  oLowF <- order(chibound[!MDglobalTF,sortlow])
  selout.low <- oLowF[1:nlineslow]
  idxout.low <- idxout[selout.low]
  for (i in 1:length(selout.low)){
    lines(propneighb*100,chibound[idxout.low[i],],col=1,lwd=1.2,lty=2, ...)
  }
}


if (!is.null(indices)){ # indices were provided to plot
  ret <- list(indices.reg=indices.reg,indices.out=indices.out)
}
else {
  ret <- list(idxreg.up=idxreg.up,idxreg.low=idxreg.low,idxout.up=idxout.up,idxout.low=idxout.low)
}
return(ret)
}

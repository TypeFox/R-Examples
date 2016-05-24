plotsom <- function(obj,grp,type=c("num","bar"),margins=c(3,2,2,2),...)
{
# plot SOM output either as barplot (type="bar") or 
# as table (type="num")
# grp is the vector of group membership

if (length(grp)!=nrow(obj$visual))
  stop("Vector grp has not same length as number of objects!")

xdim <- obj$xdim
ydim <- obj$ydim
ngrp <- length(unique(grp))  # number of groups
if (ngrp>9) stop("Too many groups for text presentation!")
gr <- factor(grp,labels=1:ngrp)
x <- obj$visual$x+1
y <- obj$visual$y+1
out <- array(0,dim=c(xdim,ydim,ngrp),
   dimnames=list(paste("x",1:xdim,sep=""),paste("y",1:ydim,sep=""),
                 levels(as.factor(grp))))
for (i in 1:length(grp)){
  out[x[i],y[i],gr[i]] <- out[x[i],y[i],gr[i]]+1
}

if (type=="num"){
  layout(1)
  par(mar=c(1,1,1,1))
  plot(0,0,xlim=c(0,xdim),ylim=c(0,ydim),type="n",xlab="",ylab="",axes=FALSE,...)
  for (i in 0:ydim){segments(0,i,xdim,i)}
  for (j in 0:xdim){segments(j,0,j,ydim)}
  for (i in 1:xdim){
    for (j in 1:ydim){
      for (k in 1:ngrp){
        if (ngrp==2) text(i-1+k*0.33,j-1+0.5,out[i,j,k])
        if (ngrp==3) text(i-1+k*0.25,j-1+0.5,out[i,j,k])
        if (ngrp==4) {xp=c(1,2,1,2);yp=c(1,1,2,2)
                     text(i-1+xp[k]*0.33,j-0.33*yp[k],out[i,j,k])}
        if (ngrp==5) {xp=c(1,2,3,1,2);yp=c(1,1,1,2,2)
                     text(i-1+xp[k]*0.25,j-0.33*yp[k],out[i,j,k])}
        if (ngrp==6) {xp=c(1,2,3,1,2,3);yp=c(1,1,1,2,2,2)
                     text(i-1+xp[k]*0.25,j-0.33*yp[k],out[i,j,k])}
        if (ngrp==7) {xp=c(1,2,3,1,2,3,1);yp=c(1,1,1,2,2,2,3)
                     text(i-1+xp[k]*0.25,j-0.25*yp[k],out[i,j,k])}
        if (ngrp==8) {xp=c(1,2,3,1,2,3,1,2);yp=c(1,1,1,2,2,2,3,3)
                     text(i-1+xp[k]*0.25,j-0.25*yp[k],out[i,j,k])}
        if (ngrp==9) {xp=c(1,2,3,1,2,3,1,2,3);yp=c(1,1,1,2,2,2,3,3,3)
                     text(i-1+xp[k]*0.25,j-0.25*yp[k],out[i,j,k])}
      }
    }
  }
}
if (type=="bar"){
  par(mar=margins)
  fieldno <- matrix(1:(xdim*ydim),nrow=ydim,ncol=xdim,byrow=TRUE)
  fieldno <- fieldno[ydim:1,]
  layout(fieldno)
  out1 <- out/(max(out))
  for (j in 1:ydim){
    for (i in 1:xdim){
      barplot(out1[i,j,],axes=FALSE,ylim=c(0,1),...)
    }
  }
}

# summary table
sumtab <- matrix(0,nrow=xdim*ydim,ncol=ngrp+4,
   dimnames=list(1:(xdim*ydim),c("cell","x","y",levels(as.factor(grp)),"sum")))
sumtab[,1] <- 1:(xdim*ydim)
sumtab[,2] <- rep(seq(1:xdim),ydim)
sumtab[,3] <- sort(rep(seq(1:ydim),xdim))
for (i in 1:xdim){
  for (j in 1:ydim){
    sumtab[(j-1)*xdim+i,4:(3+ngrp)] <- out[i,j,]
  }
}
sumtab[,ngrp+4] <- apply(sumtab[,4:(3+ngrp)],1,sum)

return(sumtab)
}


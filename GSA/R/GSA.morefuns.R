quantitative.func  <-

function(x,y,s0=0){
# regression of x on y

my=mean(y)
  yy <- y-my
  temp <- x%*%yy
mx=rowMeans(x)
syy= sum(yy^2)

  scor <- temp/syy
  b0hat <- mx-scor*my
  xhat <- matrix(b0hat,nrow=nrow(x),ncol=ncol(x))+y*matrix(scor,nrow=nrow(x),ncol=ncol(x))
  sigma <- sqrt(rowSums((x-xhat)^2)/(ncol(xhat)-2))
  sd <- sigma/sqrt(syy)
  tt <- scor/(sd+s0)

  return(list(tt=tt, numer=scor, sd=sd))

}
tCorr.func  <-
function(x,y,s0=0){

#simple correlation (Fisher trans)

corr=cor(t(x),y)
scor=.5*log((1+corr)/(1-corr))
sd=rep(1,length(scor))
tt=scor
  return(list(tt=tt, numer=scor, sd=sd))
}

taCorr.func  <-
function(x,y,s0=0){

#simple abs correlation  (Fisher trans)

corr=abs(cor(t(x),y))
scor=.5*log((1+corr)/(1-corr))
sd=rep(1,length(scor))
tt=scor
  return(list(tt=tt, numer=scor, sd=sd))
}





ttest.func <- function(x,y,s0=0, sd=NULL){

  n1 <- sum(y==1)
  n2 <- sum(y==2)
  
  p <- nrow(x)
  m1 <- rowMeans(x[,y==1,drop=F])
  m2 <- rowMeans(x[,y==2,drop=F])


if(is.null(sd)){
  sd <- sqrt( ((n2-1) * varr(x[, y==2], meanx=m2) + (n1-1) * varr(x[, y==1], meanx=m1) )*(1/n1+1/n2)/(n1+n2-2) )
} 

  numer <-  m2 - m1

  dif.obs <- (numer)/(sd + s0)
  return(list(tt=dif.obs,numer=numer, sd=sd))
}


varr <- function(x, meanx=NULL){ 
  n <- ncol(x)
  p <- nrow(x)
  Y <-matrix(1,nrow=n,ncol=1)
  if(is.null(meanx)){   meanx <- rowMeans(x)}
  ans<- rep(1, p)
  xdif <- x - meanx %*% t(Y)
  ans <- (xdif^2) %*% rep(1/(n - 1), n)
  ans <- drop(ans)
  return(ans)

} 



paired.ttest.func <- function(x,y,s0=0, sd=NULL){

  nc <- ncol(x)/2
  o <- 1:nc
  o1 <- rep(0,ncol(x)/2);o2 <- o1
  for(j in 1:nc){o1[j] <- (1:ncol(x))[y==-o[j]]}
  for(j in 1:nc){o2[j] <- (1:ncol(x))[y==o[j]]}
  d <- x[,o2,drop=F]-x[,o1,drop=F]
  su <- x[,o2,drop=F]+x[,o1,drop=F]
  if(is.matrix(d)){
    m <-  rowMeans(d)
  }
  if(!is.matrix(d)) {m <- mean(d)}
  if(is.null(sd)){
    if(is.matrix(d)){ sd <- sqrt(varr(d, meanx=m)/nc)}
    if(!is.matrix(d)){sd <- sqrt(var(d)/nc)}
  }


  dif.obs <- m/(sd+s0)
  return(list(tt=dif.obs, numer=m, sd=sd))


}


cox.func <- function(x,y,censoring.status,s0=0){
  scor <- coxscor(x,y, censoring.status)$scor
  sd <- sqrt(coxvar(x,y, censoring.status))
  tt <- scor/(sd+s0)
  return(list(tt=tt, numer=scor, sd=sd))

}

coxscor <- function(x, y, ic, offset = rep(0., length(y))) {

  ## computes cox scor function for rows of nx by n matrix  x

  ## first put everything in time order

  n <- length(y)
  nx <- nrow(x)
  yy <- y + (ic == 0.) * (1e-05)
  otag <- order(yy)
  y <- y[otag]
  ic <- ic[otag]
  x <- x[, otag, drop = F]

  ##compute  unique failure times, d=# of deaths at each failure time, 
  ##dd= expanded version of d to length n, s=sum of covariates at each
  ## failure time, nn=#obs in each risk set, nno=sum(exp(offset)) at each failure time

  offset <- offset[otag]

  a <- coxstuff(x, y, ic, offset = offset)
  nf <- a$nf
  fail.times <- a$fail.times
  s <- a$s
  d <- a$d
  dd <- a$dd
  nn <- a$nn
  nno <- a$nno
  w <- rep(0., nx)
  for(i in (1.:nf)) {
    w <- w + s[, i]
    oo<- (1.:n)[y >= fail.times[i]]
    r<-rowSums(x[, oo, drop = F] * exp(offset[oo]))
    w<- w - (d[i]/nno[i])*r
  }

  return(list(scor = w, coxstuff.obj = a))

}

coxvar <- function(x, y, ic, offset = rep(0., length(y)), coxstuff.obj = NULL){

  ## computes information elements (var) for cox
  ## x is nx by n matrix of expression  values

  nx <- nrow(x)
  n <- length(y)
  yy <- y + (ic == 0.) * (1e-06)
  otag <- order(yy)
  y <- y[otag]
  ic <- ic[otag]
  x <- x[, otag, drop = F]
  offset <- offset[otag]
  if(is.null(coxstuff.obj)) {
    coxstuff.obj <- coxstuff(x, y, ic, offset = offset)
  }

  nf <- coxstuff.obj$nf
  fail.times <- coxstuff.obj$fail.times
  s <- coxstuff.obj$s
  d <- coxstuff.obj$d
  dd <- coxstuff.obj$dd
  nn <- coxstuff.obj$nn
  nno <- coxstuff.obj$nno

  x2<- x^2
  oo <- (1.:n)[y >= fail.times[1] ]
  sx<-(1/nno[1])*rowSums(x[, oo] * exp(offset[oo]))
  s<-(1/nno[1])*rowSums(x2[, oo] * exp(offset[oo]))
  w <-  d[1] * (s - sx * sx)


  for(i in 2.:nf) {
    oo <- (1.:n)[y >= fail.times[i-1] & y < fail.times[i] ]
    sx<-(1/nno[i])*(nno[i-1]*sx-rowSums(x[, oo,drop=F] * exp(offset[oo])))
    s<-(1/nno[i])*(nno[i-1]*s-rowSums(x2[, oo,drop=F] * exp(offset[oo])))
    w <- w + d[i] * (s - sx * sx)
  }
  return(w)

}



coxstuff<- function(x, y, ic, offset = rep(0., length(y))) {

  fail.times <- unique(y[ic == 1.])
  nf <- length(fail.times)
  n <- length(y)
  nn <- rep(0., nf)
  nno <- rep(0., nf)
  for(i in 1.:nf) {
    nn[i] <- sum(y >= fail.times[i])
    nno[i] <- sum(exp(offset)[y >= fail.times[i]])
  }

  s <- matrix(0., ncol = nf, nrow = nrow(x))
  d <- rep(0., nf)

  ##expand d out to a vector of length n

  for(i in 1.:nf) {
    o <- (1.:n)[(y == fail.times[i]) & (ic == 1.)]
    d[i] <- length(o)
  }

  oo <- match(y, fail.times)
  oo[ic==0]<-NA
  oo[is.na(oo)]<- max(oo[!is.na(oo)])+1
  s<-t(rowsum(t(x),oo))
  if(ncol(s)> nf){s<-s[,-ncol(s)]}
  dd <- rep(0., n)

  for(j in 1.:nf) {
    dd[(y == fail.times[j]) & (ic == 1.)] <- d[j]
  }

  return(list(fail.times=fail.times, s=s, d=d, dd=dd, nf=nf, nn=nn, nno=nno))

}


multiclass.func <- function(x,y,s0=0){

  ##assumes y is coded 1,2...

  nn <- table(y)
  m <- matrix(0,nrow=nrow(x),ncol=length(nn))
  v <- m
  for(j in 1:length(nn)){
    m[,j] <- rowMeans(x[,y==j])
    v[,j] <- (nn[j]-1)*varr(x[,y==j], meanx=m[,j])
  }
  mbar <- rowMeans(x)
  mm <- m-matrix(mbar,nrow=length(mbar),ncol=length(nn))
  fac <- (sum(nn)/prod(nn))
  scor <- sqrt(fac*(apply(matrix(nn,nrow=nrow(m),ncol=ncol(m),byrow=TRUE)*mm*mm,1,sum)))

  sd <- sqrt(rowSums(v)*(1/sum(nn-1))*sum(1/nn))
  tt <- scor/(sd+s0)
  mm.stand=t(scale(t(mm),center=FALSE,scale=sd))
  return(list(tt=tt, numer=scor, sd=sd,stand.contrasts=mm.stand))

}

est.s0<-function(tt,sd,s0.perc=seq(0,1, by=.05)){


  ## estimate s0 (exchangeability) factor for denominator.
  ## returns the actual estimate s0 (not a percentile)

br=unique(quantile(sd,seq(0,1,len=101)))
nbr=length(br)

  a<-cut(sd,br,labels=F)
  a[is.na(a)]<-1
  cv.sd<-rep(0,length(s0.perc))
  
  for(j in 1:length(s0.perc)){
    w<-quantile(sd,s0.perc[j])
    w[j==1]<-0
    tt2<-tt*sd/(sd+w)
    tt2[tt2==Inf]=NA
    sds<-rep(0,nbr-1) 

    for(i in 1:(nbr-1)){
      sds[i]<-mad(tt2[a==i], na.rm=TRUE)
    }

    cv.sd[j]<-sqrt(var(sds))/mean(sds)
  }

  o=(1:length(s0.perc))[cv.sd==min(cv.sd)]

# we don;t allow taking s0.hat to be  0th percentile when min sd is 0
  s0.hat=quantile(sd[sd!=0],s0.perc[o])

  return(list(s0.perc=s0.perc,cv.sd=cv.sd, s0.hat= s0.hat))

}
paired.perm=function(y){
n=max(abs(y))
nn=2*n
res=1:nn
for(i in 1:n){
 o=(1:nn)[abs(y)==i]
 u=runif(1)
 if(u>.5){
    temp=res[o[1]]
    res[o[1]]=res[o[2]]
    res[o[2]]=temp
}}
return(res)
}

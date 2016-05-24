# 06.07.07 : nothing to revise because this is Rob's code



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



quantitative.func  <- function(x,y,s0=0){

    # regression of x on y

    my=mean(y)
      yy <- y-my
      temp <- x%*%yy
    mx=rowMeans(x)
    syy= sum(yy^2)

      scor <- temp/syy
      b0hat <- mx-scor*my
    ym=matrix(y,nrow=nrow(x),ncol=ncol(x),byrow=T)
      xhat <- matrix(b0hat,nrow=nrow(x),ncol=ncol(x))+ym*matrix(scor,nrow=nrow(x),ncol=ncol(x))
      sigma <- sqrt(rowSums((x-xhat)^2)/(ncol(xhat)-2))
      sd <- sigma/sqrt(syy)
      tt <- scor/(sd+s0)

      return(list(tt=tt, numer=scor, sd=sd))

  }


  



cox.func <- function(x,y,censoring.status,s0=0){
  scor <- coxscor(x,y, censoring.status)$scor
  sd <- sqrt(coxvar(x,y, censoring.status))
  tt <- scor/(sd+s0)
  return(list(tt=tt, numer=scor, sd=sd))

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



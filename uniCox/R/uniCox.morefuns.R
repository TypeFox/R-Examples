initcx=function(y,ic){
n=length(y)
iriskq<-rep(0,n)
ddq<-rep(0,n)
tq<-rep(0,n)
kq=0

#dyn.load("/home/tibs/PAPERS/coxdiag/S.so")

junk<-.Fortran("initcx",
as.integer(n),
as.single(y),
kq=as.integer(kq),
as.integer(ic),
iriskq=as.integer(iriskq),
ddq=as.integer(ddq),
tq=as.single(tq),
PACKAGE="uniCox")

kq<-junk$kq;
tq<-junk$tq[1:kq];ddq<-junk$ddq[1:kq];iriskq<-junk$iriskq[1:kq]
return(list(tq=tq,ddq=ddq,iriskq=iriskq,kq=kq))
}

calcsx=function(x,kq,iriskq,ic){
n=length(x)
sx=rep(0,n)
junk2=.Fortran("calcsx",
as.integer(n),
as.single(x),
as.integer(kq),
as.integer(iriskq),
as.integer(ic),
sx=as.single(sx),
PACKAGE="uniCox")
return(junk2$sx[1:kq])
}



findroot=function(beta1,beta2,lam,x,sx,iriskq,kq,ddq,ic,del=.001){
SMALL=1e-6
f1=func.scor(beta1,lam,x,sx,iriskq,kq,ddq,ic)
f0n=func.scor(-SMALL,lam,x,sx,iriskq,kq,ddq,ic)
f0p=func.scor(SMALL,lam,x,sx,iriskq,kq,ddq,ic)
f2=func.scor(beta2,lam,x,sx,iriskq,kq,ddq,ic)
if(sign(f1)*sign(f0n)<0){ val=binary.search(beta1,0,lam,x,sx,iriskq,kq,ddq,ic,del=del)$beta}
if(sign(f2)*sign(f0p)<0){ val=binary.search(0,beta2,lam,x,sx,iriskq,kq,ddq,ic,del=del)$beta}
if( (sign(f1)*sign(f0n)>0) & (sign(f2)*sign(f0p)>0) ){val=0}
return(val)
}

func.scor=function(beta,lam,x,sx,iriskq,kq,ddq,ic){
n=length(x)
junk3=.Fortran("calcsc",
as.single(beta),
as.single(lam),
as.integer(n),
as.single(x),
as.single(sx),
as.integer(iriskq),
as.integer(kq),
as.integer(ddq),
as.integer(ic),
val=single(1),
scrt1=single(n),
scrt2=single(n),
PACKAGE="uniCox")

return(junk3$val)
}


binary.search=function(a,b,lam,x,sx,iriskq,kq,ddq,ic,del=.001){
fa=func.scor(a,lam,x,sx,iriskq,kq,ddq,ic)
fb=func.scor(b,lam,x,sx,iriskq,kq,ddq,ic)
while(abs(a-b)>del){
cc=(a+b)/2
fc=func.scor(cc,lam,x,sx,iriskq,kq,ddq,ic)
if(sign(fa*fc)>0) {a=cc;fa=fc}
if(sign(fb*fc)>0) {b=cc;fb=fc}
# this stops the serach
if(fc==0){a=b}
}
return(list(betahat=cc,scor=fc))
}

comp.path=function(lamlist,betalims,x,sx,iriskq,kq,ddq,ic){
val=rep(0,length(lamlist))
val[1]=findroot(betalims[1],betalims[2],lamlist[1],x,sx,iriskq,kq,ddq,ic)
betahat=val[1]
go=T
ii=1
while(go & (ii<length(lamlist))){
go=F
ii=ii+1
lam=lamlist[ii]
if(betahat>0){beta1=0;beta2=betahat}
if(betahat<0){beta2=0;beta1=betahat}
val[ii]=findroot(beta1,beta2,lam,x,sx,iriskq,kq,ddq,ic)
if(val[ii]!=0){go=T}
#betahat=val[ii]
}
return(val)
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


timeorder=function(y,icensq,x){
#x is p by n
 n <- length(y)
  nx <- nrow(x)
  yy <- y + (icensq == 0.) * (1e-05)
  otag <- order(yy)
  y <- y[otag]
  icensq <- icensq[otag]
  x <- x[, otag, drop = F]
return(list(x=x,y=y,icensq=icensq))
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


coxscor2 <- function(x, y, ic, offset = rep(0., length(y))) {
# version that allows offset to be a matrix
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

  offset <- offset[,otag,drop=F]

  a <- coxstuff2(x, y, ic, offset = offset)
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
    r<-rowSums(x[, oo, drop = F] * exp(offset[,oo,drop=F]))
    w<- w - (d[i]/nno[,i])*r
  }

  return(list(scor = w, coxstuff.obj = a))

}
coxstuff2<- function(x, y, ic, offset = rep(0., length(y))) {

# version that allows offset to be a matrix
  p=nrow(x)
  fail.times <- unique(y[ic == 1.])
  nf <- length(fail.times)
  n <- length(y)
  nn <- rep(0., nf)
  nno <- matrix(0.,nrow=p,ncol= nf)
  for(i in 1.:nf) {
    nn[i] <- sum(y >= fail.times[i])
    nno[,i] <- rowSums(exp(offset)[,y >= fail.times[i],drop=F])
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

if(!is.matrix(s)){s=matrix(s,nrow=1)}
  return(list(fail.times=fail.times, s=s, d=d, dd=dd, nf=nf, nn=nn, nno=nno))

}




coxvar2 <- function(x, y, ic, offset = rep(0., length(y)), coxstuff.obj = NULL){

# version that allows offset to be a matrix

  ## computes information elements (var) for cox
  ## x is nx by n matrix of expression  values

  nx <- nrow(x)
  n <- length(y)
  yy <- y + (ic == 0.) * (1e-06)
  otag <- order(yy)
  y <- y[otag]
  ic <- ic[otag]
  x <- x[, otag, drop = F]
  offset <- offset[,otag,drop=F]
  if(is.null(coxstuff.obj)) {
    coxstuff.obj <- coxstuff2(x, y, ic, offset = offset)
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
  sx<-(1/nno[,1])*my.rowSums(x[, oo] * exp(offset[,oo]))
  s<-(1/nno[,1])*my.rowSums(x2[, oo] * exp(offset[,oo]))
 
  w <-  d[1] * (s - sx * sx)


  for(i in 2.:nf) {
    oo <- (1.:n)[y >= fail.times[i-1] & y < fail.times[i] ]
    sx<-(1/nno[,i])*(nno[,i-1]*sx-my.rowSums(x[, oo,drop=F] * exp(offset[,oo])))
    s<-(1/nno[,i])*(nno[,i-1]*s-my.rowSums(x2[, oo,drop=F] * exp(offset[,oo])))
    w <- w + d[i] * (s - sx * sx)
  }
  return(w)

}


my.rowSums=function(x){
  if(is.matrix(x)){val=rowSums(x)}
  else{val=sum(x)}
return(val)
}





 balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
   totals <- table(y)
   fmax <- max(totals)
   nfolds <- min(nfolds, fmax)     
   nfolds= max(nfolds, 2)
         # makes no sense to have more folds than the max class size
   folds <- as.list(seq(nfolds))
   yids <- split(seq(y), y) 
         # nice we to get the ids in a list, split by class
###Make a big matrix, with enough rows to get in all the folds per class
   bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
   for(i in seq(totals)) {
     if(length(yids[[i]])>1){bigmat[seq(totals[i]), i] <- sample(yids[[i]])}
     if(length(yids[[i]])==1){bigmat[seq(totals[i]), i] <- yids[[i]]}

   }
   smallmat <- matrix(bigmat, nrow = nfolds)# reshape the matrix
### Now do a clever sort to mix up the NAs
   smallmat <- permute.rows(t(smallmat))   ### Now a clever unlisting
         # the "clever" unlist doesn't work when there are no NAs
         #       apply(smallmat, 2, function(x)
         #        x[!is.na(x)])
   res <-vector("list", nfolds)
   for(j in 1:nfolds) {
     jj <- !is.na(smallmat[, j])
     res[[j]] <- smallmat[jj, j]
   }
   return(res)
 }

permute.rows <-function(x)
{
        dd <- dim(x)
        n <- dd[1]
        p <- dd[2]
        mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
        matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}


cor.func<- 

function (x, y, s0.perc ) 
{
    n <- length(y)
    xbar <- x %*% rep(1/n, n)
    sxx <- ((x - as.vector(xbar))^2) %*% rep(1, n)
    sxy <- (x - as.vector(xbar)) %*% (y - mean(y))
    syy <- sum((y - mean(y))^2)
    numer <- sxy/sxx
    sd <- sqrt((syy/sxx - numer^2)/(n - 2))

if(is.null(s0.perc)){ fudge=median(sd)}
if(!is.null(s0.perc)){
 if(s0.perc>=0){
  fudge=quantile(sd,s0.perc)
  }
 if(s0.perc<0){
   fudge=0
  }
}
    tt <- numer/(sd + fudge)

    return(list(tt = tt, numer = numer, sd = sd, fudge=fudge ))
}

coxfunc <- 
function(x, y, censoring.status, s0.perc)
{

        junk <- coxscor(x, y, censoring.status)
        scor<-junk$scor
        sd <- sqrt(coxvar(x, y, censoring.status, coxstuff.obj=junk$coxstuff.obj))

if(is.null(s0.perc)){ fudge=median(sd)}
if(!is.null(s0.perc)){
if(s0.perc>=0){
  fudge=quantile(sd,s0.perc)
  }
 if(s0.perc<0){
   fudge=0
  }


}

        tt <- scor/(sd + fudge)

        return(list(tt = tt, numer = scor, sd = sd, fudge=fudge ))
}


coxscor <- 
function(x, y, ic, offset = rep(0., length(y)))
{
        # computes cox scor function for rows of nx by n matrix  x
 # first put everything in time order
        n <- length(y)
        nx <- nrow(x)
        yy <- y + (ic == 0.) * (1e-05)
        otag <- order(yy)
        y <- y[otag]
        ic <- ic[otag]
        x <- x[, otag, drop = F]
        #compute  unique failure times, d=# of deaths at each failure time, 
        #dd= expanded version of d to length n, s=sum of covariates at each
        # failure time, nn=#obs in each risk set, nno=sum(exp(offset)) at each failure time
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



 coxvar <- 
function(x, y, ic, offset = rep(0., length(y)), coxstuff.obj = NULL)
{
        # computes information elements (var) for cox
        # x is nx by n matrix of expression  values
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




coxstuff<-
function(x, y, ic, offset = rep(0., length(y)))
{
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
        #expand d out to a vector of length n
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


 ocoxvar <-
function(x, y, ic, offset = rep(0., length(y)), coxstuff.obj = NULL)
{
        # computes information elements (var) for cox
        # x is nx by n matrix of expression  values
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
        w <- rep(0., nx)
 x2<- x^2
        for(i in 1.:nf) {
                oo <- (1.:n)[y >= fail.times[i]]
      sx<-(1/nno[i])*rowSums(x[, oo] * exp(offset[oo]))
         s<-(1/nno[i])*rowSums(x2[, oo] * exp(offset[oo]))
       w <- w + d[i] * (s - sx * sx)
        }
        return(w)
}


mysvd<-function(x,  n.components=NULL){
# finds PCs of matrix x
  p<-nrow(x)
  n<-ncol(x)

# center the observations (rows)

 feature.means<-rowMeans(x)
x<- t(scale(t(x),center=feature.means,scale=F))


  if(is.null(n.components)){n.components=min(n,p)}
  if(p>n){
    a<-eigen(t(x)%*%x)
    v<-a$vec[,1:n.components,drop=FALSE]
    d<-sqrt(a$val[1: n.components,drop=FALSE])
    
      u<-scale(x%*%v,center=FALSE,scale=d)
 
    
    return(list(u=u,d=d,v=v,  feature.means=feature.means))
  }
  else{

      junk<-svd(x,LINPACK=TRUE)
      nc=min(ncol(junk$u), n.components)
      return(list(u=junk$u[,1:nc],d=junk$d[1:nc],
                  v=junk$v[,1:nc], feature.means=feature.means))
}
}

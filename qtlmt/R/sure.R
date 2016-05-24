
#############################
#  functions for SURE model #
################################################################################

#############################
# SURE: add, drop, stepwise #
#############################
sureMove<- function(object,y,x,range=NULL,iter=10000,tol=1e-12){
# object: object from sureEst() or alike
# y: n by p matrix
# x: n by m matrix
# range: list, which x's for y_j
  yes<- TRUE
  ooo<- object
  vvv<- ooo$v
  p<- ncol(y)
  if(is.null(range)||missing(range)){
    range<- vector("list",p)
    for(i in 1:p) range[[i]]<- 1:ncol(x)
  }
  while(yes){
    yes<- FALSE
    vv0<- vvv
    oo0<- ooo
    for(j in 1:p){
      yes1<- TRUE
      while(yes1){
        yes1<- FALSE
        v0<- vv0[[j]]; vl<- length(v0)
        if(vl>0){
          for(jj in 1:vl){
            vv1<- vv0
            rng<- setdiff(range[[j]],vv0[[j]])
            if(length(rng)>0)for(ii in rng){
              vv1[[j]][jj]<- ii
              oo1<- sureEst(y,x,vv1,oo0$sigma,iter=iter,tol=tol)
              if(oo1$loglik>oo0$loglik+1e-8){
                oo0<- oo1
                vv0<- vv1
                yes1<- TRUE
              }
            }
          }
        }
      }
    }
    if(oo0$loglik>ooo$loglik+1e-8){
      vvv<- vv0
      ooo<- oo0
      yes<- TRUE
    }
  }
  for(j in 1:p) vvv[[j]]<- sort(vvv[[j]])

  list(loglik=ooo$loglik,b=ooo$b,sigma=ooo$sigma,v=vvv)
}

#############################
# SURE: estimate b and sigma #
##############################

sureEst<-
   function(y,
            x,
            v,
            sigma,
            iter=10000,
            tol=1e-12)
{
   UseMethod("sureEst")
}

sureEst.default<- function(y, x, v, sigma, iter=10000, tol=1e-12){
# y: n by p matrix
# x: n by m matrix
# v: list of x's to start with
# sigma: can be missing
  if(!is.matrix(y)) stop("\a  y should be a matrix...")
  if(!is.matrix(x)) stop("\a  x should be a matrix...")
  if(!is.list(v)) stop("\a  v should be a variable list...")
  n<- nrow(y)
  p<- ncol(y)
  m<- ncol(x)
  if(nrow(x)!=n) stop("\a  x should have the same number of rows as y...")
  if(length(v)!=p) stop("\a  length of v could be wrong...")
  
  nqs<- NULL
  qs<- NULL
  for(i in 1:p){
    nqs<- c(nqs,length(v[[i]]))
    qs<- c(qs,v[[i]])
  }
  b<- rep(0.0, sum(nqs)+p)
  if(missing(sigma))sigma<- cov(y)*(n-1)/n

  loglik<- 1e-308  
  vv<- .C("sureEstc",
    as.double(t(y)),
    as.integer(n),
    as.integer(p),
    as.double(t(x)),
    as.integer(m),
    as.integer(nqs),
    as.integer(qs),
    b=as.double(b),
    sigma=as.double(sigma),
    loglik=as.double(loglik),
    as.integer(1),
    as.integer(iter),
    as.double(tol))
  oo<- list(loglik=vv$loglik,b=vv$b,sigma=matrix(vv$sigma,p,p),v=v)
  fit<- array(NA, dim=dim(y))
  cnt<- 0
  for(j in 1:length(oo$v)){
    cnt<- cnt+1
    fit[,j]<- oo$b[cnt]
    k<- length(oo$v[[j]])
    if(k>1){
      fit[,j]<- fit[,j] + x[,oo$v[[j]]]%*%oo$b[cnt+1:k]
    }else if(k>0){
      fit[,j]<- fit[,j] + x[,oo$v[[j]]]*oo$b[cnt+1:k]
    }
    cnt<- cnt+k
  }
  oo$fitted.values<- fit
  class(oo)<- "sure"
  oo
}

sureAdd1<-
   function(object,
            y,
            x,
            range=NULL,
            iter=10000,
            tol=1e-12,
            ext=FALSE)
{
   UseMethod("sureAdd1")
}

sureAdd1.default<- function(object,y,x,range=NULL,iter=10000,tol=1e-12,ext=FALSE){
  if(!is.matrix(y)) stop("\a  y should be a matrix...")
  if(!is.matrix(x)) stop("\a  x should be a matrix...")
  n<- nrow(y)
  p<- ncol(y)
  m<- ncol(x)
  if(nrow(x)!=n) stop("\a  x should have the same number of rows as y...")

  if(is.null(range)||missing(range)){
    range<- vector("list",p)
    for(i in 1:p) range[[i]]<- 1:ncol(x)
  }

  if(ext) o<- sureMove(object,y,x,range=range,iter=iter,tol=tol) else o<- object

  o0<- o
  add<- FALSE
  for(j in 1:p){
    v0<- setdiff(range[[j]],o$v[[j]]); vl<- length(v0)
    if(vl>0)for(vt in v0){
      vv1<- o$v
      vv1[[j]]<- sort(c(o$v[[j]],vt))
      o1<- sureEst(y,x,v=vv1,sigma=o0$sigma,iter=iter,tol=tol)
      if(ext) o1<- sureMove(o1,y,x,range=range,iter=iter,tol=tol)

      if(o1$loglik>o0$loglik+1e-8){
        o0<- o1
        add<- TRUE
      }
    }
  }
  
  oo<- list(loglik=o0$loglik,b=o0$b,sigma=o0$sigma,v=o0$v,add=add)
  class(oo)<- "sure"
  oo
}

sureDrop1<-
   function(object,
            y,
            x,
            range=NULL,
            iter=10000,
            tol=1e-12,
            ext=FALSE)
{
   UseMethod("sureDrop1")
}

sureDrop1.default<- function(object,y,x,range=NULL,
  iter=10000,tol=1e-12,ext=FALSE){
  if(!is.matrix(y)) stop("\a  y should be a matrix...")
  if(!is.matrix(x)) stop("\a  x should be a matrix...")
  n<- nrow(y)
  p<- ncol(y)
  m<- ncol(x)
  if(nrow(x)!=n) stop("\a  x should have the same number of rows as y...")

  if(is.null(range)||missing(range)){
    range<- vector("list",p)
    for(i in 1:p) range[[i]]<- 1:ncol(x)
  }

  if(ext) o<- sureMove(object,y,x,range=range,iter=iter,tol=tol) else o<- object

  o0<- o; lik<- -Inf
  drop<- FALSE
  for(j in 1:p){
    v0<- o$v[[j]]; vl<- length(v0)
    if(vl>0)for(vt in v0){
      vv1<- o$v
      vv1[[j]]<- sort(setdiff(o$v[[j]],vt))
      o1<- sureEst(y,x,v=vv1,sigma=o0$sigma,iter=iter,tol=tol)
      if(ext) o1<- sureMove(o1,y,x,range=range,iter=iter,tol=tol)

      if(o1$loglik>lik+1e-8){
        o0<- o1
        lik<- o1$loglik
        drop<- TRUE
      }
    }
  }
  
  oo<- list(loglik=o0$loglik,b=o0$b,sigma=o0$sigma,v=o0$v,drop=drop)
  class(oo)<- "sure"
  oo
}

#####################
### SURE stepwise ###
#####################
sureStep<- 
   function(object,
            y,
            x,
            cv,
            direction=c("both", "backward", "forward"),
            range=NULL,
            iter=10000,
            steps=1000,
            tol=1e-12,
            ext=FALSE)
{
   UseMethod("sureStep")
}

sureStep.default<- function(object, y, x, cv, 
  direction=c("both", "backward", "forward"), range=NULL,
  iter=10000, steps=1000, tol=1e-12, ext=FALSE){
# object: sure object to start with
# y: n by p matrix
# x: n by m matrix
# cv: threshold for -2*log(H0/H1)
# range: list, which x's for y_j
  direction<- match.arg(direction)
  if(!is.matrix(y)) stop("\a  y should be a matrix...")
  if(!is.matrix(x)) stop("\a  x should be a matrix...")
  n<- nrow(y)
  p<- ncol(y)
  m<- ncol(x)
  if(nrow(x)!=n) stop("\a  x should have the same number of rows as y...")
  
  if(is.null(range)||missing(range)){
    range<- vector("list",p)
    for(i in 1:p) range[[i]]<- 1:ncol(x)
  }

  if(ext) o<- sureMove(object,y=y,x=x,range=range,iter=iter,tol=tol) else o<- object

  if(direction=="both"){
    yes<- TRUE
    if(missing(cv)){
      stop("\a\n  cv: should be positive but missing...\n\n")
    }
    while(yes){
      yes<- FALSE
      od<- sureDrop1(object=o,y=y,x=x,range=range,iter=iter,tol=tol,ext=ext)
      if(od$drop && 2*(o$loglik-od$loglik)<cv){
        o<- od
        yes<- TRUE
      }else{
        oa<- sureAdd1(object=o,y=y,x=x,range=range,iter=iter,tol=tol,ext=ext)
        if(oa$add && 2*(oa$loglik-o$loglik)>cv){
          o<- oa
          yes<- TRUE
        }else{
          od<- sureDrop1(object=oa,y=y,x=x,range=range,iter=iter,tol=tol,ext=ext)
          if(od$drop && od$loglik>o$loglik+1e-8){
             o<- od
             yes<- TRUE
          }
        }
      }
    }
  }else if(direction=="backward"){
    if(missing(cv)){
      cv<- Inf
    }
    yes<- TRUE
    while(yes){
      yes<- FALSE
      od<- sureDrop1(object=o,y=y,x=x,range=range,iter=iter,tol=tol,ext=ext)
      if(od$drop && 2*(o$loglik-od$loglik)<cv){
        o<- od
        yes<- TRUE
      }
    }
  }else if(direction=="forward"){
    if(missing(cv)){
      cv<- 0
    }
    yes<- TRUE
    while(yes){
      yes<- FALSE
      oa<- sureAdd1(object=o,y=y,x=x,range=range,iter=iter,tol=tol,ext=ext)
      if(oa$add && 2*(oa$loglik-o$loglik)>cv){
        o<- oa
        yes<- TRUE
      }
    }
  }

  oo<- list(loglik=o$loglik,b=o$b,sigma=o$sigma,v=o$v)
  class(oo)<- "sure"
  oo
}

#******************************
# stepwise variable selection #
#******************************
# y: n by p matrix
# x: n by m matrix
# v: list of x's to start with (and end up with)
# lower: list of x's of lower scope
# upper: list of x's of upper scope
# k: penalty, 0 if missing or <0
# record(+/-y_j,x_j,loglikelihood)
#-------------------------------------------------
surStep<-
   function(y,
            x,
            v,
            lower,
            upper,
            k,
            direction=c("both", "backward", "forward"),
            iter=10000,
            max.terms=200,
            steps=1000,
            tol=1e-12)
{
   UseMethod("surStep")
}

surStep.default<- function(y, x, v, lower, upper, k,
  direction=c("both", "backward", "forward"),
  iter=10000, max.terms=200, steps=1000, tol=1e-12){
  direction<- match.arg(direction)
  if(direction=="backward") direction<- 0
    else if(direction=="forward") direction<- 1
      else if(direction=="both") direction<- 2
        else stop("stepSure: wrong direction...")
  if(!is.matrix(y)) stop("y should be a matrix...")
  if(!is.matrix(x)) stop("x should be a matrix...")
  if(!is.list(v)) stop("v should be a variable list...")
  if(!missing(lower)){
    if(!is.list(lower)) stop("lower should be a variable list...")
  }
  if(!missing(upper)){
    if(!is.list(upper)) stop("upper should be a variable list...")
  }
  n<- nrow(y)
  p<- ncol(y)
  m<- ncol(x)
  if(nrow(x)!=n) stop("x should have the same number of rows as y...")
  if(length(v)!=p) stop("length of v could be wrong...")
  if(!missing(lower)){
    if(length(lower)!=p) stop("length of lower could be wrong...")
  }
  if(!missing(upper)){
    if(length(upper)!=p) stop("length of upper could be wrong...")
  }
  
  if(missing(k)){
    k<- 0
  }else if(k>1e+308){
    k<- 1e+308
  }else if(k<0){
    k<- 0
  }
    
  vin<- matrix(0,nrow=p,ncol=m)
  nlw<- NULL
  lw<- NULL
  nup<- NULL
  up<- NULL
  for(i in 1:p){
    if(missing(lower)){
      nlw<- c(nlw,0)
      lw<- c(lw,NULL)
    }else{
      nlw<- c(nlw,length(lower[[i]]))
      lw<- c(lw,lower[[i]])
    }
    if(missing(upper)){
      nup<- c(nup,m)
      up<- c(up,1:m)
    }else{
      nup<- c(nup,length(upper[[i]]))
      up<- c(up,upper[[i]])
    }
    vin[i,v[[i]]]<- 1
  }
  rec<- rep(0.0, 3*(steps+1))
 
  vv<- .C("sureStepc",
    as.double(t(y)),
    as.integer(n),
    as.integer(p),
    as.double(t(x)),
    as.integer(m),
    as.integer(nlw),
    as.integer(lw),
    as.integer(nup),
    as.integer(up),
    as.double(k),
    as.integer(direction),
    vin=as.integer(t(vin)),
    rec=as.double(rec),
    as.integer(max.terms),
    as.integer(steps),
    as.integer(iter),
    as.double(tol))
  vin<- matrix(vv$vin,nrow=p,ncol=m,byrow=T)
  record<- matrix(vv$rec,ncol=3,byrow=T)
  mx<- c(1:(steps+1))[record[,1]==9999]
  if(length(mx)>0){
    mx<- mx[1]-1
    record<- matrix(record[1:mx,],nrow=mx)
    colnames(record)<- c("y","x","loglik")
  }
  v<- list()
  for(i in 1:p){
    v[[i]]<- c(1:m)[as.logical(vin[i,])]
  }
  fit<- sureEst(y, x, v, iter=iter, tol=tol)

  oo<- list(loglik=fit$loglik,b=fit$b,sigma=fit$sigma,v=v,trace=record)
  class(oo)<- "sure"
  oo
}

#######################
### sure: epistasis ###
#######################
# y: n by p matrix
# x: n by m matrix
# k: penalty, 0 if missing or <0
#-------------------------------------------------
sureEps<-
   function(y,
            x,
            v,
            k,
            direction=c("both", "backward", "forward"),
            iter=10000,
            max.terms=200,
            steps=1000,
            tol=1e-12)
{
   UseMethod("sureEps")
}

sureEps.default<- function(y, x, v, k, direction=c("both", "backward", "forward"),
  iter=10000, max.terms=200, steps=1000, tol=1e-12){
  direction<- match.arg(direction)
  if(!is.matrix(y)) stop("y should be a matrix...")
  if(!is.matrix(x)) stop("x should be a matrix...")
  if(!is.list(v)) stop("v should be a variable list...")
  n<- nrow(y)
  p<- ncol(y)
  m<- ncol(x)
  lower<- v
  upper<- list()
  indx<- m
  tt<- NULL
  for(jj in 1:p){
    nad<- 0
    v0<- v[[jj]]
    nx<- length(v0)
    if(nx>1){
      for(i in 1:(nx-1)){
        for(j in (i+1):nx){
          nad<- nad+1
          tt<- rbind(tt,c(jj,v0[i],v0[j]))
          x<- cbind(x,x[,v0[i]]*x[,v0[j]])
        }
      }
      upper[[jj]]<- c(v0,(indx+1):(indx+nad))
      indx<- indx+nad
    }else upper[[jj]]<- v[[jj]]
  }
  xad<- (m+1):indx
  if(direction=="backward") v<- upper
  gv<- surStep(y,x,v,lower=lower,upper=upper,k=k, direction=direction,
  iter=iter, max.terms=max.terms, steps=steps, tol=tol)
  ge<- sureEst(y, x, v=gv$v, iter=iter, tol=tol)
  nb<- length(ge$b)
  xx<- matrix(0.0,nrow=n*p,ncol=nb)
  bb<- NULL
  indx<- 0
  for(jj in 1:p){
    bb<- c(bb,c(0,gv$v[[jj]]))
    xx0<- cbind(1,x[,gv$v[[jj]]])
    xx[((jj-1)*n+1):(jj*n),(indx+1):(indx+length(gv$v[[jj]])+1)]<- xx0
    indx<- indx+length(gv$v[[jj]])+1
  }
  ss<- kronecker(solve(ge$sigma),diag(1,nrow=n,ncol=n))
  ss<- t(xx)%*%ss%*%xx
  ss<- solve(ss)
  dd<- 1:nb
  dd<- dd[bb>m]
  rst<- NULL
  for(jj in dd){
    ff<- ge$b[jj]^2/ss[jj,jj]
    p0<- 1-pf(ff,1,n*p-nb)
    jj0<- xad==bb[jj]
    rst<- rbind(rst,c(tt[jj0,],p0))
  }
  colnames(rst)<- c("trait","marker1","marker2","pvalue")
  data.frame(rst)
}

################################################################################
# the end #
###########


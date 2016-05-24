# earth.leaps.R:
#
# Copied from Thomas Lumley's leaps 2.9 package for earth 3.2-6 to avoid
# use of leaps::: in the earth code, to prevent complaints from CRAN check.


# leaps.setup is modified to handle linear dependencies in x properly.  I think
# this fix is needed only if leaps.setup is called with intercept=FALSE.
#
# The fix is needed because if there are linear dependencies in the matrix
# x passed to the original leaps.setup, it gives an incorrect error
# message: "missing value where TRUE/FALSE needed".
#
# In the earth context, this happened if you call earth.update with new
# data and the bx generated from that data has linear dependencies (which
# is actually ok).
#
# I also touched up the warning messages to be more informative.

leaps.setup<-function(x,y,wt=rep(1,length(y)),force.in=NULL,
                      force.out=NULL,intercept=TRUE,
                      nvmax=8,nbest=1,warn.dep=TRUE){
  make.names<-function(np){
    if (np<27) letters[1:np] else as.character(1:np)
  }
  np<-NCOL(x)
  nn<-NROW(x)
  if (length(y)!=nn) stop("y and x different lengths")
  if (length(wt)!=nn) stop("wt and x different lengths")
  if (is.null(colnames(x))) colnames(x)<-make.names(np)
  index<-rep(0,np)
  names(index)<-colnames(x)
  index[force.in]<--1
  if (any(index[force.out]==-1)) stop("Can't force the same variable in and out")
  index[force.out]<-1
  force.in<-(index==-1) ## make force.in, force.out logical vectors
  force.out<-(index==1)
  ii<-order(index)
  xx<-x[,ii]
  force.in<-force.in[ii]
  force.out<-force.out[ii]
  ones<-rep(1,np)
  names(ones)<-colnames(x)
  first<-1+sum(ones[force.in])
  last<-np-sum(ones[force.out])
  nvmax<-min(nvmax,np)
  if (intercept){
    np<-np+1
    xnames<-c("(Intercept)",colnames(xx))
    xx<-cbind(1,xx)
    colnames(xx)<-xnames
    first<-first+1
    last<-last+1
    nvmax<-nvmax+1
    index<-c(-1,index)
  }
  vorder<-1:np
  il<-nvmax*(nvmax+1)/2
  nrbar<-np*(np-1)/2
  qrleaps<-.Fortran("makeqr",np=as.integer(np),nn=as.integer(nn),
                    wt=as.double(wt),tx=t(xx),y=as.double(y),d=numeric(np),
                    rbar=numeric(nrbar),
                    thetab=numeric(np),sserr=numeric(1),ier=as.integer(0),
                    PACKAGE="earth")
  if (qrleaps$ier!=0)
      warning(paste("MAKEQR returned error code",qrleaps$ier))
  qrleaps$tx<-NULL
  qrleaps$wt<-NULL
  tolset<-.Fortran("tolset",as.integer(np),
                   as.integer(nrbar),qrleaps$d,qrleaps$rbar,
                   tol=numeric(np),numeric(np),ier=as.integer(0), PACKAGE="earth")
  if (tolset$ier!=0)
      warning(paste("TOLSET returned error code",tolset$ier))
  ss<-.Fortran("ssleaps",as.integer(np),qrleaps$d,
               qrleaps$thetab,qrleaps$sserr,rss=numeric(np),
               ier=as.integer(0), PACKAGE="earth")
  if (ss$ier!=0)
      warning(paste("SS returned error code",ss$ier))
  sing<-.Fortran("sing",np=as.integer(qrleaps$np),nrbar=as.integer(nrbar),
                d=qrleaps$d,rbar=qrleaps$rbar,thetab=qrleaps$thetab,
                sserr=qrleaps$sserr,tol=tolset$tol,lindep=logical(qrleaps$np),
                work=numeric(qrleaps$np),ier=as.integer(0), PACKAGE="earth")
  if (sing$ier>0)
       warning(paste("SING returned error code",sing$ier))
  sing$work<-NULL
  if(any(sing$lindep)) { # linear dependencies in x?
      if (intercept) {
          new.force.out <- sing$lindep | c(FALSE,force.out)
          reordered.col.nbrs <- order(new.force.out[-1]) # put lin dep cols at end
          try.again <- any((c(new.force.out,1) - c(0,new.force.out)) < 0) # huh?
          lindep.in.force.in <- any(sing$lindep[-1] & force.in)
          colnames.with.intercept <- c("(Intercept)", colnames(x))
      } else {
          new.force.out <- sing$lindep | force.out
          reordered.col.nbrs <- order(new.force.out)
          try.again <- any((c(new.force.out,1) - c(0,new.force.out)) < 0)
          lindep.in.force.in <- any(sing$lindep & force.in)
          colnames.with.intercept <- colnames(x)
      }
      if (warn.dep)
        warning0(if(try.again) "Trying again because " else "",
                 sum(sing$lindep),
                 " linearly dependent variable",
                 if(sum(sing$lindep) > 1) "s" else "",
                 ": ",
                 paste(colnames.with.intercept[sing$lindep], collapse=", "))
      if (lindep.in.force.in)
          stop("Linear dependency in force.in variable(s)")
      if (try.again) {
          rval<-leaps.setup(x[,ii[reordered.col.nbrs],drop=FALSE], y, wt,
                            force.in[reordered.col.nbrs], force.out[reordered.col.nbrs],
                            intercept, nvmax, nbest, warn.dep=FALSE)
          rval$reorder<-ii[reordered.col.nbrs]
          return(rval)
      }
      lastsafe<-max((1:np)[!new.force.out])
      if (lastsafe<min(nvmax,last)) {
          if (warn.dep)
              warning0("nvmax reduced from ", nvmax, " to ",
                       lastsafe-intercept, " because of linearly dependencies")
          nvmax<-lastsafe
      }
  }
  if (any(sing$lindep)){
      ss<-.Fortran("ssleaps",as.integer(np),sing$d,sing$thetab,
                   sing$sserr,rss=numeric(np),ier=as.integer(0),PACKAGE="earth")
    if (ss$ier!=0)
            warning(paste("SS returned error code",ss$ier))
  }
  initr<-.Fortran("initr",as.integer(np),as.integer(nvmax),as.integer(nbest),
                  bound=numeric(np),ress=numeric(nbest*nvmax),as.integer(nvmax),
                  lopt=integer(nbest*il),as.integer(il),vorder=as.integer(vorder),
                  ss$rss,ier=as.integer(0), PACKAGE="earth")
  if (initr$ier!=0)
      warning(paste("INITR returned error code",initr$ier))
  nullrss<-if (intercept) ss$rss[1] else sum(y^2)
  rval<-c(sing,list(nn=qrleaps$nn,rss=ss$rss,bound=initr$bound,
                    ress=matrix(initr$ress,ncol=nbest),
                    lopt=matrix(initr$lopt,ncol=nbest),
                    nvmax=nvmax,nbest=nbest,nrbar=nrbar,il=il,
                    ir=nvmax,vorder=initr$vorder,
                    first=first,last=last,xnames=colnames(xx),
                    force.in=(index==-1),force.out=(index==1),
                    intercept=intercept,nullrss=nullrss))
  class(rval)<-"regsubsets"
  invisible(rval)
}

leaps.exhaustive<-function(leaps.obj,really.big=FALSE){
    if (!inherits(leaps.obj,"regsubsets")){
        stop("Not a regsubsets object -- must run leaps.setup")
    }
    nbest<-leaps.obj$nbest
    if (!really.big & (leaps.obj$np>50 || leaps.obj$nbest>40)) {
        stop("Exhaustive search will be S L O W, must specify really.big=T")
    }
    dimwk<-3*leaps.obj$last
    dimiwk<-leaps.obj$nvmax
    rval<-.Fortran("xhaust",
                   np=as.integer(leaps.obj$np),
                   nrbar=as.integer(leaps.obj$nrbar),
                   d=leaps.obj$d,
                   rbar=leaps.obj$rbar,
                   thetab=leaps.obj$thetab,
                   first=as.integer(leaps.obj$first),
                   last=as.integer(leaps.obj$last),
                   vorder=as.integer(leaps.obj$vorder),
                   tol=leaps.obj$tol,
                   rss=leaps.obj$rss,
                   bound=leaps.obj$bound,
                   nvmax=as.integer(leaps.obj$nvmax),
                   ress=leaps.obj$ress,
                   ir=as.integer(leaps.obj$ir),
                   nbest=as.integer(leaps.obj$nbest),
                   lopt=matrix(as.integer(leaps.obj$lopt),ncol=nbest),
                   il=as.integer(leaps.obj$il),
                   wk=numeric(dimwk),
                   dimwk=as.integer(dimwk),
                   iwk=integer(dimiwk),
                   dimiwk=as.integer(dimiwk),
                   ier=as.integer(0), PACKAGE="earth")
    rval$dimwk<-rval$dimiwk<-rval$iwk<-rval$wk<-NULL
    rval$xnames<-leaps.obj$xnames
    rval$method<-c("exhaustive",leaps.obj$method)
    rval$force.in<-leaps.obj$force.in
    rval$force.out<-leaps.obj$force.out
    rval$sserr<-leaps.obj$sserr
    rval$intercept<-leaps.obj$intercept
    rval$lindep<-leaps.obj$lindep
    rval$reorder<-leaps.obj$reorder
    rval$nullrss<-leaps.obj$nullrss
    rval$nn<-leaps.obj$nn
    class(rval)<-"regsubsets"
    if(rval$ier!=0) warning(paste("XHAUST returned error code",rval$ier))
    rval
}


leaps.backward<-function(leaps.obj){
  if (!inherits(leaps.obj,"regsubsets")){
      stop("Not a regsubsets object -- must run leaps.setup")
  }
  nbest<-leaps.obj$nbest
  dimwk<-2*leaps.obj$last
  rval<-.Fortran("bakwrd",np=as.integer(leaps.obj$np),
                 nrbar=as.integer(leaps.obj$nrbar),d=leaps.obj$d,
                 rbar=leaps.obj$rbar,thetab=leaps.obj$thetab,
                 first=as.integer(leaps.obj$first),
                 last=as.integer(leaps.obj$last),
                 vorder=as.integer(leaps.obj$vorder),tol=leaps.obj$tol,
                 rss=leaps.obj$rss,bound=leaps.obj$bound,
                 nvmax=as.integer(leaps.obj$nvmax),
                 ress=leaps.obj$ress,ir=as.integer(leaps.obj$ir),
                 nbest=as.integer(leaps.obj$nbest),
                 lopt=matrix(as.integer(leaps.obj$lopt),ncol=nbest),
                 il=as.integer(leaps.obj$il),wk=numeric(dimwk),
                 dimwk=as.integer(dimwk),ier=as.integer(0), PACKAGE="earth")
  rval$dimwk<-rval$wk<-NULL
  rval$xnames<-leaps.obj$xnames
  rval$method<-c("backward",leaps.obj$method)
  rval$force.in<-leaps.obj$force.in
  rval$force.out<-leaps.obj$force.out
  rval$sserr<-leaps.obj$sserr
  rval$intercept<-leaps.obj$intercept
  rval$lindep<-leaps.obj$lindep
  rval$reorder<-leaps.obj$reorder
  rval$nullrss<-leaps.obj$nullrss
  rval$nn<-leaps.obj$nn
  class(rval)<-"regsubsets"
  if(rval$ier!=0)
      warning(paste("BAKWRD returned error code",rval$ier))
  rval
}


leaps.forward<-function(leaps.obj){
  if (!inherits(leaps.obj,"regsubsets")){
      stop("Not a regsubsets object -- must run leaps.setup")
  }
  nbest<-leaps.obj$nbest
  dimwk<-3*leaps.obj$last
  rval<-.Fortran("forwrd",np=as.integer(leaps.obj$np),
                 nrbar=as.integer(leaps.obj$nrbar),
                 d=leaps.obj$d,rbar=leaps.obj$rbar,
                 thetab=leaps.obj$thetab,first=as.integer(leaps.obj$first),
                 last=as.integer(leaps.obj$last),
                 vorder=as.integer(leaps.obj$vorder),
                 tol=leaps.obj$tol,rss=leaps.obj$rss,
                 bound=leaps.obj$bound,nvmax=as.integer(leaps.obj$nvmax),
                 ress=leaps.obj$ress,ir=as.integer(leaps.obj$ir),
                 nbest=as.integer(leaps.obj$nbest),
                 lopt=matrix(as.integer(leaps.obj$lopt),ncol=nbest),
                 il=as.integer(leaps.obj$il),wk=numeric(dimwk),
                 dimwk=as.integer(dimwk),ier=as.integer(0), PACKAGE="earth")
  rval$dimwk<-rval$wk<-NULL
  rval$xnames<-leaps.obj$xnames
  rval$method<-c("forward",leaps.obj$method)
  rval$force.in<-leaps.obj$force.in
  rval$force.out<-leaps.obj$force.out
  rval$sserr<-leaps.obj$sserr
  rval$intercept<-leaps.obj$intercept
  rval$lindep<-leaps.obj$lindep
  rval$reorder<-leaps.obj$reorder
  rval$nullrss<-leaps.obj$nullrss
  rval$nn<-leaps.obj$nn
  class(rval)<-"regsubsets"
  if(rval$ier!=0)
      warning(paste("FORWARD returned error code",rval$ier))
  rval
}


leaps.seqrep<-function(leaps.obj){
    if (!inherits(leaps.obj,"regsubsets")){
        stop("Not a regsubsets object -- must run leaps.setup")
    }
    nbest<-leaps.obj$nbest
    dimwk<-3*leaps.obj$last
    rval<-.Fortran("seqrep",np=as.integer(leaps.obj$np),
                   nrbar=as.integer(leaps.obj$nrbar),
                   d=leaps.obj$d,rbar=leaps.obj$rbar,
                   thetab=leaps.obj$thetab,
                   first=as.integer(leaps.obj$first),
                   last=as.integer(leaps.obj$last),
                   vorder=as.integer(leaps.obj$vorder),
                   tol=leaps.obj$tol,rss=leaps.obj$rss,
                   bound=leaps.obj$bound,nvmax=as.integer(leaps.obj$nvmax),
                   ress=leaps.obj$ress,ir=as.integer(leaps.obj$ir),
                   nbest=as.integer(leaps.obj$nbest),
                   lopt=matrix(as.integer(leaps.obj$lopt),
                   ncol=nbest),il=as.integer(leaps.obj$il),
                   wk=numeric(dimwk),dimwk=as.integer(dimwk),
                   ier=as.integer(0), PACKAGE="earth")
  rval$dimwk<-rval$wk<-NULL
  rval$xnames<-leaps.obj$xnames
  rval$method<-c("'sequential replacement'",leaps.obj$method)
  rval$force.in<-leaps.obj$force.in
  rval$force.out<-leaps.obj$force.out
  rval$sserr<-leaps.obj$sserr
  rval$intercept<-leaps.obj$intercept
  rval$lindep<-leaps.obj$lindep
  rval$reorder<-leaps.obj$reorder
  rval$nullrss<-leaps.obj$nullrss
  rval$nn<-leaps.obj$nn
  class(rval)<-"regsubsets"
  if(rval$ier!=0)
      warning(paste("SEQREP returned error code",rval$ier))
  rval
}

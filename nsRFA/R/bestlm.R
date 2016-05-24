bestlm <- function (dip, ind, kmax=4, nbest=3) {
  # dip = vector n x 1
  # ind = matrix n x k
  nomi <- colnames(ind)
  result <- .leaps(ind, dip, method="adjr2", nbest=nbest, names=nomi)
  quali <- result$which[result$size <= (kmax+1),]
  adjr <- result$adjr[result$size <= (kmax+1)]
  indici <- order(adjr, decreasing=TRUE)
  adjr <- adjr[indici]
  quali <- quali[indici,]
  output <- NULL
  output$dip <- dip
  output$ind <- ind
  output$kmax <- kmax
  output$nbest <- nbest
  output$R2adj <- adjr
  output$subselect <- quali
  class(output) <- "bestlm"
  return(output)
}


print.bestlm <- function (x, ...) {
  N <- dim(x$subselect)[1] 
  regr <- rep(NA, length=N)
  for (i in 1:N) {
   regr[i] <- paste(colnames(x$subselect)[x$subselect[i,]], collapse=" + ")
  }
  output <- data.frame(regr, x$R2adj); names(output) <- c("model","R2adj")
  print(output)
}


summary.bestlm <- function (object, ...) {
  N <- dim(object$subselect)[1]
  regr <- rep(NA, length=N)
  diagn <- data.frame(matrix(NA, nrow=N, ncol=8)); 
   names(diagn) <- c("RMSE","RMSEP","MAE","MAEP",
                     "RMSEjk","RMSEjkP","MAEjk","MAEjkP")
  Y <- object$dip
  X <- object$ind
  for (i in 1:N) {
   regr[i] <- paste(colnames(object$subselect)[object$subselect[i,]], collapse=" + ")
   f <- formula(paste("Y ~", regr[i]))
   rr <- lm(f, X)
   fitt <- rr$fitted.values
   crossval <- jackknife1.lm(rr)
   diagn[i,1] <- RMSE(Y, fitt)
   diagn[i,2] <- RMSEP(Y, fitt)
   diagn[i,3] <- MAE(Y, fitt)
   diagn[i,4] <- MAEP(Y, fitt)
   diagn[i,5] <- RMSE(Y, crossval)
   diagn[i,6] <- RMSEP(Y, crossval)
   diagn[i,7] <- MAE(Y, crossval)
   diagn[i,8] <- MAEP(Y, crossval)
  }
  output <- data.frame(regr, object$R2adj); names(output) <- c("model","R2adj")
  output <- data.frame(cbind(output, diagn))
  print(output)
}



# ---------------------------------------------------------------------------------- #
# 
# R code for model selection using
# Alan Miller's FORTRAN routines
# (from R package leaps)
.leaps.setup<-function(x,y,wt=rep(1,length(y)),force.in=NULL,
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
        PACKAGE="nsRFA") #,DUP=FALSE)
  if (qrleaps$ier!=0)
      warning(paste("MAKEQR returned error code",qrleaps$ier))
  qrleaps$tx<-NULL
  qrleaps$wt<-NULL
  tolset<-.Fortran("tolset",as.integer(np),
                   as.integer(nrbar),qrleaps$d,qrleaps$rbar,
                   tol=numeric(np),numeric(np),ier=as.integer(0), PACKAGE="nsRFA")
  if (tolset$ier!=0)
      warning(paste("TOLSET returned error code",tolset$ier))
  ss<-.Fortran("ssleaps",as.integer(np),qrleaps$d,
               qrleaps$thetab,qrleaps$sserr,rss=numeric(np),
               ier=as.integer(0), PACKAGE="nsRFA")
  if (ss$ier!=0)
      warning(paste("SS returned error code",ss$ier))
 sing<-.Fortran("sing",np=as.integer(qrleaps$np),nrbar=as.integer(nrbar),
                d=qrleaps$d,rbar=qrleaps$rbar,thetab=qrleaps$thetab,
                sserr=qrleaps$sserr,tol=tolset$tol,lindep=logical(qrleaps$np),
                work=numeric(qrleaps$np),ier=as.integer(0), PACKAGE="nsRFA")
   if (sing$ier>0)
       warning(paste("SING returned error code",sing$ier))
  sing$work<-NULL
  if(any(sing$lindep)) {
      if (warn.dep)
          warning(paste(sum(sing$lindep)," linear dependencies found"))
      if (any(sing$lindep[-1] & force.in)) stop("Linear dependency among variables forced in")
      rightorder<-sing$lindep | c(FALSE,force.out)
      if (any((c(rightorder,1)-c(0,rightorder))<0)) {
          if (warn.dep){
              cat("Reordering variables and trying again:\n")
          }
          reorder<-if(intercept) order(rightorder[-1]) else order(rightorder)
          rval<-.leaps.setup(x[,ii[reorder],drop=FALSE],y,wt,
                            force.in[reorder],force.out[reorder],
                            intercept,nvmax,nbest,warn.dep=FALSE)
          rval$reorder<-ii[reorder]
          return(rval)
      }
      lastsafe<-max((1:np)[!rightorder])
      if (lastsafe<min(nvmax,last)) {
          nvmax<-lastsafe
          if (warn.dep)
              warning(paste("nvmax reduced to ",nvmax-intercept))
      }
  }
  if (any(sing$lindep)){
      ss<-.Fortran("ssleaps",as.integer(np),sing$d,sing$thetab,
                   sing$sserr,rss=numeric(np),ier=as.integer(0),PACKAGE="nsRFA")
  	if (ss$ier!=0)
            warning(paste("SS returned error code",ss$ier))
  }
  initr<-.Fortran("initr",as.integer(np),as.integer(nvmax),as.integer(nbest),
                  bound=numeric(np),ress=numeric(nbest*nvmax),as.integer(nvmax),
                  lopt=integer(nbest*il),as.integer(il),vorder=as.integer(vorder),
                  ss$rss,ier=as.integer(0), PACKAGE="nsRFA")
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


.leaps.exhaustive<-function(leaps.obj,really.big=FALSE){
    if (!inherits(leaps.obj,"regsubsets")){
        stop("Not a regsubsets object -- must run .leaps.setup")
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
                   ier=as.integer(0), PACKAGE="nsRFA")
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


.leaps.backward<-function(leaps.obj){
  if (!inherits(leaps.obj,"regsubsets")){
      stop("Not a regsubsets object -- must run .leaps.setup")
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
                 dimwk=as.integer(dimwk),ier=as.integer(0), PACKAGE="nsRFA")
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


.leaps.forward<-function(leaps.obj){
  if (!inherits(leaps.obj,"regsubsets")){
      stop("Not a regsubsets object -- must run .leaps.setup")
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
                 dimwk=as.integer(dimwk),ier=as.integer(0), PACKAGE="nsRFA")
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


.leaps.seqrep<-function(leaps.obj){
    if (!inherits(leaps.obj,"regsubsets")){
        stop("Not a regsubsets object -- must run .leaps.setup")
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
                   ier=as.integer(0), PACKAGE="nsRFA")
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



#   print.regsubsets<-function(x,...){
#       ll<-x #CMD check
#       cat("Subset selection object\n")
#       if (!is.null(ll$call)) {
#           cat("Call: ")
#           print(ll$call)
#       }
#       cat(ll$np-ll$intercept)
#       cat(" Variables ")
#       if (ll$intercept) cat(" (and intercept)")
#       cat("\n")
#       fmat<-cbind(ll$force.in,ll$force.out)
#       colnames(fmat)<-c("Forced in","Forced out")
#       rownames(fmat)<-ll$xnames
#       print(fmat[-1,])
#       cat(ll$nbest)
#       cat(" subsets of each size up to ")
#       cat(ll$nvmax-ll$intercept)
#       cat("\n")
#       cat("Selection Algorithm: ")
#       if (is.null(ll$method)) cat(" not done") else cat(ll$method)
#       cat("\n")
#       invisible(NULL)
#   }
#   
#   print.summary.regsubsets<-function(x,...){
#       print(x$obj)
#       print(x$outmat)
#   }
#   
.summary.regsubsets<-function(object,all.best=TRUE,matrix=TRUE,matrix.logical=FALSE,
                             df=NULL,...){
    ll<-object #CMD check
    triangle<-function(k) {j<-k-1;1+j*(j+1)/2}
    nmodl<-ll$nbest*ll$nvmax
    if(all.best) nshow<-ll$nbest else nshow<-1
    if (!is.null(df)) n1<-df else n1<-ll$nn-ll$intercept
    outmat<-NULL
    rmat<-NULL
    rnames<-NULL	
    outnames<-NULL
    rsqvec<-NULL
    cpvec<-NULL
    adjr2vec<-NULL
    bicvec<-NULL
    rssvec<-NULL
    sigma2<-ll$sserr/(n1+ll$intercept-ll$last)
    for (i in ll$first:min(ll$last,ll$nvmax)){
        if(!matrix) outmat<-NULL
        for(j in 1:nshow){
            if (ll$ress[i,j]>=1e35) next
            if ((j>1) &
                (all(ll$lopt[triangle(i):(triangle(i+1)-1),j-1]==
                     ll$lopt[triangle(i):(triangle(i+1)-1),j])))
                next
            rline<-rep(FALSE,ll$np)
            rline[ll$lopt[triangle(i):(triangle(i+1)-1),j]]<-TRUE
            outnames<-c(outnames,paste(i-ll$intercept," (",j,")"))
            rnames<-c(rnames,as.character(i-ll$intercept))
            rmat<-rbind(rmat,rline)
            vr<-ll$ress[i,j]/ll$nullrss
            rssvec<-c(rssvec,ll$ress[i,j])
            rsqvec<-c(rsqvec,1-vr)
            adjr2vec<-c(adjr2vec,1-vr*n1/(n1+ll$intercept-i))
            cpvec<-c(cpvec,ll$ress[i,j]/sigma2-(n1+ll$intercept-2*i))
            bicvec<-c(bicvec,(n1+ll$intercept)*log(vr)+i*log(n1+ll$intercept))
        }
    }
    rownames(rmat)<-rnames
    cn<-ll$xnames
    colnames(rmat)<-cn
    reorder<-if (is.null(ll$reorder)) 1:NCOL(rmat) else c(1,1+ll$reorder)
    rmat<-rmat[,order(reorder),drop=FALSE]
    if (matrix){
        if (!matrix.logical)
            outmat<-ifelse(rmat,"*"," ")
        else
            outmat<-rmat
        rownames(outmat)<-outnames
        if (ll$intercept) outmat<-outmat[,-1,drop=FALSE]
    }
    rval<-list(which=rmat,rsq=rsqvec,rss=rssvec,adjr2=adjr2vec,
               cp=cpvec,bic=bicvec,outmat=outmat,obj=ll)
    class(rval)<-"summary.regsubsets"
    rval
}
#     
#   
#   regsubsets<-function(x,...){
#     UseMethod("regsubsets",x)
#   }
#   
#   regsubsets.default<-function(x,y,weights=rep(1,length(y)),nbest=1,
#                                nvmax=8,force.in=NULL,force.out=NULL,
#                                intercept=TRUE,
#                                method=c("exhaustive","backward","forward","seqrep"),
#                                really.big=FALSE,...)
#   {
#       
#       a<-.leaps.setup(x,y,wt=weights,nbest=nbest,nvmax=nvmax,
#                      force.in=force.in,force.out=force.out,
#                      intercept=intercept)
#       switch(1+pmatch(method[1],
#                       c("exhaustive","backward","forward","seqrep"),
#                       nomatch=0),
#              stop(paste("Ambiguous or unrecognised method name :",method)),
#              .leaps.exhaustive(a,really.big=really.big),
#              .leaps.backward(a),
#              .leaps.forward(a),
#              .leaps.seqrep(a))
#   }
#   
#   
#   regsubsets.formula<-function(x,data,weights=NULL,nbest=1,nvmax=8,force.in=NULL,
#                                force.out=NULL,intercept=TRUE,
#                                method=c("exhaustive","backward","forward","seqrep"),
#                                really.big=FALSE,...){
#     formula<-x
#     rm(x)
#     mm<-match.call()
#     mm$formula<-formula
#     mm$x<-NULL
#     mm$nbest<-mm$nvmax<-mm$force.in<-mm$force.out<-NULL
#     mm$intercept<-mm$method<-mm$really.big<-NULL
#     mm[[1]]<-as.name("model.frame")
#     mm<-eval(mm,sys.frame(sys.parent()))
#     x<-model.matrix(terms(formula,data=data),mm)[,-1]
#     y<-model.extract(mm,"response")
#     wt<-model.extract(mm,"weights")
#     if (is.null(wt))
#         wt<-rep(1,length(y))
#     else
#         wt<-weights
#     a<-.leaps.setup(x,y,wt=wt,nbest=nbest,nvmax=nvmax, force.in=force.in,
#                    force.out=force.out, intercept=intercept)
#     rval<-switch(1+pmatch(method[1],
#                           c("exhaustive","backward","forward","seqrep"),
#                           nomatch=0),
#                  stop(paste("Ambiguous or unrecognised method name :",method)),
#                  .leaps.exhaustive(a,really.big),
#                  .leaps.backward(a),
#                  .leaps.forward(a),
#                  .leaps.seqrep(a))
#     rval$call<-sys.call(sys.parent())
#     rval
#   }

  
.leaps <- function(x, y, wt=rep(1,NROW(x)), int=TRUE, method=c("Cp","adjr2","r2"),
                   nbest=10, names=NULL, df=NROW(x), strictly.compatible=TRUE){
    if (!is.logical(int))
        stop("int should be TRUE or FALSE")
    if (!is.null(names))
        colnames(x)<-names
    method<-method[1]
    if (pmatch(method,c("Cp","adjr2","r2"),nomatch=0)==0)
        stop("Ambiguous or unrecognised method name")
    if (strictly.compatible){
        if (NCOL(x)>31)
            stop("leaps does not allow more than 31 variables; use regsubsets()")
        if (is.null(names))
            colnames(x)<-c(as.character(1:9),LETTERS)[1:NCOL(x)]
    }
    a<-.leaps.setup(x,y,wt=wt,nbest=nbest,nvmax=NCOL(x)+int,
                   intercept=int,warn.dep=FALSE)
    if (strictly.compatible & any(a$lindep))
        stop("leaps requires full-rank design matrix; use regsubsets()")
    b<-.leaps.exhaustive(a)
    d<-.summary.regsubsets(b)
    rval<-list(which=d$which)
    if (int)
        rval$which<-rval$which[,-1,drop=FALSE]
    rval$label<-colnames(d$which)
    rval$size<-as.numeric(rownames(d$which))+int
    if (pmatch(method,c("Cp"),nomatch=0)==1){
        rval$Cp<-d$cp
    }
    if (pmatch(method,c("r2"),nomatch=0)==1){
        rval$r2<-d$rsq
    }
    if (pmatch(method,c("adjr2"),nomatch=0)==1){
        rval$adjr2<-d$adjr2
    }
    rval
}






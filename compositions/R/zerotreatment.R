## I propose to move all zero-treatment stuff here (missingProjectors, zeroreplace, etc). Ok done.

MARvalue <- NaN
MNARvalue<- NA
SZvalue  <- -Inf
BDLvalue <- 0.0   # Use with care, negative values specifiy the detection limit

is.BDL <- function(x,mc=attr(x,"missingClassifier")) {if( is.null(mc) ) is.finite(x)&x<=0.0  else do.call(mc,list(x))=="BDL"  }
is.SZ  <- function(x,mc=attr(x,"missingClassifier")) {if( is.null(mc) ) is.infinite(x)&x<0  else do.call(mc,list(x))=="SZ" }
is.MAR <- function(x,mc=attr(x,"missingClassifier")) {if( is.null(mc) ) is.nan(x)  else do.call(mc,list(x))=="MAR" }
is.MNAR<- function(x,mc=attr(x,"missingClassifier")) {if( is.null(mc) )  is.na(x)&!is.nan(x)  else do.call(mc,list(x))=="MNAR" }
is.NMV <- function(x,mc=attr(x,"missingClassifier")) {if( is.null(mc) )  is.finite(x)&x>0.0 else do.call(mc,list(x))=="NMV" }
is.WZERO <- function(x,mc=attr(x,"missingClassifier")) {is.BDL(x)|is.SZ(x) }
is.WMNAR <- function(x,mc=attr(x,"missingClassifier")) {is.BDL(x)|is.MNAR(x) }

has.missings <- function(x,...) UseMethod("has.missings")

has.missings.default     <- function(x,mc=attr(x,"missingClassifier"),...) {
  (!is.null(x))&&!all(is.NMV(x,mc))
}

has.missings.rmult     <- function(x,mc=attr(x,"missingClassifier"),...) {
  if( is.null(x) )
    return( FALSE )
  if( !is.null(attr(x,"missingProjector")) )
    return( TRUE )
  if( !is.null(orig <- attr(x,"orig")) )
    return( has.missings(orig) )
  else
    !all(is.finite(x))
}
     
getDetectionlimit <- function(x,dl=attr(x,"detectionlimit")) {
  if( is.null(dl) ) dl <- NA
  ifelse(is.finite(x)&x<0,-x,dl)
}

missingType <- function(x,...,mc=attr(x,"missingClassifier"),values=c("NMV","BDT","MAR","MNAR","SZ","Err")) {
  #      finite, infinit, nan, na, x>0&!na  q=(x>0&!na)|nan
  # NA     -       -       -   +    -           -     
  # NaN    -       -       +   +    -           +
  # -Inf   -       +       -   -    -           -         
  # Inf    -       +       -   -    +           +
  # 0/-1   +       -       -   -    -           -
  # +1     +       -       -   -    +           +
  # 6 Faelle -> 3 Bit noetig
  #     infinit, na , !q
  # +1    0      0    0
  # 0/-1  0      0    1
  # NaN   0      1    0
  # NA    0      1    1
  # Inf   1      0    0
  # -Inf  1      0    1
  if( is.null(mc) || identical(mc,missingType)) {
    na <- is.na(x)
    values <- values[c(1,2,3,4,6,5)]
    structure(values[2+is.infinite(x)*4+na*2-((!na&x>0)|is.nan(x))],dim=dim(x),dimnames=dimnames(x))
  }
  else
    do.call(mc,list(x,...))
}



missingSummary <- function(x,...,vlabs=colnames(x),mc=attr(x,"missingClassifier"),values=eval(formals(missingType)$values)) {
  if( is.data.frame(x) )
     x <- data.matrix(x)
  missingType <- structure(c(missingType(x,...,values=1:6)),levels=as.character(values),
            class=c("ordered","factor"))
  # meaning, only faster:
  # missingType <- factor(missingType(x,values=values),levels=values)
  if( length(dim(x)==2 )) {
    if(length(vlabs)!=ncol(x))
      vlabs <- paste("V",1:ncol(x),sep="")
    variable <- structure(rep(1:ncol(x),each=nrow(x)),
                          levels=as.character(vlabs),class="factor")

    as.missingSummary(table(variable,missingType))
  } else {
    as.missingSummary(t(table(missingType)))
  }
}

as.missingSummary <- function(x,...) {
  class(x) <- c("missingSummary",class(x))
  x
}

plot.missingSummary <-function(x,...,main="Missings",legend.text=TRUE,
                    col=c("gray","lightgray","yellow","red","white","magenta")) {
  barplot(t(x),main=main,legend.text=legend.text,...,
          col=col)
}

simulateMissings <- function(x, dl=NULL, knownlimit=FALSE, 
                             MARprob=0.0, MNARprob=0.0, mnarity=0.5, SZprob=0.0) {
  if( is.data.frame(x) )
    x <- data.matrix(x)
  at <- attributes(x)
  x <- unclass(x)
  n <- length(x)
  SZ   <- runif(n)<SZprob
  MAR  <- runif(n)<MARprob
  normal <- function(x) {
    x<- rank(ifelse(is.finite(x),x,1));
    qnorm(x/(length(x)+1))
  }
  w1 <- sqrt(1-mnarity)
  w2 <- sqrt(mnarity)
  MNAR <- pnorm(rnorm(n)*w1+normal(x)*w2) < MNARprob
  x <- clo(x,detectionlimit=dl,total=NA)
  if( !knownlimit )
    x[is.BDL(x)]<-BDLvalue
  n <- length(x)
  x[MNAR]<-MNARvalue
  x[MAR]<- MARvalue
  x[SZ]<-SZvalue
  attributes(x) <- at
  x
}

zeroreplace <- function(x,d=NULL,a=2/3){
  # ensure objects are either equivalent matrices, data frames or vectors
  W = oneOrDataset(x)
  if(!is.null(d)){
    # if the detection limit is explicitly given, give it the same shape
    Losts = oneOrDataset(is.WZERO(as.matrix(W)),W) # find wide zeroes
    d = oneOrDataset(d,W)
  }else{
    # if not, extract from the data set
    Losts = oneOrDataset(is.BDL(as.matrix(W)),W) # find values below detection limit (encoded as negative)
    d = getDetectionlimit(W)
  }
  # scale down the detection limit...
  d = a*d
  # ... and replace
  W[Losts]=d[Losts]
  return(structure(W,Losts=Losts,class=class(x)))
}

#zeroreplace <- function(x,d,a=2/3){
#  # ensure all objects are either equivalent matrices, data frames or vectors
#  W = oneOrDataset(x)
#  Losts = oneOrDataset(is.WZERO(as.matrix(W)),W) # find zeroes
#  d = oneOrDataset(d,W)
#  # scale down the detection limit
#  d = a*d
#  # replace
#  W[Losts]=d[Losts]
#  attr(W,"Losts")=Losts
#  return(W)
#}



# to compute variance of a compositional data set with zero-lost values
# TO DO: generate the equivalent gsi.covwithlosts (for Y != X)
gsi.varwithlosts <- function(x,giveCenter=FALSE){
 # x contains the cpt-transformed data set
 # dat contains the original data set
 #require("tensorA")
 dat = attributes(x)$orig
 # this index is used to remove the cases x_i=x_j from the computation
 prodId = c(outer(1:nrow(dat),1:nrow(dat),"=="))

 # index to compare two samples:
    idx = expand.grid(is=1:nrow(dat), js=1:nrow(dat))
 # for each pair of data, get the set of variables observed twice:
  prodProj = is.NMV(dat[idx$is,])*is.NMV(dat[idx$js,])  # cup
   # and its associated missing projector $P_{M_i\cup M_j}$
   ls = missingProjector(dat,has=prodProj)
   # but fix it to zero if x_i=x_j:
   ls[,,prodId]=0
   #   or if they come from a fully non-observed pair:
     ls[is.na(ls)|is.nan(ls)] = 0
   # do the kronecker product of each missing projector with itself
   ls = sapply(1:nrow(idx), function(i){
     mat = ls[,,i]
     res = kronecker(mat,mat)
     return(c(res))
   })
   dim(ls)=c(ncol(dat),ncol(dat),nrow(dat))^2

 # compute the sumMissingProjector
 A = mul.tensor(as.tensor(ls), 3, as.tensor(rep(1,nrow(dat)^2)), 1)

 # compute the log-ratios between every pairs of data $x_{ik}-x_{jk}$, 
 #       for each variable k:
   aux = unclass(x)
    attr(aux,"orig") <- NULL
   prodVal = aux[idx$is,]-aux[idx$js,]
  # eliminate 1/inf, inf, inf/inf and x_i=x_j
     prodVal[is.infinite(prodVal)|is.nan(prodVal)] = 0
     prodVal[prodId,] = 0
     dimnames(prodVal) = list(s=NULL,i=NULL)
  # compute the X·X^t for each row (containing log-ratios):
 prodVal = t( sapply(1:nrow(prodVal), 
            function(i){ c(outer(prodVal[i,],prodVal[i,],"*"))  } ) )

  # compute the variance with our formula
 vr = gsi.svdsolve(A, mul.tensor(as.tensor(ls), c(3,2), as.tensor(prodVal), c(1:2)) )
 dim(vr) = c(ncol(dat),ncol(dat))
 if( giveCenter )
   attr(vr,"center")<-mean(cdtInv(x,dat),robust=FALSE)
 return(vr)
}
 
# NA   = MAR  =  (Missing at random) = not reported measured
# -Inf = MNAR =  (Missing not at random) = danger
# NaN  = SZ   =  (Structural Zero) = It makes no sense to speak about
# 0.0    = BDL  =  (Below detection limit) = Observed Zero



gsi.recodeM2C <- function(x,y=x,BDL,SZ,MAR,MNAR,NMV) {
  if( !is.numeric(x) )x<-gsi.plain(x)
  if( !is.numeric(y) ) y<- gsi.plain(y)
  if( !missing(BDL) ) y[is.BDL(x)]<-BDL
  if( !missing(SZ)  ) y[is.SZ(x)] <-SZ
  if( !missing(MAR) ) y[is.MAR(x)]<-MAR
  if( !missing(MNAR)) y[is.MNAR(x)]<-MNAR
  if( !missing(NMV) ) y[is.NMV(x)]<-NMV
  y
}

#gsi.recodeM2Clean <- function(x,y=x,BDL=NaN,SZ=NaN,MAR=NaN,MNAR=NA) {
#  y[is.BDL(x)]<-BDL
#  y[is.SZ(x)] <-SZ
#  y[is.MAR(x)]<-MAR
#  y[is.MNAR(x)]<-MNAR
#  y[is.NMV(x)]<-NMV
#  y
#}



gsi.recodeC2M <- function(x,y=x,na,nan,ninf,inf,neg,zero,pos) {
  if( !is.numeric(x) ) x<-gsi.plain(x)
  if( !is.numeric(y) ) y<- gsi.plain(y)
  if( !missing(na) && !is.null(na) ) y[is.na(x)&!is.nan(x)]<-na
  if( !missing(nan)&& !is.null(nan) ) y[is.nan(x)]          <-nan
  if( !missing(ninf)&& !is.null(ninf)) y[is.infinite(x)&x<0] <-ninf
  if( !missing(inf)&& !is.null(inf) ) y[is.infinite(x)&x>0] <-inf
  if( !missing(neg)&& !is.null(neg) ) y[is.finite(x)&x<0]   <-neg
  if( !missing(zero)&& !is.null(zero)) y[is.finite(x)&x==0]  <-zero
  if( !missing(pos)&& !is.null(pos))  y[is.finite(x)&x>0]   <-pos
  y
}

gsi.recodeM2Clean <- function(x,y=x,BDL=NaN,SZ=NaN,MAR=NaN,MNAR=NA,NMV) {
  if( !missing(BDL) ) y[is.BDL(x)]<-BDL
  if( !missing(SZ)  ) y[is.SZ(x)] <-SZ
  if( !missing(MAR) ) y[is.MAR(x)]<-MAR
  if( !missing(MNAR)) y[is.MNAR(x)]<-MNAR
  if( !missing(NMV) ) y[is.NMV(x)]<-NMV
  y
}


gsi.cleanR <- function(x) {
  x[!is.finite(x)]<-0.0
  x
}

missingProjector <- function(x,...,by="s") UseMethod("missingProjector")

missingProjector.acomp <- function(x,has=is.NMV(x),...,by="s") {
 #require("tensorA")
  if( is.tensor(has) ) {
  } else if( is.matrix(has) ) {
    has <- as.tensor(has)
    names(has) <- c(by,"i")
  } else
    has <- to.tensor(c(has),c(i=length(has)))
  iii <- function(x) {a<-ifelse(unclass(x)==0,0,1/unclass(x));attributes(a)<-attributes(x);a}
#   reorder.tensor( diag.tensor(has,by=by,mark="'") - 
#                  einstein.tensor(has,mark(has,mark="'",by=by),
#                                  diag=iii(margin.tensor(has,by=by)),
#                                  by=by))
  reorder.tensor( diag.tensor(has,by=by) - 
                 einstein.tensor(has,mark(has,by=by),
                                 diag=iii(margin.tensor(has,by=by)),
                                 by=by))
  #  gsi.diagGenerate(as.numeric(has))-(has %o% has)/sum(has) 
}




missingProjector.rcomp <- function(x,has=!(is.MAR(x)|is.MNAR(x)),...,by="s") {
 #require("tensorA")
 warning("missingProjector.rcomp: There is no established theory available for missings in nonrelative compositional geometry. Results are experimental. ")
  if( is.tensor(has) ) {
  } else  if( is.matrix(has) ) {
    has <- as.tensor(has)
    names(has) <- c(by,"i")
  } else
    has <- to.tensor(c(has),c(i=length(has)))
  iii <- function(x) {a<-ifelse(unclass(x)==0,0,1/unclass(x));attributes(a)<-attributes(x);a}
  reorder.tensor( diag.tensor(has,by=by,mark="'") -
                 einstein.tensor(has,has[[i=~"i'"]],
                                 diag=iii(margin.tensor(has,by=by)),
                                 by=by))
                 
  #  gsi.diagGenerate(as.numeric(has))-(has %o% has)/sum(has) 
}


missingProjector.aplus <- function(x,has=is.NMV(x),...,by="s") {
 #require("tensorA")
  if( is.tensor(has) ) {
  } else if( is.matrix(has) ) {
    has <- as.tensor(has)
    names(has) <- c(by,"i")
  } else
    has <- to.tensor(c(has),c(i=length(has)))
  reorder.tensor( diag.tensor(has,by=by,mark="'"), by=by) 
  #  gsi.diagGenerate(as.numeric(has)) 
}

missingProjector.rplus <- function(x,has=!(is.MAR(x)|is.MNAR(x)),...,by="s") {
 #require("tensorA")
 warning("missingProjector.rplus: There is no established theory available for missings in nonrelative positive geometry. Results are experimental. ")
  if( is.tensor(has) ) {
  } else if( is.matrix(has) ) {
    has <- as.tensor(has)
    names(has) <- c(by,"i")
  } else
    has <- to.tensor(c(has),c(i=length(has)))
  reorder.tensor( diag.tensor(has,by=by,mark="'"), by=by) 
  #  gsi.diagGenerate(as.numeric(has)) 
}
  
missingProjector.rmult <- function(x,has=is.finite(x),...,by="s") {
 #require("tensorA")
  if( !is.null(attr(x,"missingProjector") ) )
     return(attr(x,"missingProjector"))
  if( !is.null(attr(x,"orig") ) )
     return(missingProjector(attr(x,"orig") ))
  if( is.tensor(has) ) {
  } else if( is.matrix(has) ) {
    has <- as.tensor(has)
    names(has) <- c(by,"i")
  } else
    has <- to.tensor(c(has),c(i=length(has)))
  reorder.tensor( diag.tensor(has,by=by,mark="'"), by=by) 
  #  gsi.diagGenerate(as.numeric(has)) 
}


sumMissingProjector <- function(x,...) UseMethod("sumMissingProjector")

sumMissingProjector.acomp <- function(x,has=is.NMV(x),...) {
  n <- nrow(has)
  m <- ncol(has)
  canz <- rep(1,n) %*% has 
  ranz <- c(has %*% rep(1,m))
  erg <- diag(c(canz)) -  t(has/ifelse(ranz>0,ranz,1)) %*% has 
  erg
}

sumMissingProjector.aplus <- function(x,has=is.NMV(x),...) {
  n <- nrow(has)
  canz <- rep(1,n) %*% has 
  erg <- diag(c(canz))
  erg
}

sumMissingProjector.rcomp <- function(x,has=!(is.MAR(x)|is.MNAR(x)),...) {
 warning("sumMissingProjector.rcomp: There is no established theory available for missings in nonrelative compositional geometry. Results are experimental. ")
  n <- nrow(has)
  m <- ncol(has)
  canz <- rep(1,n) %*% has 
  ranz <- c(has %*% rep(1,m))
  erg <- diag(c(canz)) -  t( has / ifelse(ranz>0,ranz,1) ) %*% has 
  erg
}

sumMissingProjector.rplus <- function(x,has=!(is.MAR(x)|is.MNAR(x)),...) {
 warning("sumMissingProjector.rplus: There is no established theory available for missings in nonrelative positive geometry. Results are experimental. ")
  n <- nrow(has)
  canz <- rep(1,n) %*% has 
  erg <- diag(c(canz))
  erg
}

sumMissingProjector.rmult <- function(x,has=is.finite(x),...) {
  if( !is.null(mp<- attr(x,"missingProjector"))  )
    return( mp %e% one.tensor(c(s=dim(mp)["s"])))
  if( !is.null(orig<-attr(x,"orig")))
     return(sumMissingProjector(orig))
  n <- nrow(has)
  canz <- rep(1,n) %*% has 
  erg <- diag(canz)
  erg
}

observeWithAdditiveError <- function(x, sigma=dl/dlf, dl=sigma*dlf, dlf=3,
                                      keepObs=FALSE, digits=NA, obsScale=1,
                                      class="acomp") {
  if( missing(sigma) && missing(dl) )
    stop("You must specify at least the error standard deviation sigma or the detection limit dl")
  Yo <- x 
  Y <- gsi.recodeM2C(oneOrDataset(x),BDL=0.0,SZ=0.0,MAR=0,MNAR=0)
  if( !isTRUE(all.equal(dim(sigma),dim(Y))) ) {
    sigma <- rep(1,nrow(Y)) %o% rep(sigma,length.out=ncol(Y)) 
    stopifnot(all.equal(dim(sigma),dim(Y)))
  }
  if( !isTRUE(all.equal(dim(dl),dim(Y))) ) {
    dl <- rep(1,nrow(Y)) %o% rep(dl,length.out=ncol(Y)) 
    stopifnot(all.equal(dim(dl),dim(Y)))
  }
  if( !all(is.na(digits)) ) {
    if( !isTRUE(all.equal(dim(digits),dim(Y)))) {
      digits <- rep(1,nrow(Y)) %o% rep(digits,length.out=ncol(Y)) 
      stopifnot(all.equal(dim(digits),dim(Y)))
    }
    if( !isTRUE(all.equal(dim(obsScale),dim(Y)))) {
      obsScale <- rep(1,nrow(Y)) %o% rep(obsScale,length.out=ncol(Y)) 
      stopifnot(all.equal(dim(obsScale),dim(Y)))
    }
  }
  Z <- unclass(Y)+structure(rnorm(length(sigma),0,sigma),dim=dim(sigma))
  if( ! all(is.na(digits)) ) {
    Z <- round(Z/obsScale,digits)*obsScale
  }
  Zc <- ifelse(Z<dl,0,Z)
  to <- totals(rplus(Zc))
  Zc <- gsi.recodeM2C(Zc,
                      Yo,BDL=BDLvalue,MAR=MARvalue,MNAR=MNARvalue,SZ=SZvalue)
  switch(match.arg(class,c("acomp","aplus","rcomp","rplus")),
         "acomp"=structure(Zc/to,class="acomp",detectionlimit=dl/to,sigma=sigma/to,obs=if(keepObs) Z/to),
         "aplus"=structure(Zc,class="acomp",detectionlimit=dl,sigma=sigma,obs=if(keepObs) Z),
         "rcomp"=structure(Zc/to,class="rcomp",detectionlimit=dl/to,sigma=sigma/to,obs=if(keepObs) Z/to),
         "rplus"=structure(Zc,class="rplus",detectionlimit=dl,sigma=sigma,obs=if(keepObs) Z)
         )
}


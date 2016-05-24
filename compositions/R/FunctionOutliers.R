# FunctionOutlier.R
# (C) 2008  K.Gerald van den Boogaart, Raimon Tolosana-Delgado
# Pulished under the GNU Public License Version 2
# part of the compositions package
#
# The functions in this file are concerned with outlier detection
# based on robust mean and variance estimation provided by the 
# robustbase package.



# MahalanobisDist
# This function computes the Mahalanobis distances of the individual
# datasets to the center. It can be be used with given or automatically
# estimated center and covariances. There is also a robust parameter,
# which should be set to true for usefull results.


# 
MahalanobisDist <- function(x,center=NULL,cov=NULL,inverted=FALSE,...) UseMethod("MahalanobisDist",x)


#MahalanobisDist.default <- function(x,center=NULL,cov=NULL,inverted=FALSE,...,pairwise=FALSE,robust=FALSE) MahalanobisDist(rmult(x),center,cov,inverted=inverted,...,pairwise=pairwise,robust=robust)


MahalanobisDist.rmult <- function(x,center=NULL,cov=NULL,inverted=FALSE,...,goodOnly=NULL,pairwise=FALSE,pow=1, robust=FALSE,giveGeometry=FALSE) {
  if( is.null(goodOnly) ) 
    xx <- x
  else 
    xx<-rmult(oneOrDataset(x)[goodOnly,,drop=FALSE])
  if( is.null(cov)) {
    cov <- var(xx,robust=robust,giveCenter=TRUE)
    inverted <- FALSE
  }
  if( is.null(center) ) {
    if( is.null(attr(cov,"center")))
      center <- mean(xx,robust=robust)
    else
      center <- attr(cov,"center")
  }
  if( pairwise ) {
    z    <- unclass(idt(x))
    power <- if(inverted) 0.5 else -0.5
    SQRT <- with(svd(cov), u %*% gsi.diagGenerate(d^power) %*% t(v))
    erg <- dist(z %*% SQRT,...)^pow  
  } else {
    more <- t(unclass(oneOrDataset(x-rmult(center))))
    if( ! is.matrix(cov) )
      cov <- matrix(cov,nrow=1,ncol=1)
    ss <- if(inverted) cov %*% more else solve(cov,more,...)
    erg <- sqrt(c(rep(1,nrow(more)) %*% (ss*more)))^pow
  }
  if( giveGeometry ) {
    attr(erg,"center") <- center
    attr(erg,"cov") <-cov
  }
#      if( any( rob$mah != erg[goodOnly] ) )
#        warning("Imprecision detected:",(rob$mah - erg[goodOnly]))
  return(erg)
}


MahalanobisDist.rcomp <- MahalanobisDist.acomp <- MahalanobisDist.rplus <- MahalanobisDist.aplus <- function(x,center=NULL,cov=NULL,inverted=FALSE,...,goodOnly=NULL, pairwise=FALSE,pow=1,robust=FALSE,giveGeometry=FALSE) {
  if( any( class(x) %in% c("acomp","rcomp")) )
    if( !is.null(cov) && gsi.getD(x) == ncol(cov))
      cov <- clrvar2ilr(cov)
  erg <- MahalanobisDist(idt(x),idt(center),cov,inverted,...,goodOnly=goodOnly,pairwise=pairwise,pow=pow,robust=robust,giveGeometry=giveGeometry)
  if( giveGeometry ) {
    attr(erg,"center") <- idtInv(attr(erg,"center"),x)
    attr(erg,"cov")<-cov
  }
  return(erg)
}

# covMcd.acomp is like covMcd from rrcov but applicable to acomp objects
# and extreme situations.
#covMcd.acomp <- function(X,...) {
#  if( nrow(X) < 2*ncol(X) ) {
#    return(list(center=mean(X),cov=var(X)))
#  }
#  rob <- covMcd(ilr(X),...)
#  rob$center <- ilr.inv(rob$center)
#  rob$cov    <- ilrvar2clr(rob$cov)
#  rob
#}

#covMcd.rcomp <- function(X,...) {
#  if( nrow(X) < 2*ncol(X) ) {
#    return(list(center=mean(X),cov=var(X)))
#  }
#  rob <- covMcd(ipt(X),...)
#  rob$center <- ipt.inv(rob$center)
#  rob$cov    <- ilrvar2clr(rob$cov)
#  rob
#}


## functions to codify the ourliers in several ways
# X:   the data set (with column names!!!)
# mat: contains a matrix like X, with logicals
# numcode: contains a number, a sum of bits
# alphacode: contains a string, the concatenation of column names
#gsi.numcode <- function(mat){
#  aux = mat%*% 2^(0:(ncol(mat)-1)) 
#  return(c(aux))
#}


                                        #gsi.alphacode <- function(mat){
#  sapply(1:nrow(mat),function(i){
#   paste(colnames(mat)[mat[i,]],collapse="|")
#  })
#}
# extracts the bitcode from a numcode plus a dimension (number of columns)
#gsi.bits <- function(x,D){
# if(max(x)>(2^D-1)){stop("dimension provided is smaller than minimal")}
# output = ""
# x = x+1
# for(k in 1:ceiling(log2(x))){
#  y = x/2
#  output=paste(output,ifelse(y==floor(y),1,0),sep="")
#  x = x - floor(y)
#  if(x==1){return(paste(output,paste(rep("0",D-k),collapse=""),sep=""))}
# }
# return(output)
#}
# ??
#num2mat = function(num,D){
# aux = sapply(num,function(x){as.logical(as.integer(strsplit(gsi.bits(x,D),split="")[[1]]))})
# return(t(aux))
#}
# ??
#fac2num = function(fak){
# aux = as.integer(levels(fak)[as.integer(fak)])
# aux[is.na(aux)]=0
# return(aux)
#}
# ??
#alpha2mat = function(alpha,colnames){
# mat = matrix(FALSE,ncol=length(colnames),nrow=length(alpha))
# colnames(mat)=colnames
# for(i in 1:length(alpha)){
#  aux = strsplit(alpha[i],split="|",extended=FALSE)[[1]]
#  mat[i,aux]= TRUE
# }
# return(mat)
#}




###### MAIN FUNCTION ############
# There are 3 classifiers doing roughly similar things

######  basic classifier function, producing a factor with the kind of outlier (or "ok" if typical)
# Old Name : myClassifier.acomp
#OutlierClassifier.acomp <- function(X,...,cutlevel=4,goodOnly=1:nrow(X)) {
#  coord <- idt(X)
#  cv <- covMcd(coord[goodOnly,],...)
#  Rcenter <- cv$center
#  Cmean   <- ilr.inv(Rcenter)
#  myIvar  <- cv$cov
#  myCvar  <- ilrvar2clr(myIvar)
#  aux <- (rmult(coord)-rmult(Rcenter))
#  genOutlier <- scalar(aux,t(solve(myIvar,t(unclass(aux)))))/gsi.getD(X)
#  redOutlier <- sapply(1:NCOL(X),function(i) {
#    xIvar <- clrvar2ilr(ilrvar2clr(myIvar)[-i,-i])
#    ax <- clr2ilr(ilr2clr(aux)[,-i])
#    scalar(ax,t(solve(xIvar,t(unclass(ax)))))/(gsi.getD(X)-1)
#  })
#  colnames(redOutlier) = colnames(X)
#  #singleOut<-rmult(clr(X-Cmean)^2)/rmult(gsi.diagExtract(myCvar))
#  erg = cbind(genOutlier,redOutlier)
#  if( missing(cutlevel) )
#    return(erg)
#  genOutlier <- genOutlier > cutlevel
#  redID <- gsi.numcode(redOutlier > cutlevel)
#  Fak = ifelse(genOutlier,binary(redID),"ok")
#  factor(Fak)
#}

#myClassifier2.acomp <- function(X,...,alpha=0.01,type=c("best","all","type","outlier"),cutlevel,goodOnly=1:nrow(X)) {
#  takeIf <- function(c,x,y) {
#    y <- ifelse(c,y,y)
#    y[c]<-x
#    y
#  }
#  coord <- idt(X)
#  cv <- covMcd(coord[goodOnly,],...)
#  Rcenter <- cv$center
#  Cmean   <- ilr.inv(Rcenter)
#  myIvar  <- cv$cov
#  myCvar  <- ilrvar2clr(myIvar)
#  aux <- (rmult(coord)-rmult(Rcenter))
#  genOutlier <- scalar(aux,t(solve(myIvar,t(unclass(aux)))))#/gsi.getD(X)
#  isOutlier <- pchisq(genOutlier,df=gsi.getD(X)-1,lower.tail=FALSE)<=alpha
#  if( any( isOutlier ) ) {
#    aux2 <- aux[isOutlier,]
#    checkRedOutlier <- function(i) {
#      xIvar <- clrvar2ilr(myCvar[-i,-i])
#      ax <- clr2ilr(ilr2clr(aux2)[,-i])
#      pchisq(scalar(ax,t(solve(xIvar,t(unclass(ax))))),df=nrow(xIvar),lower.tail=FALSE)
#    }
#    redOutlier <- sapply(1:NCOL(X),checkRedOutlier)
#    colnames(redOutlier) = colnames(X)
#    oneExplains <- c(((redOutlier > alpha)+0) %*% rep(1,ncol(X))>0)
#    bestExplain <- apply(redOutlier,1,which.max)
#  } else oneExplains <- c()
##  multiOutlier <- isOutlier & ! oneExplains
#  type = match.arg(type)
#  #singleOut<-rmult(clr(X-Cmean)^2)/rmult(gsi.diagExtract(myCvar))
##  genOutlier <- genOutlier > cutlevel
##  redID <- numcode(redOutlier > cutlevel)
#  if( type == "best" )
#    Fak = takeIf(isOutlier,ifelse(oneExplains,colnames(X)[bestExplain],"?"),"ok")
#  if( type == "all" )
#    Fak = takeIf(isOutlier,ifelse(oneExplains,binary(gsi.numcode(redOutlier > alpha)),"?"),"ok")
#  if( type == "type" )
#    Fak = takeIf(isOutlier,ifelse(oneExplains,"1","?"),"ok")
#  if( type == "outlier" )
#    Fak = takeIf(isOutlier,ifelse(oneExplains,"outlier","outlier"),"ok")    
#  factor(Fak)
#}


# The pStore environment caches simulations and results from the slow
# distribution routines for the MaxMahalanobis and Mahalanobis distributions
if( ! exists("gsi.pStore") )
  gsi.pStore <- new.env(hash=TRUE,emptyenv())

# The cdf of the distribution of the Maxium Robust Mahalanobis distance
# N = Number of rows in the dataset
# d = numer of independent!!! columns (id est d=D-1 for acomp)
# replicates = the number of simulations to establish the distribution
# resample = Do a new simulation rather than using the old one.
# ...=further arguments to the covMcd-routine
pMaxMahalanobis <- function(q,N,d,...,pow=1,replicates=998,resample=FALSE,robust=TRUE) {
  id <- paste(deparse(list(N=N,d=d,...,replicates=replicates,robust=robust,type="MaxMahalanobis")),collapse="\n")
  if( resample || is.null(s<-mget(id,gsi.pStore,ifnotfound=list(NULL))[[1]])   ) {
    s <- rMaxMahalanobis(replicates,N,d,...,pow=1,robust=robust)
    assign(id,sort(s),envir=gsi.pStore)
  }
  sapply(q,function(x) mean(c(0,(s^pow)<=x,1)))
}

# The quantile function of the distribution of the Maxium Robust Mahalanobis distance
# N = Number of rows in the dataset
# d = numer of independent!!! columns (id est d=D-1 for acomp)
# replicates = the number of simulations to establish the distribution
# resample = Do a new simulation rather than using the old one.
# ...=further arguments to the covMcd-routine
qMaxMahalanobis <- function(p,N,d,...,pow=1,replicates=998,resample=FALSE,robust=TRUE) {
  id <- paste(deparse(list(N=N,d=d,...,pow=pow,replicates=replicates,robust=robust,type="MaxMahalanobis")),collapse="\n")
  if( resample || is.null(s<-mget(id,gsi.pStore,ifnotfound=list(NULL))[[1]])   ) {
    s <- rMaxMahalanobis(replicates,N,d,...,pow=1,robust=robust) 
    assign(id,sort(s),envir=gsi.pStore)
  }
  s[1+round(p*(length(s)-1))]^pow
}

# A random number generator of the distribution of the Maxium Robust Mahalanobis distance
# N = Number of rows in the dataset
# d = numer of independent!!! columns (id est d=D-1 for acomp)
# replicates = the number of simulations to establish the distribution
# ...=further arguments to the covMcd-routine
rMaxMahalanobis <- function(n,N,d,...,pow=1,robust=TRUE) {
  return( rEmpiricalMahalanobis(n,N,d,...,pow=pow,robust=robust,sorted=max) )
  # Old Code
#  rr <- robust
#  control <- attr(robust,"control")
#  if( is.logical(robust) ) rr <- if(robust) "mcd" else "pearson"
#  if( rr == "mcd" ) require(robustbase)
#  s<-switch(rr,
#            pearson={
#              if( d==1 )
#                replicate(n,{
#                  x <- rnorm(n,mean=rnorm(1))
#                  max(c(((c(x)-mean(x))^2)/c(var(x))))
#                })
#              else {
#                quickMah <- function(x) {
#                  v <- var(x)
#                  xc <- t(x)-meanCol(x)
#                  max(rep(1,ncol(v)) %*% solve(v,xc,...) * xc )
#                }
#                replicate(ntrial,quickMah(matrix(rnorm(n*d,mean=rnorm(1)),ncol=d)))
#              }
#            },
#            mcd={
#              if( d > 1 )
#                replicate(n,sqrt(max(covMcd(structure(rnorm(N*d,mean=rnorm(1)),dim=c(N,d)),control=control)$mah)))
#              else {
#                replicate(n,{
#                  x   <- structure(rnorm(N*d),dim=c(N,d))
#                  erg <- covMcd(x,control=control)
#                  max( (x-erg$center)^2/c(erg$cov) )
#                })
#              }
#              s
#            },
#            stop("rMaxMahalanobis: Unkown robustness option:",robust)
#            )
#  s
}
  
# The cdf of the distribution of the Robust Mahalanobis distance
# of a typical data point
# N = Number of rows in the dataset
# d = numer of independent!!! columns (id est d=D-1 for acomp)
# replicates = the number of simulations to establish the distribution
# resample = Do a new simulation rather than using the old one.
# ...=further arguments to the covMcd-routine
pEmpiricalMahalanobis <- function(q,N,d,...,pow=1,replicates=100,resample=FALSE,robust=TRUE) {
  #require(robustbase)
  id <- paste(deparse(list(N=N,d=d,...,replicates=replicates,robust=robust,type="EmpiricalMahalanobis")),collapse="\n")
  if( resample || is.null(s<-mget(id,gsi.pStore,ifnotfound=list(NULL))[[1]])   ) {
    s <- rEmpiricalMahalanobis(replicates,N,d,...,pow=1,robust=robust)
    assign(id,sort(s),envir=gsi.pStore)
  }
  sapply(q,function(x) mean(c(0,(s^pow)<=x,1)))
}


# The quatile function of the distribution of the Robust Mahalanobis distance
# of a typical data point
# N = Number of rows in the dataset
# d = numer of independent!!! columns (id est d=D-1 for acomp)
# replicates = the number of simulations to establish the distribution
# resample = Do a new simulation rather than using the old one.
# ...=further arguments to the covMcd-routine

qEmpiricalMahalanobis <- function(p,N,d,...,pow=1,replicates=100,resample=FALSE,robust=TRUE) {
  #require(robustbase)
  id <- paste(deparse(list(N=N,d=d,...,replicates=replicates,robust=robust,type="EmpiricalMahalanobis")),collapse="\n")
  if( resample || is.null(s<-mget(id,gsi.pStore,ifnotfound=list(NULL))[[1]])   ) {
    s <- rEmpiricalMahalanobis(replicates,N,d,...,pow=1,robust=robust) 
    assign(id,sort(s),envir=gsi.pStore)
  }
  s[1+round(p*(length(s)-1))]^pow
}


# The random number generator of the distribution of the Robust Mahalanobis distance
# of a typical data point
# N = Number of rows in the dataset
# d = numer of independent!!! columns (id est d=D-1 for acomp)
# replicates = the number of simulations to establish the distribution
# resample = Do a new simulation rather than using the old one.
# ...=further arguments to the covMcd-routine

rEmpiricalMahalanobis <- function(n,N,d,...,sorted=FALSE,pow=1,robust=TRUE) {
  if( is.logical(robust) ) rr <- if(robust) "mcd" else "pearson"
  if( is.logical(sorted) )
    if( sorted )
      sorted <- function(x) sort(c(x))
    else
      sorted <- function(x) x
  else if( is.numeric(sorted) ) {
    ndx <- sorted
    sorted <- function(x) sort(c(x))[ndx]
  }
  params <- list(...)
  s <- switch(rr,
              pearson={
                if( d > 1 ) {
                  quickMah <- function(x) {
                    v <- var(x)
                    xc <- t(x)-meanCol(x)
                    sorted(sqrt(rep(1,ncol(v)) %*% solve(v,xc,...) * xc)^pow)
                  }
                  replicate(n,quickMah(structure(rnorm(N*d),dim=c(N,d))))
                } else {
                  replicate(n,{
                    x <- rnorm(n,mean=rnorm(1))
                    od<-covMcd(x)
                    sorted(sqrt(((c(x)-mean(x))^2)/c(var(x)))^pow)
                  })
                }
              },
              mcd={
                #require("robustbase")
                if( d > 1 )
                  replicate(n,sorted(sqrt(do.call("covMcd",c(list(structure(rnorm(N*d),dim=c(N,d))),params))$mah))^pow)
                else {
                  replicate(n,{
                    x   <- structure(rnorm(N*d),dim=c(N,d))
                    erg <- do.call("covMcd",c(list(x),params))
                    sorted(sqrt( (x-erg$center)^2/c(erg$cov) )^pow)
                  })
                }
              },
              stop("rEmpiricalMahalanobis: Unkown robustness type:", robust)
              )
  if(is.matrix(s)) t(s) else s
}

# The cdf of the distribution of the portion of Mahalanobis distances
# over a given cutoff
# of a typical data point
# N = Number of rows in the dataset
# d = numer of independent!!! columns (id est d=D-1 for acomp)
# cut= the cutoff
# replicates = the number of simulations to establish the distribution
# resample = Do a new simulation rather than using the old one.
# ...=further arguments to the covMcd-routine
pPortionMahalanobis <- function(q,N,d,cut,...,replicates=1000,resample=FALSE,pow=1,robust=TRUE) {
  #require(robustbase)
  id <- paste(deparse(list(N=N,d=d,...,replicates=replicates,robust=robust,type="PortionMahalanobis")),collapse="\n")
  if( resample || is.null(s<-mget(id,gsi.pStore,ifnotfound=list(NULL))[[1]])   ) {
    s <- rEmpiricalMahalanobis(replicates,N,d,...,sorted=TRUE,robust=robust)
    assign(id,s,envir=gsi.pStore)
  }
  sapply(q,function(x) mean(c(0,(s^pow)<=x,1)))
}


rPortionMahalanobis <- function(n,N,d,cut,...,pow=1,robust=TRUE) {
  rEmpiricalMahalanobis(n,N,d,...,pow=pow,robust=robust,
                        sorted=function(x) sum(x<=cut))
}

# The quantile function of the distribution of the portion of Mahalanobis distances
# over a given cutoff
# of a typical data point
# N = Number of rows in the dataset
# d = numer of independent!!! columns (id est d=D-1 for acomp)
# cut= the cutoff
# replicates = the number of simulations to establish the distribution
# resample = Do a new simulation rather than using the old one.
# ...=further arguments to the covMcd-routine
qPortionMahalanobis <- function(p,N,d,cut,...,replicates=1000,resample=FALSE,pow=1,robust=TRUE) {
  #require(robustbase)
  id <- paste(deparse(list(N=N,d=d,...,replicates=replicates,type="PortionMahalanobis")),collapse="\n")
  if( resample || is.null(s<-mget(id,gsi.pStore,ifnotfound=list(NULL))[[1]])   ) {
    s <- rEmpiricalMahalanobis(replicates,N,d,...,sorted=TRUE,robust=robust) 
    assign(id,s,envir=gsi.pStore)
  }
  s[1+round(p*(length(s)-1))]^pow
}




# A random number generator of the sorted list of Mahalanobis distances
# of a typical data point
# N = Number of rows in the dataset
# d = numer of independent!!! columns (id est d=D-1 for acomp)
# cut= the cutoff
# replicates = the number of simulations to establish the distribution
# resample = Do a new simulation rather than using the old one.
# ...=further arguments to the covMcd-routine
#rSortedMahalanobis <- function(n,N,d,...) {
#  require(robustbase)
#  params <- list(...)
#  if( d > 1 )
#    s <- replicate(n,sqrt(sort(do.call("covMcd",c(list(structure(rnorm(N*d),dim=c(N,d))),params))$mah)))
#  else {
#    s <- replicate(n,sort({
#      x   <- structure(rnorm(N*d),dim=c(N,d))
#      erg <- do.call("covMcd",c(list(x),params))
#      c( (x-erg$center)^2/c(erg$cov))
#    }))
#  }
#  t(s)
#}

# A routine to display a matrix as an image (used for debugging only)
#showMat <- function(X) {
#  dev <- dev.cur()
#  x11()
#  image(X)
 # dev.set(dev)
#}



# The cdf of the distribution of the empirical p-quantile of Mahalanobis distances
# of a typical data point
# N = Number of rows in the dataset
# d = numer of independent!!! columns (id est d=D-1 for acomp)
# p= the probability for qunatile
# replicates = the number of simulations to establish the distribution
# resample = Do a new simulation rather than using the old one.
# ...=further arguments to the covMcd-routine

pQuantileMahalanobis <- function(q,N,d,p,...,replicates=1000,resample=FALSE,ulimit=TRUE,pow=1,robust=TRUE) {
  #require(robustbase)
  id <- paste(deparse(list(N=N,d=d,...,pow=pow,robust=robust,replicates=replicates,type="pQuantileMahalanobis")),collapse="\n")
  if( resample || is.null(PreSort<-mget(id,gsi.pStore,ifnotfound=list(NULL))[[1]])   ) {
    s <- rEmpiricalMahalanobis(replicates,N,d,...,pow=pow,robust=robust,sorted=TRUE)
#    attr(s,"palphaEnv") <- new.env(hash=TRUE,emptyenv())
    PreSort <- apply(s,2,sort)
    st <- t(s)
#    attr(s,"PreSort") <- PreSort
    realP <- sapply(1:nrow(PreSort),function(i) mean(apply(st <= PreSort[i,],2,all)))
    attr(PreSort,"realP") <- realP
    assign(id,PreSort,envir=gsi.pStore)
  } else {
    realP <- attr(PreSort,"realP")
  }
  line <- sapply(p,
                 function(pp) {
                   if( !ulimit ) {
                     line <- sum(realP < 1-pp)+1
                     if( line > length(realP) ) {
                       line <- length(realP)
                       warning("Extreme limits might be incorrect in joint confidence band")
                     }
                   }
                   else {
                     line <- sum(realP <= 1-pp)
                     if( line < 1 ) {
                       line <- 1
                       warning("Extremely small limits might be incorrect in joint confidence band")
                     }
                     
                   }
                   line
                 }
                 )
  sapply(line,function(xline) sapply(q,function(xx) mean(PreSort[xline,]<xx)))
}

# The central decision routine checking for outliers in compositional data
# X= the dataset
# ... further arguments to solve
# further arguments to mcd can be given through the attribute of robust
# alpha= The alpha level of the outlier detection
# Old Name mahOutliers.acomp
IsMahalanobisOutlier <- function(X,...,alpha=0.05,goodOnly=NULL,replicates=1000,corrected=TRUE,robust=TRUE,crit=NULL) {
  if( is.null(goodOnly) )
    N <- nrow(X)
  else
    N <- nrow(X[goodOnly,])
  nk<-0
  alphaExt <- min(alpha,1-alpha)
  if( missing(replicates) ) replicates <- max(ceiling(50/alphaExt),replicates)
  while( is.null(crit) ) {
    crit <- if(corrected)
      qMaxMahalanobis(1-alpha,N,ncol(idt(X)),...,replicates=replicates,robust=robust)
    else
      qEmpiricalMahalanobis(1-alpha,N,ncol(idt(X)),...,replicates=
                            if( missing(replicates) )
                            1000
                   else
                            replicates,robust=robust
                            )
    if( !is.finite(crit) ) {
      crit<- NULL;
      warning("IsMahalanobisOutlier: Problems computing Quantiles")
      nk<-nk+1
      if( nk>6 ) stop("IsMahalanobisOutlier: Problems computing Quantiles")
    }
  }
  if( !is.finite(crit) ) stop("IsMahalanobisOutlier: Error in crit")
  md <- MahalanobisDist(X,...,goodOnly=goodOnly,robust=robust)
  if( any(!is.finite(crit)) ) stop("IsMahalanobisOutlier: Error in MahalanobisDist")
  (md>crit)
}

# The myClassifier3 is the routine doning the classification
# concerning single component outliers.
# X The dataset
# alpha the alpha-level to derive the cutoff
# type The type of classification to be reportet:
# best : ok ? or the component giving the best explanation
# all  : ok or a bitcode reporting all explaining components
# type : ok ? 1: wheter a outlier is single Component outlier
# outlier: ok outlier: whether or not its an outlier
# grade : ok, extrem, outlier : grade of outlyingness
# cutlevel (unfortunatly not supported at this moment)
#  intendet to provide an alternate cutlevel in case a ncase Analysis
# goodOnly: The subset of X to be used for the covMcd
# ... : Further parameters to covMcd
# corrected: Should the detection be corrected for multiple testing
# RedCorrected: Should the acceptance of reduced outliers be corrected for
#   multiple testing
# The result is a factor
# Old Name: myClassifier3
OutlierClassifier1 <- function(X,...) UseMethod("OutlierClassifier1",X)
OutlierClassifier1.acomp <- function(X,...,alpha=0.05,type=c("best","all","type","outlier","grade"),goodOnly=NULL,corrected=TRUE,RedCorrected=FALSE,robust=TRUE) {
  type<-match.arg(type)
  cls <- class(X)
  reclass <- function(x) structure(x,class=cls)
  isOutlier <- IsMahalanobisOutlier(X,...,alpha=alpha,goodOnly=goodOnly,corrected=corrected,robust=robust)
  if( type=="grade" )
    isOutlier2 <- IsMahalanobisOutlier(X,...,alpha=alpha,goodOnly=goodOnly,corrected=FALSE,robust=robust)
  if( any( isOutlier ) ) {
    checkRedOutlier <- function(i) {
      IsMahalanobisOutlier(reclass(X[,-i,drop=FALSE]),alpha=alpha,corrected=RedCorrected,robust=robust)
    }
    redOutlier <- sapply(1:NCOL(X),checkRedOutlier)
    colnames(redOutlier) = colnames(X)
    oneExplains <- c(((1-redOutlier) %*% rep(1,ncol(X)))>0)
    checkRed2Outlier <- function(i) {
      MahalanobisDist(reclass(X[,-i,drop=FALSE]),robust=robust)
    }
    red2Outlier <- sapply(1:NCOL(X),checkRed2Outlier)
    bestExplain <- apply(red2Outlier,1,which.min)
  } else oneExplains <- c()
  type = match.arg(type)
  if( type == "best" ) {
    Fak = ifelse(isOutlier,ifelse(oneExplains,colnames(X)[bestExplain],"?"),"ok")
    lev = c("ok","?",colnames(X))
  }
  if( type == "all" ) {
    gsi.numcode <- function(mat){
      aux = mat%*% 2^(0:(ncol(mat)-1)) 
      return(c(aux))
    }
    Fak = ifelse(isOutlier,binary(gsi.numcode(1-redOutlier[,rev(1:ncol(redOutlier)),drop=FALSE])),"ok")
    lev = unique(c("ok",sort(Fak)))
  }
  if( type == "type" ) {
    Fak = ifelse(isOutlier,ifelse(oneExplains,"1","?"),"ok")
    lev = c("ok","?","1")
  }
  if( type == "outlier" ) {
    Fak = ifelse(isOutlier,ifelse(oneExplains,"outlier","outlier"),"ok")
    lev = c("ok","outlier")
  }
  if( type == "grade" ) {
    Fak = ifelse(isOutlier,"outlier",ifelse(isOutlier2,"extreme","ok"))
    lev = c("ok","extreme","outlier")
  }
  factor(Fak,levels=lev)
}


# A robust mean calculator
# Now implemented as mean(X,robust=TRUE)
#rmean.acomp <- function(X,...,rmethod=covMcd) {
#  coord <- ilr(X)
 # re <- rmethod(coord,...)
#  ilr.inv(re$center)
#}

# Computes robustly estimated Mahalanobis distances between the
# observations (see myMahalanobis.acomp for distances to the center)
# Now impemented as MahalanobisDist ( pairwise=TRUE)
#robustMahalanobisDist.acomp <- function(X,...,rmethod=covMcd,fullData=X) {
#  z    <- unclass(idt(X))
 # SQRT <- with(svd(rmethod(unclass(idt(fullData)))$cov), u %*% gsi.diagGenerate(d^-0.5) %*% t(v))
 # dist(z %*% SQRT,...)  
#}

# Computes distances of directional differences to the robustly estimated center
# 
#dirDist.acomp <- dirDist.rcomp <- function(X,center,...) {
#  dist(idt(normalize(X-center)),...)  
#}


# The maximum density clustering as described in the article
# X the dataset
# ... Further arguments to covMcd
# sigma The bandwidht of the kernel in normalized space
# radius the radius in which neighbouring maxima are searched in normalized space
# asig The relatively factor of assumed spread for correcting the central density of the groups.
# minGrp The minimum size of a group
ClusterFinder1 <- function(X,...) UseMethod("ClusterFinder1",X)
ClusterFinder1.acomp <- function(X,...,sigma=0.3,radius=1,asig=1,minGrp=3,robust=TRUE) {
  Dim    <- gsi.getD(X)-1
  N      <- nrow(X)
  disQ   <- as.matrix(MahalanobisDist(X,...,pow=2,robust=robust,pairwise=TRUE))
  kernel <- exp(-disQ/(2*sigma^2))/sqrt((2*pi*sigma^2)^Dim)
  density<- kernel %*% rep(1/N,N)
  neigh  <- apply(disQ,1,order)
  locMax <- sapply(1:N,function(i) max(density[disQ[i,]<radius^2*Dim]))
  isMax  <- c(density == locMax)
  alloc  <- (1:N)[isMax]
  or     <- order(density[alloc],decreasing=TRUE)
  alloc  <- alloc[or]
  while(TRUE) {
    likeli <- sapply(alloc,function(i) log(density[i])-disQ[i,]/(2*asig^2))
    grpnr  <- if( length(alloc) > 1 )
      apply(likeli,1,which.max)
    else
      rep(1,length(likeli))
    assign <- factor(grpnr)
    if( length(levels(assign))< length(alloc) ) {
      alloc <- alloc[1:length(alloc) %in% unique(grpnr)]
    } else
      break;
  }
  prob   <- acomp(exp(likeli))
  members<- table(assign)
  realGroups <- c(members>=minGrp)
  types  <- factor(ifelse( realGroups[grpnr] , grpnr , "single" ))
  list(groups=assign,isMax=isMax,prob=prob,nmembers=members,density=density,likeli=likeli,types=types,typeTbl=table(types))  
}


# plots the centers of the clusters 
#showCenters <- function(X,...,cl= ClusterFinder.acomp(X,...),add=FALSE) {
#  #plot(ilr(X),col=ifelse(cl$isMax,"red","black"),pch=3)
#  #points(rbind(ilr(X)[cl$isMax,],ilr(X)[cl$isMax,]),col="red",pch=20)
#  plot(X,col=as.numeric(cl$groups),pch=3,add=add)
#}


## function to provide a coloured biplot (admits SVD, princomp and prcomp)
#coloredBiplot <- function(xrf, scale=1, choice=c(1,2), pc.biplot=FALSE, 
#        xcol="black", ycol="red", xpch=4, ypch=1, cex=1, 
#        xarrows=FALSE, yarrows=!xarrows, xnames=NULL, ynames=NULL,...){
# # X : cases (points)
# # Y : variables (arrows)
#   if("princomp" %in% class(xrf)){
#      X = xrf$scores
#      Y = xrf$loadings
#      if(is.null(xnames)){xnames=rownames(X)}
#      if(is.null(ynames)){ynames=rownames(Y)}
#      l = xrf$sdev
#      fac = ifelse(pc.biplot,sqrt(nrow(X)),1)
#   }
#   if("prcomp" %in% class(xrf)){
#      X = xrf$x
#      Y = xrf$rotation
#      if(is.null(xnames)){xnames=rownames(X)}
#      if(is.null(ynames)){ynames=rownames(Y)}
#      l = xrf$sdev
#      fac = ifelse(pc.biplot,sqrt(nrow(X)),1)
#   }
#   if("list" %in% class(xrf)){
#     warning("Attention: biplot tries to interpret xrf as result of an svd")
#     if(pc.biplot){
#      X = as.matrix(xrf$u)*sqrt(nrow(X))
#      Y = xrf$v/sqrt(nrow(X))
#      l = xrf$d
#     }
#     if(!pc.biplot){
#      X = as.matrix(xrf$u)  %*% diag(xrf$d)
#      Y = xrf$v
#      l = xrf$d /sqrt(nrow(X))
#     }
#      warning("recall that singular value decomposition does not center the data set")
#      fac = 1
#   }
#   Xx = X[,choice[1]]*(l[choice[1]]^(1-scale))*fac
#   Xy = X[,choice[2]]*(l[choice[2]]^(1-scale))*fac
#   Yx = Y[,choice[1]]*(l[choice[1]]^(scale))/fac
#   Yy = Y[,choice[2]]*(l[choice[2]]^(scale))/fac
## abandoned attempts to scale X and Y in the same way as R does:
##      span = function(x){ c(range(x)%*%c(-1,1)) }
##      escala = 0.8*max(span(Xx)/span(Yx),span(Xy)/span(Yy))
##      Yx = Yx #* span(Xx)/span(Yx) # escala
##      Yy = Yy #* span(Xy)/span(Yy) # escala
## more abandoned attempts to scale X and Y in the same way as R does:
##    escala =0.8*as.double(xlm%*%c(-1,1)/(xlml%*%c(-1,1)))
##    xlm=c(min(xlml*escala,xlm),max(xlml*escala,xlm))
##      xlm = c(-1,1)*max(abs(xlm))
##      ylm = xlm
##    ylm=c(min(ylml*escala,ylm),max(ylml*escala,ylm))
##      ylm = c(-1,1)*max(abs(ylm))
##    ratio=as.double(xlm%*%c(-1,1))/(ylm%*%c(-1,1))
##    ylm=ratio*ylm
##      print(c(xlm,ylm))

#   par(pty="s")
#   plot(c(Xx,Yx)*1.05,c(Xy,Yy)*1.05,type="n",col=xcol,ann=FALSE,pch=xpch,cex=cex,...)
#   if( xarrows ) 
#     arrows(x0=0,y0=0, x1=Xx, y1=Xy,length=0.1,col=xcol)
#   else
#     points(x=Xx, y=Xy,col=xcol,pch=xpch)
#   if( yarrows ){
#     arrows(x0=0,y0=0, x1=Yx, y1=Yy,length=0.1,col=ycol)
#   }else{
#     points(x=Yx, y=Yy,col=ycol,pch=ypch)
#   }
#   text(x=Xx*1.05,y=Xy*1.05,labels=xnames,...)
#   text(x=Yx*1.05,y=Yy*1.05,labels=ynames,...)
#}
coloredBiplot <- function(x, ...) UseMethod("coloredBiplot")

coloredBiplot.default <-
    function(x, y, var.axes = TRUE, col, cex = rep(par("cex"), 2),
	     xlabs = NULL, ylabs = NULL, expand=1, xlim = NULL, ylim = NULL,
	     arrow.len = 0.1,
             main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
             xlabs.col = NULL, xlabs.pc=NULL, ...)
{
    n <- nrow(x)
    p <- nrow(y)
    if(missing(xlabs)) {
	xlabs <- dimnames(x)[[1]]
	if(is.null(xlabs)) xlabs <- 1:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2]])
    if(missing(ylabs)) {
	ylabs <- dimnames(y)[[1]]
	if(is.null(ylabs)) ylabs <- paste("Var", 1:p)
    }
    ylabs <- as.character(ylabs)
    dimnames(y) <- list(ylabs, dimnames(y)[[2]])

    if(length(cex) == 1) cex <- c(cex, cex)
    if(missing(col)) {
	col <- par("col")
	if (!is.numeric(col)) col <- match(col, palette(), nomatch=1)
	col <- c(col, col + 1)
    }
    else if(length(col) == 1) col <- c(col, col)

    unsigned.range <- function(x)
        c(-abs(min(x, na.rm=TRUE)), abs(max(x, na.rm=TRUE)))
    rangx1 <- unsigned.range(x[, 1])
    rangx2 <- unsigned.range(x[, 2])
    rangy1 <- unsigned.range(y[, 1])
    rangy2 <- unsigned.range(y[, 2])

    if(missing(xlim) && missing(ylim))
	xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if(missing(xlim)) xlim <- rangx1
    else if(missing(ylim)) ylim <- rangx2
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    on.exit(par(op))
    op <- par(pty = "s")
    if(!is.null(main))
        op <- c(op, par(mar = par("mar")+c(0,0,1,0)))
    if(missing(xlabs.col)){
     col1 = col[1]
    }else{
     col1 = xlabs.col
    }
    plot(x, type = "n", xlim = xlim, ylim = ylim, col = col1,
         xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
    if(missing(xlabs.pc)){ 
       text(x, xlabs, cex = cex[1], col = col1, ...)
    }else{
       points(x, cex = cex[1], col = col1, pch=xlabs.pc, ...)
    }
    par(new = TRUE)
    plot(y, axes = FALSE, type = "n", xlim = xlim*ratio, ylim = ylim*ratio,
	 xlab = "", ylab = "", col = col[1], ...)
    axis(3, col = col[2], ...)
    axis(4, col = col[2], ...)
    box(col = col[1])
    text(y, labels=ylabs, cex = cex[2], col = col[2], ...)
    if(var.axes)
	arrows(0, 0, y[,1] * 0.8, y[,2] * 0.8, col = col[2], length=arrow.len)
    invisible()
}

coloredBiplot.princomp <- function(x, choices = 1:2, scale = 1, pc.biplot=FALSE, ...)
{
    if(length(choices) != 2) stop("length of choices must be 2")
    if(!length(scores <- x$scores))
	stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
             domain = NA)
    lam <- x$sdev[choices]
    if(is.null(n <- x$n.obs)) n <- 1
    lam <- lam * sqrt(n)
    if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
    if(scale != 0) lam <- lam^scale else lam <- 1
    if(pc.biplot) lam <- lam / sqrt(n)
    coloredBiplot.default(t(t(scores[, choices]) / lam),
		   t(t(x$loadings[, choices]) * lam), ...)
    invisible()
}

coloredBiplot.prcomp <- function(x, choices = 1:2, scale = 1, pc.biplot=FALSE, ...)
{
    if(length(choices) != 2) stop("length of choices must be 2")
    if(!length(scores <- x$x))
	stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
             domain = NA)
    if(is.complex(scores))
        stop("biplots are not defined for complex PCA")
    lam <- x$sdev[choices]
    n <- NROW(scores)
    lam <- lam * sqrt(n)
    if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
    if(scale != 0) lam <- lam^scale else lam <- 1
    if(pc.biplot) lam <- lam / sqrt(n)
    coloredBiplot.default(t(t(scores[, choices]) / lam),
		   t(t(x$rotation[, choices]) * lam), ...)
    invisible()
}



## a colour code generator for outlierfactors
# outfac the factor
# family a colorpalette to generate colors for general groups
colorsForOutliers1 <- function(outfac, family=rainbow,extreme="cyan",outlier="red",ok="gray40",unknown="blue"){
  lev <- levels(outfac)
  cols = do.call(family, args=list(n=length(lev[lev!="ok"])))
  cols = ifelse(lev!="ok",cols,ok)
  names(cols) <- lev
  if( "extreme" %in% lev ) {
    cols["extreme"]<-extreme
  }
  if( "outlier" %in% lev ) 
    cols["outlier"]<-outlier
  if( "?" %in% lev ) 
    cols["?"]<-unknown
  return(cols)
}

#colorsForOutliers2 
# A colorcode generator for bitcodes factors
# outfac: the factor as a result for myClassifier3.acomp( type="all")
# use   : the bits to be encoded
# codes : color components to be used
colorsForOutliers2 <- function(outfac,use=whichBits(gsi.orSum(levels(outfac))),
                        codes=c(2^outer(c(24,16,8),1:7,"-")),ok="yellow"
                        ) {
  lev<-levels(outfac)
  oks <- lev=="ok"
  if( any(oks) ) lev[oks]<-"0"
  dat <- c((bit(lev,use)+0)%*%codes[1:length(use)])
  r <- (dat %/% 2^16)
  g <- (dat %/% 2^8) %% 2^8
  b <- dat %% 2^8
  erg <- rgb(r/255,g/255,b/255)
  erg[oks]<-ok
  erg
}

## a colour code generator for outlierfactors
# outfac : the factor
# ok : pch for good measurements
# outlier : pch for class outlier
# extreme : pch for class extreme
# unkown  : pch for class ?
# ...     : pairs class=pch for assigning pch to another class
# other   : possible pch for unmentioned classes
pchForOutliers1 = function(outfac,ok='.',outlier='\004',extreme='\003',unknown='\004',...,
  other=c('\001','\002','\026','\027','\010','\011','\012','\013','\014','\015',
    '\016',strsplit("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ","")[[1]])
    ) {
  tbl <- c(ok=ok,outlier=outlier,extreme=extreme,"?"=unknown,...)
  lev <- levels(outfac)
  nrm <- ! (lev %in% names(tbl))
  if( any(nrm) ) 
    tbl <- c(tbl,structure(other[1:sum(nrm)],names=lev[nrm]))
  return(tbl[as.character(lev)])
}



outlierplot <- function(X,...) UseMethod("outlierplot",X)

## several diagnostic plots and their auxiliary functions
## general function to obtain exploratory plots on the outlier character of the sample
outlierplot.acomp <- function(X,colcode=colorsForOutliers1,pchcode=pchForOutliers1,
  type=c("scatter","biplot","dendrogram","ecdf","portion","nout","distdist"),
  legend.position,pch=19,...,clusterMethod="ward",
  myCls=classifier(X,alpha=alpha,type=class.type,corrected=corrected),
  classifier=OutlierClassifier1,
  alpha=0.05,
  class.type="best",
  Legend,pow=1,
  main=paste(deparse(substitute(X))),
  corrected=TRUE,robust=TRUE,princomp.robust=FALSE,
                              mahRange=exp(c(-5,5))^pow,
                              flagColor="red",
                              meanColor="blue",
                              grayColor="gray40",
                              goodColor="green",
                              mahalanobisLabel="Mahalanobis Distance"
                              ){
  cl=match.call()
  what = match.arg(type)
  if(what=="biplot"){
    if( is.function(colcode) )
      colcode <- colcode(myCls)
    if( is.function(pchcode) ) {
      pchcode <- pchcode(myCls)
    }
    if( length(pchcode) > 1 )
      pch <- pchcode[as.integer(myCls)]
    if( is.logical(princomp.robust) )
      pr <- princomp(cdt(X),robust=princomp.robust)
    else
      pr <- princomp.robust
    coloredBiplot(x=pr,xlabs.col=colcode[as.integer(myCls)],pc.biplot=FALSE,xlabs.pc=pch,...)
    title(main=main)
    if(!missing(legend.position))
      legend(legend.position,legend=levels(myCls),col=colcode,pch=pch)
    erg<-list(call=cl,colcode=colcode,classes=myCls,princomp=pr)
 }

##########################################
  ### corrected influence mahalanobis
# infl.mahalanobis.acomp = function(X,m=cov.mcd(idt(X))$center, v=cov.mcd(idt(X))$cov ){
#  aux = svd(v)
#  st = clr(acomp(X)-ilr.inv(m)) %*% t( ilr2clr(t(aux$v))) %*% diag(1/sqrt(aux$d))
#  st = ilr2clr(st)
#  colnames(st)=colnames(X)
#  return(st)
# }
  infl.mahalanobis = function(X,center=attr(cov,"center"), cov=var(X,robust=robust,giveCenter=TRUE), cond = 1e-10 ,robust=TRUE){
    aux = svd(cov)  # decompose Sigma = L * Lt in clrs
    nd = sum(aux$d>cond*max(aux$d))
    Linv = diag(sqrt(1/aux$d[1:nd])) %*% t(aux$v[,1:nd])
    st = cdt(X-center)
    st = unclass(st) %*% t(Linv) %*% t(aux$v[,1:nd])
 colnames(st)=colnames(X)
 return(st)
}
#####################################
# if(what=="barplot"){ 
#  st = infl.mahalanobis(acomp(X))
#  if( is.function(colcode) )
#    colcode <- colcode(rownames(st))
#  erg <- t(st)^2/gsi.getD(st)
#  barplot(erg,col=colcode,names.arg=1:nrow(X),...,main=main)
#  abline(h=cutlevel)
#  if(!missing(legend.position))
#    legend(legend.position,legend=colnames(X),fill=colcode)
#  erg<-invisible(list(call=cl,colcode=colcode,erg=erg,outside=erg>=cutlevel,cutlevel=cutlevel))
# }
 if(what=="scatter"){ 
#   myCls = classifier(X, cutlevel=cutlevel,alpha=alpha)
    if( is.function(colcode) )
      colcode <- colcode(myCls)
    if( is.function(pchcode) ) {
      pchcode <- pchcode(myCls)
    }
    if( length(pchcode) > 1 )
      pch <- pchcode[as.integer(myCls)]
   plot(X,col=colcode[as.numeric(myCls)],pch=pch,...,main=main)
   erg<-invisible(list(call=cl,colcode=colcode,pchcode=pchcode,classes=myCls))
 }
  if(what=="dendrogram"){
#  myCls = classifier(X, cutlevel=cutlevel,alpha=alpha)
    hc = hclust(MahalanobisDist(X,pairwise=TRUE,robust=robust,pow=pow),method=clusterMethod)
    plot(hc,labels=myCls,...,main=main)
    erg<-invisible(list(call=cl,clust=hc,classes=myCls))
  }
  if(what=="ecdf") {
    K <- sort(mahalanobis<-MahalanobisDist(X,robust=robust,pow=pow))
    Ks <- exp(seq(log(min(mahRange)),log(max(mahRange)),length.out=100))
    cdf <- ecdf(K)
    plot(1,1,type="n",main=main,log="x",xlim=mahRange,ylim=c(0,1),
         xlab=mahalanobisLabel,ylab="",...)
    plot(cdf,main=main,add=TRUE,xlim=mahRange,xlab="Mahalanobis (log)")
    lines(Ks,pEmpiricalMahalanobis(Ks,N=nrow(X),d=gsi.getD(X)-1,pow=pow),col=meanColor,lwd=2)
    KpQ <- pQuantileMahalanobis(Ks,N=nrow(X),d=gsi.getD(X)-1,p=alpha,pow=pow)
    lines(Ks,KpQ,col=goodColor,lwd=2)
                                        #  print(ifelse(cdf(K) < KpQ,1-KpQ/(1-cdf(K)),0))
    lines(Ks,ifelse(cdf(Ks) < KpQ,1-(1-KpQ)/(1-cdf(Ks)),0),col=grayColor) 
    if( max(K) > max(Ks) ) {
      warning("outlierplot: Extrem outliers found, extend mahRange!")
      segments(max(Ks),cdf(max(Ks)),max(Ks),1.0,col=flagColor)
    }
    erg<-invisible(list(call=cl,mahalanobis=mahalanobis,Ks=Ks,cdf=cdf,KpQ=KpQ,K=K))
  }
  if(what=="portion") {
      K <- sort(mahalanobis<-MahalanobisDist(X,robust=robust,pow=pow))
      Ks <- exp(seq(log(min(mahRange)),log(max(mahRange)),length.out=100))
      n <- nrow(X)
      cdf <- ecdf(K) 
      Kp  <- pEmpiricalMahalanobis(Ks,N=nrow(X),d=gsi.getD(X)-1,robust=robust)
      KpQ  <- pQuantileMahalanobis(Ks,N=nrow(X),d=gsi.getD(X)-1,p=alpha,robust=robust)
      plot(1,1,type="n",xlim=mahRange,ylim=n*range(c(Kp-cdf(Ks),Kp-KpQ),0),log="x",main=main,xlab="Mahalanobis distance",ylab="Number of Outliers",...)
      lines(Ks,n*(Kp-cdf(Ks)),col=grayColor,...)
      lines(Ks,n*(Kp-KpQ),col=goodColor,lty=1,...)
      stripchart(K,method="stack",add=TRUE,...)    
      abline(h=0)
      if( max(K) > max(Ks) ) {
        warning("outlierplot: Extrem outliers found, extend mahRange!")
        segments(max(Ks),0,max(Ks),sum(K>max(Ks)),col="red")
      }
      erg<-invisible(list(call=cl,mahalanobis=mahalanobis,Ks=Ks,cdf=cdf,Kp=Kp,KpQ=KpQ,K=K))
    }
  if(what=="nout") {
###############
# MinOutlierBound computes the an alpha confidence bound of outliers above a given limit
# it returns matrix with the columns:
# x    = The Mahalanobisdistances of the limit
# nmin = The Minium number of outliers above that limit
# por  = The Minimum portion of outliers above that limit
# i    = the id of the observation with that distance
# invOr= the sequece position of the given observation
# if bestOnly is TRUE 
    MinOutlierBound <- function(X,alpha=0.05,...,bestOnly=TRUE,robust=TRUE,mah=MahalanobisDist(X,...,robust=robust),pow=1) {
      oriMah <- mah
      or  <- order(mah)
      K  <-  mah[or]
      Ks  <- mah[or]
      cdf <- ecdf(K) 
      KpQ <- pQuantileMahalanobis(c(Ks[-1],rev(Ks)[1]),N=nrow(X),d=gsi.getD(X)-1,p=alpha,...,robust=robust)
      por <- ifelse(cdf(Ks) < KpQ,1-(1-KpQ)/(1-cdf(Ks)),0)
      invOr <- 1:nrow(X)
      invOr[or] <- 1:nrow(X)
      n     <- (nrow(X)-1):0
      nmin  <- por*n
      falseClass <- max(nmin)-nmin+(n-nmin)
      data.frame(x=Ks,mahalanobis=mah,nmin=nmin,por=por,i=or,invOr=invOr,falseClass=falseClass)
    }

    ################
    bnd <- MinOutlierBound(X,alpha=alpha,robust=robust,pow=pow)
    plot(bnd$x,bnd$nmin,type="s",log="x",...,lwd=2,xlab="Mahalanobis Distance",ylab="n outliers")
    points(bnd$x,rep(0,length(bnd$x)),pch="|",col="red")
    #lines(bnd$x,bnd$por*max(bnd$nmin),col="gray40",type="s")
    erg <- bnd
  }
  if(what == "distdist" ) {
    if( is.function(colcode) )
      colcode <- colcode(myCls)
    if( is.function(pchcode) ) {
      pchcode <- pchcode(myCls)
    }
    if( length(pchcode) > 1 )
      pch <- pchcode[as.integer(myCls)]
   NormalMahalanobisDist <- MahalanobisDist(X,robust=FALSE)
   RobustMahalanobisDist <- MahalanobisDist(X,robust=TRUE)
   plot(NormalMahalanobisDist,RobustMahalanobisDist,
        col=colcode[as.numeric(myCls)],pch=pch,...,main=main)
   crit1<-qMaxMahalanobis(0.95,nrow(X),ncol(idt(X)),...,replicates=1000,robust=TRUE)
   if(is.finite(crit1)) abline(h=crit1,lty=3)
   erg<-invisible(list(call=cl,colcode=colcode,pchcode=pchcode,classes=myCls,NormalMahalanobisDist=NormalMahalanobisDist,RobustMahalanobisDist=RobustMahalanobisDist,crit=crit1))
  }
  if( !missing(Legend) ) {
    es <- substitute(Legend)
    frame <- sys.frame(sys.nframe())
    erg$legend <- eval(es,frame)
  }
  invisible(erg)
}

## influence of each part on the global mahalanobis distance
# if we extract R parts from a composition, the rest of the parts have 
#    clr(x*)=clr(x)+sum(clr(x_R))/(D-R)
#    and the sum of its squares is:
#    clr(x*)^2=clr(x)^2- [ sum(clr(x_R)) ]^2/(D-R)
# so, the Aitchison norms (without division by D) before and after are related
#   A(x*) = A(x|x_R = 0) - [ sum(clr(x_R)) ]^2
# ... and dividing by the number of parts
#   A'(x*) = D/(D-R) * A'(x|x_R = 0) -  [ sum(clr(x_R)) ]^2/(D-R)
# this would translate to Mahalanobis distances iff the st(clr(x)[,-parts]) =_A st(clr(x[,-parts]))
#   but it is not true. (say, inversion of Sigma is not subcompositionally coherent, 
#            as it is not marginally coherent in R).
#infl.mahalanobis.acomp = function(x,m=covMcd(idt(x))$center, v=covMcd(idt(x))$cov ){
# aux = svd(v)
# st = clr(acomp(x)-ilr.inv(m)) %*% ilrvar2clr(aux$v %*% diag(1/sqrt(aux$d)))
 #colnames(st)=colnames(x)
# return(st)
#}



## routine to extract the biggest typical subcomposition(s) for each individual
#Problems: cutlevel not defined, combinat necessary
#typicalsubcomp <- function(X,robust=robust){ 
# # global mahalanobis distance of the data set
# mah = MahalanobisDist(X,robust=robust)
# # encode all subcompositions
# require("combinat")
# parts = colnames(X)
# partlist = sapply(2:length(parts),function(i){combn(x=length(parts),m=i, simplify=TRUE) })
# name = paste(paste(parts[partlist[[length(partlist)]] ],collapse="+"))
# # run through all subcompositions to obtain the mahalanobis distance 
# #     on each subcomposition between each individual and the average
# k = length(partlist[[length(partlist)]])
# for(i in (length(partlist)-1):1){
#   for(j in 1:ncol(partlist[[i]]) ){
#     aux = MahalanobisDist.acomp(acomp(X[, partlist[[i]][,j] ]),robust=robust)
#     mah = cbind(mah,aux)
#     name = c(name, paste(paste(parts[partlist[[i]][,j] ],collapse="+")) )
#     k = c(k, length( partlist[[i]][,j]) )
#   }
# }
# colnames(mah)=name
# # get the number of parts of the biggest typical subcompositions
# biggerD = sapply(1:nrow(mah),function(i){max(k[mah[i,]<cutlevel])})
# # get the biggest typical subcompositions (not necessarily one)
# biggersubs = lapply(1:nrow(mah),function(i){
#    name[(mah[i,]<cutlevel)&(k==biggerD[i])]
# })
# output = biggersubs
#   attr(output,"parts")=biggerD
# return(output)
#}





# The inner product where the * and the sum may be replaced
# works much like outer()
# x: the first matrix
# y: the second matrix
# binary: the binary operation replacing *
# simp  : the summary operation replacing the sum
#inner <- function(x,y,binary="*",simp=sum) {
#  print("Inner used")
#  if( length(dim(y))> 1) {
#    sapply(1:nrow(y),function(i) Recall(x,y[,j],binary,simp))
#  } else if( length(dim(x))>1 ) {
#    c(sapply(1:nrol(x),function(i) Recall(x[i,],y,binary,simp)))
#  }  else {
#    do.call(simp,do.call(binary,list(x,y)))
#  }
#}

# compositional random normals (what is the problem with the function in the
# package ????
#rnorm.acomp <-function (n, mean, var){
#    D <- gsi.getD(mean)-1
#    perturbe(ilr.inv(matrix(rnorm(n*D), ncol = D) %*% chol(clrvar2ilr(var))), mean)
#}
# THIS IS THE ORIGINAL VERSION:
#rnorm.acomp -> function (n, mean, var)
#{
#    D <- NCOL(oneOrDataset(mean))
#    perturbe(ilr.inv(matrix(rnorm(n * length(mean)), ncol = D -
#        1) %*% chol(clrvar2ilr(var))), mean)
#}



## Dendrogram of outliering character. Ideas:
## 0.- global dendrogram: normalize the data set via Sigma^{-1/2}, cluster it with complete linkeage

# X=aux
#   aux2 = infl.mahalanobis.acomp(X)
#   hc = hclust(dist(aux2),method="complete")
#   plot(hc,labels=myCls)


## 1.- compute, for each lr, the Aitchison-Mahalanobis norm of the datum; use as distance matrix

#myDendro = function(X){
#onelr = function(x){
# D = length(x)
# aux = outer(1:D,1:D,function(i,j){log(x[i]/x[j])})
# return(aux)
#}

#D = gsi.getD(X)
#mm = mapply(function(i,j,X){median(log(X[,i]/X[,j]))},
#      i=rep(1:D,times=D), j=rep(1:D,each=D), MoreArgs=list(X=X))
#  dim(mm)=c(D,D)
#vr = mapply(function(i,j,X){mad(log(X[,i]/X[,j]))},
#      i=rep(1:D,times=D), j=rep(1:D,each=D), MoreArgs=list(X=X))
#  dim(vr)=c(D,D)
#  vr = 1/vr
#  vr[is.nan(vr)]=0

#output = list()
#for(k in c(1:gsi.getN(X))){
# dd = (onelr(X[k,])-mm)^2*vr
#  colnames(dd)=colnames(X)
#  rownames(dd)=colnames(X)
#  dd =  hclust(as.dist(dd),method="complete")
#  output[[k]]=dd
#}
#return(output)
#}



### robust svd (svd based )
#svd.mcd = function(X,covmcd=NULL){
# noms = colnames(X) # keep names to track them afterwards
# svd1 = svd(cdt(X)) # pattern of svd for clr (it will be mostly replaced)
# # cov.mcd cannot handle singular covariance matrices
# X = idt(X)
# if(is.null(covmcd)){covmcd = covMcd(X)}
# # represent the matrix of the covariance in ilrs as in clrs
# cv = ilrvar2clr(covmcd$cov)
# # extract the svd
# svd2 = svd(cv)
# # eigenvalues should be squared
# d = sqrt(svd2$d)
# # right eigenvectors are equal
# v = svd2$v
#    rownames(v)=noms
# # compute left eigenvectors of the original data set by projection/inversion
# mm = ilr2clr(covmcd$center) # robust mean, as clr
# X = ilr2clr(X)
# X = X-mm # double centered clrs
# u = X %*% v %*% diag(c(1/d))
#    u = t(unclass(normalize(rmult(t(u)))))
#    u[is.nan(u)]=1/sqrt(nrow(u))
# st = list(d=d, u=u, v=v)
# return(st)
#}





# Computes the critical Mahalanobis distances values for outliers based on
# simulations
# q: the quantile
# ntrials: the number of simulations
# n: The number of cases
# d: the number of degrees of freedom (i.e. D-1 in the simplex)
### Probabily outdated due to qMaxMahalanobis
#computeOutlierKrit <- function(p=0.95,ntrial=10000,n,d,robust=TRUE,...) {
#  rr <- robust
#  if( is.logical(robust) ) rr <- if(robust) "mcd" else "pearson"
#  s<-switch(rr,
#            pearson={
#              if( d==1 )
#                replicate(ntrial,max({
#                  x <- rnorm(n,mean=rnorm(1))
#                  od<-covMcd(x)
#                  max(c(((c(x)-mean(x))^2)/c(var(x))))
#                }))
#              else {
#                quickMah <- function(x) {
#                  v <- var(x)
#                  xc <- t(x)-meanCol(x)
#                  max(rep(1,ncol(v)) %*% solve(v,xc,...) * xc )
#                  }
#                replicate(ntrial,quickMah(matrix(rnorm(n*d,mean=rnorm(1)),ncol=d)))
#              }
#            }
#            ,
#            mcd={
#              require("robustbase")
#              if( d==1 )
#                replicate(ntrial,max({
#                  x <- rnorm(n,mean=rnorm(1))
#                  od<-covMcd(x)
#                  c(((c(x)-od$center)^2)/c(od$cov))
#                }))
#              else 
#                replicate(ntrial,sqrt(max(covMcd(matrix(rnorm(n*d,mean=rnorm(1)),ncol=d))$mah)))
#            },
#  stop("computeOutlierKrit: Unkown robustness type: \"",robust,"\"")
#  )
#quantile(s,p)
#}
# The same as above but in 1 dimension
# --> Above version generalized
#computeOutlierKrit1 <- function(q=0.95,ntrial=10000,n) {
#  s <- replicate(ntrial,max({
#    x <- rnorm(n,mean=rnorm(1))
#    od<-covMcd(x)
#    c(((c(x)-od$center)^2)/c(od$cov))
#  }))
#    quantile(s,q)
#}



#dendroSimp <- function(dendro,labels) {
#  
#  erg <- structure(lapply(dendro,dendroSimp,labels=labels))
#  leaves <- sapply(erg,function(x) attr(x,"leaf"))
#  if( all(leaves) ) {
#    labs <- sapply(erg,function(x) label
 # }
#
#
#}#

# Function should extract only big groups from a factor (now working)
#bigGroupsOnly <- function(x,nmin=5,name="single") {
#  n <- c(table(x)[x])
#  factor(ifelse(n<nmin,name,as.character(x)))
#}

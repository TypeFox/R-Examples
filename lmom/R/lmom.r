
#***********************************************************************
#*                                                                     *
#*  R code for R package "lmom"                                        *
#*                                                                     *
#*  J. R. M. HOSKING <jrmhosking@gmail.com>                                                  *
#*                                                                     *
#*  Version 2.5    February 2015                                       *
#*                                                                     *
#***********************************************************************

#-------------------------------------------------------------------------------
#  Specifications for distributions, used by cdf...(), qua...(), lmr...(), pel...()
#-------------------------------------------------------------------------------

# 'lmom.dist' contains specifications for each distribution that the "lmom" package knows about
#
# name        - 3-character abbreviation for the distribution
# npara       - number of parameters
# parnames    - names of parameters
# pardefaults - default values of parameters
# himom       - highest order L-moment needed for parameter estimation
# maxmom      - highest order L-moment that can be computed by the LMR... Fortran routine
#
lmom.dist<-list(
  exp=list(name="exp", npara=2L, parnames=c("xi","alpha"),                        pardefaults=c(0,1),       himom=2L, maxmom=20L),
  gam=list(name="gam", npara=2L, parnames=c("alpha","beta"),                      pardefaults=c(1,1),       himom=2L, maxmom= 4L),
  gev=list(name="gev", npara=3L, parnames=c("xi","alpha","k"),                    pardefaults=c(0,1,0),     himom=3L, maxmom=20L),
  glo=list(name="glo", npara=3L, parnames=c("xi","alpha","k"),                    pardefaults=c(0,1,0),     himom=3L, maxmom=20L),
  gno=list(name="gno", npara=3L, parnames=c("xi","alpha","k"),                    pardefaults=c(0,1,0),     himom=3L, maxmom=20L),
  gpa=list(name="gpa", npara=3L, parnames=c("xi","alpha","k"),                    pardefaults=c(0,1,0),     himom=3L, maxmom=20L),
  gum=list(name="gum", npara=2L, parnames=c("xi","alpha"),                        pardefaults=c(0,1),       himom=2L, maxmom=20L),
  kap=list(name="kap", npara=4L, parnames=c("xi","alpha","k","h"),                pardefaults=c(0,1,0,0),   himom=4L, maxmom=20L),
  ln3=list(name="ln3", npara=3L, parnames=c("zeta","mu","sigma"),                 pardefaults=c(0,0,1),     himom=3L, maxmom=20L),
  nor=list(name="nor", npara=2L, parnames=c("mu","sigma"),                        pardefaults=c(0,1),       himom=2L, maxmom=20L),
  pe3=list(name="pe3", npara=3L, parnames=c("mu","sigma","gamma"),                pardefaults=c(0,1,0),     himom=3L, maxmom= 4L),
  wak=list(name="wak", npara=5L, parnames=c("xi","alpha","beta","gamma","delta"), pardefaults=c(0,1,0,0,0), himom=5L, maxmom=20L),
  wa0=list(name="wa0", npara=5L, parnames=c("xi","alpha","beta","gamma","delta"), pardefaults=c(0,1,0,0,0), himom=4L, maxmom=20L),
  wei=list(name="wei", npara=3L, parnames=c("zeta","beta","delta"),               pardefaults=c(0,1,0),     himom=3L, maxmom=20L)
)

#-------------------------------------------------------------------------------
#  Distributions: cdf...(), qua...(), lmr...(), pel...()
#-------------------------------------------------------------------------------

cdfexp<-function(x,para=c(0,1)){
  if (length(para)!=2) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  return(1-exp(-(pmax(0,x-para[1]))/para[2]))
}

cdfgam<-function(x,para=c(1,1)){
  if (length(para)!=2) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (any(para<=0)) stop("distribution parameters invalid")
  result<-pgamma(x/para[2],para[1])
  return(result)
}

cdfgev<-function(x,para=c(0,1,0)){
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (para[3]==0) y<-(x-para[1])/para[2]
#                        pmax so that values outside the range
#                             generate 0s or 1s rather than NAs
  else y<--1/para[3]*log(pmax(0,1-para[3]*(x-para[1])/para[2]))
  return(exp(-exp(-y)))
}

cdfglo<-function(x,para=c(0,1,0)){
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (para[3]==0) y<-(x-para[1])/para[2]
#                        pmax so that values outside the range
#                             generate 0s or 1s rather than NAs
  else y<--1/para[3]*log(pmax(0,1-para[3]*(x-para[1])/para[2]))
  return(1/(1+exp(-y)))
}

cdfgno<-function(x,para=c(0,1,0)){
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (para[3]==0) y<-(x-para[1])/para[2]
#                        pmax so that values outside the range
#                             generate 0s or 1s rather than NAs
  else y<--1/para[3]*log(pmax(0,1-para[3]*(x-para[1])/para[2]))
  return(pnorm(y))
}

cdfgpa<-function(x,para=c(0,1,0)){
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (para[3]==0) y<-(x-para[1])/para[2]
  else y<--1/para[3]*log(pmax(0,1-para[3]*(x-para[1])/para[2]))
  return(1-exp(-pmax(y,0)))
}

cdfgum<-function(x,para=c(0,1)){
  if (length(para)!=2) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  return(exp(-exp(-(x-para[1])/para[2])))
}

cdfkap<-function(x,para=c(0,1,0,0)){
  if (length(para)!=4) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  y<-(x-para[1])/para[2]
  if (para[3]!=0) y<--1/para[3]*log(pmax(0,1-para[3]*y))
  y<-exp(-y)
  if (para[4]!=0) y<--1/para[4]*log(pmax(0,1-para[4]*y))
  y<-exp(-y)
  return(y)
}

cdfnor<-function(x,para=c(0,1)){
  if (length(para)!=2) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  return(pnorm(x,para[1],para[2]))
}

cdfpe3<-function(x,para=c(0,1,0)){
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (abs(para[3])<=1e-6) return(pnorm(x,para[1],para[2]))
  alpha<-4/para[3]^2
  z<-2*(x-para[1])/(para[2]*para[3])+alpha
  result<-pgamma(pmax(0,z),alpha)
  if (para[3]<0) result<-1-result
  return(result)
}

cdfwak<-function(x,para=c(0,1,0,0,0)){
  if (length(para)!=5) stop("parameter vector has wrong length - should be 5")
  if (any(is.na(para))) stop("missing values in parameter vector")
  ok<-is.finite(x)
  xok<-x[ok]
  if (length(xok)==0) xok<-para[1]
  nxok<-length(xok)
  fort<-.Fortran("cdfwak",PACKAGE="lmom",
        as.double(xok),
        as.integer(nxok),
        as.double(para),
        cdf=double(nxok),
        ifail=integer(nxok))
  if (fort$ifail[1]==7000) stop("distribution parameters invalid")
  if (any(fort$ifail==7010)) warning("iteration did not converge - some cdf values not computed")
  result<-rep(as.numeric(NA),length(x))
  result[x==-Inf]<-0
  result[x==Inf]<-1
  result[ok]<-fort$cdf
  return(result)
}

quaexp<-function(f,para=c(0,1)){
  if (length(para)!=2) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  result<-para[1]+para[2]*(-log(1-f))
  return(result)
}

quagam<-function(f,para=c(1,1)){
  if (length(para)!=2) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (any(para<=0)) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  result<-para[2]*pmax(0,qgamma(f,para[1]))
  return(result)
}

quagev<-function(f,para=c(0,1,0)){
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  result<-
    if (para[3]==0) para[1]-para[2]*log(-log(f))
    else para[1]+para[2]/para[3]*(1-(-log(f))^para[3])
  return(result)
}

quaglo<-function(f,para=c(0,1,0)){
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  result<-
    if (para[3]==0) para[1]+para[2]*(log(f/(1-f)))
    else para[1]+para[2]/para[3]*(1-((1-f)/f)^para[3])
  return(result)
}

quagno<-function(f,para=c(0,1,0)){
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  result<-
    if (para[3]==0) para[1]+para[2]*qnorm(f)
    else para[1]+para[2]/para[3]*(1-exp(-qnorm(f)*para[3]))
  return(result)
}

quagpa<-function(f,para=c(0,1,0)){
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  result<-
    if (para[3]==0) para[1]+para[2]*(-log(1-f))
    else para[1]+para[2]/para[3]*(1-(1-f)^(para[3]))
  return(result)
}

quagum<-function(f,para=c(0,1)){
  if (length(para)!=2) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  result<-para[1]-para[2]*log(log(1/f))
  return(result)
}

quakap<-function(f,para=c(0,1,0,0)){
  if (length(para)!=4) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  f <- if (para[4]==0) (-log(f)) else (1-f^para[4])/para[4]
  f <- if (para[3]==0) (-log(f)) else (1-f^para[3])/para[3]
  result<-para[1]+para[2]*f
  if (para[4]<=0 & para[3]<=0) result[f==1]<-Inf # not -Inf!
  return(result)
}

quanor<-function(f,para=c(0,1)){
  if (length(para)!=2) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  result<-qnorm(f,para[1],para[2])
  return(result)
}

quape3<-function(f,para=c(0,1,0)){
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  if (abs(para[3])<=1e-8) return(qnorm(f,para[1],para[2]))
  alpha<-4/para[3]^2
  beta<-abs(0.5*para[2]*para[3])
  result<-
    if (para[3]>0) para[1]-alpha*beta+beta*pmax(0,qgamma(f,alpha))
    else           para[1]+alpha*beta-beta*pmax(0,qgamma(1-f,alpha))
  return(result)
}

quawak<-function(f,para=c(0,1,0,0,0)){
  if (length(para)!=5) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (lmom.parok.wak(para)==FALSE) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  flog<- -log(1-f)
  result<-para[1]+
    (if (para[2]==0) 0 else para[2]*(if (para[3]==0) flog else  (1-exp(-para[3]*flog))/para[3]))+
    (if (para[4]==0) 0 else para[4]*(if (para[5]==0) flog else -(1-exp( para[5]*flog))/para[5]))
  return(result)
}

# lmom.parok.wak -- test validity of Wakeby distribution parameters
lmom.parok.wak<-function(para) {
  alpha<-para[2]
  beta<-para[3]
  gamma<-para[4]
  delta<-para[5]
# first a quick check for the "usual" case (does this help?)
  if (gamma>0 && alpha+gamma>0 && beta+delta>0) return(TRUE)
  if (gamma<0) return(FALSE)
  if (alpha+gamma<0) return(FALSE)
  if (beta+delta<=0 && any(c(beta,gamma,delta)!=0)) return(FALSE)
  if (alpha==0 && beta!=0) return(FALSE)
  if (gamma==0 && delta!=0) return(FALSE)
  return(TRUE)
}

# 'para' default values are from lmom.dist$...$pardefaults
# 'nmom' default values are from lmom.dist$...$himom
lmrexp<-function(para=c(0,1),      nmom=2) lmrxxx("exp",para,nmom)$lmom
lmrgam<-function(para=c(1,1),      nmom=2) lmrxxx("gam",para,nmom)$lmom
lmrgev<-function(para=c(0,1,0),    nmom=3) lmrxxx("gev",para,nmom)$lmom
lmrglo<-function(para=c(0,1,0),    nmom=3) lmrxxx("glo",para,nmom)$lmom
lmrgno<-function(para=c(0,1,0),    nmom=3) lmrxxx("gno",para,nmom)$lmom
lmrgpa<-function(para=c(0,1,0),    nmom=3) lmrxxx("gpa",para,nmom)$lmom
lmrgum<-function(para=c(0,1),      nmom=2) lmrxxx("gum",para,nmom)$lmom
lmrkap<-function(para=c(0,1,0,0),  nmom=4) lmrxxx("kap",para,nmom)$lmom
lmrnor<-function(para=c(0,1),      nmom=2) lmrxxx("nor",para,nmom)$lmom
lmrpe3<-function(para=c(0,1,0),    nmom=3) lmrxxx("pe3",para,nmom)$lmom
lmrwak<-function(para=c(0,1,0,0,0),nmom=5) lmrxxx("wak",para,nmom)$lmom

lmrxxx<-function(xxx,para,nmom){
  ddata<-lmom.dist[[xxx]]
  if (missing(para)) para<-ddata$pardefaults
  if (missing(nmom)) para<-ddata$npara
  if (length(para)!=ddata$npara) stop("lmr",ddata$name,": parameter vector has wrong length - should be ",ddata$npara)
  if (any(is.na(para))) stop("lmr",ddata$name,": missing values in parameter vector")
  nnmom<-as.integer(nmom)
  if (nnmom<=0L) return(numeric(0))
  if (nnmom<=2L) rnames <- paste("lambda",1:nnmom,sep="_")
  else rnames <- c("lambda_1","lambda_2",paste("tau",3:nnmom,sep="_"))
  result <- rep(as.numeric(NA),nnmom)
  names(result) <- rnames
  if(nnmom>ddata$maxmom) {
    warning("lmr",ddata$name,": too many L-moments requested - max. is ",ddata$maxmom)
    nnmom<-min(nnmom,ddata$maxmom)
  }
  para<-as.double(para)
  xmom<-double(nnmom)
  nnmom<-as.integer(nnmom)
  ifail<-integer(1)
  fort<-switch(xxx,
    exp=.Fortran("lmrexp", PACKAGE="lmom", para, xmom=xmom, nnmom, ifail=ifail),
    gam=.Fortran("lmrgam", PACKAGE="lmom", para, xmom=xmom, nnmom, ifail=ifail),
    gev=.Fortran("lmrgev", PACKAGE="lmom", para, xmom=xmom, nnmom, ifail=ifail),
    glo=.Fortran("lmrglo", PACKAGE="lmom", para, xmom=xmom, nnmom, ifail=ifail),
    gno=.Fortran("lmrgno", PACKAGE="lmom", para, xmom=xmom, nnmom, ifail=ifail),
    gpa=.Fortran("lmrgpa", PACKAGE="lmom", para, xmom=xmom, nnmom, ifail=ifail),
    gum=.Fortran("lmrgum", PACKAGE="lmom", para, xmom=xmom, nnmom, ifail=ifail),
    kap=.Fortran("lmrkap", PACKAGE="lmom", para, xmom=xmom, nnmom, ifail=ifail),
    nor=.Fortran("lmrnor", PACKAGE="lmom", para, xmom=xmom, nnmom, ifail=ifail),
    pe3=.Fortran("lmrpe3", PACKAGE="lmom", para, xmom=xmom, nnmom, ifail=ifail),
    wak=.Fortran("lmrwak", PACKAGE="lmom", para, xmom=xmom, nnmom, ifail=ifail)
  )
  if (fort$ifail==7000) stop("lmr",ddata$name,": distribution parameters invalid")
  if (fort$ifail==7020) stop("lmr",ddata$name,": calculations of L-moments have broken down")
  if (fort$ifail>=7100) {
    nnmom<-fort$ifail-7100
    warning("lmr",ddata$name,": iterations have not converged; only ",nnmom," L-moments calculated")
  }
  result[1:nnmom]<-fort$xmom
  return(list(lmom=result,ifail=fort$ifail))
}

pelexp<-function(lmom) pelxxx("exp",lmom)$para
pelgam<-function(lmom) pelxxx("gam",lmom)$para
pelgev<-function(lmom) {
  pel<-pelxxx("gev",lmom)
  if (pel$ifail==7020) warning("iteration did not converge; results may be unreliable")
  return(pel$para)
}
pelglo<-function(lmom) pelxxx("glo",lmom)$para
pelgno<-function(lmom) {
  pel<-pelxxx("gno",lmom)
  if (pel$ifail==7010) stop("|tau_3| too large -- must be less than 0.95")
  return(pel$para)
}
pelgpa<-function(lmom, bound=NULL) {
  if (is.null(bound)) return(pelxxx("gpa",lmom)$para)
  if (length(lmom)<2) stop("with 'bound' specified, need at least 2 L-moments")
  if (any(is.na(lmom[1:2]))) stop("missing values in L-moment vector")
  lam1<-lmom[1]-bound
  if (lam1<=0 || lmom[2]<=0 || lmom[2]>=lam1) stop("L-moments invalid")
  k<-lam1/lmom[2]-2
  out<-c(bound, (1+k)*lam1, k)
  names(out)<-lmom.dist$gpa$parnames
  return(out)
}
pelgum<-function(lmom) pelxxx("gum",lmom)$para
pelkap<-function(lmom) {
  pel<-pelxxx("kap",lmom)
  if (pel$ifail==2) stop("L-moments not consistent with any kappa distribution")
  if (pel$ifail==3) warning("iteration did not converge; results may be unreliable")
  if (pel$ifail>=4) stop("unable to compute parameters owing to numerical problems")
  return(pel$para)
}
pelnor<-function(lmom) pelxxx("nor",lmom)$para
pelpe3<-function(lmom) pelxxx("pe3",lmom)$para
pelwak<-function(lmom,bound=NULL,verbose=FALSE) {
  if (is.null(bound)) pel<-pelxxx("wak",lmom)
  else {
    lmom[1]<-lmom[1]-bound
    pel<-pelxxx("wa0",lmom)
    pel$para[1]<-pel$para[1]+bound
  }
  if (verbose) {
    if (pel$ifail==1) warning("estimates could be obtained only by fitting a generalized Pareto distribution")
  }
  return(pel$para)
}

pelxxx<-function(xxx,lmom){
  ddata<-lmom.dist[[xxx]]
  himom<-ddata$himom
  if (length(lmom)<himom) stop("pel",ddata$name,": need at least ",himom," L-moments")
  lmom<-as.double(lmom[1:himom])
  if (any(is.na(lmom))) stop("pel",ddata$name,": missing values in L-moment vector")
  para<-double(ddata$npara)
  ifail<-integer(1)
  fort<-switch(xxx,
    exp=.Fortran("pelexp", PACKAGE="lmom", lmom, para=para, ifail=ifail),
    gam=.Fortran("pelgam", PACKAGE="lmom", lmom, para=para, ifail=ifail),
    gev=.Fortran("pelgev", PACKAGE="lmom", lmom, para=para, ifail=ifail),
    glo=.Fortran("pelglo", PACKAGE="lmom", lmom, para=para, ifail=ifail),
    gno=.Fortran("pelgno", PACKAGE="lmom", lmom, para=para, ifail=ifail),
    gpa=.Fortran("pelgpa", PACKAGE="lmom", lmom, para=para, ifail=ifail),
    gum=.Fortran("pelgum", PACKAGE="lmom", lmom, para=para, ifail=ifail),
    kap=.Fortran("pelkap", PACKAGE="lmom", lmom, para=para, ifail=ifail),
    nor=.Fortran("pelnor", PACKAGE="lmom", lmom, para=para, ifail=ifail),
    pe3=.Fortran("pelpe3", PACKAGE="lmom", lmom, para=para, ifail=ifail),
    wak=.Fortran("pelwak", PACKAGE="lmom", lmom, para=para, ifail=ifail),
    wa0=.Fortran("pelwa0", PACKAGE="lmom", lmom, para=para, ifail=ifail)
  )
  if (fort$ifail==7000) stop("pel",ddata$name,": L-moments invalid")
  para<-fort$para
  names(para)<-ddata$parnames
  return(list(para=para,ifail=fort$ifail))
}

cdfln3<-function(x,para=c(0,0,1)) {
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[3]<=0) stop("distribution parameters invalid")
  plnorm(x-para[1],meanlog=para[2],sdlog=para[3])
}

qualn3<-function(f,para=c(0,0,1)) {
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[3]<=0) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  para[1]+qlnorm(f,meanlog=para[2],sdlog=para[3])
}

lmrln3<-function(para=c(0,0,1),nmom=3) {
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[3]<=0) stop("distribution parameters invalid")
  expmu<-exp(para[2])
  xmom<-lmrgno(c(para[1]+expmu,expmu*para[3],-para[3]),nmom=nmom)
  return(xmom)
}

pelln3<-function(lmom, bound=NULL) {
  if (is.null(bound)) {
    if (length(lmom)<3) stop("need at least 3 L-moments")
    if (any(is.na(lmom[1:3]))) stop("missing values in L-moment vector")
    if (lmom[2]<=0 || lmom[3]>=1 || lmom[3]<=0) stop("L-moments invalid")
    pargno<-pelgno(c(lmom[1],lmom[2],lmom[3]))
    sigma<- -pargno[3]
    expmu<-pargno[2]/sigma
    out<-c(pargno[1]-expmu,log(expmu),sigma)
  } else {
    if (length(lmom)<2) stop("with 'bound' specified, need at least 2 L-moments")
    if (any(is.na(lmom[1:2]))) stop("missing values in L-moment vector")
    lam1<-lmom[1]-bound
    if (lam1<=0 || lmom[2]<=0 || lmom[2]>=lam1) stop("L-moments invalid")
    sigma<-sqrt(2)*qnorm((1+lmom[2]/lam1)/2)
    mu<-log(lam1)-0.5*sigma^2
    out<-c(bound,mu,sigma)
  }
  names(out)<-lmom.dist$ln3$parnames
  return(out)
}

cdfwei<-function(x,para=c(0,1,1)) {
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0 || para[3]<=0) stop("distribution parameters invalid")
  ifelse(x<=para[1],0,1-exp(-((x-para[1])/para[2])^para[3]))
}

quawei<-function(f,para=c(0,1,1)) {
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0 || para[3]<=0) stop("distribution parameters invalid")
  if (isTRUE(any(f<0 | f>1))) stop("probabilities must be between 0 and 1")
  para[1]+para[2]*((-log(1-f))^(1/para[3]))
}

lmrwei<-function(para=c(0,1,1),nmom=3) {
  if (length(para)!=3) stop("parameter vector has wrong length")
  if (any(is.na(para))) stop("missing values in parameter vector")
  if (para[2]<=0 || para[3]<=0) stop("distribution parameters invalid")
  xmom<-lmrgev(c(0,para[2]/para[3],1/para[3]),nmom=nmom)
  xmom[1]<- para[1]+para[2]-xmom[1]
  xmom[3]<- -xmom[3]
  return(xmom)
}

pelwei<-function(lmom,bound=NULL) {
  if (is.null(bound)) {
    if (length(lmom)<3) stop("need at least 3 L-moments")
    if (any(is.na(lmom[1:3]))) stop("missing values in L-moment vector")
    if (lmom[2]<=0 || lmom[3]>=1 || lmom[3]<=-lmrgum(,3)[3]) stop("L-moments invalid")
    pg<-pelgev(c(-lmom[1],lmom[2],-lmom[3]))
    delta<-1/pg[3]
    beta<-pg[2]/pg[3]
    out<-c(-pg[1]-beta,beta,delta)
  } else {
    if (length(lmom)<2) stop("with 'bound' specified, need at least 2 L-moments")
    if (any(is.na(lmom[1:2]))) stop("missing values in L-moment vector")
    lam1<-lmom[1]-bound
    if (lam1<=0 || lmom[2]<=0 || lmom[2]>=lam1) stop("L-moments invalid")
    delta<- -log(2)/log(1-lmom[2]/lam1)
    beta<-lam1/gamma(1+1/delta)
    out<-c(bound,beta,delta)
  }
  names(out)<-lmom.dist$wei$parnames
  return(out)
}

#-------------------------------------------------------------------------------
#  Sample L-moments
#-------------------------------------------------------------------------------

samlmu<-function(x, nmom=4, sort.data=TRUE, ratios=sort.data, trim=0){
  if (!is.numeric(x)) stop("'x' must be numeric")
  xok<-x[!is.na(x)]
  n<-length(xok)
  if (nmom<=0) return(numeric(0))
  t1<-trim[1]
  t2<-if (length(trim)==1) t1 else trim[2]
  maxmom<-min(nmom+t1+t2,n)
  nmom.actual<-maxmom-t1-t2
  if (nmom.actual<=0) lmom<-rep(NA_real_,nmom)
    else {
    fort<-.Fortran("samlm",PACKAGE="lmom",
          as.double(xok),
          as.integer(n),
          xmom=double(maxmom),
          as.integer(maxmom),
          isort=as.integer(sort.data),
          iratio=0L)
    lmom<-fort$xmom
    # Hosking (2007), eq.(12), applied for t=1,...t2 with s=0
    for (j in seq_len(t2)) {
      maxmom<-maxmom-1
      r<-1:maxmom
      lmom<-(lmom[r]*(r+j)-lmom[r+1]*(r+1))/(2*r+j-1)
    }
    # Hosking (2007), eq.(13), applied for s=1,...t1 with t=t2
    for (j in seq_len(t1)) {
      maxmom<-maxmom-1
      r<-1:maxmom
      lmom<-(lmom[r]*(r+j+t2)+lmom[r+1]*(r+1)*(r+t2)/r)/(2*r+j+t2-1)
    }
    #
    if (all(xok[(1+t1):(n-t2)]==xok[(1+t1)])) {
      lmom[-1]<-0    # R code doesn't guarantee this
      if (ratios) warning("all data values equal", if (any(trim>0)) " (after trimming)")
    }
    if (nmom>nmom.actual) lmom<-lmom[1:nmom]
  }
  if (ratios && nmom>2) {
    lmom[3:nmom]<-lmom[3:nmom]/lmom[2]
    names(lmom)<-c("l_1","l_2",paste("t",3:nmom,sep="_"))
  }  else names(lmom)<-paste("l",1:nmom,sep="_")
  if (!missing(trim)) names(lmom)<-sub("_", paste("(",t1,",",t2,")_",sep=""),names(lmom))
  return(lmom)
}

samlmu.s<-function(x, nmom=4, sort.data=TRUE, ratios=sort.data, trim=0){
  #
  # Preliminaries
  #
  if (!is.numeric(x)) stop("'x' must be numeric")
  xok<-x[!is.na(x)]           # Remove missing values
  n<-length(xok)
  if (nmom<=0) return(numeric(0))
  t1<-trim[1]
  t2<-if (length(trim)==1) t1 else trim[2]
  maxmom<-min(nmom+t1+t2,n)   # Number of untrimmed L-moments to compute
  nmom.actual<-maxmom-t1-t2
  if (nmom.actual<=0) {
    # return(rep(as.numeric(NA),nmom)) # but should have names
    # Next line is a device to get the right names for the return value
    out<-samlmu.s(seq_len(nmom+t1+t2),nmom=nmom,trim=trim,ratios=ratios)
    is.na(out)<-TRUE
    return(out)
  }
  #
  # Sort the data if needed
  #
  if (sort.data) xok<-sort(xok)
  #
  # Compute (untrimmed) sample L-moments.  The (j in 3:maxmom) loop computes
  # discrete Legendre polynomials recursively and uses them as weights for
  # the ordered observations.
  #
  lmom<-numeric(maxmom)
  lmom[1]<-sum(xok)
  if (maxmom>1) {
    temp<-seq(1-n,n-1,by=2)
    p1<-rep(1,n)
    p<-temp/(n-1)
    lmom[2]<-sum(xok*p)
    if (maxmom>2) {
      for (j in 3:maxmom) {
        p2<-p1
        p1<-p
        p<-((2*j-3)*temp*p1-(j-2)*(n+j-2)*p2)/((j-1)*(n-j+1))
        lmom[j]<-sum(xok*p)
      }
    }
  }
  lmom<-lmom/n   # Faster than mean(xok*p) within loop (R 2.15.3 on Windows 64-bit)
  #
  # If trimmed L-moments are needed, obtain them as linear combinations
  # of untrimmed L-moments using the recursions in Hosking (2007,
  # J.Statist.Plann.Inf.), eqs. (12)-(13).
  #
  # - Hosking (2007), eq.(12), applied for t=1,...t2 with s=0
  #
  for (j in seq_len(t2)) {
    maxmom<-maxmom-1
    r<-1:(maxmom-j)
    lmom<-(lmom[r]*(r+j)-lmom[r+1]*(r+1))/(2*r+j-1)
  }
  #
  # - Hosking (2007), eq.(13), applied for s=1,...t1 with t=t2
  #
  for (j in seq_len(t1)) {
    maxmom<-maxmom-1
    r<-1:maxmom
    lmom<-(lmom[r]*(r+j+t2)+lmom[r+1]*(r+1)*(r+t2)/r)/(2*r+j+t2-1)
  }
  length(lmom)<-nmom
  if (all(xok[(1+t1):(n-t2)]==xok[(1+t1)])) {
    lmom[-1]<-0    # R code doesn't guarantee this
    if (ratios && nmom>2) warning("all data values equal", if (any(trim>0)) " (after trimming)")
  }
  #
  # Attach names to result
  #
  if (ratios && nmom>2) {
    lmom[3:nmom]<-lmom[3:nmom]/lmom[2]
    names(lmom)<-c("l_1","l_2",paste("t",3:nmom,sep="_"))
  }  else names(lmom)<-paste("l",1:nmom,sep="_")
  if (!missing(trim)) names(lmom)<-sub("_", paste("(",t1,",",t2,")_",sep=""),names(lmom))
  #
  # Finished
  #
  return(lmom)
}

# Fast computation of L-moment ratios
# - No 'sort.data' or 'ratios'
# - NAs in 'x' will cause error
# - No warning if all data values equal
# - Return value has no names
.samlmu<-function(x, nmom=4)
  .Fortran("samlm",PACKAGE="lmom",
    as.double(x),length(x),double(nmom),as.integer(nmom),1L,1L)[[3]]

#-------------------------------------------------------------------------------
#  Utility routines called by lmrp() and lmrq()
#-------------------------------------------------------------------------------

uniroot.f<-function(f, lower, upper, ftol, maxiter=1000, ...)
##
##  Univariate root finding via interval bisection.  Similar to
##  R function uniroot(), except that the tolerance (argument 'ftol')
##  is specified in terms of the function value, not the x value.
##  Function 'f' is expected to be monotonic over the specified interval.
##  Return value is as for uniroot(): a list with elements "root",
##  "f.root", "iter", and "estim.prec", each of which has the same
##  meaning as in the return value of uniroot().
##
{
  if (!is.numeric(lower) || !is.numeric(upper))
    stop("'lower' and 'upper' must be numeric")
  if (lower>=upper) stop("'lower' must be less than 'upper'")
  fl<-f(lower,...)
  fu<-f(upper,...)
  if (sign(fl)*sign(fu)>=0)
    stop("function values at end points are not of opposite sign")
  converged<-FALSE
  for (iter in seq(len=maxiter)) {
    x<-(lower+upper)/2
    fx<-f(x,...)
    if (abs(fx)<ftol) {converged<-TRUE; break}
    if (fx<0) {lower<-x; fl<-fx} else {upper<-x; fu<-fx}
  }
  return(list(root=x, f.root=fx, iter=iter, estim.prec=(upper-lower)/2))
}



#--------------------------------------------------------------------------------

# Shifted Jacobi polynomials
sjp<-function(u,m=0,a=0,b=0) {
  if (a==0 && b==0) return(slp(u,m))
  if (a==1 && b==1) return(sjp11(u,m))
  px<-rep(1,length(u))
  if (m==0) return(px)
  p<-(a+b+2)*u-(1+b)
  if (m==1) return(p)
  for (r in 1:(m-1)) {
    pxx<-px
    px<-p
    rrab<-r+r+a+b
    AA<-(rrab+1)*(rrab+2)
    BB<-(-1)*(rrab+1)*(2*r*(r+1)+(r+r+b+1)*(a+b))/rrab
    CC<-(r+a)*(r+b)*(rrab+2)/rrab
    p<-((AA*u+BB)*px-CC*pxx)/((r+1)*(r+a+b+1))
  }
  return(p)
}

# slp(u,m): Shifted Legendre polynomial of order m, evaluated at u
# Used in the integrals computed by lmrq() in the case of no trimming
slp<-function(u,m) {
  if (m==0) return(rep(1,length(u)))
  if (m==1) return(2*u-1)
  if (m==2) return(1-6*u*(1-u))
  if (m==3) return((1-10*u*(1-u))*(2*u-1))
  z<-u*(1-u)
  if (m==4) return(1-z*(20-70*z))
  if (m==5) return((1-z*(28-126*z))*(2*u-1))
  if (m==6) return(1-z*(42-z*(378-924*z)))
  x<-2*u-1
  px<-(1-z*(28-126*z))*x
  p<-1-z*(42-z*(378-924*z))
  for (r in 7:m) {
     pxx<-px
     px<-p
     p<-((2*r-1)*x*px-(r-1)*pxx)/r
  }
  return(p)
}

# sjp11(u,m): Shifted Jacobi(1,1) polynomial of order m, evaluated at u
# Used in the integrals computed by lmrp() in the case of no trimming
sjp11<-function(u,m) {
  if (m==0) return(1)
  if (m==1) return(2*(2*u-1))
  if (m==2) return(3*(1-5*u*(1-u)))
  if (m==3) return(4*((1-7*u*(1-u))*(2*u-1)))
  z<-u*(1-u)
  if (m==4) return(5*(1-z*(14-42*z)))
  if (m==5) return(6*((1-z*(18-66*z))*(2*u-1)))
  if (m==6) return(7*(1-z*(27-z*(198-429*z))))
  x<-2*u-1
  px<-(1-z*(18-66*z))*x
  p<-1-z*(27-z*(198-429*z))
  for (r in 7:m) {
    pxx<-px
    px<-p
    p<-((2*r+1)*x*px-(r-1)*pxx)/(r+2)
  }
  return(p*(m+1))
}

#-------------------------------------------------------------------------------
#  L-moment computations for general distribution by numerical integration
#-------------------------------------------------------------------------------

#  lmrp - Computes L-moments of a general distribution, specified by its
#         cumulative distribution function.  L-moments are computed by
#         numerical integration of F(u)*(1-F(u))*Z_r(F(u)) where F is the
#         cumulative distribution function and Z_r(.) is an appropriate
#         polynomial.
#
lmrp<-function(pfunc, ..., bounds=c(-Inf,Inf), symm=FALSE, order=1:4,
               ratios=TRUE, trim=0, acc=1e-6, subdiv=100, verbose=FALSE) {
  pfunc<-try(match.fun(pfunc),silent=TRUE)
  if (inherits(pfunc,"try-error"))
    stop("'pfunc' must be a character string or a function")
  t1<-trim[1]
  t2<-if (length(trim)==1) t1 else trim[2]
  maxord<-max(order)
  out.val<-rep(0,maxord)
  out.err<-rep(0,maxord)
  out.msg<-rep("OK",maxord)
  if (is.function(bounds)) {
    out<-try(bounds(...),silent=TRUE)
    if (inherits(out,"try-error")) {
      parstring<-if (length(c(...))==0) character(0)
        else paste(c(' for parameter values',...),collapse=" ")
      stop("unable to compute bounds of distribution",parstring)
    }
    bounds<-out
  }
  if (length(bounds)!=2 || !isTRUE(bounds[1]<bounds[2])) {
    parstring<-if (length(c(...))==0) character(0)
      else paste(c(' for parameter values',...),collapse=" ")
    stop("bounds of distribution are invalid",parstring)
  }
  lbound<-bounds[1]
  ubound<-bounds[2]
  if (is.numeric(symm)) {
    med<-symm
    symm<-TRUE
  } else if (!identical(symm,FALSE)) stop("'symm' must be FALSE or numeric")
  allsymm<-(symm & (t1==t2))
#
  # Sanity checks: for bounds ...
  sanity.tol<-1e-4
  if (is.finite(lbound) && isTRUE(abs(pfunc(lbound,...))>sanity.tol))
    warning("The value of the distribution function at the lower bound of the distribution is not 0")
  if (is.finite(ubound) && isTRUE(abs(1-pfunc(ubound,...))>sanity.tol))
    warning("The value of the distribution function at the upper bound of the distribution is not 1")
  # .. and for centre of symmetry
  if (allsymm) {
    if (isTRUE(abs(pfunc(med,...)-0.5)>sanity.tol))
      warning("The supplied value of 'symm' appears not to be the median of the distribution")
  }
#
  # Multipliers for the integrals (Hosking, 2007, eq. (7))
  ord<-seq.int(2,max(2,maxord))
  mult<-c(1,exp(lgamma(ord-1)+lgamma(ord+t1+t2+1)-lgamma(ord+t1)-lgamma(ord+t2))/ord)
#
  # Integral for lambda_2
#
  # The integrand is F(x)*(1-F(x)), where F(x) is the cdf.  We add a test
  # to pick up invalid values of F(x); less than 0, greater than 1,
  # or nonnumeric (which may indicate invalid parameter values).
  #   The 'pfunc.local' that occurs in the function definition is the version
  # of 'pfunc' that we use in the integration: initially it will be set to
  # 'pfunc', but if the integration fails we will try again with a rescaled
  # version -- see below.
  lam2.integrand<-function(x) {
    f<-pfunc.local(x,...)
    if(!isTRUE(all(f>=0 & f<=1))) {
      bad<-which(f<0 | f>1 | is.na(f))[1]
      parstring<-if (length(c(...))==0) character(0)
        else paste(c(", parameter(s)",...),collapse=" ")
      stop("Value of distribution function is not in the range [0,1] for argument ",
        x[bad],parstring,": function value ",f[bad])
    }
    out<-f^(1+t1)*(1-f)^(1+t2)
    return(out)
  }
  # Initial setting of 'pfunc.local'
  pfunc.local<-pfunc
#
  # For a symmetric distribution, we use its symmetry to speed up the
  # computation. The range of integration is from the centre of symmetry
  # to 'ubound' rather than from 'lbound' to 'ubound'.
  if (allsymm) lbound<-med
#
  # Compute the integral
  int<-integrate(lam2.integrand,
    lbound,ubound,rel.tol=acc/2,subdivisions=subdiv,stop.on.error=FALSE)
#
  denom   <- int$value*mult[2]
  errdenom<- int$abs.error*mult[2]
  if (denom<0 && denom>(-errdenom)) denom<-0
  out.val[2]<-denom
  out.err[2]<-errdenom
  out.msg[2]<-int$message
  offset<-0
  scale<-1
#
  # If integration failed, try again: rescale the distribution function
  # so that the semi-interquartile range corresponds to the interval (0,1),
  # and integrate over (-Inf,+Inf). This should ensure that integrate()
  # will not deem the function to be constant everywhere.
  #   Note that user-supplied pfunc() need not return sensible results
  # if given an argument outside the specified bounds.  Instead, define
  # a modified pfunc() that returns 0 (resp. 1) if given an argument
  # less than 'lbound' (resp. greater than 'ubound').
  #
  if ( int$message!='OK' || int$abs.error<=0 || int$value <= 2*int$abs.error ) {
    lbound.orig<-lbound
    ubound.orig<-ubound
    pfunc2<-function(x,...) {
      out<-x
      out[x<lbound.orig]<-0
      out[x>ubound.orig]<-1
      xok<-(x>=lbound.orig & x<=ubound.orig)
      if (any(xok)) out[xok]<-pfunc(x[xok],...)
      out
    }
    if (allsymm) {
      p1<-pfunc2(x1<-(+1),...)
      while (p1<=0.75 & is.finite(x1)) p1<-pfunc2(x1<-2*x1, ...)
      if (!is.finite(x1)) stop("Distribution function does not approach 1 for large positive values of its argument")
      uu<-uniroot.f(f=function(x) pfunc2(x,...)-0.75,med,x1,ftol=0.1)
      q3<-uu$root
      scale<-q3-med
      pfunc.local<-function(x,...) pfunc2(x*scale+med, ...)
      lbound<- 0
      ubound<- +Inf
    } else {
      p0<-pfunc2(x0<-(-1),...)
      while (p0>=0.25 & is.finite(x0)) p0<-pfunc2(x0<-2*x0, ...)
      if (!is.finite(x0)) stop("Distribution function does not approach 0 for large negative values of its argument")
      p1<-pfunc2(x1<-(+1),...)
      while (p1<=0.75 & is.finite(x1)) p1<-pfunc2(x1<-2*x1, ...)
      if (!is.finite(x1)) stop("Distribution function does not approach 1 for large positive values of its argument")
      uu<-uniroot.f(f=function(x) pfunc2(x,...)-0.25,x0,x1,ftol=0.1)
      q1<-uu$root
      uu<-uniroot.f(f=function(x) pfunc2(x,...)-0.75,x0,x1,ftol=0.1)
      q3<-uu$root
      offset<-(q3+q1)/2
      scale <-(q3-q1)/2
      pfunc.local<-function(x,...) pfunc2(x*scale+offset, ...)
      lbound<- -Inf
      ubound<- +Inf
    }
    # If 'scale' is zero (or negative, if this can happen), the rescaling failed:
    # perhaps the upper and lower quartiles of the distribution are equal.
    # We'll set 'scale' equal to 1 and hope for the best.
    if (scale<=0) scale<-1
    #
    int<-integrate(lam2.integrand,lbound,ubound,
      rel.tol=acc/2,subdivisions=subdiv,stop.on.error=FALSE)
    denom   <- int$value*mult[2]
    errdenom<- int$abs.error*mult[2]
    if (denom<0 && denom>(-errdenom)) denom<-0
    out.val[2]<-denom*scale
    out.err[2]<-errdenom*scale
    out.msg[2]<-int$message
  }
#
  # If lambda_2 is negative, the integration must have failed
  if (denom<0) {
    out.val[]<-NA
    out.err[]<-NA
    out.msg[]<-"computed value of lambda_2 is negative"
    if (int$message!="OK")
      out.msg[2]<-paste(int$message,"; ",out.msg[2],sep="")
#
  } else {
#
    # For a symmetric distribution and symmetric trimming, lambda_2 is twice
    # the integral and lambda_1 (when it exists) is the centre of symmetry
    if (allsymm) {
      out.val[2]<-2*out.val[2]
      out.err[2]<-2*out.err[2]
      if (out.msg[2]=="OK") out.val[1]<-med
      else {
        out.val[1]<-NA
        out.err[1]<-NA
        out.msg[1]<-"Unable to compute lambda_2"
      }
    }
#
    # For a general distribution, compute lambda_1 if required. Use the
    # value of lambda_2 as a guide to the accuracy that can be achieved.
    if (!allsymm && is.element(1,order)) {
      abstol<-max(denom*acc/2,errdenom)
      if (abstol==0) abstol<-sqrt(.Machine$double.eps)
      if (is.finite(lbound)) {
        int<-integrate( function(x) pbeta(pfunc(x,...),1+t1,1+t2,lower.tail=FALSE),lbound,ubound,
             rel.tol=0,abs.tol=abstol,subdivisions=subdiv,stop.on.error=FALSE)
        out.val[1]<-int$value+lbound
        out.err[1]<-int$abs.error
        out.msg[1]<-int$message
      }
      else if (is.finite(ubound)) {
        int<-integrate( function(x) pbeta(pfunc(x,...),1+t1,1+t2),lbound,ubound,
             rel.tol=0,abs.tol=abstol,subdivisions=subdiv,stop.on.error=FALSE)
        out.val[1]<-ubound-int$value
        out.err[1]<-int$abs.error
        out.msg[1]<-int$message
      } else {
#       int1<-integrate( function(x) pbeta(pfunc.local(x,...),1+t1,1+t2),
#            -Inf, 0,  rel.tol=0,abs.tol=abstol/2,subdivisions=subdiv,stop.on.error=FALSE)
#       int2<-integrate( function(x) pbeta(pfunc.local(x,...),1+t1,1+t2,lower.tail=FALSE),
#             0, +Inf, rel.tol=0,abs.tol=abstol/2,subdivisions=subdiv,stop.on.error=FALSE)
#       out.val[1]<-(int2$value-int1$value)*scale+offset
#       out.err[1]<-(int1$abs.error+int2$abs.error)*scale
#       out.msg[1]<-if (int1$message=="OK") int2$message else int1$message

        # mean = integral_0^Inf { 1 - F(x) - F(-x) } dx, generalized to include trimming
        integrand<- function(x) pbeta(pfunc.local(x,...),1+t1,1+t2,lower.tail=FALSE) -
                                pbeta(pfunc.local(-x,...),1+t1,1+t2)
        int<-integrate( integrand, 0, Inf,
          rel.tol=0,abs.tol=abstol,subdivisions=subdiv,stop.on.error=FALSE)

        out.val[1]<-int$value*scale+offset
        out.err[1]<-int$abs.error*scale
        out.msg[1]<-int$message
      }
    }
#
    # Test whether lambda_2 is zero to within the reported absolute error
    if (out.msg[2]=="OK" && abs(out.val[2])<=out.err[2]) {
#
      out.msg[2]<-"lambda_2 is effectively zero"
      if (maxord>2) {
        if (ratios) {
          out.val[3:maxord]<-NaN
          out.err[3:maxord]<-NA
          out.msg[3:maxord]<-"L-moment ratios not defined when lambda_2 is zero"
        } else {
          out.val[3:maxord]<-0
          out.err[3:maxord]<-NA
          out.msg[3:maxord]<-"Higher-order L-moments not computed when lambda_2 is zero"
        }
      }
#
    } else {
#
      # 'higher' contains the orders of the L-moments that we compute via
      # numerical integration: all orders greater than 2 in the general case,
      # only the even orders greater than 2 for a symmetric distribution
      higher<-order[order>2 & { if (allsymm) order%%2==0 else TRUE } ]
#
      # Integrals for higher-order L-moments
      for (r in higher) {
        int<-integrate(
          function(x)
            {pf<-pfunc.local(x,...); pf^(1+t1)*(1-pf)^(1+t2)*sjp(pf,r-2,1+t2,1+t1)},
          lbound, ubound, rel.tol=0, abs.tol=denom*acc/2,
          subdivisions=subdiv, stop.on.error=FALSE)
        out.val[r]<-int$value*mult[r]
        out.err[r]<-int$abs.error*mult[r]
        out.msg[r]<-int$message
      }
      if (allsymm) {
        out.val[higher]<-2*out.val[higher]
        out.err[higher]<-2*out.err[higher]
      }

#     # L-moment ratios and their accuracies
      if (ratios) {
        out.val[higher]<-out.val[higher]/out.val[2]
        out.err[higher]<-(out.err[higher]+abs(out.val[higher])*out.err[2])/out.val[2]
      }
    }
  }
#
  # Prepare and return the output
  nnames<-length(out.val) # number of names needed, equal to 'max(2,maxord)'
  rnames<-rep("lambda",nnames)
  if (ratios && nnames>2) rnames[-(1:2)]<-"tau"
  if (any(trim!=0)) rnames<-paste(rnames,"(",t1,",",t2,")",sep="")
  rnames<-paste(rnames,1:nnames, sep="_")
  if (verbose) {
    # If 'order' does not contain '2' and there was an error message
    # when calculating lambda_2, report it as a warning:
    # it won't be included in the return value.
    if (!is.element(2,order) && out.msg[2]!="OK")
      warning("in computation of L-moment of order 2: ",out.msg[2])
    dframe<-data.frame(row.names=rnames,value=out.val,abs.error=out.err,message=out.msg)
    return(dframe[order,])
  } else {
    for (r in intersect(c(2,order),which(out.msg!="OK")))
      # 'if (...)' suppresses secondary warnings (when problems computing
      # lambda_2 meant that we couldn't compute L-moments of other orders)
      if (!(r!=2 && regexpr("lambda_2",out.msg[r])>0))
        warning("in computation of L-moment of order ",r,": ",out.msg[r])
    names(out.val)<-rnames
    return(out.val[order])
  }
}

#  lmrq - Computes L-moments of a general distribution, specified by its quantile function.
#         L-moments are computed by numerical integration of Q(u)*Pstar[m](u) where
#         Q is the quantile function and Pstar[m] is a shifted Legendre polynomial.
#
lmrq<-function(qfunc, ..., symm=FALSE, order=1:4, ratios=TRUE, trim=0, acc=1e-6,
               subdiv=100, verbose=FALSE) {
  qfunc<-try(match.fun(qfunc),silent=TRUE)
  if (inherits(qfunc,"try-error"))
    stop("'qfunc' must be a character string or a function")
  t1<-trim[1]
  t2<-if (length(trim)==1) t1 else trim[2]
  maxord<-max(order)
  out.val<-rep(0,maxord)
  out.err<-rep(0,maxord)
  out.msg<-rep("OK",maxord)
  allsymm<- (symm & (t1==t2))

  # For symmetric trimming and a symmetric distribution, we use the symmetry to
  # speed up the computation. The integrals involve Q(u)-Q(0.5) rather than Q(u),
  # and the range of integration is 0 to 0.5 rather than 0 to 1.
  #   'upper' contains the upper limit of integration.
  #   'qfunc.local' is the version of the quantile function that we use.
  if (allsymm) {
    upper<-0.5
    qfunc.local<- { function(u,...) qfunc(u,...)-med }
    # (Symmetric distribution) Compute the median
    out.val[1]<-med<-qfunc(0.5,...)
  } else {
    upper<-1
    qfunc.local<-qfunc
  }

  # Multipliers for the integrals (Hosking, 2007, eq. (5))
  ord<-seq_len(max(2,maxord))
  mult<-exp(lgamma(ord)+lgamma(ord+t1+t2+1)-lgamma(ord+t1)-lgamma(ord+t2))/ord

  # Integral for lambda_2
  lam2.integrand<-function(u) {
    qval<-qfunc.local(u,...)
    if (any(is.na(qval))) {
      bad<-which(is.na(qval))[1]
      parstring<-if (length(c(...))==0) character(0)
        else paste(c(", parameter(s)",...),collapse=" ")
      stop("Quantile function returned ",qval[bad]," for argument ",u[bad],parstring)
    }
    if (!identical(order(qval[order(u)]),seq(along=qval))) {
      ou<-order(u)
      uu<-u[ou]
      qq<-qval[ou]
      bad<-which(diff(qq)<0)[1]
      parstring<-if (length(c(...))==0) character(0)
        else paste(c(" for parameter(s)",...),collapse=" ")
      stop("Quantile function is not increasing on the interval ",
        uu[bad]," to ",uu[bad+1],parstring)
    }
    # if (notrim) return(qval*(2*u-1)) else   # no speed advantage
    return(qval*u^t1*(1-u)^t2*((1+t2)*u-(1+t1)*(1-u)))
  }

  int<-integrate( lam2.integrand,0,upper,
         rel.tol=acc/2,subdivisions=subdiv,stop.on.error=FALSE)
  denom    <- int$value*mult[2]
  errdenom <- int$abs.error*mult[2]
  if (denom<0 && denom>(-errdenom)) denom<-0
  out.val[2]<-denom
  out.err[2]<-errdenom
  out.msg[2]<-int$message

  # For a symmetric distribution, lambda_2 is twice the integral
  if (allsymm) {out.val[2]<-2*denom; out.err[2]<-2*errdenom}

  # If lambda_2 is negative, the integration must have failed
  if (denom<0) {
    out.val[]<-NA
    out.err[]<-NA
    out.msg[]<-"computed value of lambda_2 is negative"
    if (int$message!="OK")
      out.msg[2]<-paste(int$message,"; ",out.msg[2],sep="")
#
  } else {
#

    # For a general distribution, compute lambda_1 if required. Use the
    # value of lambda_2 as a guide to the accuracy that can be achieved.
    if (!allsymm && is.element(1,order)) {
      int<-integrate( function(u) qfunc.local(u,...)*u^t1*(1-u)^t2,0,1,
        rel.tol=0, abs.tol=max(denom*acc/2,errdenom),
        subdivisions=subdiv, stop.on.error=FALSE)
      out.val[1]<-int$value*mult[1]
      out.err[1]<-int$abs.error*mult[1]
      out.msg[1]<-int$message
    }

    # Test whether lambda_2 is zero to within the reported absolute error
    if (out.msg[2]=="OK" && abs(denom)<=out.err[2]) {

      if (maxord>2) {
        if (ratios) {
          out.val[3:maxord]<-NaN
          out.err[3:maxord]<-NA
          out.msg[3:maxord]<-"L-moment ratios not defined when lambda_2 is zero"
        } else {
          out.val[3:maxord]<-0
          out.err[3:maxord]<-NA
          out.msg[3:maxord]<-"Higher-order L-moments not computed when lambda_2 is zero"
        }
      }

    } else {

      # 'higher' contains the orders of the L-moments that we compute via
      # numerical integration: all orders greater than 2 in the general case,
      # only the even orders greater than 2 for a symmetric distribution
      higher<-order[order>2 & { if (allsymm) order%%2==0 else TRUE } ]

      # Integrals for higher-order L-moments
      for (r in higher) {
        int<-integrate( function(u) qfunc.local(u,...)*u^t1*(1-u)^t2*sjp(u,r-1,t2,t1),
          0,upper,rel.tol=0,abs.tol=denom*acc/2,subdivisions=subdiv,stop.on.error=FALSE)
        out.val[r]<-int$value*mult[r]
        out.err[r]<-int$abs.error*mult[r]
        out.msg[r]<-int$message
      }
      if (allsymm) {
        out.val[higher]<-2*out.val[higher]
        out.err[higher]<-2*out.err[higher]
      }

#     # L-moment ratios and their accuracies
      if (ratios) {
        out.val[higher]<-out.val[higher]/out.val[2]
        out.err[higher]<-(out.err[higher]+abs(out.val[higher])*out.err[2])/out.val[2]
      }


    }
  }

  # Prepare and return the output
  nnames<-length(out.val) # number of names needed, equal to 'max(2,maxord)'
  rnames<-rep("lambda",nnames)
  if (ratios && nnames>2) rnames[-(1:2)]<-"tau"
  if (any(trim!=0)) rnames<-paste(rnames,"(",t1,",",t2,")",sep="")
  rnames<-paste(rnames,1:nnames, sep="_")
  if (verbose) {
    # If 'order' does not contain '2' and there was an error message
    # when calculating lambda_2, report it as a warning:
    # it won't be included in the return value.
    if (!is.element(2,order) && out.msg[2]!="OK")
      warning("in computation of L-moment of order 2: ",out.msg[2])
    dframe<-data.frame(row.names=rnames,value=out.val,abs.error=out.err,message=out.msg)
    return(dframe[order,])
  } else {
    for (r in intersect(c(2,order),which(out.msg!="OK")))
      # 'if (...)' suppresses secondary warnings (when problems computing
      # lambda_2 meant that we couldn't compute L-moments of other orders)
      if (!(r!=2 && regexpr("lambda_2",out.msg[r])>0))
        warning("in computation of L-moment of order ",r,": ",out.msg[r])
    names(out.val)<-rnames
    return(out.val[order])
  }
}

#-------------------------------------------------------------------------------
#  Parameter estimation for a general distribution
#-------------------------------------------------------------------------------

pelp<-function(lmom, pfunc, start, bounds=c(-Inf,Inf), type=c("n","s","ls","lss"),
               ratios=NULL, trim=NULL, method="nlm", acc=1e-5, subdiv=100, ...) {
##
  type<-match.arg(type)

  infer<-infer_trim(names(lmom))
  if (is.null(infer)) {
    if (is.null(ratios)) stop("unable to infer 'ratios' from names of 'lmom'")
    if (is.null(trim)) stop("unable to infer trim orders from names of 'lmom'")
  } else if (identical(infer,list())) {
    if (is.null(ratios)) ratios<-TRUE
    if (is.null(trim)) trim<-0
  } else {
    if (is.null(ratios)) ratios<-infer$ratios
    else if (length(lmom)>2 && !identical(as.logical(ratios),infer$ratios))
      warning("value of argument 'ratios' is incompatible with names of 'lmom'")
    if (is.null(trim)) trim<-infer$trim
    else if (!all(c(trim,trim)[1:2]==infer$trim))
      warning("value of argument 'trim' is incompatible with names of 'lmom'")
  }

  method<-try(match.arg(method,c("nlm", "uniroot", eval(formals(optim)$method))),
    silent=TRUE)
  if (class(method)=="try-error") stop(sub(".*'arg'","'method'",method))
  npara<-length(start)     # Number of parameters
  outpar<-numeric(npara)   # Vector to hold the estimated parameter values

  nlmom<-length(lmom)
  maxord<- if (type=="lss" && npara>2) 2*npara-2 else npara # Highest-order L-moment that will be used
  if (nlmom<maxord) stop("not enough L-moments, or too many parameters")
  if (nlmom>maxord) warning(nlmom," L-moments supplied, but only ",maxord," will be used")

  nlocscale<-switch(type, n=0, s=1, ls=2, lss=2) # No. of loc/scale parameters
  nshape<-npara-nlocscale                # Number of shape parameters
  shape.pos<-seq(nlocscale+1,len=nshape) # Positions of shape parameters in the parameter vector

  # If type is "n" we will equate sample and population L-moments
  # If type is "s" we will equate sample and population L-moments scaled by lambda_1
  # Convert the input (L-moment ratios, with scaling by lambda_2) if necessary
  if (is.element(type,c("n","s")) && nlmom>2) lmom[-(1:2)]<-lmom[-(1:2)]*lmom[2]
  if (type=="s" && nlmom>1) lmom[-1]<-lmom[-1]/lmom[1]

  code<-0 # Convergence indicator.  Will remain unchanged if nshape==0.

  # Are the distribution parameters supplied to qfunc as a parameter vector
  # or as separate arguments?
  parvec<-length(formals(pfunc))==2 && npara>1

  # If there are any shape parameters, estimate them by numerical optimization
  if (nshape>0) {

    # order.shape - orders of L-moments that are used to estimate the shape parameters
    order.shape<-switch(type,n=1:nshape,s=2:(nshape+1),ls=3:(nshape+2),lss=seq(4,by=2,len=nshape))

    # 'critfn' is the criterion function for optimization.  For function
    # minimization, via nlm() or optim(), it is the sum of squared
    # differences between the sample and population L-moments or L-moment
    # ratios.  For root-finding, via uniroot(), it is the difference between
    # the sample and population versions of the L-moment or L-moment ratio.
    #   critfn() makes calls to lmrp(), via 'do.call(lmrp,...)'.
    critfn<-function(para) {

      # 'fullpara' contains all parameters, including location and scale
      fullpara<-switch(type,n=para,s=c(1,para),ls=c(0,1,para),lss=c(0,1,para))

      # 'paralist' contains all parameters, as a list suitable for passing
      # to 'do.call(pfunc,...)': i.e. if pfunc() expects a single vector
      # of parameters, then 'paralist' is a list with a single element which
      # is the vector of parameters; if pfunc() expects separate arguments,
      # then 'paralist' is a list with as many elements as there are
      # parameters and one parameter in each element.
      paralist<- if (parvec) list(fullpara) else as.list(fullpara)

      # Compute L-moments
      if (type=="n") {

        # For type "n", compute L-moments (so that all elements of zmom
        # have approx. the same magnitude)
        zmom<-do.call("lmrp",c(list(pfunc,bounds=bounds,order=order.shape,ratios=ratios,trim=trim,acc=acc/10,subdiv=subdiv),paralist))
        if (nlmom>2) zmom[-(1:2)]<-zmom[-(1:2)]*zmom[2]

      } else if (type=="s") {

        # For type "s", scale the L-moments by lambda_1, not lambda_2
        zmom<-do.call("lmrp",c(list(pfunc,bounds=bounds,order=1:npara,ratios=ratios,trim=trim,acc=acc/10,subdiv=subdiv),paralist))
        if (nlmom>2) zmom[-(1:2)]<-zmom[-(1:2)]*zmom[2]
        zmom<-zmom[-1]/zmom[1]

      } else {

        # For type "ls", compute L-moment ratios
        # For type "lss", compute L-moment ratios via lmrp(...,symm=TRUE)
        symm<-if (type=="lss") 0 else FALSE
        zmom<-do.call("lmrp",c(list(pfunc,bounds=bounds,order=order.shape,ratios=ratios,trim=trim,acc=acc/10,subdiv=subdiv,symm=symm),paralist))

      }

      # Construct the criterion function value
      retval<-if (method=="uniroot") zmom-lmom[order.shape]
              else sum((zmom-lmom[order.shape])^2)

      return(retval)
    }

    # Get the user-supplied arguments for the optimization functions
    dotargs<-list(...)
    if (any(names(dotargs)=="")) stop("arguments in '...' must be named")

    # Call an optimization routine, as indicated by the 'method' argument

    if (method=="uniroot") { # Estimation of a single parameter via root-finding

      if (nshape>1) stop('method "uniroot" not available when there is more than one shape parameter')

      # If user didn't provide a 'tol' argument, set it based on the value of 'acc'
      if (is.null(dotargs$tol)) dotargs<-c(dotargs,tol=acc)

      # Call the root-finding function, uniroot()
      opt<-do.call(stats::uniroot,c(f=critfn,dotargs))

      # Copy estimated shape parameter values to output vector; set return code
      outpar[shape.pos]<-opt$root
      code<-if (opt$estim.prec<=acc) 1 else 4

    } else if (method=="nlm") {

      # If user didn't provide 'fscale' or 'ndigit' arguments,
      # set them to our own defaults
      if (is.null(dotargs$fscale)) dotargs<-c(dotargs,fscale=0)
      if (is.null(dotargs$ndigit)) dotargs<-c(dotargs,ndigit=round(-2*log10(acc)))

      # Call the optimization function, nlm()
      opt<-do.call(stats::nlm,c(f=critfn,list(p=start[shape.pos]),dotargs))

      # Copy estimated shape parameter values to output vector; set return code
      outpar[shape.pos]<-opt$estimate
      code<-opt$code

    } else { # Optimization via one of the methods of optim()

      # If user didn't provide 'ndeps' or 'abstol' elements of the 'control'
      # argument, set them to our own defaults
      if (is.null(dotargs$control)) dotargs<-c(dotargs,list(control=list()))
      if (is.null(dotargs$control$ndeps))
        dotargs$control<-c(dotargs$control,list(ndeps=rep(acc,nshape)))
      if (method!="L-BFGS-B" && is.null(dotargs$control$abstol))
        dotargs$control<-c(dotargs$control,abstol=acc^2)

      # Call the optimization function, optim()
      opt<-do.call(stats::optim,c(list(par=start[shape.pos]),fn=critfn,method=method,dotargs))

      # Copy estimated shape parameter values to output vector; set return code
      outpar[shape.pos]<-opt$par
      code<-opt$convergence
      if (code==1) code<-4
      if (code==0) code<-1
    }

    if (code==3) warning('iteration may not have converged; estimates may be unreliable')
    if (code>=4) warning('iteration has not converged; estimates may be unreliable')
  }

  # Estimate the location and scale parameters

  if (is.element(type,c("ls","lss"))) {

    outpar[1:2]<-c(0,1)
    paralist<- if (parvec) list(outpar) else as.list(outpar)
    zmom<-do.call("lmrp",c(list(pfunc,bounds=bounds,order=1:2,ratios=ratios,trim=trim,acc=acc/10,subdiv=subdiv),paralist))
    outpar[2]<-lmom[2]/zmom[2]
    outpar[1]<-lmom[1]-outpar[2]*zmom[1]

  } else if (type=="s") {

    outpar[1]<-1
    paralist<- if (parvec) list(outpar) else as.list(outpar)
    zmom<-do.call("lmrp",c(list(pfunc,bounds=bounds,order=1,ratios=ratios,trim=trim,acc=acc/10,subdiv=subdiv),paralist))
    outpar[1]<-lmom[1]/zmom[1]

  }

  if (!parvec) names(outpar)<-names(formals(pfunc))[2:(npara+1)]
  return(list(para=outpar,code=code))
}

pelq<-function(lmom, qfunc, start, type=c("n","s","ls","lss"),
                ratios=NULL, trim=NULL, method="nlm", acc=1e-5, subdiv=100, ...) {
  type<-match.arg(type)

  infer<-infer_trim(names(lmom))
  if (is.null(infer)) {
    if (is.null(ratios)) stop("unable to infer 'ratios' from names of 'lmom'")
    if (is.null(trim)) stop("unable to infer trim orders from names of 'lmom'")
  } else if (identical(infer,list())) {
    if (is.null(ratios)) ratios<-TRUE
    if (is.null(trim)) trim<-0
  } else {
    if (is.null(ratios)) ratios<-infer$ratios
    else if (length(lmom)>2 && !identical(as.logical(ratios),infer$ratios))
      warning("value of argument 'ratios' is incompatible with names of 'lmom'")
    if (is.null(trim)) trim<-infer$trim
    else if (!all(c(trim,trim)[1:2]==infer$trim))
      warning("value of argument 'trim' is incompatible with names of 'lmom'")
  }

  method<-try(match.arg(method,c("nlm", "uniroot", eval(formals(optim)$method))),
    silent=TRUE)
  if (class(method)=="try-error") stop(sub(".*'arg'","'method'",method))
  npara<-length(start)     # Number of parameters
  outpar<-numeric(npara)   # Vector to hold the estimated parameter values

  nlmom<-length(lmom)
  maxord<- if (type=="lss" && npara>2) 2*npara-2 else npara # Highest-order L-moment that will be used
  if (nlmom<maxord) stop("not enough L-moments, or too many parameters")
  if (nlmom>maxord) warning(nlmom," L-moments supplied, but only ",maxord," will be used")

  nlocscale<-switch(type, n=0, s=1, ls=2, lss=2) # No. of loc/scale parameters
  nshape<-npara-nlocscale                # Number of shape parameters
  shape.pos<-seq(nlocscale+1,len=nshape) # Positions of shape parameters in the parameter vector

  # If type is "n" we will equate sample and population L-moments
  # If type is "s" we will equate sample and population L-moments scaled by lambda_1
  # Convert the input (L-moment ratios, with scaling by lambda_2) if necessary
  if (is.element(type,c("n","s")) && nlmom>2) lmom[-(1:2)]<-lmom[-(1:2)]*lmom[2]
  if (type=="s" && nlmom>1) lmom[-1]<-lmom[-1]/lmom[1]

  code<-0 # Convergence indicator.  Will remain unchanged if nshape==0.

  # Are the distribution parameters supplied to qfunc as a parameter vector
  # or as separate arguments?
  parvec<-length(formals(qfunc))==2 && npara>1

  # If there are any shape parameters, estimate them by numerical optimization
  if (nshape>0) {

    # order.shape - orders of L-moments that are used to estimate the shape parameters
    order.shape<-switch(type,n=1:nshape,s=2:(nshape+1),ls=3:(nshape+2),lss=seq(4,by=2,len=nshape))

    # 'critfn' is the criterion function for optimization.  For function
    # minimization, via nlm() or optim(), it is the sum of squared
    # differences between the sample and population L-moments or L-moment
    # ratios.  For root-finding, via uniroot(), it is the difference between
    # the sample and population versions of the L-moment or L-moment ratio.
    #   critfn() makes calls to lmrq(), via 'do.call(lmrq,...)'.
    critfn<-function(para) {

      # 'fullpara' contains all parameters, including location and scale
      fullpara<-switch(type,n=para,s=c(1,para),ls=c(0,1,para),lss=c(0,1,para))

      # 'paralist' contains all parameters, as a list suitable for passing
      # to 'do.call(qfunc,...)': i.e. if qfunc() expects a single vector
      # of parameters, then 'paralist' is a list with a single element which
      # is the vector of parameters; if qfunc() expects separate arguments,
      # then 'paralist' is a list with as many elements as there are
      # parameters and one parameter in each element.
      paralist<- if (parvec) list(fullpara) else as.list(fullpara)

      # Compute L-moments
      if (type=="n") {

        # For type "n", compute L-moments (so that all elements of 'zmom'
        # have approx. the same magnitude)
        zmom<-do.call("lmrq",c(list(qfunc,order=order.shape,trim=trim,acc=acc/10,subdiv=subdiv),paralist))
        if (nlmom>2) zmom[-(1:2)]<-zmom[-(1:2)]*zmom[2]

      } else if (type=="s") {

        # For type "s", scale the L-moments by lambda_1, not lambda_2
        zmom<-do.call("lmrq",c(list(qfunc,order=1:npara,trim=trim,acc=acc/10,subdiv=subdiv),paralist))
        if (nlmom>2) zmom[-(1:2)]<-zmom[-(1:2)]*zmom[2]
        zmom<-zmom[-1]/zmom[1]

      } else {

        # For type "ls", compute L-moment ratios
        # For type "lss", compute L-moment ratios via lmrq(...,symm=TRUE)
        zmom<-do.call("lmrq",c(list(qfunc,order=order.shape,trim=trim,acc=acc/10,subdiv=subdiv,symm=(type=="lss")),paralist))

      }

      # Construct the criterion function value
      retval<-if (method=="uniroot") zmom-lmom[order.shape]
              else sum((zmom-lmom[order.shape])^2)

      return(retval)
    }

    # Get the user-supplied arguments for the optimization functions
    dotargs<-list(...)
    if (any(names(dotargs)=="")) stop("arguments in '...' must be named")

    if (method=="uniroot") { # Estimation of a single parameter via root-finding

      if (nshape>1) stop('method "uniroot" not available when there is more than one shape parameter')

      # If user didn't provide a 'tol' argument, set it based on the value of 'acc'
      if (is.null(dotargs$tol)) dotargs<-c(dotargs,tol=acc)

      # Call the root-finding function, uniroot()
      opt<-do.call(stats::uniroot,c(f=critfn,dotargs))

      # Copy estimated shape parameter values to output vector; set return code
      outpar[shape.pos]<-opt$root
      code<-if (opt$estim.prec<=acc) 1 else 4

    } else if (method=="nlm") {

      # If user didn't provide 'fscale' or 'ndigit' arguments,
      # set them to our own defaults
      if (is.null(dotargs$fscale)) dotargs<-c(dotargs,fscale=0)
      if (is.null(dotargs$ndigit)) dotargs<-c(dotargs,ndigit=round(-2*log10(acc)))

      # Call the optimization function, nlm()
      opt<-do.call(stats::nlm,c(f=critfn,list(p=start[shape.pos]),dotargs))

      # Copy estimated shape parameter values to output vector; set return code
      outpar[shape.pos]<-opt$estimate
      code<-opt$code

    } else {

      # If user didn't provide 'ndeps' or 'abstol' elements of the 'control'
      # argument, set them to our own defaults
      if (is.null(dotargs$control)) dotargs<-c(dotargs,list(control=list()))
      if (is.null(dotargs$control$ndeps))
        dotargs$control<-c(dotargs$control,list(ndeps=rep(acc,nshape)))
      if (method!="L-BFGS-B" && is.null(dotargs$control$abstol))
        dotargs$control<-c(dotargs$control,abstol=acc^2)

      # Call the optimization function, optim()
      opt<-do.call(stats::optim,c(list(par=start[shape.pos]),fn=critfn,method=method,dotargs))

      # Copy estimated shape parameter values to output vector; set return code
      outpar[shape.pos]<-opt$par
      code<-opt$convergence
      if (code==1) code<-4
      if (code==0) code<-1
    }
    if (code==3) warning('iteration may not have converged; estimates may be unreliable')
    if (code>=4) warning('iteration has not converged; estimates may be unreliable')
  }

  # Estimate the location and scale parameters

  if (is.element(type,c("ls","lss"))) {

    outpar[1:2]<-c(0,1)
    paralist<- if (parvec) list(outpar) else as.list(outpar)
    zmom<-do.call("lmrq",c(list(qfunc,order=1:2,trim=trim,acc=acc/10,subdiv=subdiv),paralist))
    outpar[2]<-lmom[2]/zmom[2]
    outpar[1]<-lmom[1]-outpar[2]*zmom[1]

  } else if (type=="s") {

    outpar[1]<-1
    paralist<- if (parvec) list(outpar) else as.list(outpar)
    zmom<-do.call("lmrq",c(list(qfunc,order=1,trim=trim,acc=acc/10,subdiv=subdiv),paralist))
    outpar[1]<-lmom[1]/zmom[1]

  }

  if (!parvec) names(outpar)<-names(formals(qfunc))[2:(npara+1)]
  return(list(para=outpar,code=code))
}

infer_trim<-function(nam) {
##
## Infers, from the names of a vector of L-moments or L-moment ratios
## returned by lmrp(), lmrq(), or samlmu(), whether the vector contains
## L-moments or L-moment ratios and the order of trimming of the
## L-moments or L-moment ratios.
##
## The return value is a list with components "ratios" and "trim".
##
## Such vectors of names can have the form
##
##   l(s,t)_1  l(s,t)_2  l(s,t)_3 ...
##   lambda(s,t)_1  lambda(s,t)_2  lambda(s,t)_3 ...
##
## or
##
##   l(s,t)_1  l(s,t)_2  t(s,t)_3 ...
##   lambda(s,t)_1  lambda(s,t)_2  tau(s,t)_3 ...
##
## with 's' and 't' positive integers, the orders of trimming;
## in this case the return value's "trim" component is 'c(s,t)'.
## Vectors can also have the form
##
##   l_1  l_2  l_3 ...
##   lambda_1  lambda_2  lambda_3 ...
##
## or
##
##   l_1  l_2  t_3 ...
##   lambda_1  lambda_2  tau_3 ...
##
## with no trim specification, indicating no trimming; in this case
## the return value's "trim" component is 'c(0,0)'.
##
## If the input vector has no names, the return value is an empty list.
##
## In other cases when the orders of trimming cannot be unequivocally
## identified, the return value is NULL.
##
  ratios<-TRUE
  if (is.null(nam)) return(list())
  if (length(nam)==0) return(NULL)
  opt<-options(warn=-1)
  on.exit(options(opt))
  if (all(substr(nam,1,1)=="l")) ratios<-FALSE
  else {
    if (any(substr(c(nam,"l","l")[1:2],1,1)!="l")) return(NULL)     # First character of first two elements must be "l"
    if (any(substr(nam[-(1:2)],1,1)!="t")) return(NULL)             # First character of all other elements must be "t"
    ratios<-TRUE
  }
  nam<-sub("^lambda","l",nam)                                       # Treat initial "lambda" as "l"
  nam<-sub("^tau","t",nam)                                          # ... and initial "tau" as "t"
  orders<-sub(".*_","",nam)                                         # Everything after "_" is the order of the L-moment (or ratio)
  if (!identical(as.integer(orders),seq_along(nam))) return(NULL)   # Orders must be the sequence 1,2,...
  super<-sub("^.(.*)_.*$","\\1",nam)                                # Everything between 1st character and "_" is the trim specification
  if (all(nchar(super)==0)) return(list(ratios=ratios,trim=c(0,0))) # A blank trim specification is valid and indicates no trimming
  if (!all(grepl("^\\(.*\\)$",super))) return(NULL)                 # A nonblank trim specification must begin with "(" and end with ")"
  trim1<-substr(super,2,nchar(super)-1)                             # Everything in between ...
  trim2<-strsplit(paste(trim1,",",sep=""),",")                      # ... is treated as a sequence of items separated by ","
  if (!all(sapply(trim2,length)==2)) return(NULL)                   # There must be exactly 2 items
  if (!all(sapply(trim2,identical,trim2[[1]]))) return(NULL)        # They must be the same in each element of the names vector
  trim3<-as.numeric(trim2[[1]])                                     # When converted to numbers ...
  if (!(all(is.finite(trim3)) && all(trim3>=0) && all(trim3 %%1==0))) return(NULL) # ... they must be positive integers, returned as the trim orders
  return(list(ratios=ratios,trim=trim3))
}

#-------------------------------------------------------------------------------
#  Functions and data for L-moment ratio diagram
#-------------------------------------------------------------------------------

lmrd<-function(x, y, distributions = "GLO GEV GPA GNO PE3", twopar,
               xlim, ylim, pch=3, cex, col, lty, lwd=1,
               legend.lmrd = TRUE, xlegend, ylegend,
               xlab = expression(italic(L) * "-skewness"),
               ylab = expression(italic(L) * "-kurtosis"), ...) {
## Function lmrd() -- draws an L-moment ratio diagram
# Check arguments
#
  x.missing<-missing(x)
  if (x.missing) { x<-y<-numeric(0) }
  else if (missing(y)) {
    y<-NULL
    if (inherits(x,"regdata")) {
      y<-x[[6]]; x<-x[[5]]
    } else if (length(dim(x))==2) {
      dn2<-dimnames(x)[[2]]
      mm<-match(c("t_3","t_4"),dn2)
      if (any(is.na(mm))) mm<-match(c("tau_3","tau_4"),dn2)
      if (!any(is.na(mm))) { y<-x[,mm[2] ]; x<-x[,mm[1] ] }
    } else if (is.numeric(x)) {
      nx<-names(x)
      mm<-match(c("t_3","t_4"),nx)
      if (any(is.na(mm))) mm<-match(c("tau_3","tau_4"),nx)
      if (!any(is.na(mm))) { y<-x[mm[2] ]; x<-x[mm[1] ] }
    }
    if (is.null(y))
      stop("could not find L-skewness and L-kurtosis values in 'x'")
  }
  if (identical(as.logical(distributions),FALSE)) distributions<-""
  if (length(distributions)==1) distributions<-make.words(distributions)
  matchdist<-match(toupper(distributions),lmrd.3par$distributions)
  if (any(is.na(matchdist))) stop("unknown distribution(s) ",
    paste(distributions[is.na(matchdist)],collapse=" "))
  if (missing(twopar)) twopar<-lmrd.3par$twopar[matchdist]
  else if (twopar==FALSE) twopar<-""
  twopar<-unique(make.words(twopar))
  match2<-pmatch(toupper(twopar),lmrd.2par$distributions)
  if (any(is.na(match2))) stop("unknown 2-parameter distribution(s)",
    paste(distributions[is.na(match2)],collapse=" "))
#
# Three-parameter distributions
#
  if (missing(xlim)) xlim<-range(0,0.6,x,na.rm=TRUE)
  if (missing(ylim)) ylim<-range(0,0.4,y,na.rm=TRUE)
  if (length(distributions)==0) {
    matplot(0,0,type="n",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,frame.plot=FALSE,...)
    legend.lmrd<-FALSE
  } else {
   col.lines <- if (!missing(col) && (x.missing || length(col)>1)) col else lmrd.3par$col[matchdist]
   if (missing(lty)) lty<-lmrd.3par$lty[matchdist]
   matplot(round(lmrd.data[,1],2), lmrd.data[,toupper(distributions)], type="l",
     xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, col=col.lines, lty=lty, lwd=lwd,
     frame.plot=FALSE, ...)
  }
#
# Framing elements of plot (box, lines at x=0 and y=0)
#
  fg<-list(...)$fg
  if (is.null(fg)) fg<-par("fg")
  box(col=fg)
  abline(h=0,v=0,col=fg)
#
# Two-parameter distributions
#
  if (length(twopar)>0) {
    # Restrict to points that lie within the plotting area
    pu<-par("usr")
    match2<- match2 * (
      lmrd.2par$tau3[match2]>=pu[1] & lmrd.2par$tau3[match2]<=pu[2] &
      lmrd.2par$tau4[match2]>=pu[3] & lmrd.2par$tau4[match2]<=pu[4] )
    # Plot the points, don't clip at axis box
    points(lmrd.2par$tau3[match2],lmrd.2par$tau4[match2],pch=15,col="black",xpd=TRUE)
    text(lmrd.2par$tau3[match2],lmrd.2par$tau4[match2],
         lmrd.2par$text[match2],adj=c(-0.5,-0.25),xpd=TRUE)
  }
#
# Legend
#
  if (isTRUE(legend.lmrd)) legend.lmrd<-list()
  if (is.list(legend.lmrd)) {
    parusr<-par("usr")
    if (missing(xlegend)) xlegend<-parusr[1]+0.01*(parusr[2]-parusr[1])
    if (missing(ylegend)) ylegend<-parusr[4]-0.01*(parusr[4]-parusr[3])*par("pin")[1]/par("pin")[2]
    legend.args<-list(x=xlegend,y=ylegend,legend=toupper(distributions),
      bty="n",col=col.lines,lty=lty,lwd=lwd)
    legend.args[names(legend.lmrd)]<-legend.lmrd
    do.call(legend,legend.args)
  }
#
# Data points
#
  if (!x.missing) {
    if (missing(cex)) cex<-NULL
    if (!missing(col) && length(col)==1) points(x,y,pch=pch,cex=cex,col=col)
    else points(x,y,pch=pch,cex=cex)
  }
}

make.words<-function(astring)
{
  zz<-unlist(strsplit(as.character(astring)," "))
  zz[zz!=""]
}

lmrd.data <- matrix(ncol=8, byrow=TRUE, data=c(
# tau_3    GLO      GEV      GPA      GNO      PE3    WAK.LB   ALL.LB
 -1.00,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,
 -0.99,  0.9834,  0.9761,  0.9752,  0.9775,  0.9752,  0.9752,  0.9751,
 -0.98,  0.9670,  0.9532,  0.9507,  0.9561,  0.9510,  0.9507,  0.9505,
 -0.97,  0.9508,  0.9311,  0.9267,  0.9356,  0.9271,  0.9265,  0.9261,
 -0.96,  0.9347,  0.9095,  0.9030,  0.9159,  0.9038,  0.9026,  0.9020,
 -0.95,  0.9188,  0.8885,  0.8796,  0.8968,  0.8809,  0.8791,  0.8781,
 -0.94,  0.9030,  0.8681,  0.8567,  0.8782,  0.8584,  0.8559,  0.8545,
 -0.93,  0.8874,  0.8480,  0.8340,  0.8600,  0.8364,  0.8331,  0.8311,
 -0.92,  0.8720,  0.8284,  0.8118,  0.8422,  0.8149,  0.8105,  0.8080,
 -0.91,  0.8568,  0.8092,  0.7899,  0.8247,  0.7938,  0.7883,  0.7851,
 -0.90,  0.8417,  0.7904,  0.7683,  0.8077,  0.7731,  0.7664,  0.7625,
 -0.89,  0.8268,  0.7720,  0.7471,  0.7909,  0.7529,  0.7448,  0.7401,
 -0.88,  0.8120,  0.7540,  0.7262,  0.7745,  0.7330,  0.7235,  0.7180,
 -0.87,  0.7974,  0.7363,  0.7057,  0.7583,  0.7136,  0.7026,  0.6961,
 -0.86,  0.7830,  0.7189,  0.6855,  0.7425,  0.6946,  0.6819,  0.6745,
 -0.85,  0.7688,  0.7019,  0.6657,  0.7269,  0.6761,  0.6616,  0.6531,
 -0.84,  0.7547,  0.6851,  0.6462,  0.7116,  0.6579,  0.6416,  0.6320,
 -0.83,  0.7408,  0.6687,  0.6270,  0.6966,  0.6401,  0.6218,  0.6111,
 -0.82,  0.7270,  0.6526,  0.6081,  0.6818,  0.6227,  0.6024,  0.5905,
 -0.81,  0.7134,  0.6368,  0.5896,  0.6672,  0.6057,  0.5833,  0.5701,
 -0.80,  0.7000,  0.6213,  0.5714,  0.6529,  0.5891,  0.5645,  0.5500,
 -0.79,  0.6868,  0.6061,  0.5536,  0.6389,  0.5729,  0.5460,  0.5301,
 -0.78,  0.6737,  0.5912,  0.5360,  0.6251,  0.5571,  0.5278,  0.5105,
 -0.77,  0.6608,  0.5765,  0.5188,  0.6115,  0.5416,  0.5099,  0.4911,
 -0.76,  0.6480,  0.5622,  0.5019,  0.5981,  0.5265,  0.4923,  0.4720,
 -0.75,  0.6354,  0.5480,  0.4853,  0.5850,  0.5118,  0.4750,  0.4531,
 -0.74,  0.6230,  0.5342,  0.4690,  0.5721,  0.4974,  0.4580,  0.4345,
 -0.73,  0.6108,  0.5206,  0.4530,  0.5594,  0.4834,  0.4412,  0.4161,
 -0.72,  0.5987,  0.5073,  0.4374,  0.5469,  0.4697,  0.4248,  0.3980,
 -0.71,  0.5868,  0.4942,  0.4220,  0.5347,  0.4564,  0.4087,  0.3801,
 -0.70,  0.5750,  0.4814,  0.4070,  0.5226,  0.4434,  0.3928,  0.3625,
 -0.69,  0.5634,  0.4689,  0.3922,  0.5107,  0.4308,  0.3773,  0.3451,
 -0.68,  0.5520,  0.4565,  0.3778,  0.4991,  0.4185,  0.3620,  0.3280,
 -0.67,  0.5408,  0.4445,  0.3636,  0.4877,  0.4065,  0.3470,  0.3111,
 -0.66,  0.5297,  0.4326,  0.3498,  0.4764,  0.3949,  0.3323,  0.2945,
 -0.65,  0.5188,  0.4210,  0.3362,  0.4654,  0.3836,  0.3179,  0.2781,
 -0.64,  0.5080,  0.4097,  0.3229,  0.4545,  0.3726,  0.3037,  0.2620,
 -0.63,  0.4974,  0.3986,  0.3100,  0.4439,  0.3619,  0.2899,  0.2461,
 -0.62,  0.4870,  0.3877,  0.2973,  0.4334,  0.3515,  0.2763,  0.2305,
 -0.61,  0.4768,  0.3770,  0.2849,  0.4231,  0.3414,  0.2630,  0.2151,
 -0.60,  0.4667,  0.3666,  0.2727,  0.4131,  0.3317,  0.2499,  0.2000,
 -0.59,  0.4568,  0.3564,  0.2609,  0.4032,  0.3222,  0.2372,  0.1851,
 -0.58,  0.4470,  0.3464,  0.2493,  0.3935,  0.3130,  0.2247,  0.1705,
 -0.57,  0.4374,  0.3366,  0.2380,  0.3840,  0.3041,  0.2125,  0.1561,
 -0.56,  0.4280,  0.3271,  0.2270,  0.3746,  0.2955,  0.2005,  0.1420,
 -0.55,  0.4188,  0.3177,  0.2163,  0.3655,  0.2872,  0.1888,  0.1281,
 -0.54,  0.4097,  0.3086,  0.2058,  0.3565,  0.2791,  0.1774,  0.1145,
 -0.53,  0.4008,  0.2997,  0.1956,  0.3478,  0.2713,  0.1663,  0.1011,
 -0.52,  0.3920,  0.2911,  0.1857,  0.3392,  0.2638,  0.1554,  0.0880,
 -0.51,  0.3834,  0.2826,  0.1761,  0.3307,  0.2566,  0.1448,  0.0751,
 -0.50,  0.3750,  0.2743,  0.1667,  0.3225,  0.2496,  0.1344,  0.0625,
 -0.49,  0.3668,  0.2663,  0.1575,  0.3144,  0.2428,  0.1243,  0.0501,
 -0.48,  0.3587,  0.2585,  0.1487,  0.3065,  0.2363,  0.1145,  0.0380,
 -0.47,  0.3508,  0.2508,  0.1401,  0.2988,  0.2301,  0.1049,  0.0261,
 -0.46,  0.3430,  0.2434,  0.1317,  0.2913,  0.2241,  0.0956,  0.0145,
 -0.45,  0.3354,  0.2362,  0.1236,  0.2839,  0.2183,  0.0866,  0.0031,
 -0.44,  0.3280,  0.2292,  0.1158,  0.2767,  0.2127,  0.0778, -0.0080,
 -0.43,  0.3208,  0.2224,  0.1082,  0.2697,  0.2074,  0.0692, -0.0189,
 -0.42,  0.3137,  0.2157,  0.1009,  0.2628,  0.2023,  0.0609, -0.0295,
 -0.41,  0.3068,  0.2093,  0.0938,  0.2562,  0.1974,  0.0529, -0.0399,
 -0.40,  0.3000,  0.2031,  0.0870,  0.2497,  0.1928,  0.0451, -0.0500,
 -0.39,  0.2934,  0.1971,  0.0804,  0.2433,  0.1883,  0.0375, -0.0599,
 -0.38,  0.2870,  0.1913,  0.0740,  0.2371,  0.1840,  0.0303, -0.0695,
 -0.37,  0.2808,  0.1856,  0.0679,  0.2311,  0.1800,  0.0232, -0.0789,
 -0.36,  0.2747,  0.1802,  0.0621,  0.2253,  0.1761,  0.0164, -0.0880,
 -0.35,  0.2688,  0.1750,  0.0565,  0.2196,  0.1724,  0.0098, -0.0969,
 -0.34,  0.2630,  0.1699,  0.0511,  0.2141,  0.1689,  0.0035, -0.1055,
 -0.33,  0.2574,  0.1651,  0.0459,  0.2088,  0.1656, -0.0025, -0.1139,
 -0.32,  0.2520,  0.1604,  0.0410,  0.2036,  0.1624, -0.0084, -0.1220,
 -0.31,  0.2468,  0.1559,  0.0364,  0.1986,  0.1594, -0.0139, -0.1299,
 -0.30,  0.2417,  0.1516,  0.0319,  0.1937,  0.1566, -0.0193, -0.1375,
 -0.29,  0.2368,  0.1475,  0.0277,  0.1890,  0.1539, -0.0244, -0.1449,
 -0.28,  0.2320,  0.1436,  0.0237,  0.1845,  0.1514, -0.0293, -0.1520,
 -0.27,  0.2274,  0.1399,  0.0200,  0.1802,  0.1490, -0.0339, -0.1589,
 -0.26,  0.2230,  0.1364,  0.0165,  0.1759,  0.1467, -0.0383, -0.1655,
 -0.25,  0.2188,  0.1330,  0.0132,  0.1719,  0.1446, -0.0425, -0.1719,
 -0.24,  0.2147,  0.1298,  0.0101,  0.1680,  0.1426, -0.0464, -0.1780,
 -0.23,  0.2108,  0.1269,  0.0072,  0.1643,  0.1408, -0.0501, -0.1839,
 -0.22,  0.2070,  0.1241,  0.0046,  0.1607,  0.1390, -0.0536, -0.1895,
 -0.21,  0.2034,  0.1214,  0.0022,  0.1573,  0.1374, -0.0568, -0.1949,
 -0.20,  0.2000,  0.1190,  0.0000,  0.1541,  0.1358, -0.0599, -0.2000,
 -0.19,  0.1968,  0.1167, -0.0020,  0.1510,  0.1344, -0.0626, -0.2049,
 -0.18,  0.1937,  0.1147, -0.0037,  0.1481,  0.1331, -0.0652, -0.2095,
 -0.17,  0.1908,  0.1127, -0.0053,  0.1453,  0.1319, -0.0675, -0.2139,
 -0.16,  0.1880,  0.1110, -0.0066,  0.1427,  0.1307, -0.0696, -0.2180,
 -0.15,  0.1854,  0.1095, -0.0077,  0.1403,  0.1297, -0.0715, -0.2219,
 -0.14,  0.1830,  0.1081, -0.0086,  0.1380,  0.1287, -0.0732, -0.2255,
 -0.13,  0.1808,  0.1069, -0.0093,  0.1359,  0.1278, -0.0746, -0.2289,
 -0.12,  0.1787,  0.1059, -0.0098,  0.1339,  0.1270, -0.0759, -0.2320,
 -0.11,  0.1768,  0.1051, -0.0101,  0.1321,  0.1263, -0.0769, -0.2349,
 -0.10,  0.1750,  0.1044, -0.0102,  0.1305,  0.1256, -0.0776, -0.2375,
 -0.09,  0.1734,  0.1039, -0.0101,  0.1290,  0.1250, -0.0782, -0.2399,
 -0.08,  0.1720,  0.1036, -0.0098,  0.1276,  0.1245, -0.0786, -0.2420,
 -0.07,  0.1708,  0.1034, -0.0092,  0.1265,  0.1241, -0.0787, -0.2439,
 -0.06,  0.1697,  0.1035, -0.0085,  0.1254,  0.1237, -0.0786, -0.2455,
 -0.05,  0.1688,  0.1037, -0.0076,  0.1246,  0.1233, -0.0783, -0.2469,
 -0.04,  0.1680,  0.1040, -0.0065,  0.1239,  0.1231, -0.0778, -0.2480,
 -0.03,  0.1674,  0.1046, -0.0051,  0.1233,  0.1229, -0.0771, -0.2489,
 -0.02,  0.1670,  0.1053, -0.0036,  0.1229,  0.1227, -0.0761, -0.2495,
 -0.01,  0.1668,  0.1061, -0.0019,  0.1227,  0.1226, -0.0750, -0.2499,
  0.00,  0.1667,  0.1072,  0.0000,  0.1226,  0.1226, -0.0737, -0.2500,
  0.01,  0.1667,  0.1084,  0.0021,  0.1227,  0.1226, -0.0721, -0.2499,
  0.02,  0.1670,  0.1098,  0.0044,  0.1229,  0.1227, -0.0703, -0.2495,
  0.03,  0.1674,  0.1113,  0.0069,  0.1233,  0.1229, -0.0683, -0.2489,
  0.04,  0.1680,  0.1131,  0.0095,  0.1239,  0.1231, -0.0662, -0.2480,
  0.05,  0.1687,  0.1149,  0.0124,  0.1246,  0.1233, -0.0638, -0.2469,
  0.06,  0.1697,  0.1170,  0.0154,  0.1254,  0.1237, -0.0612, -0.2455,
  0.07,  0.1707,  0.1192,  0.0186,  0.1265,  0.1241, -0.0584, -0.2439,
  0.08,  0.1720,  0.1216,  0.0220,  0.1276,  0.1245, -0.0554, -0.2420,
  0.09,  0.1734,  0.1241,  0.0256,  0.1290,  0.1250, -0.0522, -0.2399,
  0.10,  0.1750,  0.1269,  0.0294,  0.1305,  0.1256, -0.0488, -0.2375,
  0.11,  0.1767,  0.1297,  0.0334,  0.1321,  0.1263, -0.0452, -0.2349,
  0.12,  0.1787,  0.1328,  0.0375,  0.1339,  0.1270, -0.0414, -0.2320,
  0.13,  0.1807,  0.1360,  0.0418,  0.1359,  0.1278, -0.0374, -0.2289,
  0.14,  0.1830,  0.1393,  0.0463,  0.1380,  0.1287, -0.0332, -0.2255,
  0.15,  0.1854,  0.1429,  0.0510,  0.1403,  0.1297, -0.0289, -0.2219,
  0.16,  0.1880,  0.1466,  0.0558,  0.1427,  0.1307, -0.0243, -0.2180,
  0.17,  0.1907,  0.1504,  0.0608,  0.1453,  0.1319, -0.0195, -0.2139,
  0.18,  0.1937,  0.1544,  0.0660,  0.1481,  0.1331, -0.0145, -0.2095,
  0.19,  0.1967,  0.1586,  0.0714,  0.1510,  0.1344, -0.0094, -0.2049,
  0.20,  0.2000,  0.1629,  0.0769,  0.1541,  0.1358, -0.0040, -0.2000,
  0.21,  0.2034,  0.1674,  0.0826,  0.1573,  0.1374,  0.0016, -0.1949,
  0.22,  0.2070,  0.1721,  0.0885,  0.1607,  0.1390,  0.0073, -0.1895,
  0.23,  0.2107,  0.1769,  0.0946,  0.1643,  0.1408,  0.0132, -0.1839,
  0.24,  0.2147,  0.1818,  0.1008,  0.1680,  0.1426,  0.0193, -0.1780,
  0.25,  0.2188,  0.1870,  0.1071,  0.1719,  0.1446,  0.0257, -0.1719,
  0.26,  0.2230,  0.1922,  0.1137,  0.1759,  0.1467,  0.0321, -0.1655,
  0.27,  0.2274,  0.1977,  0.1204,  0.1802,  0.1490,  0.0388, -0.1589,
  0.28,  0.2320,  0.2033,  0.1273,  0.1845,  0.1514,  0.0457, -0.1520,
  0.29,  0.2367,  0.2090,  0.1343,  0.1890,  0.1539,  0.0528, -0.1449,
  0.30,  0.2417,  0.2150,  0.1415,  0.1937,  0.1566,  0.0600, -0.1375,
  0.31,  0.2467,  0.2210,  0.1489,  0.1986,  0.1594,  0.0674, -0.1299,
  0.32,  0.2520,  0.2272,  0.1564,  0.2036,  0.1624,  0.0750, -0.1220,
  0.33,  0.2574,  0.2336,  0.1641,  0.2088,  0.1656,  0.0828, -0.1139,
  0.34,  0.2630,  0.2402,  0.1719,  0.2141,  0.1689,  0.0908, -0.1055,
  0.35,  0.2687,  0.2469,  0.1799,  0.2196,  0.1724,  0.0990, -0.0969,
  0.36,  0.2747,  0.2537,  0.1881,  0.2253,  0.1761,  0.1073, -0.0880,
  0.37,  0.2807,  0.2607,  0.1964,  0.2311,  0.1800,  0.1158, -0.0789,
  0.38,  0.2870,  0.2678,  0.2048,  0.2371,  0.1840,  0.1245, -0.0695,
  0.39,  0.2934,  0.2751,  0.2135,  0.2433,  0.1883,  0.1334, -0.0599,
  0.40,  0.3000,  0.2826,  0.2222,  0.2497,  0.1928,  0.1424, -0.0500,
  0.41,  0.3067,  0.2902,  0.2311,  0.2562,  0.1974,  0.1517, -0.0399,
  0.42,  0.3137,  0.2980,  0.2402,  0.2628,  0.2023,  0.1611, -0.0295,
  0.43,  0.3207,  0.3059,  0.2494,  0.2697,  0.2074,  0.1707, -0.0189,
  0.44,  0.3280,  0.3140,  0.2588,  0.2767,  0.2127,  0.1804, -0.0080,
  0.45,  0.3354,  0.3222,  0.2683,  0.2839,  0.2183,  0.1904,  0.0031,
  0.46,  0.3430,  0.3306,  0.2780,  0.2913,  0.2241,  0.2005,  0.0145,
  0.47,  0.3507,  0.3391,  0.2878,  0.2988,  0.2301,  0.2108,  0.0261,
  0.48,  0.3587,  0.3478,  0.2978,  0.3065,  0.2363,  0.2212,  0.0380,
  0.49,  0.3667,  0.3566,  0.3079,  0.3144,  0.2428,  0.2319,  0.0501,
  0.50,  0.3750,  0.3655,  0.3182,  0.3225,  0.2496,  0.2427,  0.0625,
  0.51,  0.3834,  0.3747,  0.3286,  0.3307,  0.2566,  0.2536,  0.0751,
  0.52,  0.3920,  0.3839,  0.3391,  0.3392,  0.2638,  0.2648,  0.0880,
  0.53,  0.4007,  0.3934,  0.3498,  0.3478,  0.2713,  0.2761,  0.1011,
  0.54,  0.4097,  0.4029,  0.3606,  0.3565,  0.2791,  0.2876,  0.1145,
  0.55,  0.4187,  0.4127,  0.3716,  0.3655,  0.2872,  0.2993,  0.1281,
  0.56,  0.4280,  0.4225,  0.3827,  0.3746,  0.2955,  0.3111,  0.1420,
  0.57,  0.4374,  0.4325,  0.3940,  0.3840,  0.3041,  0.3231,  0.1561,
  0.58,  0.4470,  0.4427,  0.4054,  0.3935,  0.3130,  0.3353,  0.1705,
  0.59,  0.4567,  0.4530,  0.4169,  0.4032,  0.3222,  0.3476,  0.1851,
  0.60,  0.4667,  0.4635,  0.4286,  0.4131,  0.3317,  0.3601,  0.2000,
  0.61,  0.4767,  0.4741,  0.4404,  0.4231,  0.3414,  0.3728,  0.2151,
  0.62,  0.4870,  0.4848,  0.4523,  0.4334,  0.3515,  0.3856,  0.2305,
  0.63,  0.4974,  0.4957,  0.4644,  0.4439,  0.3619,  0.3986,  0.2461,
  0.64,  0.5080,  0.5068,  0.4766,  0.4545,  0.3726,  0.4118,  0.2620,
  0.65,  0.5187,  0.5180,  0.4889,  0.4654,  0.3836,  0.4251,  0.2781,
  0.66,  0.5297,  0.5293,  0.5014,  0.4764,  0.3949,  0.4387,  0.2945,
  0.67,  0.5407,  0.5408,  0.5140,  0.4877,  0.4065,  0.4523,  0.3111,
  0.68,  0.5520,  0.5524,  0.5268,  0.4991,  0.4185,  0.4662,  0.3280,
  0.69,  0.5634,  0.5642,  0.5396,  0.5107,  0.4308,  0.4802,  0.3451,
  0.70,  0.5750,  0.5761,  0.5526,  0.5226,  0.4434,  0.4944,  0.3625,
  0.71,  0.5867,  0.5882,  0.5658,  0.5347,  0.4564,  0.5087,  0.3801,
  0.72,  0.5987,  0.6004,  0.5790,  0.5469,  0.4697,  0.5232,  0.3980,
  0.73,  0.6107,  0.6127,  0.5924,  0.5594,  0.4834,  0.5379,  0.4161,
  0.74,  0.6230,  0.6252,  0.6059,  0.5721,  0.4974,  0.5527,  0.4345,
  0.75,  0.6354,  0.6379,  0.6196,  0.5850,  0.5118,  0.5677,  0.4531,
  0.76,  0.6480,  0.6507,  0.6333,  0.5981,  0.5265,  0.5829,  0.4720,
  0.77,  0.6607,  0.6636,  0.6472,  0.6115,  0.5416,  0.5982,  0.4911,
  0.78,  0.6737,  0.6767,  0.6612,  0.6251,  0.5571,  0.6137,  0.5105,
  0.79,  0.6867,  0.6899,  0.6754,  0.6389,  0.5729,  0.6294,  0.5301,
  0.80,  0.7000,  0.7032,  0.6897,  0.6529,  0.5891,  0.6452,  0.5500,
  0.81,  0.7134,  0.7167,  0.7040,  0.6672,  0.6057,  0.6612,  0.5701,
  0.82,  0.7270,  0.7304,  0.7186,  0.6818,  0.6227,  0.6774,  0.5905,
  0.83,  0.7407,  0.7442,  0.7332,  0.6966,  0.6401,  0.6937,  0.6111,
  0.84,  0.7547,  0.7581,  0.7479,  0.7116,  0.6579,  0.7102,  0.6320,
  0.85,  0.7687,  0.7722,  0.7628,  0.7269,  0.6761,  0.7269,  0.6531,
  0.86,  0.7830,  0.7864,  0.7778,  0.7425,  0.6946,  0.7437,  0.6745,
  0.87,  0.7974,  0.8007,  0.7929,  0.7583,  0.7136,  0.7608,  0.6961,
  0.88,  0.8120,  0.8152,  0.8082,  0.7745,  0.7330,  0.7780,  0.7180,
  0.89,  0.8267,  0.8298,  0.8235,  0.7909,  0.7529,  0.7953,  0.7401,
  0.90,  0.8417,  0.8446,  0.8390,  0.8077,  0.7731,  0.8129,  0.7625,
  0.91,  0.8567,  0.8595,  0.8546,  0.8247,  0.7938,  0.8307,  0.7851,
  0.92,  0.8720,  0.8746,  0.8703,  0.8422,  0.8149,  0.8486,  0.8080,
  0.93,  0.8874,  0.8898,  0.8861,  0.8600,  0.8364,  0.8667,  0.8311,
  0.94,  0.9030,  0.9051,  0.9020,  0.8782,  0.8584,  0.8850,  0.8545,
  0.95,  0.9187,  0.9206,  0.9181,  0.8968,  0.8809,  0.9036,  0.8781,
  0.96,  0.9347,  0.9362,  0.9342,  0.9159,  0.9038,  0.9223,  0.9020,
  0.97,  0.9507,  0.9519,  0.9505,  0.9356,  0.9271,  0.9413,  0.9261,
  0.98,  0.9670,  0.9678,  0.9669,  0.9561,  0.9510,  0.9605,  0.9505,
  0.99,  0.9834,  0.9838,  0.9834,  0.9775,  0.9752,  0.9801,  0.9751,
  1.00,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000,  1.0000),
dimnames=list(NULL,c("tau_3","GLO","GEV","GPA","GNO","PE3","WAK.LB","ALL.LB")))

lmrd.3par<-list(
  distributions=c("GLO","GEV","GPA","GNO","PE3","WAK.LB","ALL.LB"),
  col=c("blue","#00C000","red","black","#00E0E0","red","black"),
  lty=c("solid","solid","solid","solid","solid","longdash","longdash"),
  twopar=c("L","G","E U","N","E N","",""))

lmrd.2par<-list(
  distributions=c("EXP","GUM","LOG","NOR","UNI"),
  tau3=c(0.3333,0.1699,0,0,0),
  tau4=c(0.1667,0.1504,0.1667,0.1226,0),
  text=c("E","G","L","N","U")
)

lmrdpoints<-function(x, y=NULL, type="p", ...) {
## Add points to an L-moment ratio diagram
  if (is.null(y)) {
    if (inherits(x,"regdata")) {
      y<-x[[6]]; x<-x[[5]]
    } else if (length(dim(x))==2) {
      dn2<-dimnames(x)[[2]]
      mm<-match(c("t_3","t_4"),dn2)
      if (any(is.na(mm))) mm<-match(c("tau_3","tau_4"),dn2)
      if (!any(is.na(mm))) { y<-x[,mm[2] ]; x<-x[,mm[1] ] }
    } else if (is.numeric(x)) {
      nx<-names(x)
      mm<-match(c("t_3","t_4"),nx)
      if (any(is.na(mm))) mm<-match(c("tau_3","tau_4"),nx)
      if (!any(is.na(mm))) { y<-x[mm[2] ]; x<-x[mm[1] ] }
    }
    if (is.null(y))
      stop("could not find L-skewness and L-kurtosis values in 'x'")
  }
  points(x=x, y=y, type=type, ...)
}

lmrdlines<-function(x, y=NULL, type="l", ...) lmrdpoints(x=x, y=y, type=type, ...)

#-------------------------------------------------------------------------------
#  Functions for extreme-value plot
#-------------------------------------------------------------------------------

evplot<-function(y,...)
UseMethod("evplot")

evplot.default<-function(y,qfunc,para,npoints=101,plim,xlim=c(-2,5),ylim,type,
        xlab=expression("Reduced variate,  " * -log(-log(italic(F)))),
        ylab="Quantile", rp.axis=TRUE, ...) {

  if (!missing(plim)) xlim<-c(-log(-log(plim)))
  missing.ylim<-missing(ylim)
  if (missing.ylim) ylim<-c(0,1) # now missing(ylim) is TRUE in S but not in R
  if (missing(y)) {
    xx<-0; yy<-0; type<-"n"
  } else {
    yy<-sort(y[!is.na(y)])
    lyy<-length(yy)
    xx<--log(-log(((1:lyy)-0.44)/(lyy+0.12)))
    if (missing.ylim) ylim<-c(min(0,yy),max(yy))
    if (missing(plim)) xlim<-c(min(xlim[1],xx[1]),max(xlim[2],xx[lyy]))
    if (missing(type)) type<-"p"
  }

  if (!missing(qfunc) && missing.ylim ) {
    xval<-range(xlim)
    pval<-exp(-exp(-xval))
    yval <- if (missing(para)) qfunc(pval) else
      if (is.list(para)) do.call(qfunc,c(list(pval),para)) else qfunc(pval,para)
    ylim<-range(c(ylim,yval))
  }

  plot(xx,yy,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,type=type,...)
  if (!missing(qfunc)) evdistq(qfunc,para,npoints=npoints,...)

#  Return period axis:
#  Define tick marks for a fixed set of return periods.
#  Use only those tick marks that lie within the x-axis limits.
#  axis(...) draws the axis, 5% of the distance from ymin to ymax.
#  text(...) plots the axis label, centred above the midpoint of the axis
#  and as far from the axis line as the x and y axis labels are.
  if (rp.axis) {
    parusr<-par("usr")
    rp.lab<-c(2,5,10,20,50,100,200,500,1000,10000,100000,1e6,1e7,1e8)
    rp.tic<- -log(-log(1-1/rp.lab))
    crit<-(rp.tic>=parusr[1] & rp.tic<=parusr[2])
    rp.tic<-rp.tic[crit]
    rp.lab<-rp.lab[crit]
    rp.ypos<-parusr[3]+(parusr[4]-parusr[3])*0.05
    axis(side=3,at=rp.tic,labels=rp.lab,pos=rp.ypos)
    text((min(rp.tic)+max(rp.tic))*0.5,rp.ypos+par("cxy")[2]*par("mgp")[1],
      "Return period",adj=c(0.5,0))
  }

  invisible()
}

evpoints<-function(y,...) {
  yval<-sort(y[!is.na(y)])
  n<-length(yval)
  xval<--log(-log(((1:n)-0.44)/(n+0.12)))
  points(xval,yval,...)
}

evdistq<-function(qfunc,para,npoints=101,...) {
  parusr<-par("usr")
  xval<-seq(from=parusr[1],to=parusr[2],length=npoints)
  pval<-exp(-exp(-xval))
  yval <- if (missing(para)) qfunc(pval) else
    if (is.list(para)) do.call(qfunc,c(list(pval),para)) else qfunc(pval,para)
  lines(xval,yval,...)
}

evdistp<-function(pfunc,para,npoints=101,...) {
  parusr<-par("usr")
  yval<-seq(from=parusr[3],to=parusr[4],length=npoints)
  pval<-if (missing(para)) pfunc(yval) else
    if (is.list(para)) do.call(pfunc,c(list(yval),para)) else pfunc(yval,para)
  xval<- -log(-log(pval))
  lines(xval,yval,...)
}

#-------------------------------------------------------------------------------
#  lmom package's own version of integrate()
#-------------------------------------------------------------------------------

#  Differences from R's own integrate():
#
#  o Does not have the bug that affected R's own integrate() between versions
#    2.12.0 and 3.0.1 (R bug PR#15219).
#
#  o Different treatment of 'abs.tol' and 'rel.tol'.  If only one is supplied,
#    the other is set to 0. (In R, setting just 'abs.tol' leaves 'rel.tol' at
#    its default, which if 'abs.tol' is small often means that 'rel.tol' will
#    still control the test for convergence.)
#
#  o Does not have the R version's unused arguments ('keep.xy' and 'aux').

integrate<- function(f, lower, upper, ..., subdivisions=100,
  rel.tol, abs.tol, stop.on.error = TRUE)
{
  f <- match.fun(f)

  if (is.na(lower) || is.na(upper)) stop("a limit is missing")
  if (lower>upper) stop("limits incorrectly ordered")
  if (subdivisions<1) stop("'subdivisions' must be at least 1")
  if (missing(rel.tol) && missing(abs.tol)) abs.tol<-rel.tol<-.Machine$double.eps^0.25
  else if (missing(rel.tol)) rel.tol<-0
  else if (missing(abs.tol)) abs.tol<-0
  if (abs.tol<=0 && rel.tol<max(50*.Machine$double.eps, 0.5e-28))
    stop("invalid tolerance values")

  env<-new.env(parent=environment())
  assign("expr", quote(f(x, ...)), envir=env)

  a<-as.double(lower)
  b<-as.double(upper)
  limit<-as.integer(subdivisions)
  epsabs<-as.double(abs.tol)
  epsrel<-as.double(rel.tol)

  result<-numeric(1)
  abserr<-numeric(1)
  neval<-integer(1)
  ier<-integer(1)
  last<-integer(1)

  if (is.finite(lower) && is.finite(upper)) {

    cout <- .Call("cdqagse", PACKAGE="lmom",
      a, b, epsabs, epsrel, limit, result, abserr, neval, ier, last, env)

  } else {

    if (is.finite(lower)) {inf<-1L; bound<-a }        # upper bound (only) infinite
    else if (is.finite(upper)) {inf<- -1L; bound<-b } # lower bound (only) infinite
    else {inf<-2L; bound<-0 }                         # both bounds infinite

    cout <- .Call("cdqagie", PACKAGE="lmom",
      bound, inf, epsabs, epsrel, limit, result, abserr, neval, ier, last, env)
  }

  res<-list(value=result, abs.error=abserr, subdivisions=last)
  res$message <- switch(ier + 1,
    "OK",
    "maximum number of subdivisions reached",
    "roundoff error was detected",
    "extremely bad integrand behaviour",
    "roundoff error is detected in the extrapolation table",
    "the integral is probably divergent",
    "the input is invalid")
  if(ier==6 || (ier>0 && stop.on.error)) stop(res$message)
  res$call<-match.call()
  class(res)<-"integrate"
  res
}


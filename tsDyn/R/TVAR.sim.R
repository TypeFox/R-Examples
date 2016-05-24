TVAR.gen <- function(data,B,TVARobject, Thresh, nthresh=1, type=c("simul","boot", "check"), 
                     n=200, lag=1, include = c("const", "trend","none", "both"),  
                     thDelay=1,  thVar=NULL, mTh=1, starting=NULL, 
                     innov=rmnorm(n, varcov=varcov), varcov=diag(1,k), 
                     show.parMat=FALSE, round=FALSE, seed){
  
  ###check correct arguments
  if(!nthresh%in%c(0,1,2))
    stop("Arg nthresh should  be either 0, 1 or 2")
  if(!missing(n)&any(!missing(data), !missing(TVARobject)))
    stop("arg n should not be given with arg data or TVARobject")
  if(!missing(TVARobject)&any(!missing(Thresh), !missing(nthresh), !missing(lag), !missing(thDelay), !missing(mTh)))
    warning("When object TVARobject is given, only args 'type' and 'round' are relevant, others are not considered")
  if(!missing(data)&!missing(B))
    stop("You have to provide either B or y, but not both")
  
  type<-match.arg(type)
  include<-match.arg(include)
  
  ## Create some variables/parameters
  p <- switch(type, "boot"=if(!missing(TVARobject)) TVARobject$lag else lag, lag)
  ninc <- switch(include, "none"=0, "const"=1, "trend"=1, "both"=2)
  add <- switch(type, "simul"=p, 0)
  if(max(thDelay)>p)
    stop("Max of thDelay should be smaller or equal to the number of lags")
  
  
  ### possibility 1: only parameters matrix is given
  if(!missing(B)){
    if(type!="simul"){
      type<-"simul"
      warning("Type check or boot are only avalaible with pre specified data. The type simul was used")
    }
    nB<-nrow(B)
    ndig<-4
    esp<-p*nB+ninc 
    if(esp*(nthresh+1)!=ncol(B)){
      stop("Matrix B badly specified: expected ", esp*(nthresh+1), " elements ( (lag*K+ n inc)* (nthresh+1) ) but has ", ncol(B), "\n" )
    }
    y<-matrix(0,ncol=nB, nrow=n+add)
    if(!is.null(starting)){
      if(all(dim(as.matrix(starting))==c(p,nB)))
        y[seq_len(p),]<- as.matrix(starting)
      else
        stop("Bad specification of starting values. Should have nrow = lag and ncol = number of variables. But is: ", dim(as.matrix(starting)), sep="")
    }
    Bmat<-B
    temp<-TVAR_thresh(mTh=mTh,thDelay=thDelay,thVar=thVar,y=y, p=p) #Stored in misc.R
    trans<-temp$trans
    combin<-temp$combin
    
    if(nthresh==1){
      if(missing(Thresh))
        Thresh<-mean(trans)
      if(length(Thresh)!=1){
        warning("Please only one Thresh value if you choose nthresh=1. First one was chosen")
        Thresh<-Thresh[1]}
    } else  if(nthresh==2){
      if(missing(Thresh))
        Thresh<-quantile(trans, probs=c(0.25, 0.75))
      if(length(Thresh)!=2)
        stop("please give two Thresh values if you choose nthresh=2")
    }
    if(nthresh%in%c(1,2))
      Thresh<-round(Thresh,ndig)
    k <- nB 		#Number of variables
    T <- nrow(y) 		#Size of start sample
    T <- n 		#Size of start sample
  }

  ### possibility 2: only data is given: compute it with linear or selectSETAR
  else if(!missing(data)){
    if(nthresh==0){
      TVARobject<-lineVar(data, lag=p, include=include, model="VAR")
    } else{ 
      if(!missing(Thresh)){
        TVARobject<-TVAR(data, lag=p, include=include, nthresh=nthresh, plot=FALSE, trace=FALSE, gamma=Thresh)
      } else{
        TVARobject<-TVAR(data, lag=p, include=include, nthresh=nthresh, plot=FALSE, trace=FALSE)
      }
    }
  }
  ### possibility 3: setarobject is given by user (or by poss 2)
  if(!missing(TVARobject)){
    k<-TVARobject$k
    T<-TVARobject$T
    p<-TVARobject$lag
    modSpe<-TVARobject$model.specific
    res<-residuals(TVARobject)
    Bmat<-coefMat(TVARobject)
    y<-as.matrix(TVARobject$model)[,1:k]
    ndig<-getndp(y[,1])
    nthresh<-modSpe$nthresh
    include <- TVARobject$include
    if(nthresh>0){
      if(modSpe$oneMatrix)
        stop("arg commoninter in TVAR currently not implemented in TVAR.sim")
      if(attr(TVARobject, "levelTransVar")=="MTAR")
        stop("arg model=MTAR in TVAR currently not implemented in TVAR.sim")
      if(attr(TVARobject, "transVar")=="external")
        stop("arg thVaR in TVAR currently not implemented in TVAR.sim")
      Thresh<-modSpe$Thresh
      thDelay<-modSpe$thDelay
      combin<-modSpe$transCombin
    }
  }

  t <- T-p+add 		#Size of end sample
  npar<-k*(p+ninc)

  ##### put coefficients vector in right form according to arg include (arg both need no modif)
  a<-NULL
  if(include=="none")
    for(i in 0:nthresh) a<-c(a, i*(p*k)+c(1,2))
  else if(include=="const")
    for(i in 0:nthresh) a<-c(a, i*(p*k+2)+c(2))
  else if(include=="trend")
    for(i in 0:nthresh) a<-c(a, i*(p*k+2)+c(1))
  #if (include=="both"): correction useless
  Bmat<-myInsertCol(Bmat, c=a ,0)
  nparBmat<-p*k+2
    
    
  ##############################
  ###Reconstitution boot/simul
  ##############################

  #initial values
  Yb<-matrix(0, nrow=T+add, ncol=k)
  Yb[1:p,]<-y[1:p,]		

  if(nthresh>0){
    z2<-vector("numeric", length=nrow(y))
    z2[1:p]<-y[1:p,]%*%combin
  }

  trend<-1:(T+add) #   trend<-c(NA, 1:(T-1+add))

  ##Innovations Specifications 
  if(type=="simul"&&dim(innov)!=c(n,k))
    stop(paste("input innov is not of right dim, should be matrix with", n,"rows and ", k, "cols\n"))
  if(!missing(seed)) set.seed(seed)
  innov<-switch(type, "boot"=if(missing(innov)) res[sample(seq_len(t), replace=TRUE),] else innov, "simul"=innov, "check"=res)
  resb<-rbind(matrix(0,nrow=p, ncol=k),innov)	

  ## MAIN loop:

  if(nthresh==0){
	  for(i in (p+1):(T+add)){
		  Yb[i,]<-rowSums(cbind(Bmat[,1], Bmat[,2]*trend[i], Bmat[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
	  }
  } else if(nthresh==1){
	  BDown<-Bmat[,seq_len(nparBmat)]
	  BUp<-Bmat[,-seq_len(nparBmat)]

	  for(i in (p+1):(T+add)){
	    if(round(z2[i-thDelay],ndig)<=Thresh) {
	      Yb[i,]<-rowSums(cbind(BDown[,1], BDown[,2]*trend[i], BDown[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
	      if(round)
	        Yb[i,]<-round(Yb[i,], ndig)
	    }
	    else{
	      Yb[i,]<-rowSums(cbind(BUp[,1], BUp[,2]*trend[i], BUp[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
	      if(round)
	        Yb[i,]<-round(Yb[i,], ndig)
	    }
	    z2[i]<-Yb[i,]%*%combin
	  }
  }  else if(nthresh==2){
    BDown <- Bmat[,seq_len(nparBmat)]
    BMiddle <- Bmat[,seq_len(nparBmat)+nparBmat]
    BUp <- Bmat[,seq_len(nparBmat)+2*nparBmat]
    for(i in (p+1):(T+add)){
      if(round(z2[i-thDelay],ndig)<=Thresh[1]) 
        Yb[i,]<-rowSums(cbind(BDown[,1], BDown[,2]*trend[i], BDown[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
      else if(round(z2[i-thDelay],ndig)>Thresh[2]) 
        Yb[i,]<-rowSums(cbind(BUp[,1], BUp[,2]*trend[i], BUp[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
      else
        Yb[i,]<-rowSums(cbind(BMiddle[,1],BMiddle[,2]*trend[i], BMiddle[,-c(1,2)]%*%matrix(t(Yb[i-c(1:p),]), ncol=1),resb[i,]))
      z2[i]<-Yb[i,]%*%combin
    }
  }
  
  
  if(FALSE){
	  new<-TVAR_thresh(y=Yb,mTh=mTh,thDelay=thDelay,thVar=thVar, p=p)$trans
	  if(mean(ifelse(new[,thDelay]<Thresh, 1,0))>0.05)
		  list(B=Bmat,Yb=Yb)
	  else
		  Recall(B=Bmat, n=n, lag=lag, nthresh=nthresh, thDelay=thDelay, Thresh, thVar=thVar, mTh=mTh, type=type, starting=starting)
  }
  # print(cbind(y, round(Yb,3)))

  if(type=="simul"){
    Yb <- tail(Yb, -p)
    rownames(Yb) <- seq_len(nrow(Yb))
  }

  if(show.parMat)
    print(Bmat)
  res<-if(round) round(Yb, ndig) else Yb
  return(res)
}





#'Simulation of a multivariate Threshold Autoregressive model (TVAR)
#'
#'Simulate a multivariate Threshold VAR (TVAR)
#'
#'This function offers the possibility to generate series following a TVAR.
#'
#'By giving a matrix of coefficients, on can only simulate a VAR (nthresh=0) or
#'TVAR (nthresh=1 or 2). One can have a specification with constant (default),
#'trend, both or none (see arg include). Order in parameters is include/lags
#'(VECM) and include/lags/include/lags for TVECM, hence, a matrix for a TVECM
#'with 3 regimes, a const and a 2 lags would have 2 lines and 2*(1+4) columns.
#'The innovations can be given by the user (a matrix of dim nxk, here n does
#'not include the starting values!), by default it uses a multivariate normal
#'distribution, with covariance matrix specified by varcov. The starting values
#'(of dim lags x k) can be given. The user should take care for their choice,
#'since it is not sure that the simulated values will cross the threshold even
#'once.
#'
#'The matrix \sQuote{B} has to be in the form: constant, trend, lags, then
#'repeated if many regimes. In case of uncertainty, using
#'\code{who.parMat=TRUE} will print the matrix as interpreted by the function,
#'helping the user to feed the right input.
#'
#'@aliases TVAR.sim VAR.sim
#'@param B Matrix of coefficients to simulate
#'@param Thresh The threshold value(s). Vector of length nthresh
#'@param nthresh number of threshold (see details)
#'@param n Number of observations to create when type="simul"
#'@param lag Number of lags to include in each regime
#'@param include Type of deterministic regressors to include. NOT WORKING
#'PROPERLY CURRENTLY if not const
#'@param thDelay 'time delay' for the threshold variable (as multiple of
#'embedding time delay d) PLEASE NOTE that the notation is currently different
#'to univariate models in tsDyn. The left side variable is taken at time t, and
#'not t+1 as in univariate cases.
#'@param thVar external transition variable
#'@param mTh combination of variables with same lag order for the transition
#'variable. Either a single value (indicating which variable to take) or a
#'combination
#'@param starting Starting values (matrix of dimension lag x k). If not given,
#'set to zero.
#'@param innov Innovations used for simulation. Should be matrix of dim n x k.
#'By default multivariate normal.
#'@param varcov Variance-covariance matrix for the innovations. By default
#'identity matrix.
#'@param show.parMat Logical. Should the parameter matrix be shown? Useful to
#'understand how to give right input
#'@param round Rounds the series created to have the same digits (hopefully) as
#'original series.
#'@param seed Optional. Seed for the random number generation.
#'@param \dots Further arguments passed to the underlying (un-exported)
#'\code{TVAR.gen} function
#'@return A matrix with the simulated/bootstraped series.
#'@author Matthieu Stigler
#'@seealso \code{\link{TVAR}} to estimate a TVAR, \code{\link{VAR.sim}} to
#'simulate/bootstrap a VAR, \code{\link{TVECM.sim}} to simulate/bootstrap a
#'TVECM.
#'@keywords ts bootstrap
#'@export
#'@examples
#'
#'##simulate VAR as in Enders 2004, p 268
#'B1<-matrix(c(0.7, 0.2, 0.2, 0.7), 2)
#'var1<-TVAR.sim(B=B1,nthresh=0,n=100, include="none")
#'ts.plot(var1, type="l", col=c(1,2))
#'
#'
#'B2<-rbind(c(0.5, 0.5, 0.5), c(0, 0.5, 0.5))
#'varcov<-matrix(c(1,0.2, 0.3, 1),2)
#'var2<-TVAR.sim(B=B2,nthresh=0,n=100, include="const", varcov=varcov)
#'ts.plot(var2, type="l", col=c(1,2))
#'
#'
#'##Simulation of a TVAR with 1 threshold
#'B<-rbind(c(0.11928245, 1.00880447, -0.009974585, -0.089316, 0.95425564, 0.02592617),
#'         c(0.25283578, 0.09182279,  0.914763741, -0.0530613, 0.02248586, 0.94309347))
#'colnames(B) <- paste(rep(c("Const", "Lag_1_var1", "Lag_1_var2"), 2), c("Low", "High"), sep="_")
#'sim<-TVAR.sim(B=B,nthresh=1,n=500, mTh=1, Thresh=5, starting=matrix(c(5.2, 5.5), nrow=1))
#'
#'#estimate the new serie
#'TVAR(sim, lag=1, dummyToBothRegimes=TRUE)
#'
#'
#'
TVAR.sim <- function(B, Thresh, nthresh=1, n=200, lag=1, include = c("const", "trend","none", "both"),  thDelay=1,  
thVar=NULL, mTh=1, starting=NULL, innov=rmnorm(n, varcov=varcov), varcov=diag(1,nrow(B)), show.parMat=FALSE, round=FALSE, seed, ...){

	    TVAR.gen(B=B, Thresh=Thresh, nthresh=nthresh, type="simul", n=n, lag=lag, include = include,  thDelay=thDelay,
		    thVar=thVar, mTh=mTh, starting=starting, innov=innov, varcov=varcov, show.parMat=show.parMat, round=round, seed=seed, ...)
}



#'Bootstrap a multivariate Threshold Autoregressive (TVAR) model
#'
#'Recursive bootstrap of a multivariate Threshold VAR (TVAR)
#'
#'The function bootstraps a given model. This is done on a object generated by
#'TVECM (or VECM). A simple residual recursive bootstrap is done.
#'
#'@aliases TVAR.boot
#'@param TVARobject Object of class \code{TVAR} generated by function
#'\code{\link{TVAR}}
#'@param innov Innovations used for bootstrap. If missing, residuals are
#'resampled.
#'@param seed Optional. Seed for the random resampling function.
#'@param \dots Further arguments passed to the underlying (un-exported)
#'\code{TVAR.gen} function
#'@return A matrix with the bootstraped series.
#'@author Matthieu Stigler
#'@seealso \code{\link{TVAR}} to estimate a TVAR, \code{\link{VAR.sim}} to
#'simulate/bootstrap a VAR, \code{\link{TVECM.sim}} to simulate/bootstrap a
#'TVECM.
#'@keywords ts bootstrap
#'@export
#'@examples
#'
#'
#'
#'##Bootstrap a TVAR with two threshold (three regimes)
#'data(zeroyld)
#'serie<-zeroyld
#'mod <- TVAR(data=serie,lag=1, nthresh=1)
#'TVAR.boot(mod)
#'
#'
TVAR.boot <- function(TVARobject, innov, seed, ...){
  if(!inherits(TVARobject, "TVAR")) stop("Please provide an object of class 'TVAR' generated by TVAR()")
  TVAR.gen(TVARobject=TVARobject,type= "boot", innov = innov, seed=seed, ...)
}

VAR.sim2 <- function(B,  n=200, lag=1, include = c("const", "trend","none", "both"),  starting=NULL, 
innov=rmnorm(n,  varcov=varcov), varcov=diag(1,nrow(B)), show.parMat=FALSE, round=FALSE, seed, ...){
  TVAR.gen(B=B, nthresh=0, type="simul", n=n, lag=lag, include = include,   starting=starting, innov=innov, varcov=varcov, show.parMat=show.parMat, round=round, seed=seed, ...)
}

VAR.boot2 <- function(VARobject, innov, seed, ...){
  if(!inherits(VARobject, "VAR")) stop("Please provide an object of class 'VAR' generated by lineVar()")
  if(attr(VARobject, "varsLevel")!="level") stop("Does not work with VAR in diff or ADf specification")

  TVAR.gen(TVARobject=VARobject,type= "boot", innov = innov, seed=seed, ...)
}


if(FALSE){
library(tsDyn)
environment(TVAR.gen)<-environment(star)
environment(TVAR.sim)<-environment(star)
environment(TVAR.boot)<-environment(star)

##Simulation of a TVAR with 1 threshold
B<-rbind(c(0.11928245, 1.00880447, -0.009974585, -0.089316, 0.95425564, 0.02592617),c(0.25283578, 0.09182279,  0.914763741, -0.0530613, 0.02248586, 0.94309347))
sim<-TVAR.sim(B=B,nthresh=1,n=500,mTh=1, Thresh=5, starting=matrix(c(5.2, 5.5), nrow=1))

#estimate the new serie
TVAR(sim, lag=1, dummyToBothRegimes=TRUE)


##Bootstrap a TVAR 
#data(zeroyld)
serie<-zeroyld

# TVAR.sim(data=serie,nthresh=0)
# all(TVAR.sim(data=serie,nthresh=0, type="check", lag=1)==serie)

##with two threshold (three regimes)
# TVAR.sim(data=serie,nthresh=2,type="boot",mTh=1, Thresh=c(7,9))

environment(TVAR.sim)<-environment(star)

##Check case 1 (B is given) ok!
ns<-nrow(serie)

comp_tvar_sim <- function(mod, serie){
  ns <- nrow(serie)
  sim_mod <- TVAR.sim(B=coef(mod), lag=mod$lag, include=mod$include,nthresh=0, n=ns-mod$lag, innov=residuals(mod), starting=serie[1:mod$lag,,drop=FALSE])
  all.equal(sim_mod, as.matrix(serie)[-c(1:mod$lag),], check.attributes=FALSE)
}


#data(barry)
var_l1_co <-lineVar(barry, lag=1, include="const")
var_l1_tr <-lineVar(barry, lag=1, include="trend")
var_l1_bo <-lineVar(barry, lag=1, include="both")
var_l1_no <-lineVar(barry, lag=1, include="none")

var_l3_co <-lineVar(barry, lag=3, include="const")
var_l3_tr <-lineVar(barry, lag=3, include="trend")
var_l3_bo <-lineVar(barry, lag=3, include="both")
var_l3_no <-lineVar(barry, lag=3, include="none")

var_all <- list(
		var_l1_co, var_l1_tr, var_l1_bo, var_l1_no, 
		var_l3_co, var_l3_tr, var_l3_bo, var_l3_no)

lapply(var_all , comp_tvar_sim, serie=barry)
comp_tvar_sim(var_l3_tr, serie=barry)




##Check case 2 (data is given) ok!
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="const", round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="const", lag=3, round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="trend", round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="trend", lag=3, round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="both", round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="both", lag=3, round=FALSE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=0, type="check", include="both", lag=2, round=FALSE),as.matrix(serie), check.attributes=FALSE)

all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=1, type="check",mTh=1, include="const", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=1, type="check",mTh=1, include="trend", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=1, type="check",mTh=2, include="const", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=1, type="check",mTh=2, include="trend", round=TRUE),as.matrix(serie), check.attributes=FALSE)

all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=2, type="check",mTh=1, include="const", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=2, type="check",mTh=1, include="trend", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=2, type="check",mTh=2, include="const", round=TRUE),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(data=serie,nthresh=2, type="check",mTh=2, include="trend", round=TRUE),as.matrix(serie), check.attributes=FALSE)



### Check case 3: TVARobject is given
all.equal(tsDyn:::TVAR.gen(TVARobject=lineVar(serie, lag=1),type="check"),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(TVARobject=lineVar(serie, lag=1, include="trend"),type="check"),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(TVARobject=lineVar(serie, lag=1, include="both"),type="check"),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(TVARobject=lineVar(serie, lag=1, include="none"),type="check"),as.matrix(serie), check.attributes=FALSE)


all.equal(tsDyn:::TVAR.gen(TVARobject=TVAR(serie, nthresh=1, lag=1, trace=FALSE),type="check"),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(TVARobject=TVAR(serie, nthresh=1, lag=2, trace=FALSE),type="check"),as.matrix(serie), check.attributes=FALSE)

all.equal(tsDyn:::TVAR.gen(TVARobject=TVAR(serie, nthresh=2, lag=1, trace=FALSE),type="check"),as.matrix(serie), check.attributes=FALSE)
all.equal(tsDyn:::TVAR.gen(TVARobject=TVAR(serie, nthresh=2, lag=2, trace=FALSE),type="check"),as.matrix(serie), check.attributes=FALSE)

a1 <- VAR.boot(VARobject=lineVar(serie, lag=1))
a2 <- VAR.boot(VARobject=lineVar(serie, lag=2))
a3 <- TVAR.boot(TVARobject=TVAR(serie, nthresh=1, lag=1, trace=FALSE))

}

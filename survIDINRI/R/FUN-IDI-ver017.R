## =================== ##
## "FUN-IDI-ver016.R"  ##
## =================== ##

# ver014 -- ver015 
# add parameter for label size --- cex.lab=NULL

# ver015 -- ver016 
# as.real() --> as.double()

# ver016 -- ver017 
# IDI.INF() --> provide p-value

## setwd("~/Dropbox/R/R-Extension/survIDINRI/distribution/ver014-survIDINRI-v1.0-1")

####################################
#--- R- wrapper -- (used in IDI.FUN)
####################################
# dyn.load("unoecdfv1.so")

unoecdf <- function(cc, pdiff, wt) {
	#----------------------
	#-- cc: vector (cutoff) (nc x 1)
	#-- pdiff: vector (difference of the risk score) (nx1)
	#-- wt: vector (weight or weigth*PTB) (nx1)
	#-- retrun -- empirical distribution fucntion
	#----------------------
	    out <- .Fortran("unoecdf",PACKAGE = "survIDINRI",
                n=as.integer(length(pdiff)),
                nc=as.integer(length(cc)),
                score=as.double(pdiff),
                cutoff=as.double(cc),
                weight=as.double(wt),
                OUTECDF=as.double(rep(0,length(cc))))
        return(out)
   }





####################################
## Function: kmcens: Keplan-Meier ## <-- survC1
####################################


###############################################
## KM for censoring
###############################################
Ghat.FUN <- function(tt, Ti, Di,type='fl')
  {
    tmpind <- rank(tt)
    summary(survfit(Surv(Ti,1-Di)~1, se.fit=F, type=type), sort(tt))$surv[tmpind]
  }


###############################################
## Calculating the weights for the IPW approach
###############################################
WGT.FUN <- function(data, t0)
  {
    Ti <- data[,1]; Di <- data[,2]
    Wi <- rep(0,length(Ti))
    Wi[Ti <= t0] <- Di[Ti<=t0]/Ghat.FUN(Ti[Ti<=t0],Ti,Di)
    Wi[Ti >  t0] <- 1/Ghat.FUN(t0,Ti,Di)
    Wi
  }

###############################################
## Vector to Matrix
###############################################
VTM<-function(vc, dm){
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }


###############################################
## cumsum2 gets cumsum of each column = apply(X, 2, cumsum)
###############################################
cumsum2 <- function(mydat)     #cumsum by row, col remains the same
  {
    if(is.null(dim(mydat))) return(cumsum(mydat))
    else{
      out <- matrix(cumsum(mydat), nrow=nrow(mydat))
      out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
      return(out)
    }
  }



##############################################
# SUM.I
##############################################
SUM.I <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
  {  
    if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
    if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
    pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)

    if(is.null(Vi)){return(pos)}else{
      Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
      out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
      out[pos!=0,] <- Vi[pos,]
      if(is.null(dim(Vi))) out <- c(out)
      return(out) ## n.y x p
    }
  } 



##############################################
# PI.k.FUN
##############################################
PI.k.FUN <- function(tt,ebzi,xi,zi,k0=0,vi=NULL)
  {
    out = ebzi; pz = ncol(zi); nn=length(out)
    if(!is.null(vi))
      { 
        vi=as.matrix(vi); pv=ncol(vi); 
        if(k0>0){zi=zi[,rep(1:pz,pv)]; out=out*vi[,rep(1:pv,rep(pz,pv))] # v1*z1, v1*z2, ...
                }else{ out=out*vi }
      }
    if(k0==1){out = out*zi}
    if(k0==2){out = out*zi[,rep(1:pz,pz)]*zi[,rep(1:pz,rep(pz,pz))]}
    as.matrix(SUM.I(tt,"<=",xi,Vi=out)/nn)
  }



##############################################
# PTB.Gamma.FUN
##############################################
PTB.Gamma.FUN <- function(Vi,data,t0,gammahat)
  {
    logLamhat = gammahat[1]; betahat = gammahat[-1]
    xi = data[,1]; di = data[,2]; zi = data[,-(1:2),drop=F]; 
    ebzi = c(exp(zi%*%betahat)); pz = ncol(zi)
    tmpind = di==1; tj = xi[tmpind]; 
##  Vj = Vi[tmpind,]; 
    Vj = Vi[tmpind]
    nn = length(xi); pv = ncol(Vi)
    pi0.tj   = c(PI.k.FUN(tj,ebzi,xi,zi,k0=0))
    pi1.tj   = PI.k.FUN(tj,ebzi,xi,zi,k0=1)
    Ahat = PI.k.FUN(tj,ebzi,xi,zi,k0=2)*pi0.tj - pi1.tj[,rep(1:pz,pz)]*pi1.tj[,rep(1:pz,rep(pz,pz))]
    Ahat = matrix(apply(Ahat/pi0.tj^2,2,sum),ncol=pz)

    term1 = t(Vj)%*%zi[tmpind,]
    term2 = t(Vj)%*%(pi1.tj/pi0.tj)
    
    pi1.tj.V = PI.k.FUN(tj,ebzi,xi,zi,k0=1,vi=Vi);
    pi0.tj.V = PI.k.FUN(tj,ebzi,xi,zi,k0=0,vi=Vi); 
    term3 = t(rep(1,length(tj)))%*%((pi1.tj.V*pi0.tj-pi0.tj.V[,rep(1:pv,rep(pz,pv))]*pi1.tj[,rep(1:pz,pv)])/pi0.tj^2)
    term3 = t(matrix(term3,nrow=pz))
    betastar = betahat + solve(Ahat)%*%t(term1 - term2 - term3)
    
    tmpind2 = tj <= t0
##  term1 = c(t(Vj[tmpind2,])%*%(1/pi0.tj[tmpind2]))
    term1 = c(t(Vj[tmpind2 ])%*%(1/pi0.tj[tmpind2]))
### term2 = apply(pi0.tj.V[tmpind2,]/pi0.tj[tmpind2]^2,2,sum)
    term2 = sum(pi0.tj.V[tmpind2]/pi0.tj[tmpind2]^2)
    term3 = c(t(betastar-betahat) %*% apply(pi1.tj[tmpind2,,drop=F]/pi0.tj[tmpind2]^2,2,sum))
    logLamstar = logLamhat + (term1-term2-term3)/nn/exp(logLamhat)
##  rbind(logLamstar,betastar)
    c(logLamstar,betastar)
  }




## ================================================= ##
## === RISK SCORE (in prob scale )BY COX GIVEN T0 == ##
## ================================================= ##
## INPUT PARAMETER                                   ##
## ------------------------------------------------- ##
## data --- Data Structure should be:                ##
##           1st       column : min(T,C)             ##
##           2nd       column : I(T <= C)            ## 
##           remaining columns: covariates Z         ##
## t0   --- time point of interest                   ##
## Gi   --- Weigths                                  ##
## ================================================= ##
  Est.HR.FUN <- function(data, t0, Gi=NULL)
  {
    xi = data[,1]; di <- data[,2]; zi <- data[,-(1:2),drop=F]; nn <- length(xi)
    if(is.null(Gi)){
        betahat <- coxph(Surv(xi,di)~zi);   
        LamCox.t  <- basehaz(betahat,centered=FALSE); 
        LamCox.t0 <- approx(LamCox.t$time,LamCox.t$hazard,t0)$y
        #==== ADD from ver008 ====#
        if(is.na(LamCox.t0)){LamCox.t0=LamCox.t$hazard[nrow(LamCox.t)]}
    }else{
        betahat <- coxph(Surv(xi[Gi!=0],di[Gi!=0])~zi[Gi!=0,],weights=Gi[Gi!=0])
        LamCox.t  <- basehaz(betahat,centered=FALSE); 
        #==== ADD from ver008 ====#
        LamCox.t0 <- approx(LamCox.t$time,LamCox.t$hazard,t0)$y
        if(is.na(LamCox.t0)){LamCox.t0=LamCox.t$hazard[nrow(LamCox.t)]}
    }
        
    bzi <- zi%*%betahat$coef
    rs  <- 1-exp(-LamCox.t0*exp(bzi))
    
    #CoxRisk.u = exp(jump.u)*phat.km$LamCox.t0; 
    
    return(list("risk"=rs, "bz"=bzi, "beta"=betahat$coef, "LamCox.t0"=LamCox.t0,"LamCox.t"=LamCox.t, "ft"=betahat))

  }


## ================================================== ##
## === PERTURBATION OF RISK SCORE ESTIAMTES ========= ##
## ================================================== ##
## INPUT PARAMETER                                    ##
## -------------------------------------------------- ##
## data ------- Data Structure should be:             ##
##              1st       column : min(T,C)           ##
##              2nd       column : I(T <= C)          ## 
##              remaining columns: covariates Z       ##
## Gi   ------- Weigths (used for perturbation/ boot) ##
## estim ------ Fit reulst from Est.HR.FUN            ##
## ================================================== ##
  PTB.Est.HR.FUN <- function(data, t0, Gi, estim)
  {
     zi = data[,-(1:2),drop=F]
 
     Vi=as.matrix(Gi-1)

     gammahat<-c(log(estim$LamCox.t0), estim$beta)
     
     ptb<-PTB.Gamma.FUN(Vi, data, t0, gammahat)

     betastar <- ptb[-1]
     LamCox.t0.star<-exp(ptb[1])
        
     bzi <- zi%*%betastar
     rs  <- 1-exp(-LamCox.t0.star*exp(bzi))
        
     return(list("risk"=rs, ptb.gam=ptb))

  }





################################################################
## IDI.FUN (point estimate)
################################################################
##14##  IDI.FUN <- function(data, t0, pnew, pold, Wi=NULL, PTB.Vi=NULL, point=NULL)
  IDI.FUN <- function(data, t0, pnew, pold, Wi=NULL, PTB.Vi=NULL)
  {
  	
    xi = data[,1]; di <- data[,2]; nn <- length(xi)

    #===== weight ====#
    if(is.null(Wi)){
        Wi <- WGT.FUN(cbind(xi,di),t0); 
       }
       
    #===== perturb ====# when pertub !is.null weigth, point should not be null
    if(is.null(PTB.Vi)){
       PTB.Vi<-rep(1,nn)
       }

    #==== case & control at t0 =====#   
       case<- xi<t0  & di==1 
       cont<- xi>=t0


    #===== F(x) and G(x) ======#
    pdiff=pnew-pold ; 

       #--------------------
       # proc.time()-ptm ; ptm <- proc.time()
       #--------------------
       cc=sort(unique(c(-1, 0, 1, seq(-1,1,length=2000)))) ;

       ### jj=length(cc) ; FX<-rep(0,jj) ; GX<-rep(0,jj)
       ### for (i in 1:jj){
       ###     FX[i]<-as.numeric(pdiff[case] > cc[i]) %*% (Wi[case]*PTB.Vi[case]) /sum(Wi[case]*PTB.Vi[case]) ;
       ###     GX[i]<-as.numeric(pdiff[cont] > cc[i]) %*% (         PTB.Vi[cont]) /sum(         PTB.Vi[cont]) ; 
       ###     }            

       #--------------------
       # proc.time()-ptm ; ptm <- proc.time()
       #--------------------

       FX=unoecdf(cc, pdiff[case], Wi[case]*PTB.Vi[case])$OUTECDF 
       GX=unoecdf(cc, pdiff[cont],          PTB.Vi[cont])$OUTECDF 

       #--------------------
       # proc.time()-ptm ; ptm <- proc.time()
       #--------------------

       #--- IDI from F & G area ---# 
       IDI_D1 = diff(cc)%*%(FX[-1]+FX[-length(cc)])/2 - 1
       IDI_D0 = diff(cc)%*%(GX[-1]+GX[-length(cc)])/2 - 1
       IDI=IDI_D1 - IDI_D0 

       #--- Continuous NRI from F & G ---# 
       TXI_D1=FX[cc==0]
       TXI_D0=GX[cc==0]
       TXI=TXI_D1 -TXI_D0

       #--- Median ---# 
       med.fx.idx1<-max(which(FX == min(FX[FX>=0.5]))) ; 
       med.fx.idx2<-min(which(FX == max(FX[FX<=0.5]))) ; MED_D1<-mean(c(cc[med.fx.idx1], cc[med.fx.idx2]))
       med.gx.idx1<-max(which(GX == min(GX[GX>=0.5]))) ; 
       med.gx.idx2<-min(which(GX == max(GX[GX<=0.5]))) ; MED_D0<-mean(c(cc[med.gx.idx1], cc[med.gx.idx2]))
       MED=MED_D1 - MED_D0


       #--- MRD ---# 
##14##       MRD_new_D1 = pnew[case] %*% (Wi[case]*PTB.Vi[case]) /sum(Wi[case]*PTB.Vi[case])  
##14##       MRD_new_D0 = pnew[cont] %*% (         PTB.Vi[cont]) /sum(         PTB.Vi[cont])
##14##       MRD_old_D1 = pold[case] %*% (Wi[case]*PTB.Vi[case]) /sum(Wi[case]*PTB.Vi[case])  
##14##       MRD_old_D0 = pold[cont] %*% (         PTB.Vi[cont]) /sum(         PTB.Vi[cont])
##14##       MRD_new=MRD_new_D1 - MRD_new_D0
##14##       MRD_old=MRD_old_D1 - MRD_old_D0

##14##      return(list(IDI=IDI, IDI_D1=IDI_D1, IDI_D0=IDI_D0, TXI=TXI, TXI_D1=TXI_D1, TXI_D0=TXI_D0, FX=FX, GX=GX, cc=cc, MED_D1=MED_D1, MED_D0=MED_D0, MED=MED, MRD_old=MRD_old,MRD_new=MRD_new))

      return(list(IDI=IDI, IDI_D1=IDI_D1, IDI_D0=IDI_D0, TXI=TXI, TXI_D1=TXI_D1, TXI_D0=TXI_D0, FX=FX, GX=GX, cc=cc, MED_D1=MED_D1, MED_D0=MED_D0, MED=MED))

  }


################################################################
## PERTURBATION OF M1, M2, M3, F(s) and G(u)
################################################################
## I: data      (nx2)     failure time, status(1=fail, 0=censor)
##    covs0     (p2x1)    covariate of old model
##    covs1     (p1x1)    covariate of new model
##    t0        (1x1)     fixing time point of interest [year]
##    estims0   <-------- object coming form estim2
##    estims1   <-------- object coming form estim2
##    kmc       <-------- object coming form kmcens
##    est.hat   <-------- IDI.FUN object (point estimates)
##    seed      (1x1)     seed for pertubation
##    npert     (1x1)     itration for perturbation
##==============================================================   
## O: PTB.IDI      (Jx1)  petrubation result for M1 (IDI)
##    PTB.TXI      (Jx1)  petrubation result for M2 
##    PTB.MED      (Jx1)  petrubation result for M3 (Median difference)
################################################################
PTB.IDI.FUN <- function(data, covs0, covs1, t0, estim0, estim1, kmc, est.hat, seed1=NULL, npert=300, npert.rand=NULL)
  {

   #-- random numbers --#
   if(!is.null(seed1)){set.seed(seed1)}

   if(!is.null(npert.rand)){
   	  npert=ncol(npert.rand)
   	  ptb.rand=npert.rand
   	}else{
      ptb.rand=c()
       for (i in 1:npert){
	  ptb.rand=cbind(ptb.rand,rexp(nrow(data)))
	  }
	}

   #-- initialize --#
   time<-data[,1]
   status<-data[,2]
   
   distinct<-unique(sort(time))
   t<-length(distinct)
   n<-length(time)
   z0<-cbind(covs0) ; p0<-ncol(as.matrix(z0))
   z1<-cbind(covs1) ; p1<-ncol(as.matrix(z1))

   PTB.IDI<-rep(0, npert)
   PTB.TXI<-rep(0, npert)
   PTB.MED<-rep(0, npert)

   PTB.IDI.D1<-rep(0, npert)
   PTB.TXI.D1<-rep(0, npert)
   PTB.MED.D1<-rep(0, npert)

   PTB.IDI.D0<-rep(0, npert)
   PTB.TXI.D0<-rep(0, npert)
   PTB.MED.D0<-rep(0, npert)

##14##   PTB.MRD_new<-rep(0, npert)
##14##   PTB.MRD_old<-rep(0, npert)


##14##   PTB.IDI.CHK<-rep(0, npert)
##14##   PTB.TXI.CHK<-rep(0, npert)
##14##   PTB.MED.CHK<-rep(0, npert)



   #-- censoring --#
   cens.surv<-kmc$surv
   cens.psii<-kmc$psii

   #------------------#
   #-- perturbation --#
   #------------------#

##14##   PTB.FX=c() ; PTB.GX=c()

   j=1
   while (j <= npert){ 

#--------------------
ptm <- proc.time()
#--------------------

   #-- ramdom nubmber --#
     rn<-ptb.rand[,j]

#--------------------
proc.time()-ptm ; ptm <- proc.time()
#--------------------

     
   #-- G star --#
   if(sum(time < t0 & status==0)==0){
    weight=rep(1,n)
   	}else{
     pghat<- cens.surv - (t(cens.psii) %*% (rn-1)/n)

     #== censoring & weight ==#
      ghat<-pghat
      #--ghat(t0) -=# 
      ghat.t0<-0
      for (i in 1:t){if(distinct[i]<=t0){ghat.t0<-ghat[i]}}
      #-- append (ghat.X) --# 
      ghat.x=ghat[match(time, distinct)]
       #ghat.x<-rep(0,n)
       #for (i in 1:n){
       #  ghat.x[i]<-ghat[distinct==time[i]]
       #  if (is.na(ghat.x[i])){ghat.x[i]<-0}
       #}

     #-- append (pi(t), vi(t), weight(t)) --# --- OK
     pit=rep(0, n)
     pit[status==1]=ghat.x[status==1]
     pit[time>=t0]=ghat.t0
      # for (i in 1:n){
      #    if(time[i]>=t0)                {pit[i]<-ghat.t0}
      #    if(time[i]<=t0 & status[i]==1) {pit[i]<-ghat.x[i]}
      # }
    
    vit<-as.numeric(time<=t0)*status + as.numeric(time>t0)
    weight<-rep(0,n)
    weight[vit==1]=1/pit[vit==1]
      #  for (i in 1:n){
      #     if (vit[i]==0) {weight[i]<-0}
      #     if (vit[i]!=0) {weight[i]<-vit[i]/pit[i]}
      #  }
    }

#--------------------
proc.time()-ptm ; ptm <- proc.time()
#--------------------

     j=j+1
     pnew<-PTB.Est.HR.FUN(as.matrix(cbind(time, status, covs1)), t0, Gi=rn, estim1)$risk
     pold<-PTB.Est.HR.FUN(as.matrix(cbind(time, status, covs0)), t0, Gi=rn, estim0)$risk

#--------------------
proc.time()-ptm ; ptm <- proc.time()
#--------------------
    
   #--------------------------# 
###14###    PTB<-IDI.FUN(cbind(time, status), t0, pnew, pold, Wi=weight, PTB.Vi=rn, point=est.hat)
    PTB<-IDI.FUN(cbind(time, status), t0, pnew, pold, Wi=weight, PTB.Vi=rn)
    
#--------------------
proc.time()-ptm ; ptm <- proc.time()
#--------------------

    #-- IDI, NQI, MED --#
    PTB.IDI[j]<-PTB$IDI
    PTB.TXI[j]<-PTB$TXI
    PTB.MED[j]<-PTB$MED

    PTB.IDI.D1[j]<-PTB$IDI_D1
    PTB.TXI.D1[j]<-PTB$TXI_D1
    PTB.MED.D1[j]<-PTB$MED_D1

    PTB.IDI.D0[j]<-PTB$IDI_D0
    PTB.TXI.D0[j]<-PTB$TXI_D0
    PTB.MED.D0[j]<-PTB$MED_D0

##14##    PTB.MRD_new[j]<-PTB$MRD_new
##14##    PTB.MRD_old[j]<-PTB$MRD_old

#--------------------
proc.time()-ptm
#--------------------

    
    #--- F and G ---#
###14###    PTB.FX<-rbind(PTB.FX, PTB$FX)
###14###    PTB.GX<-rbind(PTB.GX, PTB$GX)
  }
   
###14### return(list(PTB.IDI=PTB.IDI, PTB.TXI=PTB.TXI, PTB.MED=PTB.MED, PTB.FX=PTB.FX, PTB.GX=PTB.GX,PTB.IDI.D0=PTB.IDI.D0, PTB.TXI.D0=PTB.TXI.D0, PTB.MED.D0=PTB.MED.D0, PTB.IDI.D1=PTB.IDI.D1, PTB.TXI.D1=PTB.TXI.D1, PTB.MED.D1=PTB.MED.D1, PTB.MRD_new=PTB.MRD_new, PTB.MRD_old=PTB.MRD_old))

return(list(PTB.IDI=PTB.IDI, PTB.TXI=PTB.TXI, PTB.MED=PTB.MED,PTB.IDI.D0=PTB.IDI.D0, PTB.TXI.D0=PTB.TXI.D0, PTB.MED.D0=PTB.MED.D0, PTB.IDI.D1=PTB.IDI.D1, PTB.TXI.D1=PTB.TXI.D1, PTB.MED.D1=PTB.MED.D1))

}


##=============================================================================================



####################################################################
## IDI.INF (input data and get point est and se and (1-alpha)CI)  
####################################################################
## I: indata     (nx2)     failure time, status(1=fail, 0=censor)
##    covs0      (nxp0)    covariate of old model
##    covs1      (nxp1)    covariate of new model
##    t0         (1x1)     fixing time point of interest [year]
##    seed1      (1x1)     seed for pertubation
##    npert      (1x1)     number of perturbation to get se (default=300)
##    npert.rand (nxM)     random numbers for perturbation if already generated (default=NULL)
##    alpha      (1x1)     gives (1-alpha) confidence interval (default=0.05)
#####################################################################
IDI.INF <- function(indata, covs0, covs1, t0, npert=300, npert.rand=NULL, seed1=NULL, alpha=0.05)
  {

  #-- data ---
  indata0<-cbind(indata, covs0)
  indata1<-cbind(indata, covs1)
    	
  #--- fit Cox and get score ---#
  ft1<-Est.HR.FUN(as.matrix(indata1),t0)
  ft0<-Est.HR.FUN(as.matrix(indata0),t0)

  #--- point estimate ---#
###14###  pest<-IDI.FUN(indata, t0, ft1$risk, ft0$risk, Wi=NULL, PTB.Vi=NULL, point=NULL)
  pest<-IDI.FUN(indata, t0, ft1$risk, ft0$risk, Wi=NULL, PTB.Vi=NULL)


  #--- confidence interval ----#
  kmc<-kmcens(indata[,1], indata[,2])



  if(!is.null(npert) | !is.null(npert.rand)){
  #============================#
  #--- perturbation -----------#
  #============================#

#--------------------
ptm <- proc.time()
#--------------------

  if(!is.null(npert.rand)){
    ptb<-PTB.IDI.FUN(indata[,1:2], covs0, covs1, t0, ft0, ft1, kmc, pest, seed1=seed1, npert.rand=npert.rand)
  	}else{
    ptb<-PTB.IDI.FUN(indata[,1:2], covs0, covs1, t0, ft0, ft1, kmc, pest, seed1=seed1, npert=npert)
  		}
  }

#--------------------
proc.time()-ptm
#--------------------

  #--- point ests ---
  m1.est=c(pest$IDI, pest$IDI_D1, pest$IDI_D0)
  m2.est=c(pest$TXI, pest$TXI_D1, pest$TXI_D0)
  m3.est=c(pest$MED, pest$MED_D1, pest$MED_D0)
##14##  mrd.est=c(pest$MRD_new, pest$MRD_old)

  #=================
  # OUTPUT
  #=================
  if(!is.null(npert) | !is.null(npert.rand)){

   #--- percentile-based p-value (add ver017) ---
    m1.p=if(pest$IDI>0){mean(ptb$PTB.IDI < 0) * 2
    	          }else{mean(ptb$PTB.IDI > 0) * 2}
    m2.p=if(pest$TXI>0){mean(ptb$PTB.TXI < 0) * 2
    	          }else{mean(ptb$PTB.TXI > 0) * 2}
    m3.p=if(pest$MED>0){mean(ptb$PTB.MED < 0) * 2
    	          }else{mean(ptb$PTB.MED > 0) * 2}
   #--- CI ---
    m1=c(pest$IDI, quantile(ptb$PTB.IDI, prob=c(alpha/2, 1-alpha/2)), m1.p) 
    m2=c(pest$TXI, quantile(ptb$PTB.TXI, prob=c(alpha/2, 1-alpha/2)), m2.p) 
    m3=c(pest$MED, quantile(ptb$PTB.MED, prob=c(alpha/2, 1-alpha/2)), m3.p) 

   

    z=list(m1=m1, m2=m2, m3=m3, m1.est=m1.est, m2.est=m2.est, m3.est=m3.est, point=pest)
    }else{
    z=list(m1.est=m1.est, m2.est=m2.est, m3.est=m3.est, point=pest)
  	}
#  class(z)="IDI"
   z
}


## ====================================== ##
## === PRINT RESULTS ==================== ##
## ====================================== ##
   IDI.INF.OUT<-function(x){
   	if(!is.null(x$m1)){
   	tmp=rbind(x$m1,x$m2,x$m3)
    colnames(tmp)=c("Est.","Lower", "Upper","p-value")
    rownames(tmp)=c("M1","M2","M3")
    print(round(tmp, digits=3))
   		}else{
   	tmp=as.matrix(rbind(x$m1.est[1],x$m2.est[1],x$m3.est[1]))
    colnames(tmp)=c("Est.")
    rownames(tmp)=c("M1","M2","M3")
    print(round(tmp, digits=3))
   			}
   	}


## ====================================== ##
## === GRAPH ============================ ##
## ====================================== ##
 shade <- function(x, yl, yu, ...){ 
  	nl <- length(x)
    for(i in 1:(nl-1))
      {polygon(c(x[i], x[i], x[i+1], x[i+1]), c(yl[i], yu[i], yu[i+1], yl[i+1]), ...)}
  }

#--- ver015 -- add arguments ---
 IDI.INF.GRAPH<-function(x, main=NULL, xlab = NULL, ylab = NULL, cex.main=NULL, cex.lab=NULL,...){

   cc.diff=x$point$cc
   diffs1=1-x$point$FX
   diffs0=1-x$point$GX
   
   Xlim<-c(min(cc.diff[diffs1>0.001 | diffs0>0.001]),max(cc.diff[diffs1<0.999 | diffs0<0.999])) ; Ylim<-c(0, 1)
   plot(cc.diff, diffs1, type="l", lty=1, lwd=3, xlim=Xlim, ylim=Ylim, xlab="", ylab="",...)
    diffs0_a=diffs0; diffs0_a[diffs0>diffs1]=diffs1[diffs0>diffs1]
    diffs1_a=diffs1; diffs1_a[diffs1>diffs0]=diffs0[diffs1>diffs0]
    shade(cc.diff, diffs1, diffs0_a,col="blue",border=NA) 
    shade(cc.diff, diffs0, diffs1_a,col="red",border=NA) 

   lines(cc.diff, diffs1, lty=1, lwd=3) 
   lines(cc.diff, diffs0, lty=1, lwd=1) 
   abline(v=0.0, col="gray", lty=1)
   abline(h=0.5, col="gray", lty=1)

   c0.idx<-which(abs(cc.diff)==min(abs(cc.diff)))
   points(0,diffs1[c0.idx],pch=19, col="black", cex=1.5) ;
   points(0,diffs0[c0.idx],pch=19, col="black", cex=1.5) ; 
   m1.idx<-which(abs(diffs1-0.5)==min(abs(diffs1-0.5)))
   m0.idx<-which(abs(diffs0-0.5)==min(abs(diffs0-0.5)))
   points(cc.diff[m1.idx[1]], 0.5, col="darkgray", pch=19, cex=1.5) ;
   points(cc.diff[m0.idx[1]], 0.5, col="darkgray", pch=19, cex=1.5) ;
   if(is.null(xlab)){xlab="s"}
   if(is.null(ylab)){ylab=expression(paste("pr(",hat(D)<=s,")"))}
   title(main=main, xlab=xlab, ylab=ylab,cex.lab=cex.lab, cex.main=cex.main)
#  title(xlab="s", ylab=expression(paste("pr(",hat(D)<=s,")")), cex.lab=cex.lab)

}





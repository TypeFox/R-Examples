# DLM_output methods
TACfilter<-function(TAC) {
  TAC[TAC<0]<-NA    # Have to robustify due to R optmization problems.. work in progress.
  TAC[TAC>(mean(TAC,na.rm=T)+5*sd(TAC,na.rm=T))]<-NA  # remove very large TAC samples
  return(TAC)
}

prodPTF<-function(depletion,n,MSY){   # Pella-Tomlinson production function required for DB-SRA
   y<-(n^(n/(n-1)))/(n-1)
   MSY*y*depletion-MSY*y*depletion^n
}

fn<-function(n,BMSY_K){               # optimizer to find parameter n according to sampled BMSY/B0 (theta)
   thetapred<-n^(-1/(n-1))
   (BMSY_K-thetapred)^2
}

getn<-function(BMSY_K){               # wrapper for n finder
   optimize(fn,c(0.01,6),BMSY_K=BMSY_K)$minimum #get the optimum
}

gety<-function(n)  (n^(n/(n-1)))/(n-1) # More DBSRA code: get the y parameter for n

FMSYref<-function(x,DLM_data,reps=100)trlnorm(reps,DLM_data@OM$A[x]*(1-exp(-DLM_data@OM$FMSY[x])),0.01)
class(FMSYref)<-"DLM_output"

FMSYref50<-function(x,DLM_data,reps=100)trlnorm(reps,DLM_data@OM$A[x]*(1-exp(-DLM_data@OM$FMSY[x]))*0.5,0.01)
class(FMSYref50)<-"DLM_output"

FMSYref75<-function(x,DLM_data,reps=100)trlnorm(reps,DLM_data@OM$A[x]*(1-exp(-DLM_data@OM$FMSY[x]))*0.75,0.01)
class(FMSYref75)<-"DLM_output"

DynF<-function(x,DLM_data,yrsmth=10,gg=2,reps=100){
    
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@Ind, DLM_data@Abun, DLM_data@Mort, DLM_data@FMSY_M"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  
  C_dat<-log(DLM_data@Cat[x,ind])
  B_dat<-log(DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x])
  C_hist<-exp(predict(loess(C_dat~ind,degree=1)))
  B_hist<-exp(predict(loess(B_dat~ind,degree=1)))
  
  ind<-2:yrsmth
  ind1<-1:(yrsmth-1)
  SP_hist<-B_hist[ind]-B_hist[ind1]+C_hist[ind1]
  
  Frat<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
  Flim <- matrix(NA, nrow=2, ncol=reps)
  Flim[1,] <- Frat * 0.5 
  Flim[2,] <- Frat * 2

  yind<-1:length(SP_hist)
  SP_mu<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1))
  SP_se<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1),se=T)$se.fit
  SP_new<-rnorm(reps,SP_mu,SP_se/2)
  Glm<-summary(lm(SP_hist~B_hist[ind1]))$coefficients[2,1:2] # plot(B_hist[ind1],SP_hist) # points(B_hist[ind1],SP_hist,col='green')
  G_new<-rnorm(reps,Glm[1],Glm[2]/2)
  #G_new[G_new>2*Frat]<-2*Frat[G_new<(2*Frat)]
  #G_new[G_new<(-2*Frat)]<--2*Frat[G_new<(-2*Frat)]
  G_new[G_new>0]<-G_new[G_new>0]*3
  newF<-Frat*exp(-G_new*gg)
  newF[newF<Flim[1]]<-Flim[1]
  newF[newF>Flim[2]]<-Flim[2]
  
  TAC<-newF*B_hist[yrsmth]
  TACfilter(TAC)
  
}
class(DynF)<-"DLM_output"

Fadapt<-function(x,DLM_data,reps=100,yrsmth=7,gg=1){
  
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@Ind, DLM_data@Abun, DLM_data@Mort, DLM_data@FMSY_M"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  
  C_dat<-log(DLM_data@Cat[x,ind])
  B_dat<-log(DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x])
  C_hist<-exp(predict(loess(C_dat~ind,degree=1)))
  B_hist<-exp(predict(loess(B_dat~ind,degree=1)))
  
  ind<-2:yrsmth
  ind1<-1:(yrsmth-1)
  SP_hist<-B_hist[ind]-B_hist[ind1]+C_hist[ind1]
  
  Frat<-DLM_data@Mort[x]*DLM_data@FMSY_M[x]
  Flim<-Frat*c(0.5,2)
  Flimr<-Flim[2]-Flim[1]
  
  yind<-1:length(SP_hist)
  SP_mu<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1))
  SP_se<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1),se=T)$se.fit
  SP_new<-rnorm(reps,SP_mu,SP_se/2)
  Glm<-summary(lm(SP_hist~B_hist[ind1]))$coefficients[2,1:2] # plot(B_hist[ind1],SP_hist) # points(B_hist[ind1],SP_hist,col='green')
  G_new<-rnorm(reps,Glm[1],Glm[2])
  
  Fold<-mean(C_hist/B_hist)
  
  if(Fold<Flim[1])Fmod1<-(-2)
  if(Fold>Flim[2])Fmod1<-2
  if(Fold>Flim[1]&Fold<Flim[2]){
    Ffrac<-(Fold-Flim[1])/Flimr
    Fmod1<-log(Ffrac/(1-Ffrac))
  }
  Fmod2<-Fmod1+gg*-G_new
  newF<-Flim[1]+(exp(Fmod2)/(1+exp(Fmod2)))*Flimr
  TAC<-newF*B_hist[yrsmth]
  TACfilter(TAC)
}
class(Fadapt)<-"DLM_output"

DepF<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Year, DLM_data@Dep, DLM_data@Mort, DLM_data@FMSY_M, DLM_data@BMSY_B0"
  Frat<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
  if (is.na(DLM_data@Dep[x]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  Bt_K<-rbeta(reps*100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K<-Bt_K[Bt_K>=0.01&Bt_K<=0.99][1:reps] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  adj<-Bt_K*(1-Bt_K)*4
  adj[Bt_K>0.5]<-1
  TAC<-Frat*DLM_data@Abun[x]*adj
  TACfilter(TAC)
}
class(DepF)<-"DLM_output"

Gcontrol<-function(x,DLM_data,reps=100,yrsmth=10,gg=2,glim=c(0.5,2)){
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@Ind, DLM_data@Abun"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  C_dat<-log(DLM_data@Cat[x,ind])
  B_dat<-log(DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x])
  C_hist<-exp(predict(loess(C_dat~ind,degree=1)))
  B_hist<-exp(predict(loess(B_dat~ind,degree=1)))
  ind<-2:yrsmth
  ind1<-1:(yrsmth-1)
  SP_hist<-B_hist[ind]-B_hist[ind1]+C_hist[ind1]
  yind<-1:length(SP_hist)
  SP_mu<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1))
  SP_se<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1),se=T)$se.fit
  SP_new<-rnorm(reps,SP_mu,SP_se/2)
  Glm<-summary(lm(SP_hist~B_hist[ind1]))$coefficients[2,1:2]
  G_new<-rnorm(reps,Glm[1],Glm[2]/2)

  TAC<-SP_new*(1-gg*G_new)
  TAC[TAC<glim[1]*C_hist[yrsmth]]<-glim[1]*C_hist[yrsmth]
  TAC[TAC>glim[2]*C_hist[yrsmth]]<-glim[2]*C_hist[yrsmth]
  
  #Carr<-cbind(array(rep(DLM_data@Cat[x,],each=reps),c(reps,length(DLM_data@Cat[x,]))),TAC)
  #Warr<-(DLM_data@Mort[x]*exp(-DLM_data@Mort[x]*(1:ncol(Carr))))[ncol(Carr):1]
  #Warr<-Warr/sum(Warr)
  #TAC<-apply(t(matrix(Warr,nrow=ncol(Carr),ncol=reps))*Carr,1,sum)
  TACfilter(TAC)
}
class(Gcontrol)<-"DLM_output"

Rcontrol<-function(x,DLM_data,reps=100,yrsmth=10,gg=2,glim=c(0.5,2)){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0DLM_data@steep, DLM_data@CV_steep, DLM_data@MaxAge, DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Cat, DLM_data@Ind"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kvec<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfvec<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0vec<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  rsamp<-getr(x,DLM_data,Mvec,Kvec,Linfvec,t0vec,hvec,maxage=DLM_data@MaxAge,r_reps=reps)

  depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  Bt_K<-rbeta(100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K<-Bt_K[Bt_K>0.01&Bt_K<0.99][1] # interval censor (0.01,0.99)  as in Dick and MacCall 2011

  G_new<-rsamp*(1-2*Bt_K)          # here is a big difference from SPHCR

  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  C_dat<-log(DLM_data@Cat[x,ind])
  B_dat<-log(DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x])
  C_hist<-exp(predict(loess(C_dat~ind,degree=1)))
  B_hist<-exp(predict(loess(B_dat~ind,degree=1)))
  ind<-2:yrsmth
  ind1<-1:(yrsmth-1)
  SP_hist<-B_hist[ind]-B_hist[ind1]+C_hist[ind1]
  yind<-1:length(SP_hist)
  SP_mu<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1))
  SP_se<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1),se=T)$se.fit
  SP_new<-rnorm(reps,SP_mu,SP_se/2)

  TAC<-SP_new*(1-gg*G_new)
  TAC[TAC<glim[1]*C_hist[yrsmth]]<-glim[1]*C_hist[yrsmth]
  TAC[TAC>glim[2]*C_hist[yrsmth]]<-glim[2]*C_hist[yrsmth]
  
  #Carr<-cbind(array(rep(DLM_data@Cat[x,],each=reps),c(reps,length(DLM_data@Cat[x,]))),TAC)
  #Warr<-(DLM_data@Mort[x]*exp(-DLM_data@Mort[x]*(1:ncol(Carr))))[ncol(Carr):1]
  #Warr<-Warr/sum(Warr)
  #TAC<-apply(t(matrix(Warr,nrow=ncol(Carr),ncol=reps))*Carr,1,sum)
  TACfilter(TAC)

}
class(Rcontrol)<-"DLM_output"

Rcontrol2<-function(x,DLM_data,reps=100,yrsmth=10,gg=2,glim=c(0.5,2)){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@steep, DLM_data@CV_steep, DLM_data@MaxAge, DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Cat, DLM_data@Ind"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kvec<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfvec<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0vec<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  rsamp<-getr(x,DLM_data,Mvec,Kvec,Linfvec,t0vec,hvec,maxage=DLM_data@MaxAge,r_reps=reps)

  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  C_dat<-log(DLM_data@Cat[x,ind])
  B_dat<-log(DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x])
  C_hist<-exp(predict(loess(C_dat~ind,degree=1)))
  B_hist<-exp(predict(loess(B_dat~ind,degree=1)))
  ind<-2:yrsmth
  ind1<-1:(yrsmth-1)
  SP_hist<-B_hist[ind]-B_hist[ind1]+C_hist[ind1]
  yind<-1:length(SP_hist)
  SP_mu<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1))
  SP_se<-predict(lm(SP_hist~yind),newdat=list(yind=length(SP_hist)+1),se=T)$se.fit
  SP_new<-rnorm(reps,SP_mu,SP_se/2)
  SParr<-array(rep(SP_hist,each=reps),dim=c(reps,yrsmth-1))
  Barr<-array(rep(B_hist[ind],each=reps),dim=c(reps,yrsmth-1))
  rarr<-array(rep(rsamp,yrsmth-1),dim=c(reps,yrsmth-1))
  b2<-apply(SParr/Barr-rarr,1,sum)*apply(Barr,1,sum)/apply(Barr^2,1,sum)
  G_new<-rsamp-2*b2*B_hist[yrsmth]

  TAC<-SP_new*(1-gg*G_new)
  TAC[TAC<glim[1]*C_hist[yrsmth]]<-glim[1]*C_hist[yrsmth]
  TAC[TAC>glim[2]*C_hist[yrsmth]]<-glim[2]*C_hist[yrsmth]
  #Carr<-cbind(array(rep(DLM_data@Cat[x,],each=reps),c(reps,length(DLM_data@Cat[x,]))),TAC)
  #Warr<-(DLM_data@Mort[x]*exp(-DLM_data@Mort[x]*(1:ncol(Carr))))[ncol(Carr):1]
  #Warr<-Warr/sum(Warr)
  #TAC<-apply(t(matrix(Warr,nrow=ncol(Carr),ncol=reps))*Carr,1,sum)
  TACfilter(TAC)
}
class(Rcontrol2)<-"DLM_output"

GB_CC<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Cref,DLM_data@Cat"
  Catrec<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  TAC<-trlnorm(reps,DLM_data@Cref[x],DLM_data@CV_Cref)
  TAC[TAC>(1.2*Catrec)]<-1.2*Catrec
  TAC[TAC<(0.8*Catrec)]<-0.8*Catrec
  TACfilter(TAC)
}
class(GB_CC)<-"DLM_output"

GB_slope<-function(x,DLM_data,reps=100,yrsmth=5,lambda=1){
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@CV_Cat, DLM_data@Ind"
  Catrec<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  I_hist<-DLM_data@Ind[x,ind]
  yind<-1:yrsmth
  slppar<-summary(lm(I_hist~yind))$coefficients[2,1:2]
  Islp <-rnorm(reps,slppar[1],slppar[2])
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-rlnorm(reps,mconv(MuC,DLM_data@CV_Cat[x]*MuC),sdconv(MuC,DLM_data@CV_Cat[x]*MuC))
  TAC<-Cc*(1+lambda*Islp)
  TAC[TAC>(1.2*Catrec)]<-1.2*Catrec
  TAC[TAC<(0.8*Catrec)]<-0.8*Catrec
  TACfilter(TAC)
}
class(GB_slope)<-"DLM_output"

GB_target<-function(x,DLM_data,reps=100,w=0.5){
  dependencies="DLM_data@Cat, DLM_data@Cref, DLM_data@Iref, DLM_data@Ind"
  Catrec<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  TACtarg<-trlnorm(reps,DLM_data@Cref[x],DLM_data@CV_Cref)
  Itarg<-trlnorm(reps,DLM_data@Iref[x],DLM_data@CV_Iref)
  Iav<-mean(DLM_data@Ind[x,(length(DLM_data@Ind[x,])-4):length(DLM_data@Ind[x,])],na.rm=T)
  Irec<-mean(DLM_data@Ind[x,(length(DLM_data@Ind[x,])-3):length(DLM_data@Ind[x,])],na.rm=T)
  I0<-0.2*Iav
  TAC<-rep(NA,reps)
  if(Irec>I0)TAC<-TACtarg*(w+(1-w)*((Irec-I0)/(Itarg-I0)))
  if(Irec<I0)TAC<-TACtarg*w*(Irec/I0)^2
  TAC[TAC>(1.2*Catrec)]<-1.2*Catrec
  TAC[TAC<(0.8*Catrec)]<-0.8*Catrec
  TACfilter(TAC)
}
class(GB_target)<-"DLM_output"

CC1<-function(x,DLM_data,reps=100,yrsmth=5,xx=0){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat"
  C_dat<-DLM_data@Cat[x,(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)]
  TAC<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5)) # mean catches over the interval
  TACfilter(TAC)
}  
class(CC1)<-"DLM_output"

CC4<-function(x,DLM_data,reps=100,yrsmth=5,xx=0.3){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat"
  C_dat<-DLM_data@Cat[x,(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)]
  TAC<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5)) # mean catches over the interval
  TACfilter(TAC)
}  
class(CC4)<-"DLM_output"


LstepCC1<-function(x,DLM_data,reps=100,yrsmth=5,xx=0,stepsz=0.05,llim=c(0.96,0.98,1.05)){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat, DLM_data@ML"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2]
  if(length(DLM_data@Year)==ylast+1) {TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  }else{TACstar<-rep(DLM_data@MPrec[x],reps)}
  step<-stepsz*TACstar
  Lrecent<-mean(DLM_data@ML[ind])
  Lave<-mean(DLM_data@ML[ind3])
  rat<-Lrecent/Lave
  if(rat<llim[1]){TAC<-TACstar-2*step
  }else if(rat<llim[2]){TAC<-TACstar-step
  }else if(rat>llim[3]){TAC<-TACstar+step
  }else{TAC<-TACstar
  }
  TACfilter(TAC)
}  
class(LstepCC1)<-"DLM_output"

LstepCC4<-function(x,DLM_data,reps=100,yrsmth=5,xx=0.3,stepsz=0.05,llim=c(0.96,0.98,1.05)){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat, DLM_data@ML"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2]
  if(length(DLM_data@Year)==ylast+1) {TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  }else{TACstar<-rep(DLM_data@MPrec[x],reps)}
  step<-stepsz*TACstar
  Lrecent<-mean(DLM_data@ML[ind])
  Lave<-mean(DLM_data@ML[ind3])
  rat<-Lrecent/Lave
  if(rat<llim[1]){TAC<-TACstar-2*step
  }else if(rat<llim[2]){TAC<-TACstar-step
  }else if(rat>llim[3]){TAC<-TACstar+step
  }else{TAC<-TACstar
  }
  TACfilter(TAC)
}  
class(LstepCC4)<-"DLM_output"

Ltarget1<-function(x,DLM_data,reps=100,yrsmth=5,xx=0,xL=1.05){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat, DLM_data@ML"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2]
  TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  Lrecent<-mean(DLM_data@ML[ind])
  Lave<-mean(DLM_data@ML[ind3])
  L0<-0.9*Lave
  Ltarget<-xL*Lave
  if(Lrecent>L0){TAC<-0.5*TACstar*(1+((Lrecent-L0)/(Ltarget-L0)))
  }else{TAC<-0.5*TACstar*(Lrecent/L0)^2                  
  }
  TACfilter(TAC)
}  
class(Ltarget1)<-"DLM_output"

Ltarget4<-function(x,DLM_data,reps=100,yrsmth=5,xx=0.2,xL=1.15){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat, DLM_data@ML"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2]
  TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  Lrecent<-mean(DLM_data@ML[ind])
  Lave<-mean(DLM_data@ML[ind3])
  L0<-0.9*Lave
  Ltarget<-xL*Lave
  if(Lrecent>L0){TAC<-0.5*TACstar*(1+((Lrecent-L0)/(Ltarget-L0)))
  }else{TAC<-0.5*TACstar*(Lrecent/L0)^2                  
  }
  TACfilter(TAC)
}  
class(Ltarget4)<-"DLM_output"


Itarget1<-function(x,DLM_data,reps=100,yrsmth=5,xx=0,Imulti=1.5){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2] 
  TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  Irecent<-mean(DLM_data@Ind[x,ind])
  Iave<-mean(DLM_data@Ind[x,ind3])
  Itarget<-Iave*Imulti
  I0<-0.8*Iave
  if(Irecent>I0){TAC<-0.5*TACstar*(1+((Irecent-I0)/(Itarget-I0)))
  }else{TAC<-0.5*TACstar*(Irecent/I0)^2}
  TACfilter(TAC)  
}  
class(Itarget1)<-"DLM_output"

Itarget4<-function(x,DLM_data,reps=100,yrsmth=5,xx=0.3,Imulti=2.5){
  dependencies="DLM_data@Cat, DLM_data@CV_Cat"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year) # recent 5 years
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  ind2<-((ylast-(yrsmth-1)):ylast) # historical 5 pre-projection years
  ind3<-((ylast-(yrsmth*2-1)):ylast) # historical 10 pre-projection years
  C_dat<-DLM_data@Cat[x,ind2] 
  TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  Irecent<-mean(DLM_data@Ind[x,ind])
  Iave<-mean(DLM_data@Ind[x,ind3])
  Itarget<-Iave*Imulti
  I0<-0.8*Iave
  if(Irecent>I0){TAC<-0.5*TACstar*(1+((Irecent-I0)/(Itarget-I0)))
  }else{TAC<-0.5*TACstar*(Irecent/I0)^2}
  TACfilter(TAC)  
}  
class(Itarget4)<-"DLM_output"


Islope1<-function(x,DLM_data,reps=100,yrsmth=5,lambda=0.4,xx=0.2){
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@CV_Cat, DLM_data@Ind"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  C_dat<-DLM_data@Cat[x,ind]
  if(length(DLM_data@Year)==ylast+1) {TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  }else{TACstar<-rep(DLM_data@MPrec[x],reps)}
    
  I_hist<-DLM_data@Ind[x,ind]
  yind<-1:yrsmth
  slppar<-summary(lm(I_hist~yind))$coefficients[2,1:2]
  Islp <-rnorm(reps,slppar[1],slppar[2])
  TAC<-TACstar*(1+lambda*Islp)
  TACfilter(TAC)
}
class(Islope1)<-"DLM_output"

Islope4<-function(x,DLM_data,reps=100,yrsmth=5,lambda=0.2,xx=0.4){
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@CV_Cat, DLM_data@Ind"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  ylast<-(DLM_data@LHYear-DLM_data@Year[1])+1 #last historical year
  C_dat<-DLM_data@Cat[x,ind]
  if(length(DLM_data@Year)==ylast+1) {TACstar<-(1-xx)*trlnorm(reps,mean(C_dat),DLM_data@CV_Cat/(yrsmth^0.5))
  }else{TACstar<-rep(DLM_data@MPrec[x],reps)}
  I_hist<-DLM_data@Ind[x,ind]
  yind<-1:yrsmth
  slppar<-summary(lm(I_hist~yind))$coefficients[2,1:2]
  Islp <-rnorm(reps,slppar[1],slppar[2])
  TAC<-TACstar*(1+lambda*Islp)
  TACfilter(TAC)
}
class(Islope4)<-"DLM_output"




IT10<-function(x,DLM_data,reps=100,yrsmth=5,mc=0.1){
  
  dependencies="DLM_data@Ind, DLM_data@Cat, DLMdata@CV_Ind, DLMdata@Iref"
  ind<-max(1,(length(DLM_data@Year)-yrsmth+1)):length(DLM_data@Year)
  
  
  
  deltaI<-mean(DLM_data@Ind[x,ind])/DLM_data@Iref[x]
  if(deltaI<(1-mc))deltaI<-1-mc
  if(deltaI>(1+mc))deltaI<-1+mc
  
  TAC<-DLM_data@MPrec[x]*deltaI*trlnorm(reps,1,DLM_data@CV_Ind[x])
  TAC
}
class(IT10)<-"DLM_output"

IT5<-function(x,DLM_data,reps=100,yrsmth=5,mc=0.05){
  
  dependencies="DLM_data@Ind, DLM_data@Cat, DLMdata@CV_Ind, DLMdata@Iref"
  ind<-max(1,(length(DLM_data@Year)-yrsmth+1)):length(DLM_data@Year)
  deltaI<-mean(DLM_data@Ind[x,ind])/DLM_data@Iref[x]
  if(deltaI<(1-mc))deltaI<-1-mc
  if(deltaI>(1+mc))deltaI<-1+mc
  
  TAC<-DLM_data@MPrec[x]*deltaI*trlnorm(reps,1,DLM_data@CV_Ind[x])
  TAC
}
class(IT5)<-"DLM_output"

ITM<-function(x,DLM_data,reps=100){
  
  dependencies="DLM_data@Ind, DLM_data@Cat, DLMdata@CV_Ind, DLMdata@Iref, DLMdata@Mort"
  mc<-(5+DLM_data@Mort[x]*25)/100
  if(mc>0.2)mc<-0.2
  yrsmth<-floor(4*(1/DLM_data@Mort[x])^(1/4))
  ind<-max(1,(length(DLM_data@Year)-yrsmth+1)):length(DLM_data@Year)
  
  deltaI<-mean(DLM_data@Ind[x,ind])/DLM_data@Iref[x]
  if(deltaI<(1-mc))deltaI<-1-mc
  if(deltaI>(1+mc))deltaI<-1+mc
  
  TAC<-DLM_data@MPrec[x]*deltaI*trlnorm(reps,1,DLM_data@CV_Ind[x])
  TAC
}
class(ITM)<-"DLM_output"

SPmod<-function(x,DLM_data,reps=100,alp=c(0.8,1.2),bet=c(0.8,1.2)){
  dependencies="DLM_data@Cat, DLM_data@Ind, DLM_data@Abun"
  Ir<-length(DLM_data@Ind[x,])
  Cr<-length(DLM_data@Cat[x,])
  rat<-trlnorm(reps,DLM_data@Ind[x,Ir],DLM_data@CV_Ind)/trlnorm(reps,DLM_data@Ind[x,Ir-1],DLM_data@CV_Ind)
  cct<-trlnorm(reps,DLM_data@Cat[x,Cr],DLM_data@CV_Cat)
  Abun<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun)
  TAC<-rep(NA,reps)
  TAC[rat<alp[1]]<-cct[rat<alp[1]]*bet[1]
  TAC[rat>alp[1]& rat<alp[2]]<-cct[rat>alp[1]& rat<alp[2]]
 
  cond<-rat>alp[2]
  reps2<-sum(cond)
  if(reps2>0){
    qq1<-trlnorm(reps2,DLM_data@Ind[x,Ir]/Abun,DLM_data@CV_Ind)
    bio1<-DLM_data@Ind[x,Ir-1]/qq1
    bio2<-DLM_data@Ind[x,Ir]/qq1
    cct1<-trlnorm(reps2,DLM_data@Cat[x,Cr-1],DLM_data@CV_Cat)
    PP<-bio2-bio1+cct1
    TAC[cond]<-bet[2]*PP
  }
  TACfilter(TAC)
}
class(SPmod)<-"DLM_output"

SPslope<-function(x,DLM_data,reps=100,yrsmth=4,alp=c(0.9,1.1),bet=c(1.5,0.9)){
  
  dependencies="DLM_data@Year, DLM_data@Cat, DLM_data@Ind, DLM_data@Abun"
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  yind<-1:yrsmth
  C_dat<-DLM_data@Cat[x,ind]
  B_dat<-DLM_data@Ind[x,ind]/DLM_data@Ind[x,ind[yrsmth]]*DLM_data@Abun[x]
  Pt_mu<-max(B_dat[yrsmth]-B_dat[yrsmth-1]+C_dat[yrsmth-1],tiny)
  Pt_1<-trlnorm(reps,Pt_mu,DLM_data@CV_Cat)
  It<-exp(predict(lm(log(B_dat)~yind),newdat=list(yind=yrsmth+1)))
  Ilast<-B_dat[yrsmth]
  MC <- max(mean(C_dat), tiny)
  Ct_1<-trlnorm(reps,MC,DLM_data@CV_Cat/(yrsmth^0.5)) # mean catches over the interval
 
  rat<-It/Ilast
  
  mult <- max((1-bet[1]*(Ilast-It)/Ilast), tiny)
  if(rat<alp[1])TAC<-mult*Ct_1
  if(rat>alp[1]&rat<alp[2])TAC<-Ct_1
  if(rat>alp[2])TAC<-bet[2]*Pt_1
  TACfilter(TAC)
}
class(SPslope)<-"DLM_output"

SBT1<-function(x,DLM_data,reps=100,yrsmth=10,k1=1.5,k2=3,gamma=1){
  dependencies="DLM_data@Cat, DLM_data@Year, DLM_data@Ind"
  Cr<-length(DLM_data@Cat[x,])
  cct<-trlnorm(reps,DLM_data@Cat[x,Cr],DLM_data@CV_Cat)
  ind<-(length(DLM_data@Year)-(yrsmth-1)):length(DLM_data@Year)
  I_hist<-DLM_data@Ind[x,ind]
  test<-summary(lm(I_hist~ind))$coefficients[2,1:2]
  lambda<-rnorm(reps,test[1],test[2])
  TAC<-cct*1+k2*lambda
  cond<-lambda<0
  TAC[cond]<-cct[cond]*1-k1*-lambda[cond]^gamma
  TACfilter(TAC) 
}
class(SBT1)<-"DLM_output"

SBT2<-function(x,DLM_data,reps=100,epsB=0.25,epsR=0.75,tauR=5,tauB=7,gamma=1){
  dependencies="DLM_data@Cref, DLM_data@Rec, DLM_data@Cat"
  #Bnow<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun)
  #testrat<-Bnow/DLM_data@Bref
  #Ctarg<-rep(NA,reps)
  #Ctarg[testrat>1]<-delta*testrat[testrat>1]^(1-epsB)
  #Ctarg[testrat<1]<-detla*testrat[testrat<1]^(1+epsB)
  Ctarg<-trlnorm(reps,DLM_data@Cref[x],DLM_data@CV_Cref)
  muR<-mean(DLM_data@Rec[x,(length(DLM_data@Rec[x,])-tauR+1):length(DLM_data@Rec[x,])])
  phi<-mean(DLM_data@Rec[x,(length(DLM_data@Rec[x,])-9):length(DLM_data@Rec[x,])])
  Rrat<-muR/phi
  deltaR<-rep(NA,reps)
  deltaR[Rrat>1]<-Rrat[Rrat>1]^(1-epsR)
  deltaR[Rrat<1]<-Rrat[Rrat<1]^(1+epsR)
  TAC<-0.5*(DLM_data@Cat[x,length(DLM_data@Cat[x,])]+Ctarg*deltaR)
  TACfilter(TAC) 
}
class(SBT2)<-"DLM_output"

DD<-function(x,DLM_data,reps=100){
#for(x in 1:nsim){
  dependencies="DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@Mort, DLM_data@CV_Mort. DLM_data@wla, DLM_data@ wlb"
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  a<-DLM_data@wla[x]
  b<-DLM_data@wlb[x]

  Winf=DLM_data@wla[x]*DLM_data@vbLinf[x]^DLM_data@wlb[x]
	age<-1:DLM_data@MaxAge
  la<-DLM_data@vbLinf[x]*(1-exp(-DLM_data@vbK[x]*((age-DLM_data@vbt0[x]))))
  wa<-DLM_data@wla[x]*la^DLM_data@wlb[x]
  a50V<-iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])
  yind<-(1:length(DLM_data@Cat[x,]))[!is.na(DLM_data@Cat[x,]+DLM_data@Ind[x,])]
  C_hist<-DLM_data@Cat[x,yind]
  E_hist<-C_hist/DLM_data@Ind[x,yind]
  E_hist<-E_hist/mean(E_hist)
  ny_DD<-length(C_hist)
  params<-log(c(DLM_data@Mort[x],mean(C_hist,na.rm=T),DLM_data@Mort[x]))
  k_DD<-ceiling(a50V)   # get age nearest to 50% vulnerability (ascending limb)  -------------
  k_DD[k_DD>DLM_data@MaxAge/2]<-ceiling(DLM_data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD<-(wa[k_DD+2]-Winf)/(wa[k_DD+1]-Winf)
  Alpha_DD<-Winf*(1-Rho_DD)
  So_DD<-exp(-DLM_data@Mort[x]) # get So survival rate
  wa_DD<-wa[k_DD]
  UMSYprior<-c(1-exp(-DLM_data@Mort[x]*0.5),0.3)
  opt<-optim(params,DD_R,opty=1,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,
                    ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,
                    C_hist=C_hist,UMSYprior=UMSYprior,
                    method="L-BFGS-B",
                    lower=log(exp(params)/20),upper=log(exp(params)*20),
                    hessian=TRUE)

  #Catfit<-DD_R(opt$par,opty=3,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  #plot(Catfit[,1],ylim=c(0,max(Catfit)))
  #lines(Catfit[,2],col="red")

  TAC<-rep(NA,reps)
  #samps<-rmvnorm(reps,opt$par,solve(opt$hessian)) # assuming log parameters are multivariate normal hessian approximation
  samps<-cbind(rnorm(reps,opt$par[1],((opt$par[1])^2)^0.5*0.1),rnorm(reps,opt$par[2],((opt$par[2])^2)^0.5*0.1),rnorm(reps,opt$par[3],((opt$par[3])^2)^0.5*0.1))
  if(reps==1)samps<-matrix(c(opt$par[1],opt$par[2],opt$par[3]),nrow=1)
  for(i in 1:reps)TAC[i]<-DD_R(samps[i,],opty=2,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  TACfilter(TAC)
}
class(DD)<-"DLM_output"


DD4010<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@Mort, DLM_data@CV_Mort. DLM_data@wla, DLM_data@ wlb"
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  a<-DLM_data@wla[x]
  b<-DLM_data@wlb[x]
  
  Winf=DLM_data@wla[x]*DLM_data@vbLinf[x]^DLM_data@wlb[x]
  age<-1:DLM_data@MaxAge
  la<-DLM_data@vbLinf[x]*(1-exp(-DLM_data@vbK[x]*((age-DLM_data@vbt0[x]))))
  wa<-DLM_data@wla[x]*la^DLM_data@wlb[x]
  a50V<-iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])
  yind<-(1:length(DLM_data@Cat[x,]))[!is.na(DLM_data@Cat[x,]+DLM_data@Ind[x,])]
  C_hist<-DLM_data@Cat[x,yind]
  E_hist<-DLM_data@Ind[x,yind]
  E_hist<-C_hist/E_hist
  E_hist<-E_hist/mean(E_hist)
  ny_DD<-length(C_hist)
  params<-log(c(DLM_data@Mort[x],mean(C_hist,na.rm=T),DLM_data@Mort[x]))
  k_DD<-ceiling(a50V)   # get age nearest to 50% vulnerability (ascending limb)  -------------
  k_DD[k_DD>DLM_data@MaxAge/2]<-ceiling(DLM_data@MaxAge/2)  # to stop stupidly high estimates of age at 50% vulnerability
  Rho_DD<-(wa[k_DD+2]-Winf)/(wa[k_DD+1]-Winf)
  Alpha_DD<-Winf*(1-Rho_DD)
  So_DD<-exp(-DLM_data@Mort[x]) # get So survival rate
  wa_DD<-wa[k_DD]
  UMSYprior<-c(1-exp(-DLM_data@Mort*0.5),0.3)
  opt<-optim(params,DD_R,opty=1,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,
             ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,
             C_hist=C_hist,UMSYprior=UMSYprior,
             method="L-BFGS-B",
             lower=log(exp(params)/20),upper=log(exp(params)*20),
             hessian=TRUE)
  
  #Catfit<-DD_R(opt$par,opty=3,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  #plot(Catfit[,1],ylim=c(0,max(Catfit)))
  #lines(Catfit[,2],col="red")
  
  TAC<-rep(NA,reps)
  dep<-rep(NA,reps)
  #samps<-rmvnorm(reps,opt$par,solve(opt$hessian)) # assuming log parameters are multivariate normal hessian approximation
  samps<-cbind(rnorm(reps,opt$par[1],((opt$par[1])^2)^0.5*0.1),rnorm(reps,opt$par[2],((opt$par[2])^2)^0.5*0.1),rnorm(reps,opt$par[3],((opt$par[3])^2)^0.5*0.1))
  if(reps==1)samps<-matrix(c(opt$par[1],opt$par[2],opt$par[3]),nrow=1)
  for(i in 1:reps)TAC[i]<-DD_R(samps[i,],opty=2,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  for(i in 1:reps)dep[i]<-DD_R(samps[i,],opty=3,So_DD=So_DD,Alpha_DD=Alpha_DD,Rho_DD=Rho_DD,ny_DD=ny_DD,k_DD=k_DD,wa_DD=wa_DD,E_hist=E_hist,C_hist=C_hist,UMSYprior=UMSYprior)
  cond1<-!is.na(dep) & dep<0.4 & dep>0.1
  cond2<-!is.na(dep) & dep<0.1
  TAC[cond1]<-TAC[cond1]*(dep[cond1]-0.1)/0.3
  TAC[cond2]<-TAC[cond2]*tiny # this has to still be stochastic albeit very small
  TACfilter(TAC)
}
class(DD4010)<-"DLM_output"

DD_R<-function(params,opty,So_DD,Alpha_DD,Rho_DD,ny_DD,k_DD,wa_DD,E_hist,C_hist,UMSYprior){
  UMSY_DD=exp(params[1])
  MSY_DD=exp(params[2])
  q_DD=exp(params[3])
  SS_DD=So_DD*(1-UMSY_DD)    # Initialise for UMSY, MSY and q leading.
  Spr_DD=(SS_DD*Alpha_DD/(1-SS_DD)+wa_DD)/(1-Rho_DD*SS_DD)
  DsprDu_DD=-So_DD*(Rho_DD/(1-Rho_DD*SS_DD)*Spr_DD+1/(1-Rho_DD*SS_DD)*(Alpha_DD/(1-SS_DD)+SS_DD*Alpha_DD/(1-SS_DD)^2))
  Arec_DD=1/(((1-UMSY_DD)^2)*(Spr_DD+UMSY_DD*DsprDu_DD))
  Brec_DD=UMSY_DD*(Arec_DD*Spr_DD-1/(1-UMSY_DD))/MSY_DD
  Spr0_DD=(So_DD*Alpha_DD/(1-So_DD)+wa_DD)/(1-Rho_DD*So_DD)
  Ro_DD=(Arec_DD*Spr0_DD-1)/(Brec_DD*Spr0_DD)
  Bo_DD=Ro_DD*Spr0_DD
  No_DD=Ro_DD/(1-So_DD)

  B_DD<-rep(NA,ny_DD+1)
  N_DD<-rep(NA,ny_DD+1)
  R_DD<-rep(NA,ny_DD+k_DD)
  Cpred_DD<-rep(NA,ny_DD)

  B_DD[1]=Bo_DD
  N_DD[1]=No_DD
  R_DD[1:k_DD]=Ro_DD

  for(tt in 1:ny_DD){

    Surv_DD=So_DD*exp(-q_DD*E_hist[tt])
    Cpred_DD[tt]=B_DD[tt]*(1-exp(-q_DD*E_hist[tt]))
    Sp_DD=B_DD[tt]-Cpred_DD[tt]
    R_DD[tt+k_DD]=Arec_DD*Sp_DD/(1+Brec_DD*Sp_DD);
    B_DD[tt+1]=Surv_DD*(Alpha_DD*N_DD[tt]+Rho_DD*B_DD[tt])+wa_DD*R_DD[tt+1]
    N_DD[tt+1]=Surv_DD*N_DD[tt]+R_DD[tt+1]

  }
  Cpred_DD[Cpred_DD<tiny]<-tiny

  if(opty==1){
    test<-dnorm(log(Cpred_DD),log(C_hist),0.25,log=T)
    test2<-dlnorm(UMSY_DD,log(UMSYprior[1]),UMSYprior[2],log=T)
    test[is.na(test)]<--1000
    test[test==(-Inf)]<--1000
    if(is.na(test2)|test2==-Inf|test2==Inf)test2<-1000
    return(-sum(test,test2))      # return objective function
  }else if(opty==2){                                  # return MLE TAC estimate
    UMSY_DD*B_DD[ny_DD]
  }else if(opty==3){ 
    B_DD[tt+1]/Bo_DD
  }else{
    cbind(C_hist,Cpred_DD)                           # return observations vs predictions
  }
}

DBSRA<-function(x,DLM_data,reps=100){  # returns a vector of DBSRA estimates of the TAC for a particular simulation x
  # for(x in 1:nsim){
  dependencies="DLM_data@Cat, DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M,DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0, DLM_data@L50"
  C_hist<-DLM_data@Cat[x,]
  TAC<-rep(NA,reps)
  DBSRAcount<-1
  if (is.na(DLM_data@Dep[x]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  while(DBSRAcount<(reps+1)){
    depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
    Bt_K<-rbeta(100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
    Bt_K<-Bt_K[Bt_K>0.00999&Bt_K<0.99001][1] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
    Mdb<-trlnorm(100,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Mdb<-Mdb[Mdb<0.9][1]    # !!!! maximum M is 0.9   interval censor
    if(is.na(Mdb))Mdb<-0.9  # !!!! maximum M is 0.9   absolute limit
    FMSY_M<-trlnorm(1,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
    BMSY_K<-rbeta(100,alphaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    tryBMSY_K <- BMSY_K[BMSY_K > 0.05 & BMSY_K < 0.95][1] # interval censor (0.05,0.95) as in Dick and MacCall, 2011
	if (is.na(tryBMSY_K)) {
	  Min <- min(BMSY_K, na.rm=TRUE)
	  Max <- max(BMSY_K, na.rm=TRUE)
	  if (Max <= 0.05) BMSY_K <- 0.05
	  if (Min >= 0.95) BMSY_K <- 0.95
	}
	if (!is.na(tryBMSY_K)) BMSY_K <- tryBMSY_K
	
    adelay<-max(floor(iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])),1)
    opt<-optimize(DBSRAopt,log(c(0.01*mean(C_hist),1000*mean(C_hist))),C_hist=C_hist,nys=length(C_hist),Mdb=Mdb,
              FMSY_M=FMSY_M,BMSY_K=BMSY_K,Bt_K=Bt_K,adelay=adelay,tol=0.01)
    # if(opt$objective<0.1){
    Kc<-exp(opt$minimum)
    BMSYc<-Kc*BMSY_K
    FMSYc<-Mdb*FMSY_M
    UMSYc<-(FMSYc/(FMSYc+Mdb))*(1-exp(-(FMSYc+Mdb)))
    MSYc<-Kc*BMSY_K*UMSYc
    TAC[DBSRAcount]<-UMSYc*Kc*Bt_K
    DBSRAcount<-DBSRAcount+1
    # }
	 
	# print(DBSRAcount)
  } # end of reps
  TACfilter(TAC)
   # print(x)
  # }
}  # end of DBSRA_apply
class(DBSRA)<-"DLM_output"

DBSRA_40<-function(x,DLM_data,reps=100){  # returns a vector of DBSRA estimates of the TAC for a particular simulation x
  dependencies="DLM_data@Cat, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0, DLM_data@L50"
  C_hist<-DLM_data@Cat[x,]
  TAC<-rep(NA,reps)
  DBSRAcount<-1
  if (is.na(DLM_data@Dep[x]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  while(DBSRAcount<(reps+1)){
    depo<-0.4
    Bt_K<-rbeta(100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
    Bt_K<-Bt_K[Bt_K>0.00999&Bt_K<0.99001][1] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
    Mdb<-rlnorm(100,mconv(DLM_data@Mort[x],DLM_data@CV_Mort[x]*DLM_data@Mort[x]),sdconv(DLM_data@Mort[x],DLM_data@CV_Mort[x]*DLM_data@Mort[x]))   # log space stdev 0.4 as in Dick and MacCall 2011
    Mdb<-Mdb[Mdb<0.9][1]    # !!!! maximum M is 0.9   interval censor
    if(is.na(Mdb))Mdb<-0.9  # !!!! maximum M is 0.9   absolute limit
    FMSY_M<-trlnorm(1,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
    BMSY_K<-rbeta(100,alphaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    BMSY_K<-BMSY_K[BMSY_K>0.05&BMSY_K<0.95][1] # interval censor (0.05,0.95) as in Dick and MacCall, 2011
    adelay<-max(floor(iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])),1)
    opt<-optimize(DBSRAopt,log(c(0.1*mean(C_hist),1000*mean(C_hist))),C_hist=C_hist,nys=length(C_hist),Mdb=Mdb,
              FMSY_M=FMSY_M,BMSY_K=BMSY_K,Bt_K=Bt_K,adelay=adelay,tol=0.01)
    #if(opt$objective<0.1){
      Kc<-exp(opt$minimum)
      BMSYc<-Kc*BMSY_K
      FMSYc<-Mdb*FMSY_M
      UMSYc<-(FMSYc/(FMSYc+Mdb))*(1-exp(-(FMSYc+Mdb)))
      MSYc<-Kc*BMSY_K*UMSYc
      TAC[DBSRAcount]<-UMSYc*Kc*Bt_K
      DBSRAcount<-DBSRAcount+1
    #}
  } # end of reps
  TACfilter(TAC)
}  # end of DBSRA_apply
class(DBSRA_40)<-"DLM_output"

DBSRA_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Cat, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0, DLM_data@L50, DLM_data@CAL, DLM_data@Year, DLM_data@Cat"
  C_hist<-DLM_data@Cat[x,]
  TAC<-rep(NA,reps)
  DBSRAcount<-1
  if (is.na(DLM_data@Dep[x]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  while(DBSRAcount<(reps+1)){
    Linfc<-trlnorm(1,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
    Kc<-trlnorm(1,DLM_data@vbK[x],DLM_data@CV_vbK[x])
    Mdb<-trlnorm(100,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Mdb<-Mdb[Mdb<0.9][1]    # !!!! maximum M is 0.9   interval censor
    if(is.na(Mdb))Mdb<-0.9  # !!!! maximum M is 0.9   absolute limit
    Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=1,MLtype="dep")
    FM<-Z-Mdb
    FM[FM<0]<-0.01
    nyears<-length(DLM_data@Year)
    Ct1<-mean(DLM_data@Cat[x,1:3])
    Ct2<-mean(DLM_data@Cat[x,(nyears-2):nyears])
    dep<-c(Ct1,Ct2)/(1-exp(-FM[,c(1,2)]))
    Bt_K<-dep[2]/dep[1]
    if(Bt_K<0.01)Bt_K<-0.01       # interval censor / temporary hack to avoid doing multiple depletion estimates that would take far too long
    if(Bt_K>0.99)Bt_K<-0.99       # interval censor / temporary hack to avoid doing multiple depletion estimates that would take far too long

    FMSY_M<-trlnorm(1,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
    BMSY_K<-rbeta(100,alphaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    BMSY_K<-BMSY_K[BMSY_K>0.05&BMSY_K<0.95][1] # interval censor (0.05,0.95) as in Dick and MacCall, 2011
    adelay<-max(floor(iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])),1)
    opt<-optimize(DBSRAopt,log(c(0.1*mean(C_hist),1000*mean(C_hist))),C_hist=C_hist,nys=length(C_hist),Mdb=Mdb,
              FMSY_M=FMSY_M,BMSY_K=BMSY_K,Bt_K=Bt_K,adelay=adelay,tol=0.01)
    if(opt$objective<0.1){
      Kc<-exp(opt$minimum)
      BMSYc<-Kc*BMSY_K
      FMSYc<-Mdb*FMSY_M
      UMSYc<-(FMSYc/(FMSYc+Mdb))*(1-exp(-(FMSYc+Mdb)))
      MSYc<-Kc*BMSY_K*UMSYc
      TAC[DBSRAcount]<-UMSYc*Kc*Bt_K
      DBSRAcount<-DBSRAcount+1
    }
  } # end of reps
  TACfilter(TAC)
}
class(DBSRA_ML)<-"DLM_output"

DBSRA4010<-function(x,DLM_data,reps=100){  # returns a vector of DBSRA estimates of the TAC for a particular simulation x
  #for(x in 1:nsim){
  dependencies="DLM_data@Cat, DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@FMSY_M, DLM_data@CV_FMSY_M,DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0, DLM_data@L50"
  C_hist<-DLM_data@Cat[x,]
  TAC<-rep(NA,reps)
  DBSRAcount<-1
  if (is.na(DLM_data@Dep[x]) | is.na(DLM_data@CV_Dep[x])) return(NA)
  while(DBSRAcount<(reps+1)){
    depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
    Bt_K<-rbeta(100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
    Bt_K<-Bt_K[Bt_K>0.00999&Bt_K<0.99001][1] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
    Mdb<-trlnorm(100,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Mdb<-Mdb[Mdb<0.9][1]    # !!!! maximum M is 0.9   interval censor
    if(is.na(Mdb))Mdb<-0.9  # !!!! maximum M is 0.9   absolute limit
    FMSY_M<-trlnorm(1,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
    BMSY_K<-rbeta(100,alphaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]*DLM_data@BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
    BMSY_K<-BMSY_K[BMSY_K>0.05&BMSY_K<0.95][1] # interval censor (0.05,0.95) as in Dick and MacCall, 2011
    adelay<-max(floor(iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])),1)
    opt<-optimize(DBSRAopt,log(c(0.01*mean(C_hist),1000*mean(C_hist))),C_hist=C_hist,nys=length(C_hist),Mdb=Mdb,
                  FMSY_M=FMSY_M,BMSY_K=BMSY_K,Bt_K=Bt_K,adelay=adelay,tol=0.01)
    #if(opt$objective<0.1){
    Kc<-exp(opt$minimum)
    BMSYc<-Kc*BMSY_K
    FMSYc<-Mdb*FMSY_M
    UMSYc<-(FMSYc/(FMSYc+Mdb))*(1-exp(-(FMSYc+Mdb)))
    MSYc<-Kc*BMSY_K*UMSYc
    TAC[DBSRAcount]<-UMSYc*Kc*Bt_K
    # 40-10 rule
    if(Bt_K<0.4 & Bt_K>0.1)TAC[DBSRAcount]<-TAC[DBSRAcount]*(Bt_K-0.1)/0.3
    if(Bt_K<0.1)TAC[DBSRAcount]<-TAC[DBSRAcount]*tiny # this has to still be stochastic albeit very small
    DBSRAcount<-DBSRAcount+1
    #}
  } # end of reps
  TACfilter(TAC)
  #}
}  # end of DBSRA_apply
class(DBSRA4010)<-"DLM_output"

DBSRAopt<-function(lnK,C_hist,nys,Mdb,FMSY_M,BMSY_K,Bt_K,adelay){         # the optimization for B0 given DBSRA assumptions
  Kc<-exp(lnK)
  n<-getn(BMSY_K)
  g<-gety(n)
  FMSY<-FMSY_M*Mdb
  UMSY<-(FMSY/(FMSY+Mdb))*(1-exp(-(FMSY+Mdb)))
  MSY<-Kc*BMSY_K*UMSY
  # Bjoin rules from Dick & MacCall 2011  ---------
  Bjoin_K<-0.5
  if(BMSY_K<0.3)Bjoin_K<-0.5*BMSY_K
  if(BMSY_K>0.3&BMSY_K<0.5)Bjoin_K<-0.75*BMSY_K-0.075
  Bjoin<-Bjoin_K*Kc
  PBjoin<-prodPTF(Bjoin_K,n,MSY)
  cp<-(1-n)*g*MSY*(Bjoin^(n-2))*Kc^-n
  Bc<-rep(NA,nys)
  Bc[1]<-Kc
  obj<-0
  for(yr in 2:nys){
    yref<-max(1,yr-adelay)
    if(Bc[yref]>Bjoin|BMSY_K>0.5){
      Bc[yr]<-Bc[yr-1]+g*MSY*(Bc[yref]/Kc)-g*MSY*(Bc[yref]/Kc)^n-C_hist[yr-1]
    }else{
      Bc[yr]<-Bc[yr-1]+Bc[yref]*((PBjoin/Bjoin)+cp*(Bc[yref]-Bjoin))-C_hist[yr-1]
    }
    if(Bc[yr]<0)obj<-obj+log(-Bc[yr])
    Bc[yr]<-max(0.000001,Bc[yr])
  }
  obj+((Bc[nys]/Kc)-Bt_K)^2
}  # end of DBSRA optimization function

# DCAC =====================================================================================================================
# Variables for parallel processing of DCAC_apply function using sfSapply()
C_tot<-nyearsDCAC<-NULL

DCAC<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@Dt, DLM_data@CV_Dt, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0"
  if (is.na(DLM_data@BMSY_B0[x]) | is.na(DLM_data@CV_BMSY_B0[x])) return(NA)
  C_tot<-DLM_data@AvC[x]*DLM_data@t[x]
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  FMSY_M<-trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x]) # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  Bt_K<-trlnorm(reps,DLM_data@Dt[x],DLM_data@CV_Dt[x])
  BMSY_K<-rbeta(reps,alphaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TACfilter(C_tot/(DLM_data@t[x]+((1-Bt_K)/(BMSY_K*FMSY_M*Mdb))))
} # end of DCAC
class(DCAC)<-"DLM_output"

DCAC4010<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@Dt, DLM_data@CV_Dt, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0"
  C_tot<-DLM_data@AvC[x]*DLM_data@t[x]
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  FMSY_M<-trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x]) # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  Bt_K<-trlnorm(reps,DLM_data@Dt[x],DLM_data@CV_Dt[x])
  BMSY_K<-rbeta(reps,alphaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TAC<-C_tot/(DLM_data@t[x]+((1-Bt_K)/(BMSY_K*FMSY_M*Mdb)))
  # 40-10 rule
  cond1<-Bt_K<0.4 & Bt_K>0.1
  cond2<-Bt_K<0.1
  if (length(cond1)>0) TAC[cond1]<-TAC[cond1]*(Bt_K[cond1]-0.1)/0.3
  if (length(cond2)>0) TAC[cond2]<-TAC[cond2]*tiny # this has to still be stochastic albeit very small
  if (length(cond1) <1 & length(cond2) < 1) return(NA)
  TACfilter(TAC)
  
} # end of DCAC
class(DCAC4010)<-"DLM_output"

DCAC_40<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0"
  C_tot<-DLM_data@AvC[x]*DLM_data@t[x]
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  FMSY_M<-trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
  Bt_K<-0.4
  BMSY_K<-rbeta(reps,alphaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TACfilter(C_tot/(DLM_data@t[x]+((1-Bt_K)/(BMSY_K*FMSY_M*Mdb))))
} # end of DCAC40
class(DCAC_40)<-"DLM_output"

DCAC_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0, DLM_data@Year, DLM_data@CAL"
  C_tot<-DLM_data@AvC[x]*DLM_data@t[x]
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # default CV of 0.5 as in MacCall 2009
  FMSY_M<-trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x]) # standard deviation of 0.2 - referred to as 'standard error' in MacCall 2009
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=reps,MLtype="dep")
  FM<-Z-Mdb
  FM[FM<0]<-0.01
  nyears<-length(DLM_data@Year)
  Ct1<-mean(DLM_data@Cat[x,1:3])
  Ct2<-mean(DLM_data@Cat[x,(nyears-2):nyears])
  dep<-rep(c(Ct1,Ct2),each=reps)/(1-exp(-FM[,c(1,2)]))
  Bt_K<-dep[,2]/dep[,1]
  BMSY_K<-rbeta(reps,alphaconv(DLM_data@BMSY_B0[x],DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x]),betaconv(DLM_data@BMSY_B0[x],DLM_data@BMSY_B0[x]*DLM_data@CV_BMSY_B0[x])) #0.045 corresponds with mu=0.4 and quantile(BMSY_K,c(0.025,0.975)) =c(0.31,0.49) as in Dick and MacCall 2011
  TAC<-C_tot/(DLM_data@t[x]+((1-Bt_K)/(BMSY_K*FMSY_M*Mdb)))
  TACfilter(TAC)
} # end of DCAC_ML
class(DCAC_ML)<-"DLM_output"

BK<-function(x,DLM_data,reps=100){   # Beddington and Kirkwood life-history analysis ==============================================
  dependencies="DLM_data@LFC, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@Abun, DLM_data@CV_Abun, DLM_data@vbK, DLM_data@CV_vbK"
  Lc<-trlnorm(reps*10,DLM_data@LFC[x],0.2)
  Linfc<-trlnorm(reps*10,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Ac<-trlnorm(reps*10,DLM_data@Abun[x],DLM_data@CV_Abun[x])
  Kc<-trlnorm(reps*10,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  TAC<-Ac*(0.6*Kc)/(0.67-(Lc/Linfc))         # robustifying for use in MSE
  TACfilter(TAC[TAC>0][1:reps])              # Interval censor only those positive catch recommendations

}  # end of BK
class(BK)<-"DLM_output"


BK_CC<-function(x,DLM_data,reps=100,Fmin=0.005){
  dependencies="DLM_data@LFC, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@CAA, DLM_data@Mort"
  Lc<-trlnorm(reps,DLM_data@LFC[x],0.2)
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Mdb<-trlnorm(reps*10,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Zdb<-CC(x,DLM_data,reps=reps*10)
  Fdb<-Zdb-Mdb
  ind<-(1:(reps*10))[Fdb>Fmin][1:reps]
  Fdb<-Fdb[ind]
  Mdb<-Mdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0 ) {
    Mdb[is.na(ind)] <- trlnorm(SM,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Fdb[is.na(ind)] <- Fmin
  }	
  
  Ac<-Cc/(1-exp(-Fdb))
  TAC<-Ac*(0.6*Kc)/(0.67-(Lc/Linfc))  
  TACfilter(TAC) 
    
}  # end of BK_CC
class(BK_CC)<-"DLM_output"


BK_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@LFC, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@CAL, DLM_data@Mort"
  Lc<-trlnorm(reps*10,DLM_data@LFC[x],0.2)
  Linfc<-trlnorm(reps*10,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps*10,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Mdb<-trlnorm(reps*10,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=reps*2,MLtype="F")
  FM<-Z-Mdb
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Ac<-Cc/(1-exp(-FM))
  FMSY<-(0.6*Kc)/(0.67-(Lc/Linfc))  # robustifying for use in MSETAC<-Ac*FMSY
  TAC<-Ac*FMSY
  TAC[TAC>0&TAC<(mean(TAC,na.rm=T)+3*sd(TAC,na.rm=T))][1:reps]
}
class(BK_ML)<-"DLM_output"

Fratio<-function(x,DLM_data,reps=100){  # FMSY / M ratio method e.g. Gulland ===============================================================================
  depends="DLM_data@Abun,DLM_data@CV_Abun,DLM_data@FMSY_M, DLM_data@CV_FMSY_M,DLM_data@Mort,DLM_data@CV_Mort"
  Ac<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun[x])
  TACfilter(Ac*trlnorm(reps,DLM_data@Mort[x],
            DLM_data@CV_Mort[x])*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x]))
} # end of Fratio
class(Fratio)<-"DLM_output"

Fratio4010<-function(x,DLM_data,reps=100){  # FMSY / M ratio method e.g. Gulland ===============================================================================
  dependencies="DLM_data@Abun, DLM_data@CV_Abun, DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@Dep"
  Ac<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun[x])
  TAC<-Ac*trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])
  Bt_K<-trlnorm(reps,DLM_data@Dt[x],DLM_data@CV_Dt[x])
  # 40-10 rule
  cond1<-Bt_K<0.4 & Bt_K>0.1
  cond2<-Bt_K<0.1
  TAC[cond1]<-TAC[cond1]*(Bt_K[cond1]-0.1)/0.3
  TAC[cond2]<-TAC[cond2]*tiny # this has to still be stochastic albeit very small
  TACfilter(TAC)  
} # end of Fratio
class(Fratio4010)<-"DLM_output"

Fratio_CC<-function(x,DLM_data,reps=100,Fmin=0.005){ # FMSY / M ratio method using catch curve analysis to determine current abundance ==================================
# for (x in 1:nsim) {
  dependencies=" DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@Cat, DLM_data@CV_Cat, DLM_data@CAA"
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Mdb<-trlnorm(reps*10,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  Zdb<-CC(x,DLM_data,reps=reps*10)
  Fdb<-Zdb-Mdb
  ind <-(1:(reps*10))[Fdb>0.005][1:reps]
  
  Fdb<-Fdb[ind]
  Mdb<-Mdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0 ) {
    Mdb[is.na(ind)] <- trlnorm(SM,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Fdb[is.na(ind)] <- Fmin
  }	
  
  Ac<-Cc/(1-exp(-Fdb))
  TAC<-Ac*Mdb*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])

  TACfilter(TAC)
# }
} # end of Fratio_CC
class(Fratio_CC)<-"DLM_output"


Fratio_ML<-function(x,DLM_data,reps=100){
  dependencies=" DLM_data@FMSY_M, DLM_data@CV_FMSY_M, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@Cat, DLM_data@CV_Cat, DLM_data@CAL"
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])   # CV of 0.5 as in MacCall 2009
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=reps,MLtype="F")
  FM<-Z-Mdb
  Ac<-Cc/(1-exp(-FM))
  TAC<-Ac*trlnorm(reps,DLM_data@FMSY_M[x],DLM_data@CV_FMSY_M[x])*trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  TACfilter(TAC)
}
class(Fratio_ML)<-"DLM_output"

SPMSY<-function(x,DLM_data,reps=100){  # Martell and Froese 2012 Schaefer SP estimate of MSY given priors on r, k and depletion
  #for(x in 1:100){
  dependencies="DLM_data@MaxAge, DLM_data@vbK, DLM_data@L50, DLM_data@Cat"
  maxage<-DLM_data@MaxAge
  nsamp<-reps*200

  # Froese 2012 http://www.fishbase.de/rfroese/Catch-MSY_SummaryFinal.doc
  rule<-rep(4,3)

  if(DLM_data@vbK[x]>0.3){   # K rules
   rule[1]<-1
  }else if(DLM_data@vbK[x]<0.3 & DLM_data@vbK[x]>0.16){
   rule[1]<-2
  }else if(DLM_data@vbK[x]<0.16 & DLM_data@vbK[x]>0.05){
   rule[1]<-3
  }
  AM<-iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x])
  if(AM<1.5){  # Age at maturity rules
   rule[2]<-1
  }else if(AM<4.5 & AM>1.5){
   rule[2]<-2
  }else if(AM<10 & AM>4.5){
   rule[2]<-3
  }

  if(DLM_data@MaxAge<4){   # Maximum age rules
   rule[1]<-1
  }else if(DLM_data@MaxAge<11 & DLM_data@MaxAge>3){
   rule[1]<-2
  }else if(DLM_data@MaxAge<31 & DLM_data@MaxAge>10){
   rule[1]<-3
  }

  if(mean(rule)<1.5) rsamp<-runif(nsamp,0.6,1.5)
  if(mean(rule)>1.5&mean(rule)<2.5)rsamp<-runif(nsamp,0.2,1)
  if(mean(rule)>2.5&mean(rule)<3.5)rsamp<-runif(nsamp,0.05,0.5)
  if(mean(rule)>3.5) rsamp<-runif(nsamp,0.015,0.1)

  Ksamp<-runif(nsamp,mean(DLM_data@Cat[x,])/rsamp,(10*mean(DLM_data@Cat[x,]))/rsamp)
  nyears<-length(DLM_data@Cat[x,])
  B<-array(NA,dim=c(nsamp,nyears))

  if(DLM_data@Cat[x,1]<(0.5*max(DLM_data@Cat[x,]))){    # Martell and Froese decision rules (makes absolutely no sense to me!)
    B[,1]<-Ksamp*runif(nsamp,0.5,0.9)
  }else{
    B[,1]<-Ksamp*runif(nsamp,0.3,0.6)
  }

  if(DLM_data@Cat[x,nyears]<(0.5*max(DLM_data@Cat[x,]))){    # Martell and Froese decision rules (makes absolutely no sense to me!)
    LB<-0.01
    UB<-0.4
  }else{
    LB<-0.3
    UB<-0.7
  }

  for(i in 2:nyears){
   B[,i]<-B[,i-1]-DLM_data@Cat[x,i-1]
   B[,i]<-B[,i]+rsamp*B[,i]*(1-B[,i]/Ksamp)
  }
  B<-B/rep(Ksamp,nyears)
  cond<-(B[,nyears]>=LB)&(B[,nyears]<=UB)
  if (sum(cond) < 1) {
	B[B[,nyears]>=UB ,nyears]<- UB
	cond<-(B[,nyears]>=LB)&(B[,nyears]<=UB)
  }	
  dep<-B[cond,nyears][1:reps]
  MSY<-rsamp[cond][1:reps]*Ksamp[cond][1:reps]/4
  Kc<-Ksamp[cond][1:reps]
  rc<-rsamp[cond][1:reps]
  TAC<-Kc*dep*rc/2

  if(sum(!is.na(TAC))<ceiling(reps/10)){ # a fudge of the original method that widens current depletion to the lowest and higest bounds to get an TAC sample
    cond<-(B[,nyears]>=0.01)&(B[,nyears]<=0.7)
    dep<-B[cond,nyears][1:reps]
    MSY<-rsamp[cond][1:reps]*Ksamp[cond][1:reps]/4
    Kc<-Ksamp[cond][1:reps]
    rc<-rsamp[cond][1:reps]
    TAC<-Kc*dep*rc/2
	
  }
  #}
  
  TACfilter(TAC)
} # end of SPMSY
class(SPMSY)<-"DLM_output"

MCD<-function(x,DLM_data,reps=100){  # Daft method to demonstrate the relative value of information of current depletion
  dependencies="DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Cat"
  depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  Bt_K<-rbeta(reps*100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K<-Bt_K[Bt_K>0.00999&Bt_K<0.99001][1:reps] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  AvC<-rlnorm(reps,log(mean(DLM_data@Cat[x,],na.rm=T)),DLM_data@CV_Cat[x])
  TAC<-AvC*2*Bt_K
  TACfilter(TAC)
}
class(MCD)<-"DLM_output"

MCD4010<-function(x,DLM_data,reps=100){  # Daft method to demonstrate the relative value of information of current depletion
  dependencies="DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Cat"
  depo<-max(0.01,min(0.99,DLM_data@Dep[x]))  # known depletion is between 1% and 99% - needed to generalise the Dick and MacCall method to extreme depletion scenarios
  Bt_K<-rbeta(reps*100,alphaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])),betaconv(depo,min(depo*DLM_data@CV_Dep[x],(1-depo)*DLM_data@CV_Dep[x])))  # CV 0.25 is the default for Dick and MacCall mu=0.4, sd =0.1
  Bt_K<-Bt_K[Bt_K>0.00999&Bt_K<0.99001][1:reps] # interval censor (0.01,0.99)  as in Dick and MacCall 2011
  AvC<-rlnorm(reps,log(mean(DLM_data@Cat[x,],na.rm=T)),DLM_data@CV_Cat[x])
  TAC<-AvC*2*Bt_K
  
  #40-10 HCR
  cond1<-Bt_K<0.4 & Bt_K>0.1
  cond2<-Bt_K<0.1
  TAC[cond1]<-TAC[cond1]*(Bt_K[cond1]-0.1)/0.3
  TAC[cond2]<-TAC[cond2]*tiny # this has to still be stochastic albeit very small
  
  TACfilter(TAC)
}
class(MCD4010)<-"DLM_output"

SPSRA<-function(x,DLM_data,reps=100){  # Surplus productin stock reduction analysis T.Carruthers - basically an SP version of DBSRA
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@Dep, DLM_data@CV_Dep, DLM_data@Cat, DLM_data@steep"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kvec<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfvec=trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0vec<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  rsamp<-getr(x,DLM_data,Mvec,Kvec,Linfvec,t0vec,hvec,maxage=DLM_data@MaxAge,r_reps=reps)
  dep<-trlnorm(reps,DLM_data@Dep[x],DLM_data@CV_Dep[x])
  Ct<-DLM_data@Cat[x,]
  Csamp<-array(rep(Ct,each=reps)*trlnorm(length(Ct)*reps,1,DLM_data@CV_Cat[x]),dim=c(reps,length(Ct)))
  Psamp<-array(trlnorm(length(Ct)*reps,1,0.1),dim=c(reps,length(Ct)))
  Ksamp<-rep(NA,reps)
  for(i in 1:reps)Ksamp[i]<-exp(optimize(SPSRAopt,log(c(mean(Csamp[i,]),1000*mean(Csamp[i,]))),dep=dep[i],r=rsamp[i],Ct=Csamp[i,],PE=Psamp[i,])$minimum)
  MSY<-Ksamp*rsamp/4
  TAC<-Ksamp*dep*rsamp/2
  TACfilter(TAC)
}
class(SPSRA)<-"DLM_output"

SPSRA_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@CAL, DLM_data@Cat, DLM_data@steep"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kvec<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfvec=trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0vec<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  rsamp<-getr(x,DLM_data,Mvec,Kvec,Linfvec,t0vec,hvec,maxage=DLM_data@MaxAge,r_reps=reps)
  Z<-MLne(x,DLM_data,Linfc=Linfvec,Kc=Kvec,ML_reps=reps,MLtype="dep")
  FM<-Z-Mvec
  FM[FM<0]<-0.01
  nyears<-length(DLM_data@Year)
  Ct1<-mean(DLM_data@Cat[x,1:3])
  Ct2<-mean(DLM_data@Cat[x,(nyears-2):nyears])
  dep<-rep(c(Ct1,Ct2),each=reps)/(1-exp(-FM[,c(1,2)]))
  dep<-dep[,2]/dep[,1]
  Ksamp<-rep(NA,reps)
  Ct<-DLM_data@Cat[x,]
  Csamp<-array(rep(Ct,each=reps)*trlnorm(length(Ct)*reps,1,DLM_data@CV_Cat[x]),dim=c(reps,length(Ct)))
  Psamp<-array(trlnorm(length(Ct)*reps,1,0.1),dim=c(reps,length(Ct)))
  for(i in 1:reps)Ksamp[i]<-exp(optimize(SPSRAopt,log(c(mean(Csamp[i,]),1000*mean(Csamp[i,]))),dep=dep[i],r=rsamp[i],Ct=Csamp[i,],PE=Psamp[i,])$minimum)
  MSY<-Ksamp*rsamp/4
  TAC<-Ksamp*dep*rsamp/2
  TACfilter(TAC)
}
class(SPSRA_ML)<-"DLM_output"

SPSRAopt<-function(lnK,dep,r,Ct,PE){
  nyears<-length(Ct)
  B<-rep(NA,nyears)
  B[1]<-exp(lnK)
  OBJ<-0
  for(y in 2:nyears){
    if((B[y-1]-Ct[y-1])<0)OBJ<-OBJ+(B[y-1]-Ct[y-1])^2
    B[y]<-max(0.01,B[y-1]-Ct[y-1])
    B[y]<-B[y]+r*B[y]*(1-B[y]/B[1])*PE[y]
  }
  return(OBJ+((B[nyears]/B[1])-dep)^2)
}


YPR<-function(x,DLM_data,reps=100){   # Yield per recruit analysis F01 - Meaghan Bryan
  #for(x in 1:10){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@Abun, DLM_data@CV_Abun, DLM_data@wla, DLM_data@wlb"
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  Mdb<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  LFS<-trlnorm(reps,DLM_data@LFS[x],DLM_data@CV_LFS[x])
  a<-DLM_data@wla[x]
  b<-DLM_data@wlb[x]
  Ac<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun[x])
  FMSY<-YPRopt(Linfc,Kc,t0c,Mdb,a,b,LFS,DLM_data@MaxAge,reps)
  TAC<-Ac*FMSY
  TACfilter(TAC)
  #}
}   # end of YPR
class(YPR)<-"DLM_output"

YPR_CC<-function(x,DLM_data,reps=100,Fmin=0.005){
  #for(x in 1:16){
    dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@wlb, DLM_data@CAA, DLM_data@Cat"
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  LFS<-trlnorm(reps,DLM_data@LFS[x],DLM_data@CV_LFS[x])
  a<-DLM_data@wla[x]
  b<-DLM_data@wlb[x]
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  
  Mdb<-trlnorm(reps*10,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Zdb<-CC(x,DLM_data,reps=reps*10)
  Fdb<-Zdb-Mdb
  ind<-(1:(reps*10))[Fdb>Fmin][1:reps]
  
  Fdb<-Fdb[ind]
  Mdb<-Mdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0 ) {
    Mdb[is.na(ind)] <- trlnorm(SM,DLM_data@Mort[x],DLM_data@CV_Mort[x])
    Fdb[is.na(ind)] <- Fmin
  }	
  
  Ac<-Cc/(1-exp(-Fdb))
  FMSY<-YPRopt(Linfc,Kc,t0c,Mdb,a,b,LFS,DLM_data@MaxAge,reps)
  TAC<-Ac*FMSY
 # }
  TACfilter(TAC)
}
class(YPR_CC)<-"DLM_output"

YPR_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@wlb, DLM_data@CAL, DLM_data@Cat"
  Linfc<-trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  LFS<-trlnorm(reps,DLM_data@LFS[x],DLM_data@CV_LFS[x])
  a<-DLM_data@wla[x]
  b<-DLM_data@wlb[x]
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=reps,MLtype="F")
  FM<-Z-Mdb
  Ac<-Cc/(1-exp(-FM))
  FMSY<-YPRopt(Linfc,Kc,t0c,Mdb,a,b,LFS,DLM_data@MaxAge,reps)
  TAC<-Ac*FMSY
  TACfilter(TAC)
}
class(YPR_ML)<-"DLM_output"

Fdem<-function(x,DLM_data,reps=100){   # Demographic FMSY estimate (FMSY=r/2)
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@wlb, DLM_data@Abun, DLM_data@CV_Abun, DLM_data@steep, DLM_data@CV_steep"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfc=trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  Ac<-trlnorm(reps,DLM_data@Abun[x],DLM_data@CV_Abun[x])
  FMSY<-getr(x,DLM_data,Mvec,Kc,Linfc,t0c,hvec,maxage=DLM_data@MaxAge,r_reps=reps)/2
  TAC<-FMSY*Ac
  TACfilter(TAC)
}
class(Fdem)<-"DLM_output"

Fdem_CC<-function(x,DLM_data,reps=100,Fmin=0.005){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@wlb, DLM_data@CAA, DLM_data@steep, DLM_data@CV_steep"
  Mvec<-trlnorm(reps*10,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfc=trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Zdb<-CC(x,DLM_data,reps=reps*10)
  Fdb<-Zdb-Mvec
  ind<-(1:(reps*10))[Fdb>Fmin][1:reps]
  
  Fdb<-Fdb[ind]
  SM <- sum(is.na(ind))
  if (SM > 0 ) {
    Fdb[is.na(ind)] <- Fmin
  }	

  Ac<-Cc/(1-exp(-Fdb))
  FMSY<-getr(x,DLM_data,Mvec,Kc,Linfc,t0c,hvec,maxage=DLM_data@MaxAge,r_reps=reps)/2
  TAC<-FMSY*Ac
  
  TACfilter(TAC)
}
class(Fdem_CC)<-"DLM_output"

# Catch curve estimate of recent F (naive) ========================================================================================
CC<-function(x,DLM_data,reps=100){
  ny<-dim(DLM_data@CAA)[2]
  CAA<-apply(DLM_data@CAA[x,max(ny-2,1):ny,],2,sum) # takes last two years as the sample (or last year if there is only one
  maxageobs<-length(CAA)
  AFS<-which.max(CAA)
  AFS[AFS>(maxageobs-3)]<-maxageobs-3   # provides at least three datapoints
  
  nS<-ceiling(sum(CAA)/2)
  y <-log(CAA[AFS:maxageobs]/sum(CAA[AFS:maxageobs],na.rm=T))
  xc<-1:length(y)
  y[y=='-Inf']<-NA
  mod <- lm(y~xc)
  chk <- sum(is.na(coef(mod))) # check if model failed
  if (chk) {
    return(NA)
  } else {
    coefs <-summary(mod,weights=CAA[AFS:maxageobs])$coefficients[2,1:2]
	coefs[is.nan(coefs)] <- tiny
   return(-rnorm(reps,coefs[1],coefs[2]))
  } 
}
# class(CC)<-"DLM_output"

Fdem_ML<-function(x,DLM_data,reps=100){
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@wlb, DLM_data@CAL, DLM_data@steep, DLM_data@CV_steep"
  Mvec<-trlnorm(reps,DLM_data@Mort[x],DLM_data@CV_Mort[x])
  Kc<-trlnorm(reps,DLM_data@vbK[x],DLM_data@CV_vbK[x])
  Linfc=trlnorm(reps,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
  t0c<--trlnorm(reps,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
  hvec<-trlnorm(reps,DLM_data@steep[x],DLM_data@CV_steep[x])
  MuC<-DLM_data@Cat[x,length(DLM_data@Cat[x,])]
  Cc<-trlnorm(reps,MuC,DLM_data@CV_Cat[x])
  Z<-MLne(x,DLM_data,Linfc=Linfc,Kc=Kc,ML_reps=reps,MLtype="F")
  FM<-Z-Mvec
  Ac<-Cc/(1-exp(-FM))
  FMSY<-getr(x,DLM_data,Mvec,Kc,Linfc,t0c,hvec,maxage=DLM_data@MaxAge,r_reps=reps)/2
  TAC<-FMSY*Ac
  TACfilter(TAC)
}
class(Fdem_ML)<-"DLM_output"

CompSRA<-function(x,DLM_data,reps=100){    # optimize for fixed F to get you to current depletion C/Fcur = abundance
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@CV_wla, DLM_data@wlb, DLM_data@CV_wlb, DLM_data@L50, DLM_data@CV_L50, DLM_data@CAA, DLM_data@steep, DLM_data@CV_steep, DLM_data@LFS, DLM_data@CV_LFS, DLM_data@LFC, DLM_data@CV_LFC, DLM_data@Cat"
  maxage<-DLM_data@MaxAge
  TAC<-rep(NA,reps)
  for(i in 1:reps){
    Mc<-trlnorm(1,DLM_data@Mort[x],DLM_data@CV_Mort)
    hc<-trlnorm(1,DLM_data@steep[x],DLM_data@CV_steep[x])
    Linfc<-trlnorm(1,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
    Kc<-trlnorm(1,DLM_data@vbK[x],DLM_data@CV_vbK[x])
    t0c<--trlnorm(1,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
    LFSc<-trlnorm(1,DLM_data@LFS[x],DLM_data@CV_LFS[x])
    LFCc<-trlnorm(1,DLM_data@LFC[x],DLM_data@CV_LFC[x])
    AMc<-trlnorm(1,iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x]),DLM_data@CV_L50[x])
    ac<-trlnorm(1,DLM_data@wla[x],DLM_data@CV_wla[x])
    bc<-trlnorm(1,DLM_data@wlb[x],DLM_data@CV_wlb[x])
    Catch<-DLM_data@Cat[x,]
    ny<-length(Catch)
    nyCAA<-dim(DLM_data@CAA)[2]
    CAA<-DLM_data@CAA[x,max(nyCAA-2,1):nyCAA,] # takes last two years as the sample (or last year if there is only one

    Nac<-exp(-Mc*((1:maxage)-1)) # put a rough range on estimate of R0 assuming a mean harvest rate of 10%
    Lac<-Linfc*(1-exp(-Kc*((1:maxage)-t0c)))
    Wac<-ac*Lac^bc
    AFC<-log(1-min(0.99,LFCc/Linfc))/-Kc+t0c
    AFS<-log(1-min(0.99,LFSc/Linfc))/-Kc+t0c
    if (AFC >= 0.7 * maxage) AFC <- 0.7 * maxage
    if (AFS >= 0.9 * maxage) AFS <- 0.9 * maxage
    KES<-max(2,ceiling(mean(c(AFC,AFS))))
    pred<-Nac*Wac
    pred[1:(KES-1)]<-0
    pred<-pred/sum(pred)
    pred<-((mean(Catch)/0.1)*pred/Wac)/exp(-(1:maxage)*Mc)
    pred<-pred[pred>0]
    R0range<-c(mean(pred)/1000,mean(pred)*1000)

    fit<-optimize(SRAfunc,log(R0range),Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA)
    Ac<-SRAfunc(fit$minimum,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA,opt=2)
    fit2<-optimize(SRAFMSY,log(c(0.0001,3)),Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc)
    # FMSY<-SRAFMSY(fit2$minimum,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,opt=F)
	FMSY <- exp(fit2$minimum)
    if ((FMSY / Mc) > 3) FMSY <- 3 * Mc
    TAC[i]<-Ac*FMSY
	# message(i, " of ", reps)
    # flush.console()
  }
  TACfilter(TAC)
}
class(CompSRA)<-"DLM_output"

CompSRA4010<-function(x,DLM_data,reps=100){    # optimize for fixed F to get you to current depletion C/Fcur = abundance
 # for (x in 1:nsim) {
  dependencies="DLM_data@Mort, DLM_data@CV_Mort, DLM_data@vbK, DLM_data@CV_vbK, DLM_data@vbLinf, DLM_data@CV_vbLinf, DLM_data@vbt0, DLM_data@CV_vbt0, DLM_data@MaxAge, DLM_data@wla, DLM_data@CV_wla, DLM_data@wlb, DLM_data@CV_wlb, DLM_data@L50, DLM_data@CV_L50, DLM_data@CAA, DLM_data@steep, DLM_data@CV_steep, DLM_data@LFS, DLM_data@CV_LFS, DLM_data@LFC, DLM_data@CV_LFC, DLM_data@Cat"
  maxage<-DLM_data@MaxAge
  TAC<-Bt_K<-rep(NA,reps)
  for(i in 1:reps){
    Mc<-trlnorm(1,DLM_data@Mort[x],DLM_data@CV_Mort)
    hc<-trlnorm(1,DLM_data@steep[x],DLM_data@CV_steep[x])
    Linfc<-trlnorm(1,DLM_data@vbLinf[x],DLM_data@CV_vbLinf[x])
    Kc<-trlnorm(1,DLM_data@vbK[x],DLM_data@CV_vbK[x])
    t0c<--trlnorm(1,-DLM_data@vbt0[x],DLM_data@CV_vbt0[x])
    LFSc<-trlnorm(1,DLM_data@LFS[x],DLM_data@CV_LFS[x])
    LFCc<-trlnorm(1,DLM_data@LFC[x],DLM_data@CV_LFC[x])
    AMc<-trlnorm(1,iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x]),DLM_data@CV_L50[x])
    ac<-trlnorm(1,DLM_data@wla[x],DLM_data@CV_wla[x])
    bc<-trlnorm(1,DLM_data@wlb[x],DLM_data@CV_wlb[x])
    Catch<-DLM_data@Cat[x,]
    ny<-length(Catch)
    nyCAA<-dim(DLM_data@CAA)[2]
    CAA<-DLM_data@CAA[x,max(nyCAA-2,1):nyCAA,] # takes last two years as the sample (or last year if there is only one
    
    Nac<-exp(-Mc*((1:maxage)-1)) # put a rough range on estimate of R0 assuming a mean harvest rate of 10%
    Lac<-Linfc*(1-exp(-Kc*((1:maxage)-t0c)))
    Wac<-ac*Lac^bc
	
	AFC<-log(1-min(0.99,LFCc/Linfc))/-Kc+t0c
	AFS<-log(1-min(0.99,LFSc/Linfc))/-Kc+t0c
	if (AFC >= 0.7 * maxage) AFC <- 0.7 * maxage
	if (AFS >= 0.9 * maxage) AFS <- 0.9 * maxage
	
    KES<-max(2,ceiling(mean(c(AFC,AFS))))
    pred<-Nac*Wac
    pred[1:(KES-1)]<-0
    pred<-pred/sum(pred)
    pred<-((mean(Catch)/0.1)*pred/Wac)/exp(-(1:maxage)*Mc)
    pred<-pred[pred>0]
    R0range<-c(mean(pred)/1000,mean(pred)*1000)
    
    fit<-optimize(SRAfunc,log(R0range),Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA)
    Ac<-SRAfunc(fit$minimum,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA,opt=2)
    Bt_K[i]<-SRAfunc(fit$minimum,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA,opt=3)
    fit2<-optimize(SRAFMSY,log(c(0.0001,3)),Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc)
    FMSY<-SRAFMSY(fit2$minimum,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,opt=F)
    if ((FMSY / Mc) > 3) FMSY <- 3 * Mc
    TAC[i]<-Ac*FMSY
  }   
  # 40-10 rule
  cond1<-Bt_K<0.4 & Bt_K>0.1
  cond2<-Bt_K<0.1
  TAC[cond1]<-TAC[cond1]*(Bt_K[cond1]-0.1)/0.3
  TAC[cond2]<-TAC[cond2]*tiny # this has to still be stochastic albeit very small

  TACfilter(TAC)
  
  # message(x, " of ", nsim)
  # flush.console()
  # }
}
class(CompSRA4010)<-"DLM_output"

# options(warn=2)
# options(warn=1)

SRAfunc<-function(lnR0c,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,Catch,CAA,opt=1){

  ny<-length(Catch)
  AFC<-log(1-min(0.99,LFCc/Linfc))/-Kc+t0c
  AFS<-log(1-min(0.99,LFSc/Linfc))/-Kc+t0c
  if (AFC >= 0.7 * maxage) AFC <- 0.7 * maxage
  if (AFS >= 0.9 * maxage) AFS <- 0.9 * maxage
  KES<-max(2,ceiling(mean(c(AFC,AFS))))
  vul<-rep(1,maxage)
  vul[1:(KES-1)]<-0
  Mac<-rep(1,maxage)
  Mac[1:max(1,floor(AMc))]<-0
  Lac<-Linfc*(1-exp(-Kc*((1:maxage)-t0c)))
  Wac<-ac*Lac^bc
  R0c<-exp(lnR0c)
  N<-exp(-Mc*((1:maxage)-1))*R0c
  SSN<-Mac*N                                 # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  SSB<-SSN*Wac                               # Calculate spawning stock biomass

  B0<-sum(Biomass)
  SSB0<-sum(SSB)
  SSN0<-SSN
  SSBpR<-sum(SSB)/R0c                              # Calculate spawning stock biomass per recruit
  SSNpR<-SSN/R0c

  CN<-array(NA,dim=c(ny,maxage))
  HR<-rep(0,maxage)
  pen<-0
  for(y in 1:ny){
    # set up some indices for indexed calculation
    VB<-Biomass[KES:maxage]*exp(-Mc)
    CB<-Catch[y]*VB/sum(VB)
    testHR<-CB[1]/VB[1]
    if(testHR>0.8)pen<-pen+(testHR-0.8)^2
    HR[KES:maxage]<-min(testHR,0.8)
    FMc<--log(1-HR)                                        # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc<-FMc+Mc

    CN[y,]<-N * (1-exp(-Zc))*(FMc/Zc)
    N[2:maxage]<-N[1:(maxage-1)]*exp(-Zc[1:(maxage-1)])         # Total mortality
    N[1]<-(0.8*R0c*hc*sum(SSB))/(0.2*SSBpR*R0c*(1-hc)+(hc-0.2)*sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    #print(N[1])
    Biomass<-N*Wac
    SSN<-N*Mac
    SSB<-SSN*Wac

  } # end of year
  
  CN[CN<0] <- 0 # stop any negative catches
  syear<-ny-dim(CAA)[1]+1
  pred<-CN[syear:ny,]
  pred<-pred/array(apply(pred,1,sum),dim=c(dim(CAA)[1],maxage))
 
  fobj<-pen-sum(log(pred+tiny)*CAA,na.rm=T)
  if(opt==1){return(fobj)
  }else if(opt==2){return(sum(Biomass))
  }else if(opt==3){sum(SSB)/sum(SSB0)
  }
  #CBc<-sum(CB)
}

SRAFMSY<-function(lnFMc,Mc,hc,maxage,LFSc,LFCc,Linfc,Kc,t0c,AMc,ac,bc,opt=T){

  FMc<-exp(lnFMc)
  ny<-100
  AFC<-log(1-min(0.99,LFCc/Linfc))/-Kc+t0c
  AFS<-log(1-min(0.99,LFSc/Linfc))/-Kc+t0c
  if (AFC >= 0.7 * maxage) AFC <- 0.7 * maxage
  if (AFS >= 0.9 * maxage) AFS <- 0.9 * maxage
  
  KES<-max(2,ceiling(mean(c(AFC,AFS))))
  vul<-rep(1,maxage)
  vul[1:(KES-1)]<-0
  Mac<-rep(1,maxage)
  Mac[1:max(1,floor(AMc))]<-0
  Lac<-Linfc*(1-exp(-Kc*((1:maxage)-t0c)))
  Wac<-ac*Lac^bc
  R0c<-1
  N<-exp(-Mc*((1:maxage)-1))*R0c
  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  SSB<-SSN*Wac                               # Calculate spawning stock biomass

  B0<-sum(Biomass)
  SSB0<-sum(SSB)
  SSN0<-SSN
  SSBpR<-sum(SSB)/R0c                              # Calculate spawning stock biomass per recruit
  SSNpR<-SSN/R0c

  N<-N/2
  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  SSB<-SSN*Wac

  for(y in 1:ny){
    # set up some indices for indexed calculation
                                           # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc<-FMc*vul+Mc
    CN<-N*(1-exp(-Zc))*(FMc/Zc)
    CB<-CN*Wac
    Biomass<-N*Wac
    N[2:maxage]<-N[1:(maxage-1)]*exp(-Zc[1:(maxage-1)])         # Total mortality
    N[1]<-(0.8*R0c*hc*sum(SSB))/(0.2*SSBpR*R0c*(1-hc)+(hc-0.2)*sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    #print(N[1])
    SSN<-N*Mac
    SSB<-SSN*Wac

  } # end of year
  
  if(opt){return(-sum(CB))
  }else{return(FMc)
  }
}


# Yield per recruit estimate of FMSY Meaghan Bryan 2013 ==================================
YPRopt=function(Linfc,Kc,t0c,Mdb,a,b,LFS,maxage,reps=100) {

  nf<-200
  frates<-seq(0,3,length.out=nf)

  Winf=a*Linfc^b
	rat<-LFS/Linfc
	rat[rat>0.8]<-0.8     # need to robustify this for occasionally very high samples of LFS
  tc=log(1-rat)/-Kc+t0c
	tc=round(tc,0)
  tc[tc<1]<-1
  tc[tc>maxage]<-maxage

	vul<-array(0,dim=c(reps,maxage))
	mat<-array(0,dim=c(reps,maxage))
	lx<-array(NA,dim=c(reps,maxage))
	lxo<-array(NA,dim=c(reps,maxage))

  ypr<-array(NA,dim=c(reps,nf))
	sbpr<-array(NA,dim=c(reps,nf))
	sbpr.ratio<-array(NA,dim=c(reps,nf))
	sbpr.dif<-array(NA,dim=c(reps,nf))

  f.max<-array(NA,dim=c(reps,maxage))

	#average weight at age - follow von Bertalanffy growth
	age<-array(rep(1:maxage,each=reps),dim=c(reps,maxage))
  la<-Linfc*(1-exp(-Kc*((age-t0c))))
  wa<-a*la^b

	#vulnerability schedule - assumes knife-edge vulnerability, where all individuals age tc to maxage are fully vulnerbale
	#all individulas less than age tc are not vulnerable
	for(i in 1:reps){
    if(tc[i]>0)vul[i,tc[i]:maxage]<-1
  	if(tc[i]>1)mat[i,max(1,tc[i]-1):maxage]<-1
  }

  lx[,1]<-1
  lxo[,1]<-1
  for(k in 1:nf){
    for(i in 2:maxage){
      lx[,i]=lx[,i-1]*exp(-(Mdb+vul[,i-1]*frates[k]))
      lxo[,i]=lx[,i]*exp(-Mdb)
    }
    phi_vb=apply(lx*wa*vul,1,sum)
		sbpro=apply(lxo*wa*mat,1,sum)

		ypr[,k]=(1-exp(-frates[k]))*phi_vb
		sbpr[,k]=apply(lx*wa*mat,1,sum)
		sbpr.ratio[,k]=sbpr[,k]/sbpro
		sbpr.dif[,k]=abs(sbpr.ratio[,k]-0.3)	#hard code comparison ratio
  }

  #frates[apply(ypr,1,which.max)]    Fmaxypr

  # More code that derived F0.1 in 'per recruit analysis.R' (Meaghan Bryan)
  slope.origin=(ypr[,2]-ypr[,1])/(frates[2]-frates[1])
	slope.10=round(0.1*slope.origin,2)

	slope=array(NA,dim=dim(ypr))#vector(length=length(ypr))
	slope[,1]=slope.origin
	for(i in 3:ncol(ypr))
	{
		slope[,i-1]=round((ypr[,i]-ypr[,i-1])/(frates[i]-frates[i-1]),2)
	}
	dif=abs(slope-slope.10)
  dif[is.na(dif)]<-10e10
	frates[apply(dif,1,which.min)]#frates[which.min(dif)]
}


MLne<-function(x,DLM_data,Linfc,Kc,ML_reps=100,MLtype="F"){
  year<-1:dim(DLM_data@CAL)[2]
  nlbin<-ncol(DLM_data@CAL[x,,])
  nlyr<-nrow(DLM_data@CAL[x,,])
  mlbin<-(DLM_data@CAL_bins[1:nlbin]+DLM_data@CAL_bins[2:(nlbin+1)])/2
  nbreaks<-1
  Z<-matrix(NA,nrow=ML_reps,ncol=nbreaks+1)
  Z2<-rep(NA,nrow=ML_reps)
  for(i in 1:ML_reps){
    mlen<-rep(NA,length(year))
     ss<-ceiling(apply(DLM_data@CAL[x,,],1,sum)/2)
    if(MLtype=="dep"){
      for(y in 1:length(year)) {
	    if (sum(DLM_data@CAL[x,y,] > 0) > 0.25 * length(DLM_data@CAL[x,y,])) {
	      mlen[y]<-mean(sample(mlbin,ceiling(sum(DLM_data@CAL[x,y,])/2),replace=T,prob=DLM_data@CAL[x,y,]), na.rm=TRUE)
		}  
	  }
      Z[i,]<-bhnoneq(year=year,mlen=mlen,ss=ss,K=Kc[i],Linf=Linfc[i],Lc=DLM_data@LFS[x],nbreaks=nbreaks,
           styrs=ceiling(length(year)*((1:nbreaks)/(nbreaks+1))),stZ=rep(0.05,nbreaks+1),stsigma=20,graph=F)
    }else{
      
      ind<-(which.min(((DLM_data@CAL_bins-DLM_data@LFS[x])^2)^0.5)-1):(length(DLM_data@CAL_bins)-1)
      for(y in 1:length(year)) {
	    if (sum(DLM_data@CAL[x,y,] > 0) > 0.25 * length(DLM_data@CAL[x,y,])) {
		  mlen[y]<-mean(sample(mlbin[ind],ceiling(sum(DLM_data@CAL[x,y,ind])/2),replace=T,prob=DLM_data@CAL[x,y,ind]), na.rm=TRUE)
		}
      }		
      mlen<-mean(mlen[(length(mlen)-2):length(mlen)], na.rm=TRUE)
      Z2<-bheq(K=Kc[i],Linf=Linfc[i],Lc=DLM_data@LFS[x],Lbar=mlen)
    }
  }
  if(MLtype=="F")return(Z2)
  if(MLtype=="dep")return(Z)
}

bheq<-function(K,Linf,Lc,Lbar){
  K*(Linf-Lbar)/(Lbar-Lc)
}

bhnoneq<-function (year = NULL, mlen = NULL, ss = NULL, K = NULL, Linf = NULL,
    Lc = NULL, nbreaks = NULL, styrs = NULL, stZ = NULL, stsigma = NULL,
    graph = TRUE)  {
    if (is.null(mlen))
        stop("mean length vector does not exist")
    if (is.null(year))
        stop("year vector does not exist")
    if (is.null(ss))
        stop("numbers vector does not exist")
    if (!is.numeric(mlen))
        stop("vector is not numeric")
    if (is.null(stZ))
        stop("Initial Z vector does not exist")
    if (is.null(stsigma))
        stop("Initial sigma value does not exist")
    if (is.null(K))
        stop("K not specified")
    if (is.null(Linf))
        stop("Linf not specified")
    if (is.null(Lc))
        stop("Lc not specified")
    if (is.null(nbreaks))
        stop("Number of mortality breaks not specified")
    if (is.null(styrs))
        stop("Starting guesses for years of mortality breaks not specified ")
    if (length(mlen) != length(year))
        stop("vectors have different lengths")
    gyr <- styrs
    xx <- as.data.frame(cbind(year, mlen, ss))
    fyr <- min(xx$year)
    lyr <- max(xx$year)
    names(xx) <- c("year", "mlen", "m")
    nyr <- length(xx[!is.na(xx[, 2]), 2])
    count <- length(xx[, 1])
    mdata <- xx
    ggyr <- gyr - fyr + 1
    tm <- array(0, dim = c(nbreaks, 1))
    if (length(stZ) == nbreaks + 1)
        parms <- c(stZ, ggyr, stsigma)
    if (length(stZ) < nbreaks + 1)
        stop("The number of stZ values does not equal nbreak+1")
    if (length(stZ) > nbreaks + 1)
        stop("Too many stZ values for match to nbreak+1")
    Lpred <- NULL
    results <- NULL
    Zest <- function(y) {
        Z <- y[1:as.numeric(nbreaks + 1)]
        sigma <- y[length(y)]
        tm[1:nbreaks, 1] <- y[as.numeric(nbreaks + 2):as.numeric(length(y) -
            1)]
        dy <- array(0, dim = c(nbreaks, count))
        for (i in 1:nbreaks) {
            for (j in 1:count) {
                dy[i, j] <- ifelse(tm[i, 1] >= j, 0, j - tm[i,
                  1])
            }
        }
        if (nbreaks > 1) {
            for (i in 1:as.numeric(nbreaks - 1)) {
                for (j in 1:count) {
                  if (j > round(tm[i + 1, 1], 0)) {
                    dy[i, j] <- dy[i, j - 1]
                  }
                }
            }
        }
        a <- array(0, dim = c(nbreaks + 1, count))
        s <- array(0, dim = c(nbreaks + 1, count))
        r <- array(0, dim = c(nbreaks + 1, count))
        denom <- rep(0, count)
        numersum <- rep(0, count)
        numerator <- rep(0, count)
        for (m in 1:count) {
            a[1, m] <- 1
            r[1, m] <- 1
            for (i in 1:as.numeric(nbreaks + 1)) {
                a[i, m] <- 1
                r[i, m] <- 1
                if (i < as.numeric(nbreaks + 1)) {
                  s[i, m] <- 1 - exp(-(Z[nbreaks + 2 - i] + K) *
                    dy[nbreaks + 1 - i, m])
                }
                if (i == as.numeric(nbreaks + 1)) {
                  s[i, m] <- 1
                }
                for (j in 1:as.numeric(i - 1)) {
                  if (i > 1) {
                    a[i, m] <- a[i, m] * exp(-Z[nbreaks + 2 -
                      j] * dy[nbreaks + 1 - j, m])
                    r[i, m] <- r[i, m] * exp(-(Z[nbreaks + 2 -
                      j] + K) * dy[nbreaks + 1 - j, m])
                  }
                }
                if (i <= nbreaks) {
                  denom[m] <- denom[m] + a[i, m] * ((1 - exp(-Z[nbreaks +
                    2 - i] * dy[nbreaks + 1 - i, m]))/Z[nbreaks +
                    2 - i])
                }
                if (i == as.numeric(nbreaks + 1)) {
                  denom[m] <- denom[m] + a[i, m]/Z[nbreaks +
                    2 - i]
                }
                numersum[m] <- numersum[m] + (-((1 - Lc/Linf) *
                  r[i, m] * s[i, m])/(Z[nbreaks + 2 - i] + K))
            }
        }
        numerator <- Linf * (denom + numersum)
        Lpred <<- numerator/denom
        LL <- -nyr * log(sigma) - sum((mdata[, 3]/(2 * sigma^2)) *
            (mdata[, 2] - Lpred)^2, na.rm = T)
        LL * -1
    }
    results <- optim(parms, Zest, gr = NULL, control = list(maxit = 1e+06,
        abstol = 1e-07), hessian = TRUE)
    return(results$par[1:(nbreaks + 1)])

}

getdep<-function(lnFF,targ,Md,Linfd,Kd,t0d,AFSd,ad,bd,maxage,opt){

   FF<-exp(lnFF)
   Z<-rep(Md,maxage)
   Z[1]<-0
   Z[AFSd:maxage]<-Z[AFSd:maxage]+FF
   for(a in 2:maxage)Z[a]<-Z[a-1]+Z[a]
   Nd<-exp(-Z)
   Nobs<-Nd
   Nobs[1:max(1,(AFSd-1))]<-0
   Ld<-Linfd*(1-exp(-Kd*(((1:maxage)-0.5)-t0d)))

   if(opt)return(((sum(Nobs*Ld)/sum(Nobs))-targ)^2)
   if(!opt){
     Nd0<-exp(-Md*((1:maxage)-0.5))
     Wd<-ad*Ld^bd
     return((sum(Nd*Wd)/sum(Nd))/(sum(Nd0*Wd)/sum(Nd0)))
   }
}

getr <- function(x,DLM_data,Mvec,Kvec,Linfvec,t0vec,hvec,maxage,r_reps=100){
  r<-rep(NA,r_reps)
  for(i in 1:r_reps){
    log.r=log(0.3)

    opt=optimize(demofn,lower=log(0.0001),upper=log(1.4),
    M=Mvec[i],
    amat=iVB(DLM_data@vbt0[x],DLM_data@vbK[x],DLM_data@vbLinf[x],DLM_data@L50[x]),
    sigma=0.2,
    K=Kvec[i],
    Linf=Kvec[i],
    to=t0vec[i],
    hR=hvec[i],
    maxage=maxage,
    a=DLM_data@wla[x],
    b=DLM_data@wlb[x])
    #demographic2(opt$minimum,M[x],ageM[x],0.2,K[x],Linf,t0,steepness[x],maxage,a,b)$r
    r[i]<-exp(opt$minimum)
  }
  r
}

iVB<-function(t0,K,Linf,L)((-log(1-L/Linf))/K+t0) # Inverse Von-B

EDCAC<-function (x, DLM_data, reps = 100) # extended depletion-corrected average catch (Harford and Carruthers 2015)
{
  dependencies = "DLM_data@AvC, DLM_data@t, DLM_data@Mort, DLM_data@CV_Mort, DLM_data@Dt, DLM_data@CV_Dt, DLM_data@BMSY_B0, DLM_data@CV_BMSY_B0"
  C_tot <- DLM_data@AvC[x] * DLM_data@t[x]
  Mdb <- trlnorm(reps, DLM_data@Mort[x], DLM_data@CV_Mort[x])
  FMSY_M <- trlnorm(reps, DLM_data@FMSY_M[x], DLM_data@CV_FMSY_M[x])
  Bt_K <- trlnorm(reps, DLM_data@Dt[x], DLM_data@CV_Dt[x])
  BMSY_K <- rbeta(reps, alphaconv(DLM_data@BMSY_B0[x], DLM_data@CV_BMSY_B0[x]), 
                  betaconv(DLM_data@BMSY_B0[x], DLM_data@CV_BMSY_B0[x]))
  dcac<-C_tot/(DLM_data@t[x] + ((1 - Bt_K)/(BMSY_K * FMSY_M * Mdb)))
  TAC<-dcac*Bt_K/BMSY_K
  TACfilter(TAC)
}
class(EDCAC)<-"DLM_output"

AvC<-function(x,DLM_data,reps=100)rlnorm(reps,log(mean(DLM_data@Cat[x,],na.rm=T)),0.2)
class(AvC)<-"DLM_output"


LBSPR_ItTAC <- function(x, DLM_data, yrsmth=1,reps=reps) {
 dependencies="DLM_data@CAL, DLM_data@CAL_bins, DLM_data@vbLinf, DLM_data@vbK, DLM_data@Mort, LM_data@vbK, DLM_data@L50, DLM_data@L95, DLM_data@wlb" 
  
  MiscList <- LBSPR(x, DLM_data, yrsmth=yrsmth,reps=reps)

  XX <- 1:4 
  YY <- MiscList[[2]][(length(MiscList[[2]]) - (max(XX)-1)):length(MiscList[[2]])]
  
  EstSPR <- YY[length(YY)]
  
  TgSPR <- 0.4
  h <- DLM_data@steep[x]
  SPRLim <- -(2*(h-1))/(3*h+1) # SPR that results in 0.5 R0
  
  phi1 <- 6
  phi2 <- 1
  
  MaxDw <- -0.3
  MaxUp <- 0.3
  
  minSlope <- 0.01
  
  Slope <- coef(lm(YY~XX))[2]  
  # if (abs(Slope) < minSlope) Slope <- 0 
  Dist <- EstSPR - TgSPR 
  
  # Control Rule #
  Mod <- 0 
  Buff <- 0.1
  Buffer <- c(TgSPR - Buff,  TgSPR + Buff)
  inBuff <- FALSE
  belowTG <- FALSE 
  aboveTG <- FALSE
  slopeUp <- FALSE
  slopeDw <- FALSE 
  belowLim <- FALSE
  if (Dist < 0) belowTG <- TRUE 
  if (Dist > 0) aboveTG <- TRUE 
  if (EstSPR > min(Buffer) & EstSPR < max(Buffer)) inBuff <- TRUE
  if (Slope <= 0) slopeDw <- TRUE
  if (Slope > 0) slopeUp <- TRUE
  if (EstSPR < SPRLim) belowLim <- TRUE
   
  # If within buffer zone - only slope
  if (inBuff) Mod <- phi1 * Slope
  if (slopeUp & aboveTG) Mod <- phi1 * Slope +  phi2 * Dist
  if (slopeUp & belowTG) Mod <- phi1 * Slope 
  
  if (slopeDw & aboveTG) Mod <- phi1 * Slope 
  if (slopeDw & belowTG) Mod <- phi1 * Slope +  phi2 * Dist
  
  if (belowLim) Mod <- MaxDw
  
  Mod[Mod > MaxUp] <- MaxUp
  Mod[Mod < MaxDw] <- MaxDw
  Mod <- Mod + 1 
 
  TAC <- DLM_data@MPrec[x] * Mod
  TAC <- TACfilter(TAC)
 
  Out <- list()
  Out[[1]] <- TAC 
  Out[[2]] <- MiscList
 
  return(Out) 
}
class(LBSPR_ItTAC)<-"DLM_output"





# DLMtool Miscellaneous Functions 
# January 2016
# Tom Carruthers UBC (t.carruthers@fisheries.ubc.ca)
# Adrian Hordyk (a.hordyk@murdoch.edu.au)

# Collection of miscellaneous functions.
# All functions have accompanying help files.

# Generic class finder
avail<-function(classy){
  chk <- "package:DLMtooldev" %in% search() 
  if (chk) { # development version
    return(unique(c(ls('package:DLMtooldev')[unlist(lapply(ls('package:DLMtooldev'),
      getclass,classy=classy))], ls(envir=.GlobalEnv)[unlist(
	  lapply(ls(envir=.GlobalEnv),getclass,classy=classy))]))) 
  } else {	  
    return(unique(c(ls('package:DLMtool')[unlist(lapply(ls('package:DLMtool'),
    getclass,classy=classy))], ls(envir=.GlobalEnv)[unlist(
	lapply(ls(envir=.GlobalEnv),getclass,classy=classy))])))
  }	
}

# A function that finds all methods in the environment and searches the function
# text for slots in the DLM data object
Required <- function(funcs=NA){
  if(is.na(funcs[1]))funcs<-c(avail("DLM_output"),avail("DLM_input"))
  slots<-slotNames('DLM_data')
  slotnams<-paste("DLM_data@",slotNames('DLM_data'),sep="")
  repp<-rep("",length(funcs))

  for(i in 1:length(funcs)){
    temp<-format(match.fun(funcs[i]))
    temp<-paste(temp[1:(length(temp))],collapse=" ")
    rec<-""
    for(j in 1:length(slotnams))if(grepl(slotnams[j],temp))rec<-c(rec,slots[j])
    if(length(rec)>1)repp[i]<-paste(rec[2:length(rec)],collapse=", ")
  }
  cbind(funcs,repp,deparse.level=0)
}

# A way of locating where the package was installed so you can find example 
# data files and code etc.
DLMDataDir<-function(stock=NA){
  chk <- "package:DLMtooldev" %in% search() 
  if (chk) { # dev version 
    if(is.na(stock)){
      return(paste(searchpaths()[match("package:DLMtooldev",search())],"/",sep=""))
    }else{
      return(paste(searchpaths()[match("package:DLMtooldev",search())],"/",stock,".csv",sep=""))
    }
  } else {
    if(is.na(stock)){
      return(paste(searchpaths()[match("package:DLMtool",search())],"/",sep=""))
    }else{
      return(paste(searchpaths()[match("package:DLMtool",search())],"/",stock,".csv",sep=""))
    }
  }  
}

# What MPs may be run (best case scenario) for various data-availability 
#  scenarios?
Fease<-function(feaseobj,outy="table"){
  
  if(class(feaseobj)!="DLM_fease")stop("Incorrect format: you need an object of class DLM_fease")
  
  sloty<-c("Cat","Ind","AvC","Dt","Rec","CAA","CAL","Mort","L50","L95","vbK",
           "vbLinf","vbt0","wla","wlb","steep","LFC","LFS","Cref","Bref","Iref","Dep","Abun")
  
  type<-c("Catch","Index","Catch","Index","Recruitment_index","Catch_at_age","Catch_at_length",
          "Natural_mortality_rate","Maturity_at_length","Maturity_at_length","Growth","Growth","Growth",
          "Length_weight_conversion","Length_weight_conversion","Stock_recruitment_relationship",
          "Fleet_selectivity","Fleet_selectivity","Target_catch","Target_biomass","Target_index",
          "Index","Abundance") 
  
  ncases<-length(feaseobj@Case)
  slots<-slotNames(feaseobj)
  ns<-length(slots)
  ftab<-array(TRUE,c(ns-2,ncases))
  for(j in 3:ns)ftab[j-2,]<-as.logical(as.numeric(slot(feaseobj,slots[j])))
  
  req<-Required()
  nMPs<-nrow(req)
  gridy<-array("",c(nMPs,ncases))
  for(i in 1:ncases){
    types<-slotNames(feaseobj)[3:17][ftab[,i]]
    slots<-sloty[type%in%types]
    for(m in 1:nMPs){
      brec<-unlist(strsplit(req[m,2],", "))
      brec<-brec[grep("CV_",brec,invert=T)] #remove CV dependencies (we think we can guess these...)
      brec<-brec[brec!="Year"&brec!="MaxAge"&brec!="FMSY_M"&brec!="BMSY_B0"&brec!="t"&brec!="OM"&brec!="MPrec"&brec!="CAL_bins"]
      nr<-length(brec) 
      if(nr==0){
        gridy[m,i]<-"Yes"
      }else{ 
        cc<-0
        for(r in 1:nr){ #loop over requirements
          if(brec[r]%in%slots)cc<-cc+1
        }
        if(cc==nr)gridy[m,i]<-"Yes"
      }
    }
  }
  gridy<-as.data.frame(gridy)
  row.names(gridy)=req[,1]
  names(gridy)=feaseobj@Case
  if(outy=="table")return(gridy)
  if(outy!="table"&class(outy)!="numeric")return(req[,1][gridy[,1]=="Yes"])
  if(class(outy)=="numeric"){
    if(outy<(ncases+1)){
      return(req[,1][gridy[,as.integer(outy)]=="Yes"])
    }else{
      return(req[,1][gridy[,1]=="Yes"])
    }  
  }
  
}

# Internal functions used in runMSE 
getmov<-function(x,Prob_staying,Frac_area_1){
  test<-optim(par=c(0,0,0),movfit,method="L-BFGS-B",lower=rep(-6,3),upper=rep(6,3),prb=Prob_staying[x],frac=Frac_area_1[x])
  mov<-array(c(test$par[1],test$par[2],0,test$par[3]),dim=c(2,2))
  mov<-exp(mov)
  mov/array(apply(mov,1,sum),dim=c(2,2))
}

movfit<-function(par,prb,frac){
  mov<-array(c(par[1],par[2],0,par[3]),dim=c(2,2))
  mov<-exp(mov)
  mov<-mov/array(apply(mov,1,sum),dim=c(2,2))
  dis<-c(frac,1-frac)
  for(i in 1:100)dis<-apply(array(dis,c(2,2))*mov,2,sum)
  (log(mov[1,1])-log(prb))^2+(log(frac)-log(dis[1]))^2
}

getq<-function(x,dep,Find,Perr,Marray,hs,Mat_age,Wt_age,R0,V,nyears,maxage,mov,Spat_targ,SRrel,aR,bR){
  opt<-optimize(qopt,log(c(0.0075,15)),depc=dep[x],Fc=Find[x,],Perrc=Perr[x,],
                     Mc=Marray[x,],hc=hs[x],Mac=Mat_age[x,],Wac=Wt_age[x,,],
                     R0c=R0,Vc=V[x,,],nyears=nyears,maxage=maxage,movc=mov[x,,],
                     Spat_targc=Spat_targ[x],SRrelc=SRrel[x],aRc=aR[x,],bRc=bR[x,])
  return(exp(opt$minimum))
}

qopt<-function(lnq,depc,Fc,Perrc,Mc,hc,Mac,Wac,R0c,Vc,nyears,maxage,movc,Spat_targc,SRrelc,aRc,bRc,opt=T){
  qc<-exp(lnq)
  nareas<-nrow(movc)
  #areasize<-c(asizec,1-asizec)
  idist<-rep(1/nareas,nareas)
  for(i in 1:300)idist<-apply(array(idist,c(2,2))*movc,2,sum)

  N<-array(exp(-Mc[1]*((1:maxage)-1))*R0c,dim=c(maxage,nareas))*array(rep(idist,each=maxage),dim=c(maxage,nareas))
  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac[,1]
  SSB<-SSN*Wac[,1]                               # Calculate spawning stock biomass

  B0<-sum(Biomass)
  R0a<-idist*R0c
  SSB0<-apply(SSB,2,sum)
  SSBpR<-SSB0/R0a                              # Calculate spawning stock biomass per recruit

  for(y in 1:nyears){
    # set up some indices for indexed calculation
    targ<-(apply(Vc[,y]*Biomass,2,sum)^Spat_targc)/mean(apply(Vc[,y]*Biomass,2,sum)^Spat_targc)
    FMc<-array(qc*Fc[y]*Vc[,y],dim=c(maxage,nareas))*array(rep(targ,each=maxage),dim=c(maxage,nareas))                                           # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc<-FMc+Mc[y]
    N[2:maxage,]<-N[1:(maxage-1),]*exp(-Zc[1:(maxage-1),])         # Total mortality
    if(SRrelc==1){
      N[1,]<-Perrc[y]*(0.8*R0a*hc*apply(SSB,2,sum))/(0.2*SSBpR*R0a*(1-hc)+(hc-0.2)*apply(SSB,2,sum))  # Recruitment assuming regional R0 and stock wide steepness
    }else{
      N[1,]<- aRc*apply(SSB,2,sum)*exp(-bRc*apply(SSB,2,sum)) 
    }
      
    #print(N[1])
    indMov<-as.matrix(expand.grid(1:nareas,1:nareas,1:maxage)[3:1])
    indMov2<-indMov[,1:2]
    indMov3<-indMov[,2:3]
    temp<-array(N[indMov2]*movc[indMov3],dim=c(nareas,nareas,maxage))
    N<-apply(temp,c(3,1),sum)
    SSN<-N*Mac
    SSB<-SSN*Wac[,y]
    Biomass<-N*Wac[,y]
    SBiomass<-SSN*Wac[,y]
    #print(sum(Biomass))
  } # end of year
  return((log(depc)-log(sum(SBiomass)/sum(SSB0)))^2)
}


## Operating Model Functions ---------------------------------------------------
# These functions are used to manually specify, choose, or estimate various 
# parameters of the Operating Model.
# The functions typically take OM object (or Stock or Fleet) and return the same
# object with the relevant parameters populated.

# A highly dubious means of getting very uncertain estimates of current stock 
# biomass and (equilibrium) fishing mortality rate from growth, natural 
# mortality rate, recruitment and fishing selectivity.
ML2D<-function(OM,ML,nsim=100,ploty=T,Dlim=c(0.05,0.6)){
  
  maxage<-OM@maxage
  M<-runif(nsim,OM@M[1],OM@M[2]) # Natural mortality rate
  h<-runif(nsim,OM@h[1],OM@h[2]) # Steepness
  Linf<-runif(nsim,OM@Linf[1],OM@Linf[2]) # Maximum length
  K<-runif(nsim,OM@K[1],OM@K[2]) # Maximum growth rate
  t0<-runif(nsim,OM@t0[1],OM@t0[2]) # Theorectical length at age zero

  LFS<-runif(nsim,OM@LFS[1],OM@LFS[2])*mean(OM@L50)
  AFS<-L2A(t0,Linf,K,LFS,maxage)

  L5<-runif(nsim,OM@L5[1],OM@L5[2])*mean(OM@L50)
  age05<-L2A(t0,Linf,K,L5,maxage)

  Vmaxage<-runif(nsim,OM@Vmaxlen[1],OM@Vmaxlen[2])#runif(BT_fleet@Vmaxage[1],BT_fleet@Vmaxage[2]) # selectivity of oldest age class

  LM<-runif(nsim,OM@L50[1],OM@L50[2])
  AM<-L2A(t0,Linf,K,LM,maxage)

  # age at maturity
  a<-OM@a # length-weight parameter a
  b<-OM@b # length-weight parameter b

  mod<-AFS          # the age at modal (or youngest max) selectivity
  deriv<-getDNvulnS(mod,age05,Vmaxage,maxage,nsim)           # The vulnerability schedule
  vuln<-deriv[[1]]

  Agearray<-array(rep(1:maxage,each=nsim),c(nsim,maxage))
  mat<-1/(1+exp((AM-(Agearray))/(AM*0.1)))  # Maturity at age array

  nyears<-100
  #bootfun<-function(dat,ind)mean(dat[ind])
  #MLo<-boot(MLt,bootfun,nsim)
  #ML<-MLo$t
  out<-CSRA(M,h,Linf,K,t0,AM,a,b,vuln,mat,ML=rep(ML,nsim),NA,NA,maxage,nyears)
  cond<-out[,1]>Dlim[1]&out[,1]<Dlim[2]&out[,2]<2.5 # Stock levels are unlikely to be above 80% unfished, F is unlikely to be above 2.5
  
  if(ploty){
    par(mfrow=c(1,2))
    plot(density(out[cond,1],from=0,adj=0.4),main="Depletion")
    plot(density(out[cond,2],from=0,adj=0.4),main="Fishing mortality rate")
    OM@D<-quantile(out[cond,1],c(0.05,0.95))
  }

  OM
}

# Composition stock reduction analysis 
CSRA<-function(M,h,Linf,K,t0,AM,a,b,vuln,mat,ML,CAL,CAA,maxage,nyears){ 
  nsim<-length(M)
  Dep<-rep(NA,nsim)
  Fm<-rep(NA,nsim)
  for(i in 1:nsim){
    fit<-optimize(CSRAfunc,log(c(0.0001,5)),Mc=M[i],hc=h[i],maxage,nyears,Linfc=Linf[i],Kc=K[i],t0c=t0[i],AMc=AM[i],
                  ac=a,bc=b,vulnc=vuln[i,],matc=mat[i,],MLc=ML[i],CAL=NA,CAA=NA,opt=T)
    
    
    out<-CSRAfunc(fit$minimum,Mc=M[i],hc=h[i],maxage,nyears,Linfc=Linf[i],Kc=K[i],t0c=t0[i],AMc=AM[i],
                  ac=a,bc=b,vulnc=vuln[i,],matc=mat[i,],MLc=ML[i],CAL=NA,CAA=NA,opt=3)
    
    Dep[i]<-out[1]
    Fm[i]<-out[2]
    
    
  }
  cbind(Dep,Fm)
}

# The function that CSRA operates on
CSRAfunc<-function(lnF,Mc,hc,maxage,nyears,AFSc,AFCc,Linfc,Kc,t0c,AMc,ac,bc,vulnc,matc,MLc,CAL,CAA,opt=T,meth="ML"){
  
  Fm<-exp(lnF)
  Fc<-vulnc*Fm
  Lac<-Linfc*(1-exp(-Kc*((1:maxage)-t0c)))
  Wac<-ac*Lac^bc
  N<-exp(-Mc*((1:maxage)-1))
  SSN<-matc*N                                 # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  SSB<-SSN*Wac                               # Calculate spawning stock biomass
  
  B0<-sum(Biomass)
  SSB0<-sum(SSB)
  SSN0<-SSN
  SSBpR<-sum(SSB)                             # Calculate spawning stock biomass per recruit
  SSNpR<-SSN
  Zc<-Fc+Mc
  CN<-array(NA,dim=c(nyears,maxage))
  HR<-rep(0,maxage)
  pen<-0
  for(y in 1:nyears){
    VB<-Biomass*vulnc*exp(-Mc)
    CN[y,]<-N*(1-exp(-Zc))*(Fc/Zc)
    N[2:maxage]<-N[1:(maxage-1)]*exp(-Zc[1:(maxage-1)])         # Total mortality
    N[1]<-(0.8*hc*sum(SSB))/(0.2*SSBpR*(1-hc)+(hc-0.2)*sum(SSB))  # Recruitment assuming regional R0 and stock wide steepness
    Biomass<-N*Wac
    SSN<-N*matc
    SSB<-SSN*Wac
  } # end of year
  
  pred<-sum((CN[nyears,]*Lac))/sum(CN[nyears,])
  fobj<-(pred-MLc)^2 # Currently a least squares estimator. Probably not worth splitting hairs WRT likelihood functions!
  if(opt==1){return(fobj)
  }else{c(sum(SSB)/sum(SSB0),Fm)
  }
}

# Stochastic inverse growth curve used to back-calculate age at first capture from length at first capture
getAFC<-function(t0c,Linfc,Kc,LFC,maxage){ 
  nsim<-length(t0c)
  agev<-c(0.0001,1:maxage)
  agearray<-matrix(rep(agev,each=nsim),nrow=nsim)
  Larray<-Linfc*(1-exp(-Kc*(agearray-t0c)))
  matplot(agev,t(Larray),type='l')
  abline(h=LFC,col="#ff000030",lwd=2)
  AFC<-(log(1-(LFC/Linfc))/-Kc)+t0c
  abline(v=AFC,col="#0000ff30",lwd=2)
  AFC
}  

L2A<-function(t0c,Linfc,Kc,Len,maxage){ 
  nsim<-length(t0c)
  agev<-c(0.0001,1:maxage)
  agearray<-matrix(rep(agev,each=nsim),nrow=nsim)
  Larray<-Linfc*(1-exp(-Kc*(agearray-t0c)))
  matplot(agev,t(Larray),type='l')
  abline(h=Len,col="#ff000030",lwd=2)
  age<-(log(1-(Len/Linfc))/-Kc)+t0c
  abline(v=age,col="#0000ff30",lwd=2)
  age
}  

# Set Cyclic Recruitment Pattern
SetRecruitCycle <- function(x=1, Period, Amplitude, TotYears, Shape=c("sin", "shift")) {
  Shape <- match.arg(Shape) 
  Npers <- ceiling(TotYears/min(Period))
  pers <- round(runif(Npers, min(Period), max(Period)),0)
  amp <- runif(Npers, min(Amplitude), max(Amplitude))
  ct <- 1; Dir <- sample(c(-1, 1),1)
  Rm <- rep(NA, Npers)
  if (Shape == "sin") {
    for (X in 1:length(pers)) {
      yrper <- pers[X]
	  Period <- pi/yrper  #round(runif(1,-1, 1),0) * pi/yrper 
	  xs <- seq(from=0, by=1, to=pers[X])
	  Dir <- ifelse(Dir >= 0, -1, 1) # change direction each cycle
	  Rm[ct:(ct+pers[X])] <- amp[X] * sin(xs*Period) * Dir
	  ct <- ct+pers[X]
    }
  }
  if (Shape == "shift") {
    for (X in 1:length(pers)) {
	 Dir <- ifelse(Dir >= 0, -1, 1) # change direction each cycle
     if (X==1) Rm[ct:(ct+pers[X])] <- 0
     if (X > 1) Rm[ct:(ct+pers[X])] <- amp[X] * Dir
	 Rm[ct:(ct+pers[X])] <- amp[X] * Dir
     ct <- ct+pers[X]
    }
  }	
  Rm <- Rm[1:TotYears] + 1
  return(Rm)
}

# Sketch trends in historical fishing mortality --------------------------------
# Takes Fleet object, runs Sketch function for user to specify historical effort
# Then returns Fleet object with Effort objects populated with output of Sketch
ChooseEffort <- function(FleetObj, Years=NULL) {
  nyears <- FleetObj@nyears
  runSketch <- SketchFun(nyears, Years)
  FleetObj@EffYears <- runSketch[,1]
  FleetObj@EffLower <- runSketch[,2]
  FleetObj@EffUpper <- runSketch[,3]
  return(FleetObj)
}

# Sketch Historical Selectivity Patterns ---------------------------------------
ChooseSelect <- function(Fleet, Stock=NULL, FstYr=NULL, SelYears=NULL) {
  
  LastYr <- as.numeric(format(Sys.Date(), format="%Y"))
  if (is.null(FstYr)) {
    message("*****************")
    message("Enter first historical year")
    message("Note: *nyears* will be specified from this value")
    FstYr <- as.numeric(readline("First Historical Year: "))
	message("\n")
  }
  if (is.null(SelYears)) {
    message("Enter each selectivity break point year, seperated by a comma")
    message("Note: break points are the years where selectivity pattern changed")
    message("Note: break points must be within year range")
    inString <- readline("Enter each selectivity break point year, seperated by a comma: ")
    options(warn=-1)
    SelYears <- as.numeric(unlist(strsplit(inString, ",")))
    if (is.na(SelYears))  SelYears <- as.numeric(unlist(strsplit(inString, " ")))
    options(warn=0)
  }	
  SelYears <- sort(SelYears)
  if (SelYears[1] != FstYr) SelYears <- c(FstYr, SelYears)
  message("Break Points Years are: ")
  print(SelYears)
  if (length(SelYears) < 2) stop("Must be more than one year")
  if (max(SelYears) > LastYr) stop("Must specify historical year")
  if (min(SelYears) < FstYr) stop("Year before first year")
  flush.console()
  Selnyears <- length(SelYears)
 
  Years <- FstYr:LastYr #SelYears[1]:SelYears[length(SelYears)]
  Fleet@nyears <- length(Years)
  ind <- round((Range(SelYears, Max=LastYr, Min=FstYr)) * Fleet@nyears,0) +1
  ind[length(ind)] <- max(ind) - 1  
  Fleet@AbsSelYears <- SelYears 
  Fleet@SelYears <- ind
  
  tempL5 <- matrix(0, nrow=Selnyears, ncol=2)
  tempLFS <- matrix(0, nrow=Selnyears, ncol=2)
  tempmaxlen <- matrix(0, nrow=Selnyears, ncol=2)
  
  # if(is.null(Stock)) Stock <- NA
  set.par <- par(no.readonly=TRUE)
  message("Select selectivity points on plot")
  flush.console()
  for (N in 1:Selnyears) {
    BlankSelPlot(Stock=Stock, Yr=SelYears[N], N=N)
    L5Out <- ChooseL5()
    tempL5[N,] <- sort(L5Out[,1])
    LFSout <- ChooseLFS(L5Out)
    tempLFS[N,] <- sort(LFSout[,1])
    Vmaxout <- ChooseVmaxlen()
    tempmaxlen[N,] <- sort(Vmaxout[,2])
    polygon(x=c(0, max(tempL5[N,]), max(tempLFS[N,]), 3, 
     rev(c(0, min(tempL5[N,]), min(tempLFS[N,]), 3))),
	 y= c(0, 0.05, 1, min(tempmaxlen[N,]),
	 rev(c(0, 0.05, 1, max(tempmaxlen[N,])))), col="grey")
    par(ask=TRUE)
  }	
  par(set.par)
  # CheckSelect(Fleet, Stock)
  Fleet@L5Lower <- tempL5[,1]
  Fleet@L5Upper <- tempL5[,2]
  Fleet@LFSLower <- tempLFS[,1]
  Fleet@LFSUpper <- tempLFS[,2]
  Fleet@VmaxLower <- tempmaxlen[,1]
  Fleet@VmaxUpper <- tempmaxlen[,2]
  Fleet
}


# Kalman filter and Rauch-Tung-Striebel smoother
KalmanFilter <- function(RawEsts, R=1, Q=0.1, Int=100) {
  # Kalman smoother and Rauch-Tung-Striebel smoother #http://read.pudn.com/downloads88/ebook/336360/Kalman%20Filtering%20Theory%20and%20Practice,%20Using%20MATLAB/CHAPTER4/RTSvsKF.m__.htm
  # R  # Variance of sampling noise
  # Q  # Variance of random walk increments
  # Int # Covariance of initial uncertainty
  Ppred <-  rep(Int, length(RawEsts))
  nNA <- sum(is.na(RawEsts))
  while(nNA > 0) { # NAs get replaced with last non-NA
    RawEsts[is.na(RawEsts)] <- RawEsts[which(is.na(RawEsts))-1]
    nNA <- sum(is.na(RawEsts))
  }
  
  Pcorr <- xcorr <- xpred <- rep(0, length(RawEsts))
  # Kalman Filter
  for (X in 1:length(Ppred)) {
    if (X !=1) {
	  Ppred[X] <- Pcorr[X-1] + Q
	  xpred[X] <- xcorr[X-1]
	}
	W <- Ppred[X]/(Ppred[X] + R)
	xcorr[X] <- xpred[X] + W * (RawEsts[X] - xpred[X]) # Kalman filter estimate
	Pcorr[X] <- Ppred[X] - W * Ppred[X]
  }
  # Smoother 
  xsmooth <- xcorr
  for (X in (length(Pcorr)-1):1) {
    A <- Pcorr[X]/Ppred[X+1]
	xsmooth[X] <- xsmooth[X] + A*(xsmooth[X+1] - xpred[X+1]) 
  }
  return(xsmooth)

}


# LBSPR - Hordyk et al ICES 2015 (slow!)
LBSPRSim <- function(StockPars, FleetPars, SizeBins=NULL, P=0.001, Nage=101) {

  MK <- StockPars$MK 
  Linf <- StockPars$Linf
  CVLinf <- StockPars$CVLinf 
  L50 <- StockPars$L50 
  L95 <- StockPars$L95 
  Beta <- StockPars$FecB 
  MaxSD <- StockPars$MaxSD
  
  # Assumed constant CV here
  SDLinf <- CVLinf * Linf # Standard Deviation of Length-at-Age 
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 5
	SizeBins$ToSize <- Linf + MaxSD * SDLinf
  }
  if (is.null(SizeBins$ToSize)) SizeBins$ToSize <- Linf + MaxSD * SDLinf
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
  
  FM <- FleetPars$FM 
  SL50 <- FleetPars$SL50 
  SL95 <- FleetPars$SL95 
  
  LenBins <- seq(from=0, by=Linc, to=ToSize)	
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=length(LenBins)-1)
  x <- seq(from=0, to=1, length.out=Nage) # relative age vector
  EL <- (1-P^(x/MK)) * Linf # length at relative age 
  rLens <- EL/Linf # relative length 
  SDL <- EL * CVLinf # standard deviation of length-at-age
  
  Nlen <- length(LenMids) 
  Prob <- matrix(NA, nrow=Nage, ncol=Nlen)
  Prob[,1] <- pnorm((LenBins[2] - EL)/SDL, 0, 1) # probablility of length-at-age
  for (i in 2:(Nlen-1)) {
    Prob[,i] <- pnorm((LenBins[i+1] - EL)/SDL, 0, 1) - 
		pnorm((LenBins[i] - EL)/SDL, 0, 1)
  }
  Prob[,Nlen] <- 1 - pnorm((LenBins[Nlen] - EL)/SDL, 0, 1)
  
  # Truncate normal dist at MaxSD 
  mat <- array(1, dim=dim(Prob))
  for (X in 1:Nage) {
    ind <- which(abs((LenMids - EL[X]) /SDL[X]) >= MaxSD)
    mat[X,ind] <- 0
  }
  
  Prob <- Prob * mat

  SL <- 1/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50))) # Selectivity at length
  Sx <- apply(t(Prob) * SL, 2, sum) # Selectivity at relative age 
  MSX <- cumsum(Sx) / seq_along(Sx) # Mean cumulative selectivity for each age 
  Ns <- (1-rLens)^(MK+(MK*FM)*MSX) # number at relative age in population
  
  Cx <- t(t(Prob) * SL) # Conditional catch length-at-age probablilities  
  Nc <- apply(Ns * Cx, 2, sum) # 
  Pop <- apply(Ns * Prob, 2, sum)
  
  Ml <- 1/(1+exp(-log(19)*(LenMids-L50)/(L95-L50))) # Maturity at length
  Ma <-  apply(t(Prob) * Ml, 2, sum) # Maturity at relative age 
  
  N0 <- (1-rLens)^MK # Unfished numbers-at-age 
  SPR <- sum(Ma * Ns * rLens^Beta)/sum(Ma * N0 * rLens^Beta)
  
  Output <- NULL 
  Output$SPR <- SPR 
  Output$LenMids <- LenMids
  Output$PropLen <- Nc/sum(Nc)
  Output$Pop <- Pop
  
  Output$LCatchFished <- Nc/sum(Nc)
  Output$LPopFished <- Pop
  Output$LCatchUnfished <- apply(N0 * Cx, 2, sum)
  return(Output)
}  

OptFun <- function(tryFleetPars, LenDat, StockPars, SizeBins=NULL, 
	mod=c("GTG", "LBSPR")) {
  Fleet <- NULL
  Fleet$SL50 <- exp(tryFleetPars[1]) * StockPars$Linf
  Fleet$SL95 <- Fleet$SL50  + (exp(tryFleetPars[2]) * StockPars$Linf)
  Fleet$MLLKnife <- NA
  Fleet$FM <- exp(tryFleetPars[3])
  
  # if (mod == "GTG") runMod <-  GTGLBSPRSim(StockPars, Fleet, SizeBins)
  if (mod == "LBSPR") runMod <- LBSPRSim(StockPars, Fleet, SizeBins)
  
  LenDat <- LenDat + 1E-15 # add tiny constant for zero catches
  LenProb <- LenDat/sum(LenDat)
  predProb <- runMod$LCatchFished 
  predProb <- predProb + 1E-15 # add tiny constant for zero catches
  NLL <- -sum(LenDat * log(predProb/LenProb))
  
  if(!is.finite(NLL)) return(1E9)
  
  # add penalty for SL50 
  trySL50 <- exp(tryFleetPars[1])
  PenVal <- NLL
  Pen <- dbeta(trySL50, shape1=5, shape2=0.01) * PenVal
  if (Pen == 0) Pen <- PenVal * trySL50
  
  # plot(xx, dbeta(xx, shape1=5, shape2=0.01) )
  
  NLL <- NLL+Pen 

  return(NLL)
}

DoOpt <- function(StockPars, LenDat, SizeBins=NULL, mod=c("GTG", "LBSPR")) {
  
  SDLinf <- StockPars$CVLinf * StockPars$Linf
  if (is.null(SizeBins)) {
    SizeBins$Linc <- 5
	SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  }
  if (is.null(SizeBins$ToSize)) 
	SizeBins$ToSize <- StockPars$Linf + StockPars$MaxSD * SDLinf
  
  Linc <- SizeBins$Linc 
  ToSize <- SizeBins$ToSize
 
  LenBins <- seq(from=0, by=Linc, to=ToSize)	
  LenMids <- seq(from=0.5*Linc, by=Linc, length.out=length(LenBins)-1)
  
  sSL50 <- LenMids[which.max(LenDat)]/StockPars$Linf # Starting guesses
  sDel <- 0.2 * LenMids[which.max(LenDat)]/StockPars$Linf
  sFM <- 0.5 
  Start <- log(c(sSL50, sDel, sFM))
  
  # opt <- nlminb(Start, OptFun, LenDat=LenDat, StockPars=StockPars, 
	# SizeBins=SizeBins, mod=mod) 
	
  opt2 <- nlm(OptFun, Start, steptol=1e-4,gradtol=1e-4, LenDat=LenDat, 
	StockPars=StockPars, SizeBins=SizeBins, mod=mod)
  opt <- NULL
  opt$objective <- opt2$minimum
  opt$par <- opt2$estimate
	
  ModFailed <- FALSE 	
  if (opt$objective	== 1E9) ModFailed <- TRUE 
	# ,control= list(iter.max=300, eval.max=400, abs.tol=1E-20))
  # barplot(LenDat, names.arg=LenMids) 
  
  newFleet <- NULL 
  newFleet$FM <- exp(opt$par[3])
  newFleet$SL50 <- exp(opt$par[1]) * StockPars$Linf 
  newFleet$SL95 <- newFleet$SL50 + exp(opt$par[2]) * StockPars$Linf

  # if (mod == "GTG") runMod <-  GTGLBSPRSim(StockPars, newFleet, SizeBins) # not used
  if (mod == "LBSPR") runMod <- LBSPRSim(StockPars, newFleet, SizeBins)
  
  Out <- NULL 
  Out$Ests <- c(FM=newFleet$FM, SL50=newFleet$SL50, SL95=newFleet$SL95, 
	SPR=runMod$SPR)
  Out$PredLen <- runMod$LCatchFished * sum(LenDat)
  Out$ModFailed <- ModFailed
  return(Out)
}

# Run LBSPR Model for time-series of catch length composition data 
LBSPR <- function(x, DLM_data, yrsmth=1,reps=reps) {
  # Save other stuff for smoothing estimates
  TotYears <- nrow(DLM_data@CAL[1,,]) # How many years of length data exist
  if (length(DLM_data@Misc[[x]]) == 0) { # Misc List is empty
    # Create Empty List Object
	MiscList <- rep(list(0), 5) # Create empty list
	MiscList[[1]] <- rep(NA, TotYears) # SPR ests
	MiscList[[2]] <- rep(NA, TotYears) # Smoothed SPR ests
	MiscList[[3]] <- rep(NA, TotYears) # FM ests
	MiscList[[4]] <- rep(NA, TotYears) # Smoothed FM ests
	MiscList[[5]] <- list()
  }
  if (length(DLM_data@Misc[[x]]) != 0) MiscList <- DLM_data@Misc[[x]]
  
  # Add Extra Row when needed 
  if (length(MiscList[[1]]) < TotYears) {
    Diff <- TotYears - length(DLM_data@Misc[[x]][[1]])
    MiscList[[1]] <- append(MiscList[[1]], rep(NA,Diff)) 
	MiscList[[2]] <- append(MiscList[[2]], rep(NA,Diff)) 
	MiscList[[3]] <- append(MiscList[[3]], rep(NA,Diff)) 
	MiscList[[4]] <- append(MiscList[[4]], rep(NA,Diff)) 
  }

  NEmpty <- sum(is.na(MiscList[[1]])) # Number of empty spots
  IsEmpty <- which(is.na(MiscList[[1]]))

 StockPars <- NULL
 StockPars$MK <- DLM_data@Mort[x] / DLM_data@vbK[x]
 StockPars$Linf <- DLM_data@vbLinf[x]
 StockPars$CVLinf <- 0.1 # NEED TO ADD THIS TO INPUT VARIABLES
 StockPars$L50 <- DLM_data@L50[x] 
 StockPars$L95 <- DLM_data@L95[x]
 StockPars$FecB <- DLM_data@wlb[x]
 StockPars$MaxSD <- 2
 
 # yrsmth not implemented here
 LenMatrix <- DLM_data@CAL[x, IsEmpty,]
 
 binWidth <- DLM_data@CAL_bins[2] - DLM_data@CAL_bins[1]
 CAL_binsmid <- seq(from=0.5*binWidth, by=binWidth, length=length(DLM_data@CAL_bins)-1)
     
 SizeBins <- NULL
 SizeBins$Linc <- binWidth
 SizeBins$ToSize <- max(CAL_binsmid) + 0.5*binWidth

 # if(sfIsRunning()){
    # AllOpt <- sfSapply(1:length(IsEmpty), function (X) 
		# DoOpt(StockPars, LenDat=LenMatrix[X,], SizeBins=SizeBins, mod="LBSPR"))
  # } else {
    AllOpt <- sapply(1:length(IsEmpty), function (X)
		DoOpt(StockPars, LenDat=LenMatrix[X,], SizeBins=SizeBins, mod="LBSPR"))
  # } 
  
 EstFM <- sapply(AllOpt[1,], "[[", 1)
 estSL50 <- sapply(AllOpt[1,], "[[", 2)
 estSL95 <- sapply(AllOpt[1,], "[[", 3)
 EstSPR <- sapply(AllOpt[1,], "[[", 4)
 EstFM[EstFM > 5] <- 5 
 
 Fails <- which(sapply(AllOpt[3,], "[[", 1))
 if (length(Fails) > 0) {
   EstFM[Fails] <- NA 
   EstSPR[Fails] <- NA
 }
 while(sum(is.na(EstFM)) > 0) { # if model failed, make same as last time 
   EstFM[is.na(EstFM)] <- EstFM[which(is.na(EstFM))-1]
   EstSPR[is.na(EstSPR)] <- EstSPR[which(is.na(EstSPR))-1]
 }
 
 MiscList[[1]][IsEmpty] <- EstSPR # Save estimate of SPR for smoothing
 MiscList[[3]][IsEmpty] <- EstFM # Save estimate of F/M for smoothing 
  
  # Smoothed estimates - SPR
  MiscList[[2]] <- KalmanFilter(RawEsts=MiscList[[1]]) 
  MiscList[[2]][MiscList[[2]] <0] <- 0.05
  MiscList[[2]][MiscList[[2]] > 1] <- 0.99
 
  # Smoothed estimates - FM
  MiscList[[4]] <- KalmanFilter(RawEsts=MiscList[[3]])
  
  return(MiscList)
}



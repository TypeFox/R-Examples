# DLMtool Internal Functions 
# January 2016
# Tom Carruthers UBC (t.carruthers@fisheries.ubc.ca)
# Adrian Hordyk (a.hordyk@murdoch.edu.au)

# Various short functions that are used internally in the DLMtool.
# Functions are only used internally, and not directly accessible or available 
# to users of the package and therefore do not appear in the help manual.

# These can (and perhaps should) be exported to the package namespace 
# and have help files written for them.


# Misc. Defintions -------------------------------------------------------------
utils::globalVariables(c("R0","Mdb","mod","i","nareas","nsim","dFfinal"))
tiny <- 1E-15 # define tiny variable
proportionMat<-TL<-Wa<-SurvWeiMat<-r<-lx<-logNormDensity<-sumlogNormDen<-NULL
proportionMat=vector()


# Various Small Helper Functions -----------------------------------------------
# generic class finder  e.g. available('data.frame') avalable('DLM')
getclass <- function(x,classy) inherits(get(x),classy)

# Return TRUE if x is NULL or NA
NAor0 <- function(x){
  if(length(x)==0) return(TRUE)
  if(length(x)>0) return(is.na(x[1]))
}  

cv <- function(x) sd(x)/mean(x)
# get log normal standard deviation from transformed space mean 
# and standard deviation
sdconv <- function(m,sd)(log(1+((sd^2)/(m^2))))^0.5  
# get log normal mean from transformed space mean and standard deviation      
mconv<-function(m,sd)log(m)-0.5*log(1+((sd^2)/(m^2)
))    
alphaconv<-function(m,sd)m*(((m*(1-m))/(sd^2))-1)

betaconv<-function(m,sd)(1-m)*(((m*(1-m))/(sd^2))-1)

trlnorm<-function(reps,mu,cv)return(rlnorm(reps,mconv(mu,mu*cv),sdconv(mu,mu*cv)))

condmet<-function(vec)TRUE%in%vec

sampy<-function(x) sample(x,1,prob=!is.na(x))

range01 <- function(x) (x-min(x))/(max(x)-min(x)) # function to standardize to minimum and maximum
Range <- function(x, Max, Min) (x-Min)/(Max-Min) # function to standardize to minimum and maximum


# Functions --------------------------------------------------------------------
# Simulate historical fishing effort trajectory
getEffhist <- function(Esd, nyears, EffYears, EffLower, EffUpper) { 
  if(length(EffLower) == length(EffUpper) & length(EffUpper) == length(EffYears)) {
    nsim <- length(Esd) # get nsim 
    refYear <- floor(range01(EffYears+0.5) * nyears) + 1 # standardize years 
    refYear[length(refYear)] <- nyears # first year is year 1 
	Effs <- mapply(runif, n=nsim, min=EffLower, max=EffUpper) # sample Effort
    if (nsim > 1) {
	  effort <- t(sapply(1:nsim, function (x) approx(x=refYear, y = Effs[x,],  method = "linear", n=nyears)$y)) # linear interpolation
	}  
	if (nsim == 1) {
	  # Effs <- Effs/max(Effs)
	  effort <- approx(x=refYear, y = Effs,  method = "linear", n=nyears)$y
 	}
    
	effort <- range01(effort)
	
    Emu <- -0.5*Esd^2
    Eerr <-array(exp(rnorm(nyears*nsim,rep(Emu,nyears),rep(Esd,nyears))),c(nsim,nyears)) # calc error
    out <- NULL
    eff <-effort * Eerr # add error 
	out[[1]] <- eff
	out[[2]] <- (effort[,nyears]-effort[,nyears-4])/5
    return(out)
  } else {
    message("Input vectors of effort years and bounds not of same length")
    return(NULL)
  }
}

# Calculate FMSY and related
getFMSY<-function(x,Marray,hs,Mat_age,Wt_age,R0,V,maxage,nyears,proyears,Spat_targ,mov,SRrel,aR,bR){
  opt<-optimize(FMSYopt,log(c(0.001,5)),
                     Mc=Marray[x,nyears],hc=hs[x],Mac=Mat_age[x,],Wac=Wt_age[x,,nyears],
                     R0c=R0,Vc=V[x,],maxage=maxage,nyears=nyears,proyears=proyears,Spat_targc=Spat_targ[x],movc=mov[x,,],SRrelc=SRrel[x],aRc=aR[x,],bRc=bR[x,],Opt=T)
				
  return(FMSYopt(opt$minimum,
                     Mc=Marray[x,nyears],hc=hs[x],Mac=Mat_age[x,],Wac=Wt_age[x,,nyears],
                     R0c=R0,Vc=V[x,],maxage=maxage,nyears=nyears,proyears=proyears,Spat_targc=Spat_targ[x],movc=mov[x,,],SRrelc=SRrel[x],aRc=aR[x,],bRc=bR[x,],Opt=F)
         )
}

FMSYopt<-function(lnF,Mc,hc,Mac,Wac,R0c,Vc,maxage,nyears,proyears,Spat_targc,movc,SRrelc,aRc,bRc,Opt=T){

  FMSYc<-exp(lnF)
  nareas<-nrow(movc)
  #areasize<-c(asizec,1-asizec)
  idist<-rep(1/nareas,nareas)
  for(i in 1:100)idist<-apply(array(idist,c(2,2))*movc,2,sum)

  N<-array(exp(-Mc*((1:maxage)-1))*R0c,dim=c(maxage,nareas))*array(rep(idist,each=maxage),dim=c(maxage,nareas))
  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  VBiomass<-Biomass*Vc
  SSB<-SSN*Wac                              # Calculate spawning stock biomass
  
  B0 <- sum(Biomass)
  VB0<-sum(VBiomass)
  R0a<-idist*R0c
  SSB0<-apply(SSB,2,sum)
  SSBpR<-SSB0/R0a

  N<-N/2                              # Calculate spawning stock biomass per recruit
  SSN<-Mac*N   # Calculate initial spawning stock numbers
  Biomass<-N*Wac
  SSB<-SSN*Wac                              # Calculate spawning stock biomass

  for(y in 1:nyears){
    # set up some indices for indexed calculation
    dis<-apply(Vc*Biomass,2,sum)/sum(Vc*Biomass)
    targ<-(dis^Spat_targc)/mean(dis^Spat_targc)
    FMc<-array(FMSYc*Vc,dim=c(maxage,nareas))*array(rep(targ,each=maxage),dim=c(maxage,nareas))                                           # Fishing mortality rate determined by effort, catchability, vulnerability and spatial preference according to biomass
    Zc<-FMc+Mc
    CN<-N*(1-exp(-Zc))*(FMc/Zc)
    CB<-CN*Wac

    N[2:maxage,]<-N[1:(maxage-1),]*exp(-Zc[1:(maxage-1),])         # Total mortality
    if(SRrelc==1){
      N[1,]<-(0.8*R0a*hc*apply(SSB,2,sum))/(0.2*SSBpR*R0a*(1-hc)+(hc-0.2)*apply(SSB,2,sum))  # Recruitment assuming regional R0 and stock wide steepness
    }else{
      N[1,]<- aRc*apply(SSB,2,sum)*exp(-bRc*apply(SSB,2,sum)) 
    }
    #print(N[1])
    N[1,]<-apply(array(N[1,],c(2,2))*movc,2,sum)
    SSN<-N*Mac
    SSB<-SSN*Wac
    Biomass<-N*Wac
    VBiomass<-Biomass*Vc
    #print(sum(Biomass))
  } # end of year
 
  CBc<-sum(CB)
  if(Opt){
    return(-CBc)
  }else{
    return(c(CBc,-log(1-(CBc/(sum(VBiomass)+CBc))),sum(SSB)/sum(SSB0)))
  }
}

getFhist<-function(nsim,Esd,nyears,dFmin,dFmax,bb){

  ne<-nsim*3                                                         # Number of simulated effort datasets
  dEfinal<-runif(ne,dFmin,dFmax)#(exp(rnorm(ne,mean=demu,sd=desd))-1)*6               # Sample the final gradient in effort
  a<-(dEfinal-bb)/nyears                                         # Derive slope to get there from intercept
  a<-array(a,dim=c(ne,nyears))                                  # Slope array
  bb<-array(bb,dim=c(ne,nyears))                                  # Intercept array
  x<-array(rep(1:nyears,each=ne),dim=c(ne,nyears))              # Year array
  dE<-a*x+bb                                                     # Change in effort
  # E<-array(NA,dim=c(ne,nyears))                                 # Define total effort array
  # E[,1]<-dE[,1]
  # for(y in 2:nyears){
    # E[,y]<-apply(dE[,1:y],1,sum)
  # }
  # E<-E/array(apply(E,1,mean),dim=c(ne,nyears))                  # Standardise Effort to average 1
  E2 <- t(apply(dE, 1,cumsum))  								# Define total effort array
  E2 <- E2/array(apply(E2,1,mean),dim=c(ne,nyears))             # Standardise Effort to average 1
  E <- E2 # 
  cond <- apply(E,1,min)>0
  pos<-(1:ne)[cond]
  pos<-pos[1:nsim]
  #environment("dEfinal")<-asNamespace('DLMtool')#assign("dFfinal",dEfinal[pos],envir=.GlobalEnv)
  
  E<-E[pos,]                                 # Sample only those without negative effort
  Emu<--0.5*Esd^2
  Eerr<-array(exp(rnorm(nyears*nsim,rep(Emu,nyears),rep(Esd,nyears))),c(nsim,nyears))
  outy<-new('list')
  outy[[1]]<-E*Eerr
  outy[[2]]<-dEfinal[pos]
  outy
}

getFref<-function(x,Marray,Wt_age,Mat_age,Perr,N_s,SSN_s, Biomass_s,VBiomass_s,SSB_s,
                  Vn,hs,R0a,nyears,proyears,nareas,maxage,mov,SSBpR,aR,bR,SRrel){
  
  opt<-optimize(doprojPI,log(c(0.001,10)),
                Mvec=Marray[x,(nyears+1):(nyears+proyears)],Wac=Wt_age[x,,(nyears+1):(nyears+proyears)],Mac=Mat_age[x,],
                Pc=Perr[x,(nyears+1):(nyears+proyears)],N_c=N_s[x,,],SSN_c=SSN_s[x,,],Biomass_c=Biomass_s[x,,],
                VBiomass_c=VBiomass_s[x,,],SSB_c=SSB_s[x,,],Vc=Vn[x,],hc=hs[x],R0ac=R0a[x,],proyears,nareas,maxage,movc=mov[x,,],SSBpRc=SSBpR[x],aRc=aR[x,],bRc=bR[x,],SRrelc=SRrel[x])
  # print(exp(opt$minimum))
  return(-opt$objective)
  
}

doprojPI<-function(lnF,Mvec,Wac,Mac,Pc,N_c,SSN_c,Biomass_c,VBiomass_c,SSB_c,Vc,hc,R0ac,proyears,nareas,maxage,movc,SSBpRc,aRc,bRc,SRrelc){
  
  FF<-exp(lnF)
  
  N_P<-array(NA,dim=c(maxage,proyears,nareas))
  Biomass_P<-array(NA,dim=c(maxage,proyears,nareas))
  VBiomass_P<-array(NA,dim=c(maxage,proyears,nareas))
  SSN_P<-array(NA,dim=c(maxage,proyears,nareas))
  SSB_P<-array(NA,dim=c(maxage,proyears,nareas))
  FM_P<-array(NA,dim=c(maxage,proyears,nareas))
  Z_P<-array(NA,dim=c(maxage,proyears,nareas))
  CB_P<-rep(NA,proyears)
  
  AYR<-as.matrix(expand.grid(1:maxage,1,1:nareas))
  YA<-as.matrix(expand.grid(1,1:maxage))         # Projection year
  Y<-YA[,1]
  A<-YA[,2]
  AY<-YA[,c(2,1)]
  
  N_P[AYR]<-N_c#[AYRL]
  SSN_P[AYR]<-SSN_c#SSN[AYRL]
  Biomass_P[AYR]<-Biomass_c#[AYRL]
  VBiomass_P[AYR]<-VBiomass_c#[AYRL]
  SSB_P[AYR]<-SSB_c#[AYRL]
  
  FM_P[AYR]<-FF*Vc[A]
  Z_P[AYR]<-FM_P[A]+Mvec[Y]
  
  for(y in 2:proyears){
    
    AY1R<-as.matrix(expand.grid(1:maxage,y-1,1:nareas))
    AYR<-as.matrix(expand.grid(1:maxage,y,1:nareas))
    Y<-AYR[,2]
    A<-AYR[,1]
    AY<-AYR[,1:2]
    R<-AYR[,3]
    A2YR<-as.matrix(expand.grid(2:maxage,y,1:nareas))
    A1YR<-as.matrix(expand.grid(1:(maxage-1),y-1,1:nareas))
    A1Y<-as.matrix(expand.grid(1:(maxage-1),y-1))
    
    indMov<-as.matrix(expand.grid(1:nareas,1:nareas,y,1:maxage)[4:1])
    indMov2<-indMov[,c(1,2,3)]
    indMov3<-indMov[,c(3,4)]
    
    N_P[A2YR]<-N_P[A1YR]*exp(-Z_P[A1Y])         # Total mortality
    
    if(SRrelc==1){
      N_P[1,y,]<-Pc[y]*(0.8*R0ac*hc*apply(SSB_P[,y-1,],2,sum))/(0.2*SSBpRc*R0ac*(1-hc)+(hc-0.2)*apply(SSB_P[,y-1,],2,sum))  # Recruitment assuming regional R0 and stock wide steepness
    }else{
      N_P[1,y,]<-Pc[y]*aRc*apply(SSB_P[,y-1,],2,sum)*exp(-bRc*apply(SSB_P[,y-1,],2,sum))  
    }
    
    temp<-array(N_P[indMov2]*movc[indMov3],dim=c(nareas,nareas,maxage))  # Move individuals
    N_P[,y,]<-apply(temp,c(3,1),sum)
    
    Biomass_P[AYR]<-N_P[AYR]*Wac[AY]                                    # Calculate biomass
    VBiomass_P[AYR]<-Biomass_P[AYR]*Vc[A]                       # Calculate vulnerable biomass
    SSN_P[AYR] <-N_P[AYR]*Mac[A]                                       # Calculate spawning stock numbers
    SSB_P[AYR]<-SSN_P[AYR]*Wac[AY] # Calculate spawning stock biomass
    FM_P[AYR]<-FF*Vc[A]
    Z_P[AYR]<-FM_P[AYR]+Mvec[Y]
    CNtemp<-N_P[,y,]*exp(Z_P[,y,])*(1-exp(-Z_P[,y,]))*(FM_P[,y,]/Z_P[,y,])
    CB_P[y]<-sum(Biomass_P[,y,]*exp(Z_P[,y,])*(1-exp(-Z_P[,y,]))*(FM_P[,y,]/Z_P[,y,]))
	
	# CB_P[y] <- sum(CNtemp*Wac[AY]) 
  } # end of year
  return(-mean(CB_P[(proyears-min(4,(proyears-1))):proyears],na.rm=T))
  
}


# Calculate initial spatial distribution
getinitdist<-function(tol,mov,indMain){
  init<-array(1/nareas,dim=c(nsim,nareas))
  ind4<-as.matrix(cbind(rep(1:nsim,each=nareas*nareas),indMain[,2]))
  i<-0
  delta<-1
  #for(i in 1:100){
  while(delta > tol){
    i<-i+1
    trial<-init
    temp<-array(init[ind4]*mov[indMain],dim=c(nareas,nareas,nsim))
    init<-apply(temp,c(3,1),sum)
    delta<-max((trial-init)^2)
  }
  print(paste("Converged in ",i," iterations"))
  init
}

# Demographic model 
demofn<-function(log.r,M,amat,sigma,K,Linf,to,hR,maxage,a,b)demographic2(log.r,M,amat,sigma,K,Linf,to,hR,maxage=maxage,a,b)$epsilon

demographic2=function(log.r,M,amat,sigma,K,Linf,to,hR,maxage,a,b){		#switch on and off to use either S or m in MC simulations
  r=exp(log.r)
  lx=exp(-M)^((1:maxage)-1)				#survivorship
  logNormDensity=(dnorm(x=log((1:maxage)),mean=log(amat),sd=sigma))/(1:maxage)	#Maturity ogive calculation
  logNormDensity[1]=0
  sumlogNormDen=sum(logNormDensity)
  NormalisedMaturity=logNormDensity/sumlogNormDen
  proportionMat[1]=NormalisedMaturity[1]
    for(i in 2:maxage)  proportionMat[i]=proportionMat[i-1]+NormalisedMaturity[i]
  TL=Linf*(1-exp(-K*((1:maxage)-to)))		#length at age
  Wa=a*TL^b					#wegith at age
  SurvWeiMat=lx*Wa*proportionMat	#survivorship X weight X maturity
  SBPR=sum(SurvWeiMat)		#Spawner biomass per recruit
  RPS=1/(SBPR*(1-hR)/(4*hR)) # Beverton Holt
  #RPS=(5*hR)^(5/4)/SBPR			# Ricker Recruitment per spawner biomass
  RPF=Wa*proportionMat*RPS		#Recruits per female
  Lotka=lx*RPF*exp(-(1:maxage)*r)
  sumLotka=sum(Lotka)
  epsilon=(1-sumLotka)^2	 				#objective function
  return(list(epsilon=epsilon,r=r))
}

# Various functions for determing selectivity curve from input parameters
densnorm<-function(sd1){   # difference in density from 0.05 given a standard deviation sd1 (sd_asc) and age at maximum vulnerability modo
  (0.05-(dnorm(0,mod[i],sd1)/dnorm(mod[i],mod[i],sd1)))^2
}

densnormasc<-function(sd1,age_05,mody){
  (0.05-(dnorm(age_05,mody,sd1)/dnorm(mody,mody,sd1)))^2
}

getsdasc<-function(sm,age05,mod){
  optimize(densnormasc,interval=c(0.5,100),age_05=age05[sm],mody=mod[sm])$minimum
}

densnormdesc<-function(sd2,V_maxage,maxy,mody){
  (V_maxage-(dnorm(maxy,mody,sd2)/dnorm(mody,mody,sd2)))^2
}

getsddesc<-function(sm,Vmaxage,maxage,mod){
  optimize(densnormdesc,interval=c(0.5,10000),V_maxage=Vmaxage[sm],maxy=maxage,mody=mod[sm])$minimum
}

getDNvulnS<-function(mod,age05,Vmaxage,maxage,nsim){
  sd_asc<-sapply(1:nsim,getsdasc,age05=age05,mod=mod)
  sd_desc<-sapply(1:nsim,getsddesc,Vmaxage=Vmaxage,maxage=maxage,mod=mod)
  V<-array(NA,dim=c(nsim,maxage))
  for(i in 1:nsim){
    V[i,1:ceiling(mod[i])]<-dnorm(1:ceiling(mod[i]),mod[i],sd_asc[i])
    V[i,(1+ceiling(mod[i])):maxage]<-dnorm((1+ceiling(mod[i])):maxage,mod[i],sd_desc[i])
    V[i,(1+ceiling(mod[i])):maxage]<-V[i,(1+ceiling(mod[i])):maxage]/V[i,1+ceiling(mod[i])]#/V[i,floor(mod[i])+1]
    V[i,1:ceiling(mod[i])]<-V[i,1:ceiling(mod[i])]/dnorm(mod[i],mod[i],sd_asc[i])#,mod[i],sd_asc[i])#V[i,floor(mod[i])]

  }
  outy<-new('list')
  outy[[1]]<-V
  outy[[2]]<-mod-1.18*sd_asc
  outy
}

# Creates a time series per simulation that has gradient grad and random normal 
# walk with sigma
gettempvar<-function(targ,targsd,targgrad,nyears,nsim){   
  mutemp<--0.5*targsd^2
  temp<-array(1,dim=c(nsim,nyears))
  for(i in 2:nyears){
    temp[,i]<-temp[,i]*exp(rnorm(nsim,mutemp,targsd))
  }
  yarray<-array(rep((1:nyears)-1,each=nsim),dim=c(nsim,nyears))
  temp<-temp*(1+targgrad/100)^yarray
  targ*temp/apply(temp,1,mean)
}

# Growth-Type-Group Functions (not currently used - Jan 2016)
GenLenFun <- function(NatAGTG, LenatAgeGTG, LenBin, LenMid) {
  Nbins <- length(LenMid)
  SizeComp <- rep(0, Nbins)
  for (L in 1:length(LenMid)) {
    temp <- NatAGTG
    ind <- LenatAgeGTG <= LenBin[L+1] & LenatAgeGTG > LenBin[L]
    temp[!ind] <- 0
    SizeComp[L] <- SizeComp[L] + sum(temp)
  }
  return(SizeComp) 
}

# Two sided selectivity curve
TwoSidedFun <- function(L1, s1, s2, Lens) {
  Sl <- rep(0, length(Lens))
  Sl[Lens < L1] <- exp( -((Lens[Lens < L1] - L1)^2)/(2*s1^2))
  Sl[Lens >=L1 ] <- exp( -((Lens[Lens >= L1] - L1)^2)/(2*s2^2))
  return(Sl)
}

getroot<-function(X,ageM,age95)uniroot(getSlopeFun, interval=c(0.0001, 5), age50=ageM[X], age95=age95[X])$root
getSlope1 <- function(tst, L1, L0.05 ) (0.05 - TwoSidedFun(L1=L1, s1=tst, s2=1000, L0.05 ))^2
getSlope2 <- function(tst, L1, s1, Linf, MaxSel) (MaxSel - TwoSidedFun(L1=L1, s1=s1, s2=tst, Linf))^2
getSlopeFun <- function(SD, age50, age95) 0.95 - (1/(1+exp((age50-age95)/(age50*SD))))

# Selectivity at length function 
SelectFun <- function(i, SL0.05, SL1, MaxSel, Linfs, Lens) {
  s1 <- optimise(getSlope1, interval=c(0, 100), L1=SL1[i], L0.05=SL0.05[i])$minimum
  s2 <- optimise(getSlope2, interval=c(0, 1000), L1=SL1[i], s1=s1, Linf=Linfs[i], MaxSel=MaxSel[i])$minimum 
  TwoSidedFun(L1=SL1[i], s1=s1, s2=s2, Lens=Lens[i,])
}
# Selectivity at length function for GTG model 
SelectFunGTG <- function(i, SL0.05, SL1, MaxSel, Linfs, LenGTG) {
  s1 <- optimise(getSlope1, interval=c(0, 100), L1=SL1[i], L0.05=SL0.05[i])$minimum
  s2 <- optimise(getSlope2, interval=c(0, 1000), L1=SL1[i], s1=s1, Linf=Linfs[i], MaxSel=MaxSel[i])$minimum 
  NGTG <- dim(LenGTG)[1]
  t(sapply(1:NGTG, function (X) TwoSidedFun(L1=SL1[i], s1=s1, s2=s2, Lens=LenGTG[X,i,])))
}

# Obj value for opt routine
FitSelect <- function(Pars, V, Linf, Lens) {
  SL0.05 <- (Pars[1])
  SL1 <- (Pars[2])
  MaxSel <- (Pars[3])
  Lens <- t(as.matrix(Lens))
  SS <- sum((V-SelectFun(1, SL0.05, SL1, MaxSel, Linf, Lens))^2)
  return(SS)
}


# BlankSelPlot
BlankSelPlot <- function(Stock=NULL, Yr=NULL, N=NULL) {
  Max <- 3 
  AxCex <- 1.3
  By <- 0.05 
  par(mfrow=c(1,1), mai=c(2, 1, .5, .3), oma=c(1,1,1,1))
  plot(c(0,3), c(0,1), type="n", xlab="", ylab="", axes=FALSE)
  mtext(side=2, line=3, "Selectivity", cex=AxCex)
  axis(side=2)
  Xax <- seq(from=0, to=Max-By, by=2*By)
  axis(side=1, at=Xax)
  axis(side=1, at=c(2, Max), labels=c("", "Lmax"), xpd=NA)
  mtext(side=1, line=3.5, "Relative Length", cex=AxCex)
  axis(side=1, at=1, line=1.5, labels="L50")
  if (!is.null(Stock) & class(Stock) == "Stock") {
    L50 <- mean(Stock@L50) # mean length at maturity
    MatAx <- L50 * Xax
    axis(side=1, line=5.5, at=Xax, labels=MatAx)
    axis(side=1, line=5.5, at=c(2, Max), labels=c("", "Lmax"), xpd=NA)
    mtext(side=1, line=8.5, "Approx. Length", cex=AxCex)
    axis(side=1, at=1, line=6.5, labels="Mean L50")
  } 
  if (N == 1) {
    title(paste("Choose selectivity points for Year", Yr, "(First Year)"))
  } else {
    title(paste("Choose selectivity points for Year", Yr))
  }	
}
# Choose L5 
ChooseL5 <- function() {
  By <- 0.05
  Xs <-seq(from=0, to=1, by=By)
  Ys <- rep(0.05, length(Xs))
  points(Xs, Ys, col="gray", cex=0.5)
  text(0.5, 0.2, "Choose two points for L5")
  L5out <- identifyPch(x=Xs, y=Ys, tolerance=0.1, n=2)
  L5out
}
# Choose LFS 
ChooseLFS <- function(L5out) {
  Max <- 3
  By <- 0.05
  Xs <-seq(from=max(L5out[,1]), to=Max, by=By)
  Ys <- rep(1, length(Xs))
  points(Xs, Ys, col="gray", cex=0.5)
  text(1.5, 0.5, "Choose two points for LFS")
  LFSout <- identifyPch(x=Xs, y=Ys, tolerance=0.1, n=2)
  LFSout
}
# Choose Vmaxlen
ChooseVmaxlen <- function(L5out) {
  Max <- 3 
  Ys <- seq(from=0, to=1, by=0.05) 
  Xs <- rep(Max, length(Ys))
  points(Xs, Ys, col="gray", cex=0.5)
  text(2, 0.8, "Choose two points for selectivity\n at maximum length")
  Vmaxout <- identifyPch(x=Xs, y=Ys, tolerance=0.1, n=2)
  Vmaxout
}
	
# Rough Plot of Historical Selectivity Patterns --------------------------------
CheckSelect <- function(Fleet, Stock=NULL) {
 # NEEDS TO BE FIXED 
 if (length(Fleet@SelYears) < 1) stop("No break points in selectivity pattern")
 n <- length(Fleet@SelYears)
 if (n < 4) par(mfrow=c(n, 1), mar=c(4,4,1,1), oma=c(2,3,1,1), bty="l")
 if(n >= 4) {
   N <- ceiling(n/2)
   par(mfrow=c(N, 2), , mar=c(4,4,1,1), oma=c(2,3,1,1), bty="l")
 }
 for (X in 1:n) {
   plot(c(0,3), c(0,1), type="n", xlab="", ylab="")
   if(length(Fleet@AbsSelYears) > 0) title(Fleet@AbsSelYears[X])
   if(length(Fleet@AbsSelYears) == 0) title(Fleet@SelYears[X])
   polygon(x=c(0, max(Fleet@L5[X,]), max(Fleet@LFS[X,]), 3, 
     rev(c(0, min(Fleet@L5[X,]), min(Fleet@LFS[X,]), 3))),
	 y= c(0, 0.05, 1, min(Fleet@Vmaxlen[X,]),
	 rev(c(0, 0.05, 1, max(Fleet@Vmaxlen[X,])))), col="grey")
  lines(c(1,1), c(0,1), lty=3)
  text(1.1, 0.2, "L50", cex=1.25)
}
  mtext(side=2, outer=TRUE, "Selectivity", cex=1.25, xpd=NA)
  mtext(side=1, outer=TRUE, "Relative Length", cex=1.25, xpd=NA)
}
  
# Helper functions for above
identifyPch <- function(x, y = NULL, n = length(x), pch = 19, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
        ans <- identify(x[!sel], y[!sel], n = 1, plot = FALSE, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch)
        sel[ans] <- TRUE
		res <- c(res, ans)
    }
    out <- cbind(x[res], y[res])
	out <- out[order(out[,1]),]
	return(out)
}

# function to allow users to pick points off a plotted grid
SketchFun <- function(nyears=NULL, Years=NULL) {
  
  if (length(Years) == 0 & length(nyears) == 0) stop()
  if (length(Years) > 0) { 
    nyears <- length(Years)
	years <- Years
  }	
  if (length(Years) == 0) years <- 1:nyears
  par(mfrow=c(1,1), mar=c(5, 4, 5, 2))
  ys <- seq(from=0, to=1, by=0.05)
  years1 <- seq(from=years[1], by=2, to=max(years))
  years2 <- seq(from=years[2], by=2, to=max(years))
  grd1 <- expand.grid(years1, ys)
  grd2 <- expand.grid(years2, ys)
  grd3 <- expand.grid(years, ys)
  Xs <- grd3[,1]
  Ys <- grd3[,2]

  plot(grd1[,1], grd1[,2], col="black", pch=19, cex=0.2, xlab="Years", ylab="Variable", cex.lab=1.5)
  points(grd2[,1], grd2[,2], col="darkgrey", pch=19, cex=0.2)
  
  mtext(side=3, "Right Click to Finish. Escape to Quit.", xpd=NA, cex=1.25)
  line1 <- "Use mouse to select points on the grid"
  line2 <- "First and last year must be selected."
  line3 <- "Select two points in a single year to represent range of uncertainty"
  par(xpd=TRUE) 
  text(years[1],par("usr")[4]+0.15*(par("usr")[4]-par("usr")[3]), line1,cex=1, pos=4) 
  text(years[1],par("usr")[4]+0.1125*(par("usr")[4]-par("usr")[3]), line2,cex=1, pos=4) 
  text(years[1],par("usr")[4]+0.075*(par("usr")[4]-par("usr")[3]),line3,cex=1, pos=4)  
  par(xpd=FALSE)
  message(line1,"\n", line2, "\n", line3, "\n")
  flush.console()
  
  par()
  out <- NULL
  out <- identifyPch(x=Xs, y=Ys, tolerance=0.1)
  while(is.null(dim(out))) {
    message("Must choose more than one point")
	flush.console()
    out <- identifyPch(x=Xs, y=Ys, tolerance=0.1)
  }
  while(min(out[,1]) != years[1]) {
    message("Choose point(s) for first year (usually 0)")
	flush.console()
	dat <- rbind(out, identifyPch(x=Xs, y=Ys, tolerance=0.1))
	out <- dat[order(dat[,1]),]
  }
  while(max(out[,1]) != years[length(years)]) {
    message("Choose point(s) for last year (nyear)")
	flush.console()
	dat <- rbind(out, identifyPch(x=Xs, y=Ys, tolerance=0.1))
	out <- dat[order(dat[,1]),]
  }
  ord <- order(out[,1], out[,2])
  out <- out[ord,]
  yrs <- unique(out[,1])
  mat <- matrix(NA, nrow=length(yrs), ncol=3)
  mat[,1] <- yrs
  ind <- which(!duplicated(out[,1]))
  mat[,2:3] <- out[ind,2]
  for (X in seq_along(yrs)) {
    chk <- out[,1] %in% yrs[X]
	ind <- range(which(chk))
    if(sum(chk) > 1) {
	  mat[X,2:3] <- out[ind,2] 
	}
  }

  lines(mat[,1], mat[,2])
  lines(mat[,1], mat[,3])

  colnames(mat) <- c("Years", "Lower", "Upper")
  return(mat)
}


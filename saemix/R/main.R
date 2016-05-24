saemix<-function(model,data,control=list()) {

# Convergence plots during fit (special function, not user-level)
  convplot.infit<-function(allpar,K1,niter=0) {
# Convergence plots for all the fixed effects, random effects and residual variability
    np<-dim(allpar)[2]
    K<-dim(allpar)[1]
    n1<-round(sqrt(np))
    n2<-ceiling(np/n1)
    if(n1>5 | n2>5) {n1<-3;n2<-4}
    if(niter==0) niter<-K
    par(mfrow=c(n1,n2))
    for(j in 1:np) {
      plot(1:niter,allpar[1:niter,j],type="l", xlab="Iteration", ylab=colnames(allpar)[j])
      abline(v=K1)
    }
  }
  if(class(model)!="SaemixModel") {
    cat("Please provide a valid model object (see the help page for SaemixModel)\n")
    return()
  }
  if(class(data)!="SaemixData") {
    cat("Please provide a valid data object (see the help page for SaemixData)\n")
    return()
  }

  saemixObject<-new(Class="SaemixObject",data=data,model=model,options=control)
#  saemixObject<-new(Class="SaemixObject",data=saemix.data, model=saemix.model,options=saemix.options)
  opt.warn<-getOption("warn")
  if(!saemixObject["options"]$warnings) options(warn=-1)

  saemix.options<-saemixObject["options"]
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
#  showall(saemixObject)

# structural model, check nb of parameters
  structural.model<-saemix.model["model"]
  nb.parameters<-saemix.model["nb.parameters"]
  N<-saemix.data["N"]

# Stepsize
  stepsize<-rep(1,saemix.options$nbiter.tot)
stepsize[(saemix.options$nbiter.saemix[1]+1):saemix.options$nbiter.tot]<-1/
(1:saemix.options$nbiter.saemix[2])
  stepsize[1:saemix.options$nbiter.burn]<-0
  alpha1.sa<-saemix.options$alpha.sa
  alpha0.sa<-10^(-3/saemix.options$nbiter.sa)

# Random generator
  OLDRAND<-TRUE
  set.seed(saemix.options$seed)

# ECO TODO: integrate all this section in the object creation ?
# Initialisation: 
# create local copies modified of omega.init and covariate.model in saemix.model
# A la fin: i1.omega2 renomme en indx.omega et ind.res en indx.res
  i0.omega2<-which((1-mydiag(saemix.model["covariance.model"]))>0) # index of parameters without IIV
  indest.omega<-which(saemix.model["covariance.model"]>0)
#  i1.omega2<-which(mydiag(saemix.model$covariance.model)>0)
  i1.omega2<-saemix.model@indx.omega # index of parameters with IIV
  ind.res<-saemix.model["indx.res"]

# Covariate model & design matrix
  id<-saemix.data["data"][,saemix.data["name.group"]]
  if(length(saemix.data["name.covariates"])==0) tab<-data.frame(id=id) else
    tab<-data.frame(id=id,saemix.data["data"][, saemix.data["name.covariates",drop=FALSE]])
   temp2<-unique(tab)
   temp<-tab[!duplicated(id),,drop=FALSE]
   if(dim(temp)[1]!=dim(temp2)[1]) cat("Some covariates have time-varying values; only the first is taken into account in the current version of the algorithm.\n")
#temp<-temp[order(temp[,1]),]
  if(length(saemix.data["name.covariates"])>0) {
    Mcovariates<-data.frame(id=rep(1,N),temp[,2:dim(temp)[2]])} else {
    Mcovariates<-data.frame(id=rep(1,N))
   }
# removing from model unused lines
  j.cov<-which(rowSums(saemix.model["betaest.model"])>0)
  betaest.model<-saemix.model["betaest.model"][j.cov,,drop=FALSE]
  Mcovariates<-Mcovariates[,j.cov,drop=FALSE] # eliminate all the unused covariates
  for(icol in dim(Mcovariates)[2])
    if(is.factor(Mcovariates[,icol])) Mcovariates[,icol]<-as.numeric(Mcovariates[,icol])-1
  
#  if(length(j.cov)==1) {
#    betaest.model<-matrix(betaest.model,nrow=1, dimnames=list(c("Fixed"),colnames(saemix.model["betaest.model"])))
#    Mcovariates<-matrix(Mcovariates)
#  }
  saemix.model["betaest.model"]<-betaest.model
  temp1<-betaest.model[-c(1),,drop=FALSE]
#  if(is.null(dim(temp1))) temp1<-matrix(temp1,nrow=1, dimnames=list(rownames(betaest.model)[-c(1)], colnames(betaest.model)))  
  saemix.model["covariate.model"]<-temp1

  fixedpsi.ini<-saemix.model["psi0"][1,] # initial fixed effects (original parametrization)
  betaI.ini<-transpsi(matrix(fixedpsi.ini,nrow=1),saemix.model["transform.par"]) #initial fixed effects (Gaussian parametrization)
  fixed.ini<-saemix.model["betaest.model"]*0
  fixed.ini[1,]<-betaI.ini
  nr.psi0<-dim(saemix.model["psi0"])[1]
  nr.cov<-dim(saemix.model["betaest.model"])[1]
  if(nr.psi0>nr.cov) {
    saemix.model["psi0"]<-saemix.model["psi0"][1:nr.cov,,drop=FALSE]
    nr.psi0<-dim(saemix.model["psi0"])[1]
  }
#t1<-NULL
  if(nr.psi0<nr.cov) {
#  t1<-t(covariate.model[(nr.psi0+1):nr.cov,])
    psi1<-saemix.model["psi0"][nr.psi0,]
    for(j in (nr.psi0+1):nr.cov)
      saemix.model["psi0"]<-rbind(saemix.model["psi0"],psi1)
    nr.psi0<-dim(saemix.model["psi0"])[1]
  }
  if(nr.psi0>1) fixed.ini[2:nr.psi0,]<-saemix.model["psi0"][2:nr.psi0,]

#covariate.estim<-matrix(c(rep(saemix.model$fixed.estim,nr.psi0),t1),byrow=TRUE, nrow=nr.cov)
  covariate.estim<-matrix(rep(saemix.model["fixed.estim"],nr.psi0),byrow=TRUE, ncol=length(saemix.model["fixed.estim"]))
#if(!is.null(dim(t1))) covariate.estim<-rbind(covariate.estim,t1) 
  covariate.estim<-covariate.estim*saemix.model["betaest.model"]

  betas.ini<-fixed.ini[which(saemix.model["betaest.model"]>0)]
  betas.ini<-matrix(betas.ini,ncol=1)

  nb.betas<-sum(saemix.model["betaest.model"])
  ind.covariate<-which(saemix.model["betaest.model"]==1)
#matrix(which(covariate.model==1),nrow=1)

# Initialising
  yobs<-saemix.data["data"][,saemix.data["name.response"]]


# # Residual Error model.
# error models are a + bf described by [a b]
# error models :
#   constant            y = f + a*e
#   proportional        y = f + b*f*e
#   combined            y = f + (a+b*f)*e
#   exponential         y = f*exp(a*e)    ( <=>  log(y) = log(f) + a*e )
  ares.ini<-saemix.model["error.init"][1] # initial residual error model : constant coefficient
  bres.ini<-saemix.model["error.init"][2] #  initial residual error model : proportional coefficient
  ares<-ares.ini
  bres<-bres.ini
  pres<-c(ares,bres)

  nb.theta<-nb.parameters+length(i1.omega2)+length(ind.res)
  parpop<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1),ncol=nb.theta)
  colnames(parpop)<-c(saemix.model["name.modpar"],saemix.model["name.random"], saemix.model["name.res"][ind.res])
  allpar<-matrix(data=0,nrow=(saemix.options$nbiter.tot+1), ncol=(nb.betas+length(i1.omega2)+length(ind.res)))
  colnames(allpar)<-c(saemix.model["name.fixed"],saemix.model["name.random"], saemix.model["name.res"][ind.res])
  var.eta<-mydiag(saemix.model["omega.init"])
  theta0<-c(fixedpsi.ini,var.eta[i1.omega2],pres[ind.res])
  parpop[1,]<-theta0

# the covariates
LCOV<-MCOV<-matrix(data=0,nrow=nb.betas,ncol=nb.parameters)
j1<-1
COV<-matrix(nrow=dim(Mcovariates)[1],ncol=0)
pfix<-matrix(data=0,nrow=1,ncol=nb.parameters)
mean.phi<-matrix(data=0,nrow=N,ncol=nb.parameters)
for(j in 1:nb.parameters) {
  jcov<-which(saemix.model["betaest.model"][,j]==1)
  lambdaj<-fixed.ini[jcov,j]
  aj<-as.matrix(Mcovariates[,jcov])
  COV<-cbind(COV,aj)
  nlj<-length(lambdaj)
  j2<-j1+nlj-1
  LCOV[j1:j2,j]<-matrix(data=1,nrow=nlj,ncol=1)
  j1<-j2+1
  if(length(jcov)<=1) mean.phi[,j]<-aj*lambdaj else mean.phi[,j]<-aj%*%lambdaj
  pfix[j]<-length(lambdaj)
}
indx.betaI<-cumsum(c(0,pfix[1:(nb.parameters-1)]))+1
idx<-1:nb.betas
indx.betaC<-idx[is.na(match(idx,indx.betaI))]
saemix.model["indx.fix"]<-indx.betaI
saemix.model["indx.cov"]<-indx.betaC

l1<-betas.ini
l1[indx.betaI]<-transphi(matrix(l1[indx.betaI],nrow=1),saemix.model["transform.par"])
allpar[1,]<-c(l1,var.eta[i1.omega2],pres[ind.res])

COV2<-t(COV)%*%COV
j.covariate<-which(LCOV==1)
MCOV[j.covariate]<-betas.ini
betas<-betas.ini

ind.fix1<-which(covariate.estim[ind.covariate]==1)
ind.fix0<-which(covariate.estim[ind.covariate]==0)
COV1<-COV[,ind.fix1]
#if(length(ind.fix0)==1) dstatphi<-matrix(COV[,ind.fix0],ncol=1)%*%MCOV[ind.fix0,] else 
dstatphi<-COV[,ind.fix0,drop=FALSE]%*%MCOV[ind.fix0,]

covariate.estim1<-covariate.estim
covariate.estim1[,i0.omega2]<-0
ind.fix11<-which(covariate.estim1[ind.covariate]==1)
covariate.estim0<-covariate.estim
covariate.estim0[,i1.omega2]<-0
ind.fix10<-which(covariate.estim0[ind.covariate]==1)
MCOV0<-MCOV[ind.fix10,i0.omega2,drop=FALSE]
#if(is.null(dim(MCOV0)) & length(MCOV0)>0) MCOV0<-matrix(MCOV0,ncol=1)
COV0<-COV[,ind.fix10]
j0.covariate<-which(LCOV[ind.fix10,i0.omega2]==1)
flag.fmin<-as.integer(sum(covariate.estim0[1,])>0)

# using several Markov chains
  chdat<-new(Class="SaemixRepData",data=saemix.data, nb.chains=saemix.options$nb.chains)
  NM<-chdat["NM"]
  IdM<-chdat["dataM"]$IdM
  yM<-chdat["dataM"]$yM
  XM<-chdat["dataM"][,saemix.data["name.predictors"],drop=FALSE]
  io<-matrix(data=0,nrow=N,ncol=max(saemix.data["nind.obs"]))
  for(i in 1:N)
    io[i,1:saemix.data["nind.obs"][i]]<-1
  ioM<-do.call(rbind,rep(list(io),saemix.options$nb.chains))
  ind.ioM <- which(t(ioM)!=0)
  DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])

# Initialisation of phiM
if(length(i0.omega2)>0) {
  xmat<-covariate.estim[,i0.omega2]
  if(is.null(dim(xmat))) xmat<-matrix(xmat,ncol=length(i0.omega2))
  i0.temp<-which(colSums(xmat)==0)
  ind0.eta<-i0.omega2[i0.temp] # ind0.eta: index of parameters without IIV
} else ind0.eta<-c()
if(length(ind0.eta)>0) { # ind.eta: index of parameters with IIV
  idx<-1:nb.parameters
  ind.eta<-idx[-c(ind0.eta)]
} else ind.eta<-1:nb.parameters
nb.etas<-length(ind.eta)

itest.phi<-1:NM
ltest.phi<-length(itest.phi)
phiM<-matrix(data=0,nrow=NM,ncol=nb.parameters)
etaM<-matrix(data=0,nrow=NM,ncol=nb.etas)

mean.phiM<-do.call(rbind,rep(list(mean.phi),saemix.options$nb.chains))
kt<-0
fixed.psi<-fixedpsi.ini
omega<-saemix.model["omega.init"]
chol.omega<-try(chol(omega[ind.eta,ind.eta]),silent=TRUE)
if(class(chol.omega)=="try-error") {
#	cat("ind.eta=",ind.eta,"\n")
#	print(saemix.model["omega.init"])
#	print(omega[ind.eta,ind.eta])
	chol.omega<-saemix.model["omega.init"][ind.eta,ind.eta]<-omega[ind.eta, ind.eta]<-mydiag(nrow=length(ind.eta),ncol=length(ind.eta))
  cat("Problem inverting covariance matrix, setting initial Omega to diagonal.\n")
}

# Find a valid set of parameters wrt to the structural.model.
# Any parameter set that does not generate NaN, inf or imaginary numbers
# will satisfy this criteria.
phiMc<-mean.phiM
while (ltest.phi>0) {
    kt<-kt+1
    if (kt==100) 
        stop("stats:fit.saemix:FailedInitialParameterGuess\nFailed to find a valid initial parameter guess\n")
    end   
    etaMc<-0.5*matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
    phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
    etaM[itest.phi,]<-etaMc[itest.phi,]
    phiM[itest.phi,]<-phiMc[itest.phi,]
    psiM<-transphi(phiM,saemix.model["transform.par"])
    fpred<-structural.model(psiM, IdM, XM)
    inan<-(is.na(fpred)+is.infinite(fpred)+(Im(fpred)!=0))
    itest.phi<-unique(IdM[inan])
    ltest.phi<-length(itest.phi)
}

# initialization of the sufficient statistics
phi<-array(data=0,dim=c(N,nb.parameters, saemix.options$nb.chains))
statphi1<-0
statphi2<-0
statphi3<-0
statrese<-0

############################################
#  The Algorithm

omega.eta<-omega[ind.eta,ind.eta] # IIV matrix for estimated parameters
diag.omega<-mydiag(omega)
domega2<-do.call(cbind,rep(list((sqrt(mydiag(omega.eta)))*saemix.options$rw.ini),nb.etas))
VK<-rep(c(1:nb.etas),2)

Uargs<-list(i0.omega2=i0.omega2,MCOV0=MCOV0,COV0=COV0,j0.covariate=j0.covariate,
nmc=saemix.options$nb.chains,IdM=IdM,XM=XM,structural.model=structural.model,
transform.par=saemix.model["transform.par"],error.model=saemix.model["error.model"], ind.ioM=ind.ioM,yM=yM)

# hw=waitbar(1,'Estimating the population parameters (SAEM). Wait...');
if(saemix.options$displayProgress) par(ask=FALSE)
cat("Running main SAEM algorithm\n")
print(date())
for (kiter in 1:saemix.options$nbiter.tot) { # Iterative portion of algorithm

  if(kiter%%saemix.options$nbdisplay==0) {
    cat(".")
    if(saemix.options$displayProgress)    
      try(convplot.infit(allpar,saemix.options$nbiter.saemix[1],niter=(kiter-2)))
  }

  if(flag.fmin && kiter==saemix.options$nbiter.sa) {
  	COV1<-COV[,ind.fix11]
  	ind.prov<-!(ind.eta %in% i0.omega2)
  	domega2<-domega2[ind.prov,ind.prov,drop=FALSE] # keep in domega2 only indices of parameters with IIV
  	ind0.eta<-i0.omega2
    ind.eta<-1:nb.parameters  	
    if(length(ind0.eta)>0) ind.eta<-ind.eta[!(ind.eta %in% ind0.eta)] # update ind.eta, now only parameters with IIV
  	nb.etas<-length(ind.eta)
    VK<-rep(c(1:nb.etas),2)
    statphi1<-0
    statphi2<-0
    statphi3<-0
  }

  nb.etas<-length(ind.eta)
  domega<-cutoff(mydiag(omega[ind.eta,ind.eta]),.Machine$double.eps)
  omega.eta<-omega[ind.eta,ind.eta,drop=FALSE]
  omega.eta<-omega.eta-mydiag(mydiag(omega[ind.eta,ind.eta]))+mydiag(domega)
#  print(omega.eta)
  chol.omega<-try(chol(omega.eta))
# "/" dans Matlab = division matricielle, selon la doc "roughly" B*INV(A) (et *= produit matriciel...)
  d1.omega<-LCOV[,ind.eta]%*%solve(omega.eta)
  d2.omega<-d1.omega%*%t(LCOV[,ind.eta])
  comega<-COV2*d2.omega


# Simulation MCMC
  mean.phiM<-do.call(rbind,rep(list(mean.phi),saemix.options$nb.chains))
  phiM[,ind0.eta]<-mean.phiM[,ind0.eta]
  psiM<-transphi(phiM,saemix.model["transform.par"])
  fpred<-structural.model(psiM, IdM, XM)
  if(saemix.model["error.model"]=="exponential")
     fpred<-log(cutoff(fpred))
  gpred<-cutoff(ares+bres*abs(fpred))
  DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)^2+log(gpred)
  U.y<-colSums(DYF)

  etaM<-phiM[,ind.eta]-mean.phiM[,ind.eta,drop=FALSE]
#  if(length(ind.eta)==1) etaM<-matrix(etaM,ncol=1)
  phiMc<-phiM
  for(u in 1:saemix.options$nbiter.mcmc[1]) { # 1er noyau
    etaMc<-matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
    phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
    psiMc<-transphi(phiMc,saemix.model["transform.par"])
    fpred<-structural.model(psiMc, IdM, XM)
    if(saemix.model["error.model"]=="exponential")
      fpred<-log(cutoff(fpred))
    gpred<-cutoff(ares+bres*abs(fpred))
    DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)^2+log(gpred)
    Uc.y<-colSums(DYF)
    deltau<-Uc.y-U.y
    ind<-which(deltau<(-1)*log(runif(NM)))
    etaM[ind,]<-etaMc[ind,]
    U.y[ind]<-Uc.y[ind]
  }

  U.eta<-0.5*rowSums(etaM*(etaM%*%solve(omega.eta)))

# Second stage

  if(saemix.options$nbiter.mcmc[2]>0) {
    nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
    nrs2<-1
    for (u in 1:saemix.options$nbiter.mcmc[2]) {
     for(vk2 in 1:nb.etas) {
       etaMc<-etaM
       etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(NM*nrs2), ncol=nrs2)%*%mydiag(domega2[vk2,nrs2],nrow=1) # 2e noyau ? ou 1er noyau+permutation?
       phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
       psiMc<-transphi(phiMc,saemix.model["transform.par"])
       fpred<-structural.model(psiMc, IdM, XM)
       if(saemix.model["error.model"]=="exponential")
         fpred<-log(cutoff(fpred))
       gpred<-cutoff(ares+bres*abs(fpred))
       DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)**2+log(gpred)
       Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
       Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%solve(omega.eta)))
       deltu<-Uc.y-U.y+Uc.eta-U.eta
       ind<-which(deltu<(-1)*log(runif(NM)))
       etaM[ind,]<-etaMc[ind,]
       U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
       U.eta[ind]<-Uc.eta[ind]
       nbc2[vk2]<-nbc2[vk2]+length(ind)
       nt2[vk2]<-nt2[vk2]+NM
     }
    }
    domega2[,nrs2]<-domega2[,nrs2]*(1+saemix.options$stepsize.rw* (nbc2/nt2-saemix.options$proba.mcmc))
  }

  if(saemix.options$nbiter.mcmc[3]>0) {
    nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
    nrs2<-kiter%%(nb.etas-1)+2
    if(is.nan(nrs2)) nrs2<-1 # to deal with case nb.etas=1
    for (u in 1:saemix.options$nbiter.mcmc[3]) {
      if(nrs2<nb.etas) {
        vk<-c(0,sample(c(1:(nb.etas-1)),nrs2-1))
        nb.iter2<-nb.etas
      } else {
        vk<-0:(nb.etas-1)
#        if(nb.etas==1) vk<-c(0)
        nb.iter2<-1
      }
      for(k2 in 1:nb.iter2) {
        vk2<-VK[k2+vk]
        etaMc<-etaM
        etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(NM*nrs2), ncol=nrs2)%*%mydiag(domega2[vk2,nrs2])
        phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
        psiMc<-transphi(phiMc,saemix.model["transform.par"])
        fpred<-structural.model(psiMc, IdM, XM)
        if(saemix.model["error.model"]=="exponential")
         fpred<-log(cutoff(fpred))
        gpred<-cutoff(ares+bres*abs(fpred))
        DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)**2+log(gpred)
        Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
        Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%solve(omega.eta)))
        deltu<-Uc.y-U.y+Uc.eta-U.eta
        ind<-which(deltu<(-log(runif(NM))))
        etaM[ind,]<-etaMc[ind,]
#        if(kiter<20 | (kiter>150 & kiter<170)) {
#        	cat("kiter=",kiter,length(ind),"  ind.eta=",ind.eta,"  nrs2=",nrs2,"\n")
#        	print(head(etaMc))
#        }
        U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
        U.eta[ind]<-Uc.eta[ind]
        nbc2[vk2]<-nbc2[vk2]+length(ind)
        nt2[vk2]<-nt2[vk2]+NM
      }
    }
    domega2[,nrs2]<-domega2[,nrs2]*(1+saemix.options$stepsize.rw* (nbc2/nt2-saemix.options$proba.mcmc))
  }
  
  phiM[,ind.eta]<-mean.phiM[,ind.eta]+etaM
  psiM<-transphi(phiM,saemix.model["transform.par"])

# Algorithme, Part 2

  if(stepsize[kiter]>0) {
############# Stochastic Approximation
    fpred<-structural.model(psiM, IdM, XM)
    if(saemix.model["error.model"]=="exponential")
      fpred<-log(cutoff(fpred))
    ff<-matrix(fpred,nrow=saemix.data["ntot.obs"],ncol=saemix.options$nb.chains)
    for(k in 1:saemix.options$nb.chains) phi[,,k]<-phiM[((k-1)*N+1):(k*N),]
# overall speed similar
#    phi<-aperm(array(phiM,c(N,saemix.options$nb.chains,3)),c(1,3,2))
    stat1<-apply(phi[,ind.eta,,drop=FALSE],c(1,2),sum) # sommer sur les composantes ind.eta de phi, ? travers la 3?me dimension
    stat2<-matrix(data=0,nrow=nb.etas,ncol=nb.etas)
    stat3<-apply(phi**2,c(1,2),sum) # somme sur phi**2, ? travers la 3?me dimension
    statr<-0
    for(k in 1:saemix.options$nb.chains) {
     phik<-phi[,ind.eta,k]
     stat2<-stat2+t(phik)%*%phik
     fk<-ff[,k]
     if(!is.na(match(saemix.model["error.model"],c("constant","exponential"))))
      resk<-sum((yobs-fk)**2) else {
      if(saemix.model["error.model"]=="proportional")
       resk<-sum((yobs-fk)**2/cutoff(fk**2,.Machine$double.eps)) else resk<-0
     }
     statr<-statr+resk
    }
# Update sufficient statistics
    statphi1<-statphi1+stepsize[kiter]*(stat1/saemix.options$nb.chains-statphi1)
    statphi2<-statphi2+stepsize[kiter]*(stat2/saemix.options$nb.chains-statphi2)
    statphi3<-statphi3+stepsize[kiter]*(stat3/saemix.options$nb.chains-statphi3)
    statrese<-statrese+stepsize[kiter]*(statr/saemix.options$nb.chains-statrese)

############# Maximisation
##### fixed effects

    if (flag.fmin && kiter>=saemix.options$nbiter.sa) {
    temp<-d1.omega[ind.fix11,]*(t(COV1)%*%(statphi1-dstatphi[,ind.eta]))
    betas[ind.fix11]<-solve(comega[ind.fix11,ind.fix11],rowSums(temp)) 
# ECO TODO: utiliser optimise dans le cas de la dimension 1
    beta0<-optim(par=betas[ind.fix10],fn=compute.Uy,phi=phiM,ares=ares,bres=bres, args=Uargs,DYF=DYF,control=list(maxit=saemix.options$maxim.maxiter))$par # else
    betas[ind.fix10]<-betas[ind.fix10]+stepsize[kiter]*(beta0-betas[ind.fix10])
    } else {
    temp<-d1.omega[ind.fix1,]*(t(COV1)%*%(statphi1-dstatphi[,ind.eta]))
    betas[ind.fix1]<-solve(comega[ind.fix1,ind.fix1],rowSums(temp)) 
    }
 
    MCOV[j.covariate]<-betas
    mean.phi<-COV%*%MCOV
    e1.phi<-mean.phi[,ind.eta,drop=FALSE]
    
# Covariance of the random effects
    omega.full<-matrix(data=0,nrow=nb.parameters,ncol=nb.parameters)
    omega.full[ind.eta,ind.eta]<-statphi2/N + t(e1.phi)%*%e1.phi/N - t(statphi1)%*%e1.phi/N - t(e1.phi)%*%statphi1/N
    omega[indest.omega]<-omega.full[indest.omega]
    
# Simulated annealing (applied to the diagonal elements of omega)
    if (kiter<=saemix.options$nbiter.sa) {
        diag.omega.full<-mydiag(omega.full)
        vec1<-diag.omega.full[i1.omega2]
	vec2<-diag.omega[i1.omega2]*alpha1.sa
	idx<-as.integer(vec1<vec2)
	diag.omega[i1.omega2]<-idx*vec2+(1-idx)*vec1
        diag.omega[i0.omega2]<-diag.omega[i0.omega2]*alpha0.sa
    } else {
        diag.omega<-mydiag(omega)
    }
    omega<-omega-mydiag(mydiag(omega))+mydiag(diag.omega)
  
# Residual error
    if (saemix.model["error.model"]=="constant" | saemix.model["error.model"]=="exponential") {
      sig2<-statrese/saemix.data["ntot.obs"]
      ares<-sqrt(sig2)
    }
    if (saemix.model["error.model"]=="proportional") {
      sig2<-statrese/saemix.data["ntot.obs"]
      bres<-sqrt(sig2)
    }
    if (saemix.model["error.model"]=="combined") {
# ECO TODO: ? v?rifier (ici fpred<0 donc NaN, et puis que faire si bres<0 ???)
      ABres<-optim(par=c(ares,bres),fn=error,y=yM,f=fpred)$par
      if (kiter<=saemix.options$nbiter.saemix[1]) {
        ares<-max(ares*alpha1.sa,ABres[1])
        bres<-max(bres*alpha1.sa,ABres[2])
      } else {
        ares<-ares+stepsize[kiter]*(ABres[1]-ares)
        bres<-bres+stepsize[kiter]*(ABres[2]-bres)
      }
    }
    pres<-c(ares,bres)
    beta.I<-betas[indx.betaI]
    fixed.psi<-transphi(matrix(beta.I,nrow=1),saemix.model["transform.par"])
    betaC<-betas[indx.betaC]
    var.eta<-mydiag(omega)
    l1<-betas.ini
    l1[indx.betaI]<-fixed.psi
    l1[indx.betaC]<-betaC
    allpar[(kiter+1),]<-c(l1,var.eta[i1.omega2],pres[ind.res])
  } else { #end of loop on if (stepsize[kiter]>0)
    allpar[(kiter+1),]<-allpar[kiter,]
  }
  theta<-c(fixed.psi,var.eta[i1.omega2],pres[ind.res])
  parpop[(kiter+1),]<-theta

# End of loop on kiter
}
cat("\n    Minimisation finished\n")
print(date())

fixed.effects<-0*betas
fixed.effects[indx.betaI]<-fixed.psi
fixed.effects[indx.betaC]<-betaC
omega[i0.omega2,]<-0
omega[,i0.omega2]<-0

############# After end of iterations
##### SAEM convergence plots
# ECO TODO
# plot.saemix(saemix.res,saemix.data,saemix.model,saemix.options, niter=saemix.options$nbiter.tot)

##### Compute the individual parameters (MAP)
phi[,i0.omega2,1:saemix.options$nb.chains]<-mean.phi[,i0.omega2]
phi.samp<-phi
phi<-apply(phi,c(1,2),mean)

##### Conditional means and variances used for the estimation of the log-likelihood via Importance Sampling

cond.mean.phi<-phi
sphi1<-phi
sphi1[,ind.eta]<-statphi1
cond.mean.phi[,i1.omega2]<-sphi1[,i1.omega2]
cond.var.phi<-array(data=0,dim=dim(phi))
cond.var.phi[,i1.omega2]<-statphi3[,i1.omega2]-cond.mean.phi[,i1.omega2]**2
cond.mean.psi<-transphi(cond.mean.phi,saemixObject["model"]["transform.par"])

cond.mean.eta<-matrix(0,nrow=dim(etaM)[1],ncol=nb.parameters)
cond.mean.eta[,ind.eta]<-etaM
cond.mean.eta<-array(t(cond.mean.eta),dim=c(nb.parameters,N, saemix.options$nb.chains))
cond.mean.eta<-t(apply(cond.mean.eta,c(1,2),mean))

# Updating objects
  saemix.model["Mcovariates"]<-Mcovariates
  saemix.model["indx.res"]<-ind.res
  saemix.model["indx.fix"]<-indx.betaI
  saemix.model["indx.cov"]<-indx.betaC
  saemix.model["indx.omega"]<-i1.omega2

# Filling in result object
    saemix.res<-new(Class="SaemixRes",name.fixed=saemix.model["name.fixed"], name.random=saemix.model["name.random"],name.res=saemix.model["name.res"], fixed.effects=c(fixed.effects),fixed.psi=c(fixed.psi),betas=betas,betaC=betaC, omega=omega,respar=pres,cond.mean.phi=cond.mean.phi,cond.var.phi=cond.var.phi, mean.phi=mean.phi, phi=phi,phi.samp=phi.samp,parpop=parpop,allpar=allpar,MCOV=MCOV)
  saemix.res["indx.res"]<-ind.res
  saemix.res["indx.fix"]<-indx.betaI
  saemix.res["indx.cov"]<-indx.betaC
  saemix.res["indx.omega"]<-i1.omega2
  nb.parest<-sum(covariate.estim)+ sum(saemix.model["covariance.model"][upper.tri(saemix.model["covariance.model"], diag=TRUE)])+1+as.integer(saemix.model["error.model"]=="combined")
  saemix.res["npar.est"]<-nb.parest
  saemix.res["cond.mean.psi"]<-cond.mean.psi
  saemix.res["cond.mean.eta"]<-cond.mean.eta

# Updating elements of saemixObject
  saemixObject["model"]<-saemix.model
  saemixObject["results"]<-saemix.res
  saemixObject["options"]<-saemix.options
  saemixObject["rep.data"]<-chdat # Utile ? maybe remove rep.data

# ECO TODO check
# a la fin: mais verifier, pe pb de distribution ??? ie allpar sur l'echelle des betas et pas parpop ? a verifier
# saemix.res["allpar"]<-allpar
# saemix.res["parpop"]<-allpar[,-c(indx.betaC)]

#### Final computations
# Compute the MAP estimates of the PSI_i's 
  if(saemix.options$algorithms[1]) saemixObject<-map.saemix(saemixObject)

# Compute the Fisher Information Matrix & update saemix.res
  if(saemix.options$algorithms[2]) saemixObject<-fim.saemix(saemixObject)

# Estimate the log-likelihood via importance Sampling
  if(saemix.options$algorithms[3]) saemixObject<-llis.saemix(saemixObject)

#### Pretty printing the results (TODO finish in particular cov2cor)
  if(saemix.options$print) print(saemixObject,digits=2)

#### Save the results to a file
  if(saemix.options$save | saemix.options$save.graphs) {
# create directory to save the results
     if(saemix.options$directory!="") xsave<-dir.create(saemix.options$directory)
     if(!xsave) {
# Check that we're not trying to create a directory with the same name as a file
       if(!file_test("-d",saemix.options$directory)) {
         cat("Unable to create directory",saemix.options$directory)
         saemix.options$directory<-"newdir"
         dir.create(saemix.options$directory)         
         xsave<-file_test("-d",saemix.options$directory)
         if(!xsave) {
           saemix.options$directory<-""
           xsave<-TRUE
           cat(", saving in current directory.\n")
         } else cat(", saving results in newdir instead.\n")
       } else {
       xsave<-TRUE
       cat("Overwriting files in directory",saemix.options$directory,"\n")
       }
     }
   }
  if(saemix.options$save) {
    namres<-ifelse(saemix.options$directory=="","pop_parameters.txt", file.path(saemix.options$directory,"pop_parameters.txt"))
    xtry<-try(sink(namres))
    if(class(xtry)!="try-error") {
    print(saemixObject)
    sink()
    namres<-ifelse(saemix.options$directory=="","indiv_parameters.txt", file.path(saemix.options$directory,"indiv_parameters.txt"))
    if(length(saemixObject["results"]["map.psi"])>0)
       write.table(saemixObject["results"]["map.psi"],namres,quote=FALSE, row.names=FALSE)
     } else {
       cat("Unable to save results, check writing permissions and/or path to directory.\n")
     }
  }

# ECO TODO finish, adding all
  if(saemix.options$save.graphs) {
    saemixObject<-predict(saemixObject)
    if(saemix.options$directory=="") namgr<-"diagnostic_graphs.ps" else
      namgr<-file.path(saemix.options$directory,"diagnostic_graphs.ps")
    xtry<-try(postscript(namgr,horizontal=TRUE))
    if(class(xtry)!="try-error") {
    par(mfrow=c(1,1))
    try(plot(saemixObject,plot.type="data"))

    try(plot(saemixObject,plot.type="convergence"))

    if(length(saemixObject["results"]["ll.is"])>0) {
      par(mfrow=c(1,1))
      try(plot(saemixObject, plot.type="likelihood"))
    }

    try(plot(saemixObject,plot.type="observations.vs.predictions"))

    try(plot(saemixObject,plot.type="random.effects"))

    try(plot(saemixObject,plot.type="correlations"))

# Note: can replace all this by:
#    default.saemix.plots(saemixObject)

    dev.off()
    
    if(saemix.options$directory=="") namgr<-"individual_fits.ps" else
      namgr<-file.path(saemix.options$directory,"individual_fits.ps")
    postscript(namgr,horizontal=FALSE)
    try(plot(saemixObject,plot.type="individual.fit"))
    dev.off()
    } else {
       cat("Unable to save results, check writing permissions and/or path to directory.\n")
     }
  }

  options(warn=opt.warn)

  return(saemixObject)
}
###########################	Likelihood by IS	#############################

llis.saemix<-function(saemixObject) {
# Estimate the log-likelihood via importance Sampling

  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  ncov<-length(saemix.data)["name.covariates"]
  npred<-length(saemix.data["name.predictors"])
  yobs<-saemix.data["data"][,saemix.data["name.response"]]
  xind<-saemix.data["data"][,saemix.data["name.predictors"],drop=FALSE]

  i1.omega2<-saemix.model["indx.omega"]
  Omega<-saemix.res["omega"]
  a.res<-saemix.res["respar"][1]
  b.res<-saemix.res["respar"][2]
  cond.var.phi<-saemix.res["cond.var.phi"]
  cond.mean.phi<-saemix.res["cond.mean.phi"]
  nphi1<-length(i1.omega2)
  IOmega.phi1<-solve(Omega[i1.omega2,i1.omega2])
  mean.phi1<-saemix.res["mean.phi"][,i1.omega2]
  
  MM<-100
  KM<-round(saemixObject["options"]$nmc.is/MM)
  log.const<-0
  if(saemix.model["error.model"]=="exponential")
    log.const<-(-sum(yobs))
  IdM<-rep(c(0:(MM-1)),each=saemix.data["ntot.obs"])*saemix.data["N"]+ rep(saemix.data["data"][,"index"],MM)
  yM<-rep(yobs,MM)
  XM<-matrix(rep(t(xind),MM),ncol=npred, byrow=TRUE)

  io<-matrix(0,nrow=saemix.data["N"],ncol=max(saemix.data["nind.obs"]))
  for(isuj in 1:saemix.data["N"])
    io[isuj,1:saemix.data["nind.obs"][isuj]]<-1
  ioM<-matrix(rep(t(io),MM),ncol=dim(io)[2],byrow=TRUE)
  ind.ioM <- which(t(ioM)!=0)
  DYF<-matrix(0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])
  mean.phiM1<-matrix(rep(t(mean.phi1),MM),byrow=TRUE,ncol=nphi1)
  mtild.phiM1<-matrix(rep(t(cond.mean.phi[,i1.omega2]),MM),byrow=TRUE,ncol=nphi1)
# ECO TODO: securiser cas i1.omega2 de longueur 1
  cond.var.phi1<-cond.var.phi[,i1.omega2,drop=FALSE]
  for(i in dim(cond.var.phi1)[2])
     cond.var.phi1[,i]<-cutoff(cond.var.phi1[,i])
  stild.phiM1<-matrix(rep(t(sqrt(cond.var.phi1)),MM),byrow=TRUE,ncol=nphi1)
  phiM<-matrix(rep(t(cond.mean.phi),MM),byrow=TRUE,ncol=dim(cond.mean.phi)[2])
  meana<-rep(0,saemix.data["N"])
  LL<-matrix(0,nrow=KM,ncol=1)
  
  c2<- log(det(Omega[i1.omega2,i1.omega2,drop=FALSE])) + nphi1*log(2*pi)
  c1<-log(2*pi)
  if(saemixObject["options"]$print.is) par(mfrow=c(1,1))
  
  tit<-"Estimation of the log-likelihood"
  kmin<-min(10,ceiling(KM/4))
  for(km in 1:KM) {
    if(saemixObject["options"]$print.is & km>kmin & (trunc(KM/5))%%km==0) {
    x1<-MM*c(kmin:km)
    y1<-(-2)*LL[kmin:km]
      if(sum(!is.na(y1))) try(plot(x1,y1,type="l",xlab="Size of the Monte-Carlo sample", ylab="'-2xLog-Likelihood",main=tit))
    }
    r<-trnd.mlx(saemixObject["options"]$nu.is,saemix.data["N"]*MM,nphi1)
    phiM1<-mtild.phiM1+stild.phiM1*r
    dphiM<-phiM1-mean.phiM1
    
    d2<-(-0.5)*(rowSums(dphiM*(dphiM%*%IOmega.phi1)) + c2)
    e2<-matrix(d2,nrow=saemix.data["N"],ncol=MM)
    pitild.phi1<-rowSums(log(tpdf.mlx(r,saemixObject["options"]$nu.is)))
    e3<-matrix(pitild.phi1,nrow=saemix.data["N"],ncol=MM)- matrix(rep(0.5*rowSums(log(cond.var.phi1)),MM),ncol=MM)
    
    phiM[,i1.omega2]<-phiM1
    psiM<-transphi(phiM,saemix.model["transform.par"])
    f<-saemix.model["model"](psiM,IdM,XM)
    if(saemix.model["error.model"]=="exponential")
      f<-log(cutoff(f))
    g<-cutoff(a.res+b.res*abs(f))
    DYF[ind.ioM] <- -0.5*((yM-f)/g)**2 - log(g) - 0.5*c1
    e1<-matrix(colSums(DYF),nrow=saemix.data["N"],ncol=MM)
    sume<-e1+e2-e3
    newa<-rowMeans(exp(sume),na.rm=TRUE)
# ECO 11/05/03: added this line to avoid LL becoming NA due to NaN predicted values
#    newa[is.na(newa)]<-meana[is.na(newa)]
# 
    meana<-meana+1/km*(newa-meana)
    LL[km]<-sum(log(cutoff(meana)))+ log.const
  }

  x1<-MM*c(kmin:KM)
  y1<-(-2)*LL[kmin:KM]
  if(sum(!is.na(y1))) try(plot(x1,y1,type="l",xlab="Size of the Monte-Carlo sample", ylab="'-2xLog-Likelihood",main=tit)) else cat("Likelihood cannot be computed by Importance Sampling.\n")
  
  saemixObject["results"]["LL"]<-c(LL)
  saemixObject["results"]["ll.is"]<-LL[KM]
  saemixObject["results"]["aic.is"]<-(-2)*saemixObject["results"]["ll.is"]+ 2*saemixObject["results"]["npar.est"]
  saemixObject["results"]["bic.is"]<-(-2)*saemixObject["results"]["ll.is"]+ log(saemixObject["data"]["N"])*saemixObject["results"]["npar.est"]

  return(saemixObject)
}

###########################	Likelihood by GQ	#############################

llgq.saemix<-function(saemixObject) {
# RES = MLXGQ(RES) Estimate the log-likelihood using Gaussian Quadrature (multidimensional grid)
  nnodes.gq<-saemixObject["options"]$nnodes.gq  # number of nodes on each 1-D grid
  nsd.gq<-saemixObject["options"]$nsd.gq  # the integral is computed on the interval [E(eta|y) +- nsd_gq*SD(eta|y)]

  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  xind<-saemix.data["data"][,saemix.data["name.predictors"],drop=FALSE]
  yobs<-saemix.data["data"][,saemix.data["name.response"]]

  i1.omega2<-saemixObject["model"]["indx.omega"]
  Omega<-saemix.res["omega"]
  a.res<-saemix.res["respar"][1]
  b.res<-saemix.res["respar"][2]
  cond.var.phi<-saemix.res["cond.var.phi"]
  cond.mean.phi<-saemix.res["cond.mean.phi"]
  nphi1<-length(i1.omega2)
  IOmega.phi1<-solve(Omega[i1.omega2,i1.omega2])
  mean.phi1<-saemix.res["mean.phi"][,i1.omega2]

  io<-matrix(0,nrow=saemix.data["N"],ncol=max(saemix.data["nind.obs"]))
  for(isuj in 1:saemix.data["N"])
    io[isuj,1:saemix.data["nind.obs"][isuj]]<-1
  ind.io <- which(t(io)!=0)
  DYF<-matrix(0,nrow=dim(io)[2],ncol=dim(io)[1])
  
  phi<-saemix.res["mean.phi"]
  y<-gqg.mlx(nphi1,nnodes.gq)
  x<-(y$nodes-0.5)*2
  w<-(y$weights)*(2**nphi1)
# ECO TODO check dimensions (unclear in matlab)
  nx<-dim(x)[1]
  condsd.eta<-sqrt(cond.var.phi[,i1.omega2])
  xmin<-cond.mean.phi[,i1.omega2]-nsd.gq*condsd.eta
  xmax<-cond.mean.phi[,i1.omega2]+nsd.gq*condsd.eta
  a<-(xmin+xmax)/2
  b<-(xmax-xmin)/2
    log.const<-0
  if(saemixObject["model"]["error.model"]=="exponential")
    log.const<-(-sum(yobs))
 
  Q<-0
  for (j in 1:nx) {
    phi[,i1.omega2] <- a+b*matrix(rep(x[j,],saemix.data["N"]),ncol=nphi1,byrow=TRUE)
    psi<-transphi(phi,saemixObject["model"]["transform.par"])
    f<-saemixObject["model"]["model"](psi, saemix.data["data"][,"index"], xind)
    if(saemixObject["model"]["error.model"]=="exponential")
      f<-log(cutoff(f))
    g<-cutoff(a.res+b.res*abs(f))
    DYF[ind.io] <- -0.5*((yobs-f)/g)**2 - log(g)
    ly<-colSums(DYF)
    dphi1<-phi[,i1.omega2]-saemix.res["mean.phi"][,i1.omega2]
    lphi1<-(-0.5)*rowSums((dphi1%*%IOmega.phi1)*dphi1)
    ltot<-ly+lphi1
    ltot[is.na(ltot)]<-(-Inf)
    Q<-Q+w[j]*exp(ltot)
  }
  S<-saemix.data["N"]*log(det(Omega[i1.omega2,i1.omega2]))+ saemix.data["N"]*nphi1*log(2*pi)+ saemix.data["ntot.obs"]*log(2*pi)
  ll<-(-S/2) + sum(log(Q)+rowSums(log(b)))+ log.const
  saemixObject["results"]["ll.gq"]<-ll
  saemixObject["results"]["aic.gq"]<-(-2)*saemixObject["results"]["ll.gq"]+ 2*saemixObject["results"]["npar.est"]
  saemixObject["results"]["bic.gq"]<-(-2)*saemixObject["results"]["ll.gq"]+ log(saemixObject["data"]["N"])*saemixObject["results"]["npar.est"]

  return(saemixObject)
}

gqg.mlx<-function(dim,nnodes.gq) {
#GQG.MLX Nodes and weights for numerical integration on grids
#(multidimensional Gaussian Quadrature)
#    dim  : dimension of the integration problem
#    nnodes.gq   : number of points on any 1-D grid
#
#    x    = matrix of nodes with dim columns
#    w    = row vector of corresponding weights
#
  if(nnodes.gq>25) {
    cat("The number of nodes for Gaussian Quadrature should be less than 25.\n")
    return(list(nodes=NULL,weights=c()))
  }
  if(nnodes.gq==1) {
    n<-c(5.0000000000000000e-001)
    w<-c(1.0000000000000000e+000) }
  if(nnodes.gq==2) {
    n<-c(7.8867513459481287e-001)
    w<-c(5.0000000000000000e-001) }
  if(nnodes.gq==3) {
    n<-c(5.0000000000000000e-001, 8.8729833462074170e-001)
    w<-c(4.4444444444444570e-001, 2.7777777777777712e-001) }
  if(nnodes.gq==4) {
    n<-c(6.6999052179242813e-001, 9.3056815579702623e-001)
    w<-c(3.2607257743127516e-001, 1.7392742256872484e-001) }
  if(nnodes.gq==5) {
    n<-c(5.0000000000000000e-001, 7.6923465505284150e-001, 9.5308992296933193e-001)
    w<-c(2.8444444444444655e-001, 2.3931433524968501e-001, 1.1846344252809174e-001) }
  if(nnodes.gq==6) {
    n<-c(6.1930959304159849e-001, 8.3060469323313235e-001, 9.6623475710157603e-001)
    w<-c(2.3395696728634746e-001, 1.8038078652407072e-001, 8.5662246189581834e-002) }
  if(nnodes.gq==7) {
    n<-c(5.0000000000000000e-001, 7.0292257568869854e-001, 8.7076559279969723e-001, 9.7455395617137919e-001)
    w<-c(2.0897959183673620e-001, 1.9091502525256090e-001, 1.3985269574463935e-001, 6.4742483084431701e-002) }
  if(nnodes.gq==8) {
    n<-c(5.9171732124782495e-001, 7.6276620495816450e-001, 8.9833323870681348e-001, 9.8014492824876809e-001)
    w<-c(1.8134189168918213e-001, 1.5685332293894469e-001, 1.1119051722668793e-001, 5.0614268145185180e-002) }
  if(nnodes.gq==9) {
    n<-c(5.0000000000000000e-001, 6.6212671170190451e-001, 8.0668571635029518e-001, 9.1801555366331788e-001, 9.8408011975381304e-001)
    w<-c(1.6511967750063075e-001, 1.5617353852000226e-001, 1.3030534820146844e-001, 9.0324080347429253e-002, 4.0637194180784583e-002) }
  if(nnodes.gq==10) {
    n<-c(5.7443716949081558e-001, 7.1669769706462361e-001, 8.3970478414951222e-001, 9.3253168334449232e-001, 9.8695326425858587e-001)
    w<-c(1.4776211235737713e-001, 1.3463335965499873e-001, 1.0954318125799158e-001, 7.4725674575290599e-002, 3.3335672154342001e-002) }
  if(nnodes.gq==11) {
    n<-c(5.0000000000000000e-001, 6.3477157797617245e-001, 7.5954806460340585e-001, 8.6507600278702468e-001, 9.4353129988404771e-001, 9.8911432907302843e-001)
    w<-c(1.3646254338895086e-001, 1.3140227225512388e-001, 1.1659688229599563e-001, 9.3145105463867520e-002, 6.2790184732452625e-002, 2.7834283558084916e-002) }
  if(nnodes.gq==12) {
    n<-c(5.6261670425573451e-001, 6.8391574949909006e-001, 7.9365897714330869e-001, 8.8495133709715235e-001, 9.5205862818523745e-001, 9.9078031712335957e-001)
    w<-c(1.2457352290670189e-001, 1.1674626826917781e-001, 1.0158371336153328e-001, 8.0039164271673444e-002, 5.3469662997659276e-002, 2.3587668193254314e-002) }
  if(nnodes.gq==13) {
    n<-c(5.0000000000000000e-001, 6.1522915797756739e-001, 7.2424637551822335e-001, 8.2117466972017006e-001, 9.0078904536665494e-001, 9.5879919961148907e-001, 9.9209152735929407e-001)
    w<-c(1.1627577661543741e-001, 1.1314159013144903e-001, 1.0390802376844462e-001, 8.9072990380973202e-002, 6.9436755109893875e-002, 4.6060749918864378e-002, 2.0242002382656228e-002) }
  if(nnodes.gq==14) {
    n<-c(5.5402747435367183e-001, 6.5955618446394482e-001, 7.5762431817907705e-001, 8.4364645240584268e-001, 9.1360065753488251e-001, 9.6421744183178681e-001, 9.9314190434840621e-001)
    w<-c(1.0763192673157916e-001, 1.0259923186064811e-001, 9.2769198738969161e-002, 7.8601583579096995e-002, 6.0759285343951711e-002, 4.0079043579880291e-002, 1.7559730165874574e-002) }
  if(nnodes.gq==15) {
    n<-c(5.0000000000000000e-001, 6.0059704699871730e-001, 6.9707567353878175e-001, 7.8548608630426942e-001, 8.6220886568008503e-001, 9.2410329170521366e-001, 9.6863669620035298e-001, 9.9399625901024269e-001)
    w<-c(1.0128912096278091e-001, 9.9215742663556039e-002, 9.3080500007781286e-002, 8.3134602908497196e-002, 6.9785338963077315e-002, 5.3579610233586157e-002, 3.5183023744054159e-002, 1.5376620998057434e-002) }
  if(nnodes.gq==16) {
    n<-c(5.4750625491881877e-001, 6.4080177538962946e-001, 7.2900838882861363e-001, 8.0893812220132189e-001, 8.7770220417750155e-001, 9.3281560119391593e-001, 9.7228751153661630e-001, 9.9470046749582497e-001)
    w<-c(9.4725305227534431e-002, 9.1301707522462000e-002, 8.4578259697501462e-002, 7.4797994408288562e-002, 6.2314485627767105e-002, 4.7579255841246545e-002, 3.1126761969323954e-002, 1.3576229705875955e-002) }
  if(nnodes.gq==17) {
    n<-c(5.0000000000000000e-001, 5.8924209074792389e-001, 6.7561588172693821e-001, 7.5634526854323847e-001, 8.2883557960834531e-001, 8.9075700194840068e-001, 9.4011957686349290e-001, 9.7533776088438384e-001, 9.9528773765720868e-001)
    w<-c(8.9723235178103419e-002, 8.8281352683496447e-002, 8.4002051078225143e-002, 7.7022880538405308e-002, 6.7568184234262890e-002, 5.5941923596702053e-002, 4.2518074158589644e-002, 2.7729764686993612e-002, 1.2074151434273140e-002) }
  if(nnodes.gq==18) {
    n<-c(5.4238750652086765e-001, 6.2594311284575277e-001, 7.0587558073142131e-001, 7.7988541553697377e-001, 8.4584352153017661e-001, 9.0185247948626157e-001, 9.4630123324877791e-001, 9.7791197478569880e-001, 9.9578258421046550e-001)
    w<-c(8.4571191481571939e-002, 8.2138241872916504e-002, 7.7342337563132801e-002, 7.0321457335325452e-002, 6.1277603355739306e-002, 5.0471022053143716e-002, 3.8212865127444665e-002, 2.4857274447484968e-002, 1.0808006763240719e-002) }
  if(nnodes.gq==19) {
    n<-c(5.0000000000000000e-001, 5.8017932282011264e-001, 6.5828204998181494e-001, 7.3228537068798050e-001, 8.0027265233084055e-001, 8.6048308866761469e-001, 9.1135732826857141e-001, 9.5157795180740901e-001, 9.8010407606741501e-001, 9.9620342192179212e-001)
    w<-c(8.0527224924391946e-002, 7.9484421696977337e-002, 7.6383021032929960e-002, 7.1303351086803413e-002, 6.4376981269668232e-002, 5.5783322773667113e-002, 4.5745010811225124e-002, 3.4522271368820669e-002, 2.2407113382849821e-002, 9.7308941148624341e-003) }
  if(nnodes.gq==20) {
    n<-c(5.3826326056674867e-001, 6.1389292557082253e-001, 6.8685304435770977e-001, 7.5543350097541362e-001, 8.1802684036325757e-001, 8.7316595323007540e-001, 9.1955848591110945e-001, 9.5611721412566297e-001, 9.8198596363895696e-001, 9.9656429959254744e-001)
    w<-c(7.6376693565363113e-002, 7.4586493236301996e-002, 7.1048054659191187e-002, 6.5844319224588346e-002, 5.9097265980759248e-002, 5.0965059908620318e-002, 4.1638370788352433e-002, 3.1336024167054569e-002, 2.0300714900193556e-002, 8.8070035695753026e-003) }
  if(nnodes.gq==21) {
    n<-c(5.0000000000000000e-001, 5.7278092708044759e-001, 6.4401065840120053e-001, 7.1217106010371944e-001, 7.7580941794360991e-001, 8.3356940209870611e-001, 8.8421998173783889e-001, 9.2668168229165859e-001, 9.6004966707520034e-001, 9.8361341928315316e-001, 9.9687608531019478e-001)
    w<-c(7.3040566824845346e-002, 7.2262201994985134e-002, 6.9943697395536658e-002, 6.6134469316668845e-002, 6.0915708026864350e-002, 5.4398649583574356e-002, 4.6722211728016994e-002, 3.8050056814189707e-002, 2.8567212713428641e-002, 1.8476894885426285e-002, 8.0086141288864491e-003) }
  if(nnodes.gq==22) {
    n<-c(5.3486963665986109e-001, 6.0393021334411068e-001, 6.7096791044604209e-001, 7.3467791899337853e-001, 7.9382020175345580e-001, 8.4724363159334137e-001, 8.9390840298960406e-001, 9.3290628886015003e-001, 9.6347838609358694e-001, 9.8503024891771429e-001, 9.9714729274119962e-001)
    w<-c(6.9625936427816129e-002, 6.8270749173007697e-002, 6.5586752393531317e-002, 6.1626188405256251e-002, 5.6466148040269712e-002, 5.0207072221440600e-002, 4.2970803108533975e-002, 3.4898234212260300e-002, 2.6146667576341692e-002, 1.6887450792407110e-002, 7.3139976491353280e-003) }
  if(nnodes.gq==23) {
    n<-c(5.0000000000000000e-001, 5.6662841214923310e-001, 6.3206784048517251e-001, 6.9515051901514546e-001, 7.5475073892300371e-001, 8.0980493788182306e-001, 8.5933068156597514e-001, 9.0244420080942001e-001, 9.3837617913522076e-001, 9.6648554341300807e-001, 9.8627123560905761e-001, 9.9738466749877608e-001)
    w<-c(6.6827286093053176e-002, 6.6231019702348404e-002, 6.4452861094041150e-002, 6.1524542153364815e-002, 5.7498320111205814e-002, 5.2446045732270824e-002, 4.6457883030017563e-002, 3.9640705888359551e-002, 3.2116210704262994e-002, 2.4018835865542369e-002, 1.5494002928489686e-002, 6.7059297435702412e-003) }
  if(nnodes.gq==24) {
    n<-c(5.3202844643130276e-001, 5.9555943373680820e-001, 6.5752133984808170e-001, 7.1689675381302254e-001, 7.7271073569441984e-001, 8.2404682596848777e-001, 8.7006209578927718e-001, 9.1000099298695147e-001, 9.4320776350220048e-001, 9.6913727600136634e-001, 9.8736427798565474e-001, 9.9759360999851066e-001)
    w<-c(6.3969097673376246e-002, 6.2918728173414318e-002, 6.0835236463901793e-002, 5.7752834026862883e-002, 5.3722135057982914e-002, 4.8809326052057039e-002, 4.3095080765976693e-002, 3.6673240705540205e-002, 2.9649292457718385e-002, 2.2138719408709880e-002, 1.4265694314466934e-002, 6.1706148999928351e-003) }
  if(nnodes.gq==25) {
    n<-c(5.0000000000000000e-001, 5.6143234630535521e-001, 6.2193344186049426e-001, 6.8058615290469393e-001, 7.3650136572285752e-001, 7.8883146512061142e-001, 8.3678318423673415e-001, 8.7962963151867890e-001, 9.1672131438041693e-001, 9.4749599893913761e-001, 9.7148728561448716e-001, 9.8833196072975871e-001, 9.9777848489524912e-001)
    w<-c(6.1588026863357799e-002, 6.1121221495155122e-002, 5.9727881767892461e-002, 5.7429129572855862e-002, 5.4259812237131867e-002, 5.0267974533525363e-002, 4.5514130991481903e-002, 4.0070350167500532e-002, 3.4019166906178545e-002, 2.7452347987917691e-002, 2.0469578350653148e-002, 1.3177493307516108e-002, 5.6968992505125535e-003)
  }
  n1<-1-n
  if(nnodes.gq%%2==0) {
    x<-c(rev(n1),n)
    w<-c(rev(w),w)
  } else {
    x<-c(rev(n1[-1]),n)
    w<-c(rev(w[-1]),w)  
  }
  mw<-nodes<-matrix(0,nrow=nnodes.gq**dim,ncol=dim)
  for(j in 1:dim) {
    nodes[,j]<-rep(rep(x,each=nnodes.gq**(dim-j)),nnodes.gq**(j-1))
    mw[,j]<-rep(rep(w,each=nnodes.gq**(dim-j)),nnodes.gq**(j-1))
  }  
  weights<-apply(mw,1,prod)
  return(list(nodes=nodes,weights=weights))
}

#######################	Conditional means estimates of PSI_i ########################


conddist.saemix<-function(saemixObject,nsamp=1,max.iter=NULL,...) {
# Estimate conditional means and estimates for the individual parameters PSI_i using the MCMC algorithm
# nsamp= number of MCMC samples
# kmax= max nb of iterations
# returns an array 
  N<-saemixObject["data"]["N"]
  nb.parameters<-saemixObject["model"]["nb.parameters"]
  if(is.null(max.iter)) kmax<-sum(saemixObject["options"]$nbiter.saemix)*2 else kmax<-max.iter
# using several Markov chains
  chdat<-new(Class="SaemixRepData",data=saemixObject["data"], nb.chains=nsamp)
  NM<-chdat["NM"]
  IdM<-chdat["dataM"]$IdM
  yM<-chdat["dataM"]$yM
  XM<-chdat["dataM"][,saemixObject["data"]["name.predictors"],drop=FALSE]
  io<-matrix(data=0,nrow=N,ncol=max(saemixObject["data"]["nind.obs"]))
  for(i in 1:N)
    io[i,1:saemixObject["data"]["nind.obs"][i]]<-1
  ioM<-do.call(rbind,rep(list(io),nsamp))
  ind.ioM <- which(t(ioM)!=0)
  DYF<-matrix(data=0,nrow=dim(ioM)[2],ncol=dim(ioM)[1])

  ind.eta<-saemixObject["model"]["indx.omega"]
  nb.etas<-length(ind.eta)
  omega.eta<-saemixObject["results"]["omega"][ind.eta,ind.eta]
  diag.omega<-mydiag(saemixObject["results"]["omega"])
  domega2<-do.call(cbind,rep(list((sqrt(mydiag(omega.eta)))* saemixObject["options"]$rw.ini), nb.etas))
  VK<-rep(c(1:nb.etas),2)
  phi<-array(data=0,dim=c(N,nb.parameters, nsamp))

  domega<-cutoff(mydiag(saemixObject["results"]["omega"][ind.eta,ind.eta]), .Machine$double.eps)
  omega.eta<-saemixObject["results"]["omega"][ind.eta,ind.eta]
  omega.eta<-omega.eta-mydiag(mydiag(saemixObject["results"]["omega"][ind.eta, ind.eta]))+mydiag(domega)
  chol.omega<-chol(omega.eta)
  ares<-saemixObject["results"]["respar"][1]
  bres<-saemixObject["results"]["respar"][2]

# Preparing plots
  if(saemixObject["options"]$displayProgress) {
    plot.opt<-saemixObject["prefs"]
    plot.opt$xlab<-"Iteration"
    plot.opt<-replace.plot.options(plot.opt,...)
    change.ylab<-FALSE
    if(plot.opt$ylab!=saemixObject["prefs"]$ylab & length(plot.opt$which.par)==1) change.ylab<-TRUE
    change.main<-FALSE
    if(plot.opt$main!=saemixObject["prefs"]$main & length(plot.opt$which.par)==1) change.main<-TRUE
    np<-nb.etas
    if(length(plot.opt$mfrow)==0) {
      n1<-round(sqrt(np))
      n2<-ceiling(np/n1)
      if(n1>5 | n2>5) {
        n1<-3
        n2<-4
#      cat("Changing the plot layout\n")
      }
      plot.opt$mfrow<-c(n1,n2)
    }
  }
# Simulation MCMC
# initialisation a phiM=estimation des parametres individuels
  mean.phiM<-do.call(rbind,rep(list(saemixObject["results"]["mean.phi"]),nsamp))
  phiM<-do.call(rbind,rep(list(saemixObject["results"]["cond.mean.phi"]),nsamp))
  etaM<-phiM[,ind.eta]-mean.phiM[,ind.eta]  
  psiM<-transphi(phiM,saemixObject["model"]["transform.par"])
  fpred<-saemixObject["model"]["model"](psiM, IdM, XM)
  if(saemixObject["model"]["error.model"]=="exponential")
     fpred<-log(cutoff(fpred))
  gpred<-cutoff(ares+bres*abs(fpred))
  DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)^2+log(gpred)
  U.y<-colSums(DYF)
  phiMc<-phiM
  
  econd<-sdcond<-array(0,dim=c(nb.parameters,N*kmax,nsamp)) # parametres individuels
  ebar<-sdbar<-array(0,dim=c(nb.parameters,kmax, nsamp)) # moyenne des parametres individuels
  cat("Estimating the conditional mean and variance of the distribution of individual parameters\n")
  k<-1
  while(k<=kmax) { # Set a maximum nb of iterations

  for(u in 1:saemixObject["options"]$nbiter.mcmc[1]) { # 1er noyau
    etaMc<-matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
    phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
    psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
    fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
    if(saemixObject["model"]["error.model"]=="exponential")
      fpred<-log(cutoff(fpred))
    gpred<-cutoff(ares+bres*abs(fpred))
    DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)^2+log(gpred)
    Uc.y<-colSums(DYF)
    deltau<-Uc.y-U.y
    ind<-which(deltau<(-1)*log(runif(NM)))
    etaM[ind,]<-etaMc[ind,]
    U.y[ind]<-Uc.y[ind]
  }

  U.eta<-0.5*rowSums(etaM*(etaM%*%solve(omega.eta)))

# Second stage

  if(saemixObject["options"]$nbiter.mcmc[2]>0) {
    nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
    nrs2<-1
    for (u in 1:saemixObject["options"]$nbiter.mcmc[2]) {
     for(vk2 in 1:nb.etas) {
       etaMc<-etaM
       etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(NM*nrs2), ncol=nrs2)%*%mydiag(domega2[vk2,nrs2],nrow=1) # 2e noyau ? ou 1er noyau+permutation?
       phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
       psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
       fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
       if(saemixObject["model"]["error.model"]=="exponential")
         fpred<-log(cutoff(fpred))
       gpred<-cutoff(ares+bres*abs(fpred))
       DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)**2+log(gpred)
       Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
       Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%solve(omega.eta)))
       deltu<-Uc.y-U.y+Uc.eta-U.eta
       ind<-which(deltu<(-1)*log(runif(NM)))
       etaM[ind,]<-etaMc[ind,]
       U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
       U.eta[ind]<-Uc.eta[ind]
       nbc2[vk2]<-nbc2[vk2]+length(ind)
       nt2[vk2]<-nt2[vk2]+NM
     }
    }
    domega2[,nrs2]<-domega2[,nrs2]*(1+saemixObject["options"]$stepsize.rw* (nbc2/nt2-saemixObject["options"]$proba.mcmc))
  }

  if(saemixObject["options"]$nbiter.mcmc[3]>0) {
    nt2<-nbc2<-matrix(data=0,nrow=nb.etas,ncol=1)
    nrs2<-k%%(nb.etas-1)+2
#    if(is.nan(nrs2)) nrs2<-1 # to deal with case nb.etas=1
    for (u in 1:saemixObject["options"]$nbiter.mcmc[3]) {
      if(nrs2<nb.etas) {
        vk<-c(0,sample(c(1:(nb.etas-1)),nrs2-1))
        nb.iter2<-nb.etas
      } else {
        vk<-0:(nb.etas-1)
        nb.iter2<-1
      }
      for(k2 in 1:nb.iter2) {
        vk2<-VK[k2+vk]
        etaMc<-etaM
        etaMc[,vk2]<-etaM[,vk2]+matrix(rnorm(NM*nrs2), ncol=nrs2)%*%mydiag(domega2[vk2,nrs2])
        phiMc[,ind.eta]<-mean.phiM[,ind.eta]+etaMc
        psiMc<-transphi(phiMc,saemixObject["model"]["transform.par"])
        fpred<-saemixObject["model"]["model"](psiMc, IdM, XM)
        if(saemixObject["model"]["error.model"]=="exponential")
         fpred<-log(cutoff(fpred))
        gpred<-cutoff(ares+bres*abs(fpred))
        DYF[ind.ioM]<-0.5*((yM-fpred)/gpred)**2+log(gpred)
        Uc.y<-colSums(DYF) # Warning: Uc.y, Uc.eta = vecteurs
        Uc.eta<-0.5*rowSums(etaMc*(etaMc%*%solve(omega.eta)))
        deltu<-Uc.y-U.y+Uc.eta-U.eta
        ind<-which(deltu<(-log(runif(NM))))
        etaM[ind,]<-etaMc[ind,]
        U.y[ind]<-Uc.y[ind] # Warning: Uc.y, Uc.eta = vecteurs
        U.eta[ind]<-Uc.eta[ind]
        nbc2[vk2]<-nbc2[vk2]+length(ind)
        nt2[vk2]<-nt2[vk2]+NM
      }
    }
    domega2[,nrs2]<-domega2[,nrs2]*(1+saemixObject["options"]$stepsize.rw* (nbc2/nt2-saemixObject["options"]$proba.mcmc))
  }
  
  phiM[,ind.eta]<-mean.phiM[,ind.eta]+etaM
    if(k==1) {
      eik<-array(t(phiM),dim=c(nb.parameters,N,nsamp))
      varik<-0*eik
      ebar[,k,]<-apply(eik,c(1,3),mean)
      ij<-(k-1)*N
      econd[,(ij+1):(ij+N),]<-eik
      sdcond[,(ij+1):(ij+N),]<-varik
    } else {
      eik1<-eik
      eik<-eik*(k-1)/k+array(t(phiM),dim=c(nb.parameters,N,nsamp))/k
      varik<-varik*(k-1)/k+(array(t(phiM),dim=c(nb.parameters,N,nsamp))**2)/k+ (eik1**2)*(k-1)/k-(eik**2)
      sdik<-sqrt(varik)
      ij<-(k-1)*N
      econd[,(ij+1):(ij+N),]<-eik
      sdcond[,(ij+1):(ij+N),]<-sdik
      ebar[,k,]<-apply(eik,c(1,3),mean)
      sdbar[,k,]<-apply(sdik,c(1,3),mean)
      if(k>=saemixObject["options"]$ipar.lmcmc) {
      ibeg<-max(1,k-saemixObject["options"]$ipar.lmcmc)
      ekmax<-(ebar[,k,]+abs(ebar[,k,])* saemixObject["options"]$ipar.rmcmc)
      ekmin<-(ebar[,k,]-abs(ebar[,k,])*saemixObject["options"]$ipar.rmcmc)
        if(nsamp==1) {
          vec<-apply(ebar[,ibeg:k,],1,max)
          vec2<-apply(ebar[,ibeg:k,],1,min)
          iek<-sum(vec>ekmax)+sum(vec2<ekmin)
          vec<-apply(sdbar[,ibeg:k,],1,max)
          vec2<-apply(sdbar[,ibeg:k,],1,min)
          sdek<-sum(vec>(sdbar[,k,]+abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))+ sum(vec2<(sdbar[,k,]-abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))
        } else {
        vec<-apply(ebar[,ibeg:k,],c(1,3),max)
	vec2<-apply(ebar[,ibeg:k,],c(1,3),min)
        iek<-sum(vec>ekmax)+sum(vec2<ekmin)
	if(ibeg==1) ibeg<-2
        vec<-apply(sdbar[,ibeg:k,],c(1,3),max)
        vec2<-apply(sdbar[,ibeg:k,],c(1,3),min)
        sdek<-sum(vec>(sdbar[,k,]+abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))+ sum(vec2<(sdbar[,k,]-abs(sdbar[,k,])* saemixObject["options"]$ipar.rmcmc))
        }
        if(is.na(sdek)) sdek<-0
        if(iek==0 & sdek==0) break
      }
    }
    if((k%%50)==50) {
      cat(".")
      if(saemixObject["options"]$displayProgress) {
#        par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
        par(mfrow=plot.opt$mfrow)
        limy<-cbind(apply(cbind(apply(ebar[,1:k,],1,min),ekmin),1,min), apply(cbind(apply(ebar[,1:k,],1,max),ekmax),1,max))
        for(ipar in ind.eta) {
        laby<-saemixObject["model"]["name.modpar"][ipar]
	plot(1:k,ebar[ipar,1:k,1],type="n",xlab=plot.opt$xlab, ylab=laby,ylim=limy[ipar,],main=plot.opt$main)
        for(isamp in 1:nsamp) 
	  lines(1:k,ebar[ipar,1:k,isamp],col=plot.opt$col,lty=plot.opt$lty, lwd=plot.opt$lwd)
        abline(h=apply(ebar[,k,],1,mean), col=plot.opt$ablinecol, lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd)
        abline(h=apply(ekmax,1,mean)[ipar], col=plot.opt$ablinecol, lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd)
        abline(h=apply(ekmin,1,mean)[ipar], col=plot.opt$ablinecol, lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd)
        }
      }
    }
    k<-k+1
  }
  cat("\n")
  if(k>=kmax)
    cat("Computing the empirical conditional mean and variance: maximum number of iterations reached without meeting convergence criterion (max.iter=",kmax,")\n") else cat("Convergence achieved in",k,"iterations\n")
  eta.cond<-matrix(0,nrow=dim(etaM)[1],ncol=nb.parameters)
  eta.cond[,ind.eta]<-etaM
  eta.cond<-array(t(eta.cond),dim=c(nb.parameters,N,nsamp))
  eta.cond<-t(apply(eta.cond,c(1,2),mean))
  cond.shrinkage<-100*(1-apply(eta.cond,2, var)/mydiag(saemixObject["results"]["omega"]))
  names(cond.shrinkage)<-paste("Sh.",names(cond.shrinkage),".%",sep="")
  resh.eik<-resh.varik<-array(0,dim=c(nrow=N,ncol=nb.parameters,nsamp))
  for(isamp in 1:nsamp) {
    resh.eik[,,isamp]<-t(eik[,,isamp])
    resh.varik[,,isamp]<-t(varik[,,isamp])
  }

  saemixObject["results"]["cond.mean.phi"]<-apply(resh.eik,c(1,2),mean)
  saemixObject["results"]["cond.var.phi"]<-apply(resh.varik,c(1,2),mean)
  saemixObject["results"]["cond.mean.eta"]<-eta.cond
  saemixObject["results"]["cond.shrinkage"]<-cond.shrinkage
  saemixObject["results"]["phi.samp"]<-resh.eik
  saemixObject["results"]["phi.samp.var"]<-resh.varik

# ECO TODO: verifier ce qui se passe quand un parametre n'est pas estime => faut-il prendre la moyenne des valeurs ou grace a chol.omega qui vaut 0, ca marche ? Verifier d'ailleurs si on a chol.omega quand il manque un parametre (normalement oui)

  return(saemixObject)
}

###########################  Individual MAP estimates 	#############################

map.saemix<-function(saemixObject) {
# Compute the MAP estimates of the individual parameters PSI_i
  i1.omega2<-saemixObject["model"]["indx.omega"]
  iomega.phi1<-solve(saemixObject["results"]["omega"][i1.omega2,i1.omega2])
  id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
  xind<-saemixObject["data"]["data"][,saemixObject["data"]["name.predictors"], drop=FALSE]
  yobs<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  id.list<-unique(id)
  phi.map<-saemixObject["results"]["phi"]
  
  cat("Estimating the individual parameters, please wait a few moments...\n")
  for(i in 1:saemixObject["data"]["N"]) {
    cat(".")
    isuj<-id.list[i]
    xi<-xind[id==isuj,,drop=FALSE]
#    if(is.null(dim(xi))) xi<-matrix(xi,ncol=1)
    yi<-yobs[id==isuj]
    idi<-rep(1,length(yi))
    mean.phi1<-saemixObject["results"]["mean.phi"][i,i1.omega2]
    phii<-saemixObject["results"]["phi"][i,]
    phi1<-phii[i1.omega2]
    phi1.opti<-optim(par=phi1, fn=conditional.distribution, phii=phii,idi=idi,xi=xi,yi=yi,mphi=mean.phi1,idx=i1.omega2,iomega=iomega.phi1, trpar=saemixObject["model"]["transform.par"], model=saemixObject["model"]["model"], pres=saemixObject["results"]["respar"], err=saemixObject["model"]["error.model"])
    phi.map[i,i1.omega2]<-phi1.opti$par
  }
  cat("\n")
  map.psi<-transphi(phi.map,saemixObject["model"]["transform.par"])
  map.psi<-data.frame(id=id.list,map.psi)
  map.phi<-data.frame(id=id.list,phi.map)
  colnames(map.psi)<-c(saemixObject["data"]["name.group"], saemixObject["model"]["name.modpar"])
  saemixObject["results"]["map.psi"]<-map.psi
  saemixObject["results"]["map.phi"]<-map.phi
  return(saemixObject)
}

compute.eta.map<-function(saemixObject) {
# Compute individual estimates of the MAP random effects from the MAP estimates of the parameters
# returns the parameters (psi), newly computed if needs be, the corresponding random effects, and the associated shrinkage
  if(length(saemixObject["results"]["map.psi"])) {
      saemixObject<-map.saemix(saemixObject)
  }
  psi<-saemixObject["results"]["map.psi"][,-c(1)]
  phi<-transpsi(as.matrix(psi),saemixObject["model"]["transform.par"])

# Computing COV again here (no need to include it in results)  
#  COV<-matrix(nrow=dim(saemix.model["Mcovariates"])[1],ncol=0)
  COV<-matrix(nrow=saemixObject["data"]["N"],ncol=0)
  for(j in 1:saemixObject["model"]["nb.parameters"]) {
    jcov<-which(saemixObject["model"]["betaest.model"][,j]==1)
    aj<-as.matrix(saemixObject["model"]["Mcovariates"][,jcov])
    COV<-cbind(COV,aj)
  }
  eta<-phi-COV%*%saemixObject["results"]["MCOV"] 
  shrinkage<-100*(1-apply(eta,2,var)/mydiag(saemixObject["results"]["omega"]))
  names(shrinkage)<-paste("Sh.",names(shrinkage),".%",sep="")
  colnames(eta)<-paste("ETA(",colnames(eta),")",sep="")
  eta<-cbind(id=saemixObject["results"]["map.psi"][,1],eta)
  
  saemixObject["results"]["map.eta"]<-eta
  saemixObject["results"]["map.shrinkage"]<-shrinkage
  
  return(saemixObject)
}

###########################  Fisher Information Matrix 	#############################

fim.saemix<-function(saemixObject) {
# Estimate the Fisher Information Matrix and the s.e. of the estimated parameters  

	saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  xind<-saemix.data["data"][,saemix.data["name.predictors"],drop=FALSE]
  yobs<-saemix.data["data"][,saemix.data["name.response"]]

  covariance.model<-0*saemix.model["covariance.model"]
  diag(covariance.model)<-mydiag(saemix.model["covariance.model"])
  omega<-0*saemix.res["omega"]
  diag(omega)<-mydiag(saemix.res["omega"])
  hat.phi<-saemix.res["cond.mean.phi"]
  nphi<-dim(hat.phi)[2]
  dphi<-cutoff(abs(colMeans(hat.phi))*1e-4,1e-10)
  coefphi<-c(0,-1,1)

  F<-array(data=0,dim=c(saemix.data["ntot.obs"],nphi,length(coefphi)))
  gs<-matrix(0,saemix.data["ntot.obs"],4)

  for (l in 1:length(coefphi)) {
    for (j in 1:nphi) {
        phi<-hat.phi
        phi[,j]<-phi[,j]+coefphi[l]*dphi[j]
        psi<-transphi(phi,saemix.model["transform.par"])
        f <- saemix.model["model"](psi, saemix.data["data"][,"index"],xind)
      	if(saemix.model["error.model"]=='exponential') 
		f<-log(cutoff(f))        
        F[,j,l]<-f
    	}
  }

  ind.covariates<-which(saemix.model["betaest.model"]>0)
  f0<-F[,1,1]
  g0<-cutoff(saemix.res["respar"][1]+saemix.res["respar"][2]*abs(f0))
  DF<-(F[,,3]-F[,,2])/matrix(rep(dphi,each=saemix.data["ntot.obs"]), ncol=length(dphi))/2 #gradient of f 
  z<-matrix(0,saemix.data["ntot.obs"],1)
  j2<-0

	for (i in 1:saemix.data["N"]) {
    j1<-j2+1
    j2<-j2+saemix.data["nind.obs"][i]
    z[j1:j2]<-yobs[j1:j2] - f0[j1:j2] + DF[j1:j2,,drop=FALSE]%*%hat.phi[i,]
  }

# ECO ici modifie car role de covariate.estim pas clair
# covariate.estim=si un parametre (et ses covariables associees) sont estimees ou non
  covariate.estim<-matrix(rep(saemix.model["fixed.estim"], dim(saemix.model["betaest.model"])[1]),byrow=TRUE, ncol=length(saemix.model["fixed.estim"]))*saemix.model["betaest.model"]
  j<-which(saemix.model["betaest.model"]>0)
  ind.fixed.est<-(covariate.estim[j]>0)

# hw=waitbar(1,'Estimating the population parameters (SAEM). Wait...')

  ll.lin<- -0.5*saemix.data["ntot.obs"]*log(2*pi)
  Fmu<-0
  FO<-0
  j2<-0
	for (i in 1:saemix.data["N"]) {
#waitbar(i/N,hw)
    ni<-saemix.data["nind.obs"][i]
    j1<-j2+1
    j2<-j2+ni
    yi<-yobs[j1:j2]
    DFi<-DF[j1:j2,,drop=FALSE]
    f0i<-f0[j1:j2]
    g0i<-g0[j1:j2]
    zi<-z[j1:j2]
    Ai<-kronecker(diag(nphi),as.matrix(saemix.model["Mcovariates"][i,]))
    Ai<-Ai[,ind.covariates,drop=FALSE]
    DFAi<-DFi%*%Ai
    Gi<-DFi%*%omega%*%t(DFi) + mydiag(g0i^2,nrow=ni)  #variance of zi
  
    Gi<-round(Gi*1e10)/1e10
    VD<-try(eigen(Gi))
    if(class(VD)=="try-error") {
      cat("Unable to compute the FIM by linearisation.\n")
      return(saemixObject)
    }
    D<-Re(VD$values)
    V<-Re(VD$vectors)
    IGi<-V%*%mydiag(1/D,nrow=length(D))%*%t(V)
    Dzi<-zi-DFAi%*%saemix.res["betas"]
    
    if (sum(ind.fixed.est)>0) {
        DFAiest<-DFAi[,ind.fixed.est,drop=FALSE]
        Fmu<-Fmu-t(DFAiest)%*%IGi%*%DFAiest
    }
    
    ###########################################################
    OP<-NULL
    for (k in 1:nphi) {
        for (l in 1:nphi) {
            if (covariance.model[k,l]==1) {
                OPkl<-DFi[,k,drop=FALSE]%*%t(DFi[,l,drop=FALSE])
                OP<-cbind(OP,c(OPkl))
        	}
        }
    }
    if (sum(saemix.res["indx.res"]==1)>0) {
        SIi<-2*g0i
        dSIi<-mydiag(SIi)
        OP<-cbind(OP,c(dSIi))
    	}
    if (sum(saemix.res["indx.res"]==2)>0) {
        SIi<-2*f0i%*%g0i
        dSIi<-mydiag(SIi)
	  OP<-cbind(OP,c(dSIi))
    	}
    kl<-0
    FG<-matrix(0,dim(OP)[2],ni*ni)
    for (k in 1:ni)
	{
        for (l in 1:ni)
	{
            FGkl<- -IGi[,k]%*%t(IGi[l,])/2
            kl<-kl+1
            FG[,kl]<-t(OP)%*%c(FGkl)
        }
    }
    FO<-FO+FG%*%OP
    
    ###############################################################
    ll.lin <- ll.lin - 0.5*log(det(Gi)) - 0.5*t(Dzi)%*%IGi%*%Dzi 
  }
#partie precedente verifiee pas a pas avec matlab

  if(saemix.model["error.model"]=='exponential')
    ll.lin<-ll.lin-sum(yobs)

  if (sum(ind.fixed.est)>0) {
    Mparam<-matrix(0,dim(saemix.model["betaest.model"])[1], dim(saemix.model["betaest.model"])[2])
    Mparam[1,]<-saemix.model["transform.par"]
    Mtp<-Mparam[saemix.model["betaest.model"]>0]    
    Mtp<-Mtp[ind.fixed.est]
    dbetas <- dtransphi(saemix.res["betas"][ind.fixed.est],Mtp)
    Mupth<-mydiag(1/dbetas,nrow=length(dbetas))
    Fth<-t(Mupth)%*%Fmu%*%Mupth
    Cth<-try(solve(-Fth))
    if(class(Cth)=="try-error") {
      cat("Error computing the Fisher Information Matrix: singular system.\n")
      Cth<-NA*Fth
    }
  } else {
    Cth<-NULL
  }

  fim<-rbind(cbind(Fth,matrix(0,dim(Fth)[1],dim(FO)[2])), cbind(matrix(0,dim(FO)[1],dim(Fth)[2]),FO)) 

  sTHest<-sqrt(mydiag(Cth))
#sTH<-matrix(0,1,length(saemix.res["betas"]))
  sTH<-rep(0,length(saemix.res["betas"]))
  sTH[ind.fixed.est]<-sTHest
  se.fixed<-sTH

  CO<-try(solve(-FO))
    if(class(CO)=="try-error") {
      CO<-NA*FO
      cat("Error computing the Fisher Information Matrix: singular system.\n")
  }
  sO<-sqrt(mydiag(CO))
  nb.omega2<-length(saemix.model["indx.omega"])
  se.omega<-matrix(0,nphi,1)
  se.omega[saemix.model["indx.omega"]]<-sO[1:nb.omega2]
  se.res<-matrix(0,2,1)
  se.res[saemix.res["indx.res"]]<-sO[(nb.omega2+1):length(sO)]    

# ECO TODO : pourquoi negatif ??? FIM = -fim calculee ici ?
  saemix.res["se.fixed"]<-se.fixed
  saemix.res["se.omega"]<-c(se.omega)
  saemix.res["se.respar"]<-c(se.res)
  saemix.res["ll.lin"]<-c(ll.lin )
  saemix.res["fim"]<-fim
  saemix.res["aic.lin"]<-(-2)*saemix.res["ll.lin"]+ 2*saemix.res["npar.est"]
  saemix.res["bic.lin"]<-(-2)*saemix.res["ll.lin"]+ log(saemix.data["N"])*saemix.res["npar.est"]

##################################
#delete(hw)
  saemixObject["results"]<-saemix.res
  return(saemixObject)
}

#######################	Model simulations and residuals ########################

simul.saemix<-function(saemixObject,nsim=saemixObject["options"]$nb.sim, predictions=TRUE,res.var=TRUE,uncertainty=FALSE) {
# Simulate individual parameters from the population distribution
# predictions: if TRUE, use the parameters to predict observations
# res.var: if TRUE, add residual error to the predictions to obtain simulated data
# uncertainty: if TRUE, add uncertainty when simulating (not implemented yet)
  saemix.model<-saemixObject["model"]
  saemix.data<-saemixObject["data"]
  saemix.res<-saemixObject["results"]
  xind<-saemix.data["data"][,saemix.data["name.predictors"],drop=FALSE]

  N<-saemix.data["N"]
  ind.eta<-saemix.model["indx.omega"]
  nb.etas<-length(ind.eta)
  NM <- N*nsim  
  domega<-cutoff(mydiag(saemix.res["omega"][ind.eta, ind.eta]),.Machine$double.eps)
  omega.eta<-saemix.res["omega"][ind.eta,ind.eta]
  omega.eta<-omega.eta-mydiag(mydiag(saemix.res["omega"][ind.eta,ind.eta]))+mydiag(domega)
  chol.omega<-chol(omega.eta)

  phiM<-mean.phiM<-do.call(rbind,rep(list(saemix.res["mean.phi"]),nsim))
  etaM<-matrix(rnorm(NM*nb.etas),ncol=nb.etas)%*%chol.omega
  phiM[,ind.eta]<-mean.phiM[,ind.eta]+etaM
  psiM<-transphi(phiM,saemix.model["transform.par"])

  if(predictions) {
    index<-rep(1:N,times=saemix.data["nind.obs"])
    IdM<-kronecker(c(0:(nsim-1)),rep(N,saemix.data["ntot.obs"]))+rep(index,nsim)
    XM<-do.call(rbind,rep(list(xind),nsim))
    ares<-saemix.res["respar"][1]
    bres<-saemix.res["respar"][2]
    sim.pred<-sim.data<-NULL
      fpred<-saemix.model["model"](psiM, IdM, XM)
      sim.pred<-fpred
      if(res.var) {
        if(saemix.model["error.model"]=="exponential")
          fpred<-log(cutoff(fpred))
        gpred<-ares+bres*abs(fpred)
        eps<-rnorm(length(fpred))
        sim.data<-fpred+gpred*eps
      }
  } else {
    sim.pred<-sim.data<-IdM<-c()
  }
  sim.psi<-data.frame(id=rep(unique(saemix.data["data"][, saemix.data["name.group"]]),nsim),psiM)
  colnames(sim.psi)<-c(saemix.data["name.group"],saemix.model["name.modpar"])
  datasim<-data.frame(idsim=rep(index,nsim),irep=rep(1:nsim, each=saemix.data["ntot.obs"]),ypred=sim.pred,ysim=sim.data)
  ysim<-new(Class="SaemixSimData",saemix.data,datasim)
  ysim["sim.psi"]<-sim.psi
  saemixObject["sim.data"]<-ysim
  
  return(saemixObject)
}

compute.sres<-function(saemixObject) {
# Compute standardised residuals (WRES, npd and npde) using simulations
# saemix.options$nb.sim simulated datasets used to compute npd, npde, and VPC
# saemix.options$nb.simpred simulated datasets used to compute ypred and WRES ? for the moment saemix.options$nb.sim used for both
  nsim<-saemixObject["options"]$nb.sim
  if(length(saemixObject["sim.data"]["N"])==0 || saemixObject["sim.data"]["nsim"]!=nsim) {
	  cat("Simulating data using nsim =",nsim,"simulated datasets\n")
	  saemixObject<-simul.saemix(saemixObject,nsim)
  }
# ECO TODO: maybe here be more clever and use simulations if available (adding some if not enough, truncating if too much ?)  
  ysim<-saemixObject["sim.data"]["datasim"]$ysim
  idsim<-saemixObject["sim.data"]["datasim"]$idsim
  idy<-saemixObject["data"]["data"][,"index"]
  yobsall<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  ypredall<-pd<-npde<-wres<-c()
#  pde<-c()
  cat("Computing WRES and npde ")
  for(isuj in 1:saemixObject["data"]["N"]) {
    if(isuj%%10==1) cat(".")
    ysimi<-matrix(ysim[idsim==isuj],ncol=nsim)
#    ysimi.pred<-ysimi[,1:saemixObject["options"]$nb.simpred]
    yobs<-yobsall[idy==isuj]
    tcomp<-apply(cbind(ysimi,yobs),2,"<",yobs)
    if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
    pdsuj<-rowMeans(tcomp)
#      pdsuj[pdsuj==0]<-1/nsim
#      pdsuj[pdsuj==1]<-1-1/nsim
    pdsuj<-pdsuj+0.5/(nsim+1)
# pdsuj=0 pour yobs<min(tabobs$yobs)
# pdsuj=1 : jamais
# dc dist va de 0 ? 1-1/(nsim+1)
# dc dist+0.5/(nsim+1) va de 0.5/(nsim+1) ? 1-0.5/(nsim+1)
# est-ce que ?a recentre ma distribution ? ECO TODO CHECK
    pd<-c(pd,pdsuj)
    ypred<-rowMeans(ysimi)
    ypredall<-c(ypredall,ypred)
    xerr<-0
    if(length(yobs)==1) {
      npde<-c(npde,qnorm(pdsuj))
      wres<-c(wres,(yobs-ypred)/sd(t(ysimi)))
    } else {
# Decorrelation
    vi<-cov(t(ysimi))
    xmat<-try(chol(vi))
    if(is.numeric(xmat)) {
      sqvi<-try(solve(xmat))
      if(!is.numeric(sqvi)) 
        xerr<-2
      } else xerr<-1
    if(xerr==0) {
    #decorrelation of the simulations
      decsim<-t(sqvi)%*%(ysimi-ypred)
      decobs<-t(sqvi)%*%(yobs-ypred)
      wres<-c(wres,decobs)
#    ydsim<-c(decsim)
#    ydobs<-decobs
    #Computing the pde
      tcomp<-apply(cbind(decsim,decobs),2,"<",decobs)
      if(!is.matrix(tcomp)) tcomp<-t(as.matrix(tcomp))
      pdesuj<-rowMeans(tcomp)
      pdesuj<-pdesuj+0.5/(nsim+1)
#      pde<-c(pde,pdesuj)
      npde<-c(npde,qnorm(pdesuj))
    } else {
      npde<-c(npde,rep(NA,length(yobs)))
      wres<-c(wres,rep(NA,length(yobs)))
    }
    }
  }
  cat("\n")
  saemixObject["results"]["npde"]<-npde
  saemixObject["results"]["wres"]<-wres
  saemixObject["results"]["ypred"]<-ypredall # = [ E_i(f(theta_i)) ]
  saemixObject["results"]["pd"]<-pd
  if(length(saemixObject["results"]["predictions"])==0) saemixObject["results"]["predictions"]<-data.frame(ypred=ypredall,wres=wres,pd=pd,npde=npde) else {
  	saemixObject["results"]["predictions"]$ypred<-ypred
  	saemixObject["results"]["predictions"]$wres<-wres
  	saemixObject["results"]["predictions"]$npde<-npde
  	saemixObject["results"]["predictions"]$pd<-pd
  }
#  return(list(ypred=ypredall,pd=pd,npd=npd,wres=wres,sim.data=ysim, sim.pred=x$sim.pred))
  return(saemixObject)
}

###########################	Computational fcts	#############################
# Redefining diag function, too many problems with the R version
mydiag <- function (x = 1, nrow, ncol) {
	if (is.matrix(x)) {
		if (nargs() > 1L) 
			stop("'nrow' or 'ncol' cannot be specified when 'x' is a matrix")
		if ((m <- min(dim(x))) == 0L) 
			return(vector(typeof(x), 0L))
		y <- c(x)[1L + 0L:(m - 1L) * (dim(x)[1L] + 1L)]
		nms <- dimnames(x)
		if (is.list(nms) && !any(sapply(nms, is.null)) && identical((nm <- nms[[1L]][seq_len(m)]), 
																																nms[[2L]][seq_len(m)])) 
			names(y) <- nm
		return(y)
	}
	if (is.array(x) && length(dim(x)) != 1L) 
		stop("'x' is an array, but not 1D.")
	if (missing(x)) 
		n <- nrow
	else n <- length(x)
	if (!missing(nrow)) 
		n <- nrow
	if (missing(ncol)) 
		ncol <- n
	p <- ncol
	y <- array(0, c(n, p))
	if ((m <- min(n, p)) > 0L) 
		y[1L + 0L:(m - 1L) * (n + 1L)] <- x
	y
}

cutoff<-function(x,seuil=.Machine$double.xmin) {x[x<seuil]<-seuil; return(x)}
cutoff.max<-function(x) max(x,.Machine$double.xmin)
cutoff.eps<-function(x) max(x,.Machine$double.eps)
cutoff.res<-function(x,ares,bres) max(ares+bres*abs(x),.Machine$double.xmin)

# Inverse of the normal cumulative distribution fct: using erfcinv from ?pnorm
norminv<-function(x,mu=0,sigma=1)  mu-sigma*qnorm(x,lower.tail=FALSE)

# Truncated gaussian distribution (verifie par rapport a definition de erf/matlab)
normcdf<-function(x,mu=0,sigma=1)
  cutoff(pnorm(-(x-mu)/sigma,lower.tail=FALSE),1e-30)

error<-function(ab,y,f) {
  g=abs(ab[1]+ab[2]*f)
  e=sum(((y-f)/g)**2+2*log(g))
  return(e)
}

transpsi<-function(psi,tr) {
  phi<-psi
#  if(is.null(dim(psi))) phi<-as.matrix(t(phi),nrow=1)
# ECO TODO: pourquoi ce test ??? Dans le cas ou psi est un vecteur ?
  i1<-which(tr==1) # log-normal
  phi[,i1]<-log(phi[,i1])
  i2<-which(tr==2) # probit
  phi[,i2]<-norminv(phi[,i2])
  i3<-which(tr==3) # logit
  phi[,i3]<-log(phi[,i3]/(1-phi[,i3]))
  if(is.null(dim(psi))) phi<-c(phi)
  return(phi)
}

transphi<-function(phi,tr) {
  psi<-phi
#  if(is.null(dim(psi))) psi<-as.matrix(t(psi),nrow=1)
  i1<-which(tr==1) # log-normal
  psi[,i1]<-exp(psi[,i1])
  i2<-which(tr==2) # probit
  psi[,i2]<-normcdf(psi[,i2])
  i3<-which(tr==3) # logit
  psi[,i3]<-1/(1+exp(-psi[,i3]))
  if(is.null(dim(phi))) psi<-c(psi)
  return(psi)
}

derivphi<-function(phi,tr) {
# Fonction calculant la derivee de h pour tracer la distribution des parametres
  psi<-phi # identite
  i1<-which(tr==1) # log-normal
  psi[,i1]<-1/exp(phi[,i1])
  i2<-which(tr==2) # probit
  psi[,i2]<-1/(sqrt(2*pi))*exp(-(phi[,i2]**2)/2)
  i3<-which(tr==3) # logit
  psi[,i3]<-2+exp(phi[,i3])+exp(-phi[,i3])
  if(is.null(dim(phi))) psi<-c(psi)
  return(psi)
}

dtransphi<-function(phi,tr) {
  psi<-phi
  if(is.null(dim(phi))) {
     dpsi<-as.matrix(t(rep(1,length(phi))))
     psi<-as.matrix(t(phi),nrow=1)
  } else 
    dpsi<-matrix(1,dim(phi)[1],dim(phi)[2])
  i1<-which(tr==1) # log-normal
  dpsi[,i1]<-exp(psi[,i1])
  i2<-which(tr==2) # probit
  dpsi[,i2]<-1/dnorm(qnorm(dpsi[,i2]))   # derivee de la fonction probit, dqnorm <- function(p) 1/dnorm(qnorm(p))
  i3<-which(tr==3) # logit
  dpsi[,i3]<-1/(2+exp(-psi[,i3])+exp(psi[,i3]))
  if(is.null(dim(phi))) dpsi<-c(dpsi)
  return(dpsi)
}

compute.Uy<-function(b0,phiM,ares,bres,args,DYF) {
# Attention, DYF variable locale non modifiee en dehors
  args$MCOV0[args$j0.covariate]<-b0
  phi0<-args$COV0 %*% args$MCOV0
  phiM[,args$i0.omega2]<-do.call(rbind,rep(list(phi0),args$nmc))
  psiM<-transphi(phiM,args$transform.par)
  fpred<-args$structural.model(psiM,args$IdM,args$XM)
  if(args$error.model=="exponential")
     fpred<-log(cutoff(fpred))
  gpred<-cutoff(ares+bres*abs(fpred))
  DYF[args$ind.ioM]<-0.5*((args$yM-fpred)/gpred)**2+log(gpred)
  U<-sum(DYF)
  return(U)
}

conditional.distribution<-function(phi1,phii,idi,xi,yi,mphi,idx,iomega,trpar,model,pres,err) {
  phii[idx]<-phi1
  psii<-transphi(matrix(phii,nrow=1),trpar)
  if(is.null(dim(psii))) psii<-matrix(psii,nrow=1)
  fi<-model(psii,idi,xi)
  if(err=="exponential")
    fi<-log(cutoff(fi))
  gi<-cutoff((pres[1]+pres[2]*abs(fi)))
  Uy<-sum(0.5*((yi-fi)/gi)**2+log(gi))
  dphi<-phi1-mphi
  Uphi<-0.5*sum(dphi*(dphi%*%iomega))
  return(Uy+Uphi)
}

trnd.mlx<-function(v,n,m) {
  r<-rnorm(n*m)*sqrt(v/2/gammarnd.mlx(v/2,n,m))
  return(r=matrix(r,nrow=n,ncol=m))
}

gammarnd.mlx<-function(a,n,m) {
  nm<-n*m
  y0 <- log(a)-1/sqrt(a)
  c <- a - exp(y0)
  b <- ceiling(nm*(1.7 + 0.6*(a<2)))
  y <- log(runif(b))*sign(runif(b)-0.5)/c + log(a)
  f <- a*y-exp(y) - (a*y0 - exp(y0))
  g <- c*(abs((y0-log(a))) - abs(y-log(a)))
  reject <- ((log(runif(b)) + g) > f)
  y<-y[!reject]
  if(length(y)>=nm) x<-exp(y[1:nm]) else 
    x<-c(exp(y),gammarnd.mlx(a,(nm-length(y)),1))
#  x<-matrix(x,nrow=n,ncol=m) # not useful ?
  return(x)
}

tpdf.mlx<-function(x,v) {
# TPDF_MLX  Probability density function for Student's T distribution

    term<-exp(lgamma((v + 1) / 2) - lgamma(v/2))
    return(term/(sqrt(v*pi)*(1+(x**2)/v)**((v+1)/2)))
}

###########################	Functions for npde	#############################

kurtosis<-function (x) 
{
#from Snedecor and Cochran, p 80
    x<-x[!is.na(x)]
    m4<-sum((x - mean(x))^4)
    m2<-sum((x - mean(x))^2)
    kurt<-m4*length(x)/(m2**2)-3
    return(kurtosis=kurt)
}
skewness<-function (x) 
{
#from Snedecor and Cochran, p 79
    x<-x[!is.na(x)]
    m3<-sum((x - mean(x))^3)
    m2<-sum((x - mean(x))^2)
    skew<-m3/(m2*sqrt(m2/length(x)))
    return(skewness=skew)
}
   
testnpde<-function(npde) 
{
    cat("---------------------------------------------\n")
    cat("Distribution of npde:\n")
    sev<-var(npde)*sqrt(2/(length(npde)-1))
    sem<-sd(npde)/sqrt(length(npde))
    cat("           mean=",format(mean(npde),digits=4),"  (SE=",format(sem,digits=2),")\n")
    cat("       variance=",format(var(npde),digits=4),"  (SE=",format(sev,digits=2),")\n")
    cat("       skewness=",format(skewness(npde),digits=4),"\n")
    cat("       kurtosis=",format(kurtosis(npde),digits=4),"\n")
    cat("---------------------------------------------\n\n")
    myres<-rep(0,4)
    y<-wilcox.test(npde)
    myres[1]<-y$p.val
    y<-shapiro.test(npde)
    myres[3]<-y$p.val

    # test de variance pour 1 ?chantillon
    # chi=s2*(n-1)/sigma0 et test de H0={s=sigma0} vs chi2 ? n-1 df
    semp<-sd(npde)
    n1<-length(npde)
    chi<-(semp**2)*(n1-1)
    y<-2*min(pchisq(chi,n1-1),1-pchisq(chi,n1-1))
    myres[2]<-y
    xcal<-3*min(myres[1:3])
    myres[4]<-min(1,xcal)
    names(myres)<-c("  Wilcoxon signed rank test ","  Fisher variance test      ",
    "  SW test of normality      ","Global adjusted p-value     ")
    cat("Statistical tests\n")
    for(i in 1:4) {
      cat(names(myres)[i],": ")
      #if (myres[i]<1) 
      cat(format(myres[i],digits=3)) 
      #else cat(myres[i])
      if(as.numeric(myres[i])<0.1 & as.numeric(myres[i])>=0.05) cat(" .")
      if(as.numeric(myres[i])<0.05) cat(" *")
      if(as.numeric(myres[i])<0.01) cat("*")
      if(as.numeric(myres[i])<0.001) cat("*")
      cat("\n")
    }
      cat("---\n")
      cat("Signif. codes: '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 \n")
    cat("---------------------------------------------\n")
    return(myres)
}

#######################	   Default plot options (list)	 ########################

saemix.plot.setoptions<-function(saemixObject) {
# setting default plot options
  plot.opt<-list(
# Options for plot types
    ilist=c(1:saemixObject["data"]["N"]),
    level=0:1,
    smooth=FALSE,
    line.smooth="s",
    indiv.par="map",			# type of individual parameters
    which.par="all",			# which parameters to plot 
    which.cov="all",			# which covariates to plot 
    which.resplot=c("res.vs.x","res.vs.pred","dist.qqplot","dist.hist"), # which type of residual plots
    which.pres=c("wres","npde"),	# which population weighted residuals
    which.poppred=c("ppred"),		# which population predictions to use (ypred=E(f(theta_i)), ppred=f(population parameters))
    indiv.histo=FALSE,			# whether to include an histogram of estimated individual parameters
    cov.value=rep(NA,length(saemixObject["model"]["name.cov"])),
# General graphical options
    new=TRUE,				# whether a new page should be called
    ask=FALSE,				# whether the program should ask before creating a new page
    interactive=FALSE,			# whether the user should be prompted before computing predictions or performing simulations for VPC, npde and wres
    mfrow=c(),				# page layout (if empty, defaults to the default layout for each graph type)
    main="",				# title
    xlab="",
    ylab="",
    col="black",
    pch=20,
    lty=1,
    lwd=1,
    xlim=c(),
    ylim=c(),
    xlog=FALSE,
    ylog=FALSE,
    type="b",
    cex=1,
    cex.axis=1,
    cex.lab=1,
    cex.main=1,
    obs.pch=20,
    pred.pch=20,
    obs.col="black",
    obs.lty=1,
    obs.lwd=1,
    obs.pch=20,
    ipred.col="black",
    ipred.lty=2,
    ipred.lwd=1,
    ipred.pch=20,
    ppred.col="black",
    ppred.lty=3,
    ppred.lwd=1,
    ppred.pch=20,
    pcol="black",
    lcol="black",
    fillcol="lightblue1",
    ablinecol="DarkRed",
    ablinelty=2,
    ablinelwd=2,
# 
    range=3,
    col.fillmed="pink",
    col.fillpi="slategray1",
    col.lmed="indianred4",
    col.lpi="slategray4",
    col.pobs="steelblue4",
    col.lobs="steelblue4",
    lty.lmed=2,
    lty.lpi=2,
    lwd.lmed=2,
    lwd.lpi=1,
    lwd.lobs=2,
    lty.lobs=1,
# Options for VPC plot
    vpc.method="equal",			# method (one of "equal"=same nb of points in each interval, "width"=equally spaced intervals (on the log-scale if xlog=TRUE), "user"=user-defined breaks, "optimal"=Marc's optimal binning algorithm); for "user", the breaks must be specified in vpc.breaks (otherwise defaults back to "equal"), while for the other methods the number of bins must be specified in vpc.bin
    vpc.breaks=NULL,			# user-defined breaks
    vpc.bin=10,				# nb of bins
    vpc.beta=0.2,			# value of beta used to compute the variance-based criterion (Jopt,beta(I)) in the clustering algorithm
    vpc.lambda=0.3,			# value of lambda used in the penalised criterion to select the number of bins (if vpc.bin=NULL)
    vpc.interval=0.95,
    vpc.pi=TRUE,
    vpc.obs=TRUE)
    
     if(is.null(plot.opt$name.X)) {
        if(length(saemixObject["data"]["name.X"])>0) plot.opt$name.X<-saemixObject["data"]["name.X"] else plot.opt$name.X<-saemixObject["data"]["name.predictors"][1]
    }
    plot.opt$xlab<-paste(plot.opt$name.X," (",saemixObject["data"]["units"]$x,")", sep="")
     if(length(saemixObject["data"]["name.response"])>0)
    plot.opt$ylab<-paste(saemixObject["data"]["name.response"]," (", saemixObject["data"]["units"]$y,")",sep="")
   return(plot.opt)
}

#################    Function to supersede default plot options	 ##################

replace.plot.options<-function(plot.opt,...) {
  args1<-match.call(expand.dots=TRUE)
  if(length(args1)>2) {
# General arguments: col, pch
    i1<-match("col",names(args1))
    if(!is.na(i1)) {
      plot.opt$col<-eval(args1[[i1]])
      plot.opt$obs.col<-eval(args1[[i1]])
      plot.opt$ipred.col<-eval(args1[[i1]])
      plot.opt$ppred.col<-eval(args1[[i1]])
      plot.opt$pcol<-eval(args1[[i1]])
      plot.opt$lcol<-eval(args1[[i1]])
    }
    i1<-match("pch",names(args1))
    if(!is.na(i1)) {
      plot.opt$pch<-eval(args1[[i1]])
      plot.opt$obs.pch<-eval(args1[[i1]])
      plot.opt$ipred.pch<-eval(args1[[i1]])
      plot.opt$ppred.pch<-eval(args1[[i1]])
    }
# Other arguments
    for(i in 3:length(args1)) {
      if(match(names(args1)[i],names(plot.opt),nomatch=0)>0)    
#    plot.opt[[names(args1)[i]]]<-args1[[i]] else {
    plot.opt[[names(args1)[i]]]<-eval(args1[[i]]) else {
      if(names(args1)[i]!="plot.type") cat("Argument",names(args1)[i],"not available, check spelling.\n")
    }
   }
  }
  return(plot.opt)
}

#####################################################################################
###########################		Plots		#############################
#####################################################################################
###############################	   Wrapper functions  #############################

saemix.plot.select<-function(saemixObject,data=FALSE,convergence=FALSE, likelihood=FALSE,individual.fit=FALSE,population.fit=FALSE,both.fit=FALSE, observations.vs.predictions=FALSE,residuals.scatter=FALSE, residuals.distribution=FALSE,random.effects=FALSE,correlations=FALSE, parameters.vs.covariates=FALSE,randeff.vs.covariates=FALSE, marginal.distribution=FALSE,vpc=FALSE,npde=FALSE,...) {
# Function selecting which plots are to be drawn
  namObj<-deparse(substitute(saemixObject))
  interactive<-saemixObject["prefs"]$interactive
  boolsim<-boolpred<-boolres<-FALSE
  if(vpc) {
    if(length(saemixObject["sim.data"]["nsim"])==0) boolsim<-TRUE
  }
  if(individual.fit | both.fit | observations.vs.predictions) {
    if(length(saemixObject["results"]["ipred"])==0) boolpred<-TRUE
  }
  if(population.fit | both.fit | observations.vs.predictions) {
    if(saemixObject["prefs"]$which.poppred=="ppred" & length(saemixObject["results"]["ppred"])==0) boolpred<-TRUE
    if(saemixObject["prefs"]$which.poppred=="ypred" & length(saemixObject["results"]["ypred"])==0) boolres<-TRUE
  }
  if(residuals.scatter | residuals.distribution) {
      if(length(saemixObject["results"]["npde"])==0) boolres<-TRUE
  }
  if(boolsim & !boolres & interactive) {
    cok<-readline(prompt="Simulations will be performed to obtain residuals, VPC and npde. This might take a while, proceed ? (y/Y) [default=yes] ")
    if(!cok %in% c("y","Y","yes","")) boolsim<-FALSE 
  }
  if(boolres & interactive) {
    cok<-readline(prompt="Simulations will be performed to obtain residuals, VPC and npde. This might take a while, proceed ? (y/Y) [default=yes] ")
    if(!cok %in% c("y","Y","yes","")) boolres<-FALSE 
  }
  if(boolpred & interactive) {
    cok<-readline(prompt="Computations will be performed to obtain model predictions, proceed ? (y/Y) [default=yes] ")
    if(!cok %in% c("y","Y","yes","")) boolpred<-FALSE 
  }
  if(boolsim & !boolres) {
    saemixObject<-simul.saemix(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  if(boolpred) {
    saemixObject<-predict(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
    if(boolres) {
    saemixObject<-compute.sres(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  if(parameters.vs.covariates) {
    if(length(saemixObject["results"]["map.psi"])==0) {
      saemixObject<-map.saemix(saemixObject)
      assign(namObj,saemixObject,envir=parent.frame())
    }
  }
# ECO TODO: replace with partial matching
  if(data) plot(saemixObject,plot.type="data",...)
  if(convergence) plot(saemixObject,plot.type="convergence",...)
  if(likelihood) plot(saemixObject,plot.type="likelihood",...)
  if(observations.vs.predictions) plot(saemixObject,plot.type="observations.vs.predictions", ...)
  if(individual.fit) plot(saemixObject,plot.type="individual.fit",...)
  if(population.fit) plot(saemixObject,plot.type="population.fit",...)
  if(both.fit) plot(saemixObject,plot.type="both.fit",...)
  if(residuals.scatter) plot(saemixObject,plot.type="residuals.scatter",...)
  if(residuals.distribution) plot(saemixObject,plot.type="residuals.distribution",...)
  if(random.effects) plot(saemixObject,plot.type="random.effects",...)
  if(correlations) plot(saemixObject,plot.type="correlations",...)
  if(parameters.vs.covariates) plot(saemixObject,plot.type="parameters.vs.covariates", ...)
  if(randeff.vs.covariates) plot(saemixObject,plot.type="randeff.vs.covariates",...)
  if(marginal.distribution) plot(saemixObject,plot.type="marginal.distribution",...)
  if(vpc) plot(saemixObject,plot.type="vpc",...)
  if(npde) plot(saemixObject,plot.type="npde",...)
}

#### Meta-niveau
default.saemix.plots<-function(saemixObject,...) {
# When plot(saemixObject) is called without plot.type  
  namObj<-deparse(substitute(saemixObject))
  if(length(saemixObject["results"]["ipred"])==0) {
    saemixObject<-predict(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  if(length(saemixObject["results"]["npde"])==0) {
    saemixObject<-compute.sres(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  saemix.plot.select(saemixObject,data=TRUE,convergence=TRUE,likelihood=TRUE, observations.vs.predictions=TRUE,residuals.scatter=TRUE, residuals.distribution=TRUE,random.effects=TRUE,correlations=TRUE, marginal.distribution=TRUE,vpc=TRUE,...)
}

basic.gof<-function(saemixObject,...) {
# Basic goodness of fit plots
  cat("Now producing basic goodness of fit plots\n")
  namObj<-deparse(substitute(saemixObject))
  if(length(saemixObject["results"]["ipred"])==0) {
    saemixObject<-predict(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  saemix.plot.select(saemixObject,convergence=TRUE,likelihood=TRUE, observations.vs.predictions=TRUE, ...)
}

advanced.gof<-function(saemixObject,...) {
# Advanced goodness of fit plots
  cat("Now producing advanced goodness of fit plots\n")
  namObj<-deparse(substitute(saemixObject))
  if(length(saemixObject["results"]["ipred"])==0) {
    saemixObject<-predict(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  if(length(saemixObject["results"]["npde"])==0) {
    saemixObject<-compute.sres(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  saemix.plot.select(saemixObject,residuals.scatter=TRUE,residuals.distribution=TRUE, vpc=TRUE,...)
}

individual.fits<-function(saemixObject,...) {
# Individual plots
  cat("Now producing plots of individual fits\n")
  namObj<-deparse(substitute(saemixObject))
  if(length(saemixObject["results"]["ipred"])==0) {
    saemixObject<-predict(saemixObject)
    assign(namObj,saemixObject,envir=parent.frame())
  }
  plot(saemixObject,plot.type="individual.fit",...)
}

covariate.fits<-function(saemixObject,which="parameters",...) {
# Parameters or random effects versus covariates
  if(which=="parameters") {
    cat("Now producing plots of parameters versus covariates\n")
    plot(saemixObject,plot.type="parameters.vs.covariates",...)
  } else {
    cat("Now producing plots of random effects versus covariates\n")
    plot(saemixObject,plot.type="randeff.vs.covariates",...)
  }
}

###############################	   	Data	 #################################

# ECO FINISH THIS ONE (redo without using data part of object)
saemix.plot.data<-function(saemixObject,...) {
# Plot of the data as spaghetti plot
# options: change data point, line type, line color, lines plotted or not, points plotted or not...
  plot.opt<-saemixObject["prefs"]
  plot.opt$new<-TRUE
  plot.opt$plot.type<-"l"
  plot.opt<-replace.plot.options(plot.opt,...)
  if(plot.opt$new) {
    mfrow=c(1,1)
    if(length(plot.opt$mfrow)>0) mfrow<-plot.opt$mfrow
  par(mfrow=mfrow,ask=plot.opt$ask)
  }
  plot(saemixObject["data"],plot.type=plot.opt$plot.type,...)
}

#######################	   Convergence plots & LL	 ########################

saemix.plot.convergence<-function(saemixObject,niter=0,...) {
# Convergence plots for all the fixed effects, random effects and residual variability
  plot.opt<-saemixObject["prefs"]
  plot.opt$xlab<-"Iteration"
  plot.opt<-replace.plot.options(plot.opt,...)
  change.ylab<-FALSE
  if(plot.opt$ylab!=saemixObject["prefs"]$ylab & length(plot.opt$which.par)==1) change.ylab<-TRUE
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main & length(plot.opt$which.par)==1) change.main<-TRUE
  K<-dim(saemixObject["results"]["allpar"])[1]
  if(niter==0) niter<-K
  if(plot.opt$which.par[1]=="all")
     np<-dim(saemixObject["results"]["allpar"])[2] else  
     np<-length(plot.opt$which.par)
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)==0) {
    n1<-round(sqrt(np))
    n2<-ceiling(np/n1)
    if(n1>5 | n2>5) {
      n1<-3
      n2<-4
#      cat("Changing the plot layout\n")
    }
    par(mfrow=c(n1,n2),ask=plot.opt$ask)
  } else par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
  }
  if(plot.opt$which.par[1]=="all") { # default convergence plot, for all parameters
    for(j in 1:np) {
      laby<-"" #colnames(saemixObject["results"]["allpar"])[j]
      maintit<-colnames(saemixObject["results"]["allpar"])[j]
      plot(1:niter,saemixObject["results"]["allpar"][1:niter,j],type="l", xlab=plot.opt$xlab,ylab=laby, main=maintit,col=plot.opt$col,lty=plot.opt$lty, lwd=plot.opt$lwd,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
      abline(v=saemixObject["options"]$nbiter.saemix[1],col=plot.opt$ablinecol, lwd=plot.opt$ablinelwd)
    }
  } else {
      for(ipar in 1:length(plot.opt$which.par)) {
      j<-as.integer(plot.opt$which.par[ipar])
      if(is.na(j)) j<-which(colnames(saemixObject["results"]["allpar"])== plot.opt$which.par[ipar])
      if(length(j)>0) {
        laby<-""
        maintit<-colnames(saemixObject["results"]["allpar"])[j]
        if(change.ylab) laby<-plot.opt$ylab
        if(change.main) maintit<-plot.opt$main
        plot(1:niter,saemixObject["results"]["allpar"][1:niter,j],type="l", xlab=plot.opt$xlab,ylab=laby,main=maintit,col=plot.opt$col,lty=plot.opt$lty, lwd=plot.opt$lwd,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
      abline(v=saemixObject["options"]$nbiter.saemix[1],col=plot.opt$ablinecol, lwd=plot.opt$ablinelwd)
    }}
  }
}

saemix.plot.llis<-function(saemixObject,...) {
# Plot of the evolution of the log-likelihood by importance sampling
    plot.opt<-saemixObject["prefs"]
    plot.opt$main<-"-2xLL by Importance Sampling"
    plot.opt$xlab<-"Iteration"
    plot.opt$ylab<-"-2 x LL"
    plot.opt<-replace.plot.options(plot.opt,...)
    MM<-100
    KM<-round(saemixObject["options"]$nmc.is/MM)
    kmin<-min(10,ceiling(KM/4))
    x1<-MM*c(kmin:KM)
    y1<-(-2)*saemixObject["results"]["LL"][kmin:KM]
    if(plot.opt$new) {
      if(length(plot.opt$mfrow)==0) mfrow=c(1,1) else mfrow<-plot.opt$mfrow
      par(mfrow=mfrow,ask=plot.opt$ask)
    }
    if(sum(!is.na(y1))) plot(x1,y1,type="l",xlab=plot.opt$xlab, ylab=plot.opt$ylab,main=plot.opt$main,col=plot.opt$col,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
}

#######################	   Basic GOF plots & residuals	 ########################

saemix.plot.obsvspred<-function(saemixObject,...) {
# Predictions versus observations
  plot.opt<-saemixObject["prefs"]
  plot.opt$ylab<-"Observations"
  plot.opt$xlab<-"Predictions"
  plot.opt$main<-""
  plot.opt<-replace.plot.options(plot.opt,...)
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  if(plot.opt$new) {
    mfrow<-c(1,length(plot.opt$level))
    if(length(plot.opt$mfrow)>0) mfrow<-plot.opt$mfrow
    par(mfrow=mfrow,ask=plot.opt$ask)
  }
  if(saemixObject["model"]["error.model"]=="exponential")
    ydat<-saemixObject["data"]["yorig"] else ydat<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  if(length(grep(0,plot.opt$level))>0) {
    if(!change.main) main<-"Population predictions" else main<-plot.opt$main
    if(plot.opt$which.poppred=="ppred") xpl<-saemixObject["results"]["ppred"] else xpl<-saemixObject["results"]["ypred"]
    if(length(xpl)==length(ydat)) {
    plot(xpl,ydat,xlab=plot.opt$xlab, ylab=plot.opt$ylab,pch=plot.opt$pch, col=plot.opt$col,main=main,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(0,1,col=plot.opt$ablinecol,lty=plot.opt$ablinelty, lwd=plot.opt$ablinelwd)
    }
     }
  if(length(grep(1,plot.opt$level))>0) {
    if(!change.main) main<-paste("Individual predictions", ifelse(plot.opt$indiv.par=="map","MAP","Cond mean"),sep=", ") else main<-plot.opt$main
    if(plot.opt$indiv.par=="map") xpl<-saemixObject["results"]["ipred"] else xpl<-saemixObject["results"]["icpred"]
    if(length(xpl)==length(ydat)) {
    plot(xpl,ydat,xlab=plot.opt$xlab, ylab=plot.opt$ylab,pch=plot.opt$pch, col=plot.opt$col,main=main,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(0,1,col=plot.opt$ablinecol,lty=plot.opt$ablinelty,lwd=plot.opt$ablinelwd)
    }
   }
}

saemix.plot.distribresiduals<-function(saemixObject,...) {
# Histogram and QQ-plot
  plot.opt<-saemixObject["prefs"]
  plot.opt$main<-""
  plot.opt$level<-0:1
  plot.opt$smooth<-TRUE
  plot.opt$which.resplot<-c("dist.qqplot","dist.hist")
  plot.opt$which.pres<-c("wres","npde")
  plot.opt<-replace.plot.options(plot.opt,...)
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  change.xlab<-FALSE
  if(plot.opt$xlab!=saemixObject["prefs"]$xlab) change.xlab<-TRUE
  change.ylab<-FALSE
  if(plot.opt$ylab!=saemixObject["prefs"]$ylab) change.ylab<-TRUE
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)>0) mfrow<-plot.opt$mfrow else {
      ncol<-as.integer(1%in%plot.opt$level)+as.integer(0%in%plot.opt$level)* length(plot.opt$which.pres)
      mfrow=c(length(plot.opt$which.resplot),ncol)
    }
    par(mfrow=mfrow,ask=plot.opt$ask)
  }
  plot.ind<-FALSE
  if(1%in%plot.opt$level) {
    if(length(saemixObject["results"]["iwres"])==0) {
      cat("Please compute individual residuals first using predict().\n")
      return()
    }
    plot.ind<-TRUE
    if(plot.opt$indiv.par=="map") {
      iwres<-saemixObject["results"]["iwres"]
    } else {
      iwres<-saemixObject["results"]["icwres"]
    }
  }
  plot.pop<-FALSE
  if(0%in%plot.opt$level) {
    if(length(saemixObject["results"]["wres"])==0 | length(saemixObject["results"]["npde"])==0) {
      cat("Please compute WRES and npde first by using compute.sres().\n")
      return()
    }
    plot.pop<-TRUE
    wres<-saemixObject["results"]["wres"]
    npde<-saemixObject["results"]["npde"]
  }
  if("dist.qqplot"%in%plot.opt$which.resplot) {
  if(plot.pop & "wres"%in%plot.opt$which.pres) {
    laby<-"Sample quantiles"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Theoretical quantiles"
    if(change.xlab) labx<-plot.opt$xlab
    main<-"Population weighted residuals"
    if(change.main) main<-plot.opt$main
    qqnorm(wres,xlab=labx,ylab=laby,main=plot.opt$main, col=plot.opt$col)
    qqline(wres,lty=plot.opt$lty,col=plot.opt$col)
  }
  if(plot.ind) {
    laby<-"Sample quantiles"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Theoretical quantiles"
    if(change.xlab) labx<-plot.opt$xlab
    main<-"Individual weighted residuals"
    if(change.main) main<-plot.opt$main
    qqnorm(iwres,xlab=labx,ylab=laby,main=plot.opt$main, col=plot.opt$col)
    qqline(iwres,lty=plot.opt$lty,col=plot.opt$col)
  }
  if(plot.pop & "npde"%in%plot.opt$which.pres) {
    laby<-"Sample quantiles"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Theoretical quantiles"
    if(change.xlab) labx<-plot.opt$xlab
    main<-"NPDE"
    if(change.main) main<-plot.opt$main
    qqnorm(npde,xlab=labx,ylab=laby,main=plot.opt$main, col=plot.opt$col)
    qqline(npde,lty=plot.opt$lty,col=plot.opt$col)
  }
  }
  if("dist.hist"%in%plot.opt$which.resplot) {
  if(plot.pop & "wres"%in%plot.opt$which.pres) {
    labx<-"Population weighted residuals"
    if(change.xlab) labx<-plot.opt$xlab
    vec<-wres
    xh<-hist(vec,nclass=10,main=plot.opt$main, xlab=labx)
    if(plot.opt$smooth) {
      xpl<-min(vec)+c(0:100)/100*(max(vec)-min(vec))
      ypl<-dnorm(xpl)
      ypl<-ypl/max(ypl)*max(xh$counts)
      lines(xpl,ypl,lwd=2)
    }
  }
  if(plot.ind) {
    labx<-"Individual weighted residuals"
    if(change.xlab) labx<-plot.opt$xlab
    vec<-iwres
    xh<-hist(vec,nclass=10,main=plot.opt$main, xlab=labx)
    if(plot.opt$smooth) {
      xpl<-min(vec)+c(0:100)/100*(max(vec)-min(vec))
      ypl<-dnorm(xpl)
      ypl<-ypl/max(ypl)*max(xh$counts)
      lines(xpl,ypl,lwd=2)
    }
  }
  if(plot.pop & "npde"%in%plot.opt$which.pres) {
    labx<-"NPDE"
    if(change.xlab) labx<-plot.opt$xlab
    vec<-npde
    xh<-hist(vec,nclass=10,main=plot.opt$main, xlab=labx)
    if(plot.opt$smooth) {
      xpl<-min(vec)+c(0:100)/100*(max(vec)-min(vec))
      ypl<-dnorm(xpl)
      ypl<-ypl/max(ypl)*max(xh$counts)
      lines(xpl,ypl,lwd=2)
    }
  }
  }
}

saemix.plot.scatterresiduals<-function(saemixObject,...) {
# Graphs of residuals versus time and predictions
  plot.opt<-saemixObject["prefs"]
  plot.opt$main<-""
  plot.opt$level<-0:1
  plot.opt$which.resplot<-c("res.vs.x","res.vs.pred")
  plot.opt$which.pres<-c("wres","npde")
  plot.opt<-replace.plot.options(plot.opt,...)
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  change.xlab<-FALSE
  if(plot.opt$xlab!=saemixObject["prefs"]$xlab) change.xlab<-TRUE
  change.ylab<-FALSE
  if(plot.opt$ylab!=saemixObject["prefs"]$ylab) change.ylab<-TRUE
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)>0) mfrow<-plot.opt$mfrow else {
      ncol<-as.integer(1%in%plot.opt$level)+as.integer(0%in%plot.opt$level)* length(plot.opt$which.pres)
      mfrow=c(length(plot.opt$which.resplot),ncol)
    }
    par(mfrow=mfrow,ask=plot.opt$ask)
  }
  plot.ind<-FALSE
  if(1%in%plot.opt$level) {
    if(length(saemixObject["results"]["iwres"])==0) {
      cat("Please compute individual residuals first using predict().\n")
      return()
    }
    plot.ind<-TRUE
    if(plot.opt$indiv.par=="map") {
      iwres<-saemixObject["results"]["iwres"]
      ipred<-saemixObject["results"]["ipred"]
    } else {
      iwres<-saemixObject["results"]["icwres"]
      ipred<-saemixObject["results"]["icpred"]
    }
  }
  plot.pop<-FALSE
  if(0%in%plot.opt$level) {
    if(length(saemixObject["results"]["wres"])==0 | length(saemixObject["results"]["npde"])==0) {
      cat("Please compute WRES and npde first by using compute.sres().\n")
      return()
    }
    plot.pop<-TRUE
    wres<-saemixObject["results"]["wres"]
    npde<-saemixObject["results"]["npde"]
    if(plot.opt$which.poppred=="ppred") ppred<-saemixObject["results"]["ppred"] else ppred<-saemixObject["results"]["ypred"]
  }
  if("res.vs.x"%in%plot.opt$which.resplot) {
  if(plot.pop & "wres"%in%plot.opt$which.pres) {
    laby<-"Population weighted residuals"
    if(change.ylab) laby<-plot.opt$ylab
    plot(saemixObject["data"]["data"][,saemixObject["data"]["name.X"]],wres, pch=plot.opt$pch, col=plot.opt$col,main=plot.opt$main,xlab=plot.opt$xlab,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  if(plot.ind) {
    laby<-"Individual weighted residuals"
    if(change.ylab) laby<-plot.opt$ylab
    plot(saemixObject["data"]["data"][,saemixObject["data"]["name.X"]],iwres, pch=plot.opt$pch,col=plot.opt$col,main=plot.opt$main,xlab=plot.opt$xlab,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  if(plot.pop & "npde"%in%plot.opt$which.pres) {
    laby<-"NPDE"
    if(change.ylab) laby<-plot.opt$ylab
    plot(saemixObject["data"]["data"][,saemixObject["data"]["name.X"]],npde, pch=plot.opt$pch,col=plot.opt$col,main=plot.opt$main,xlab=plot.opt$xlab,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  }
  if("res.vs.pred"%in%plot.opt$which.resplot) {
  if(plot.pop & "wres"%in%plot.opt$which.pres) {
    laby<-"Population weighted residuals"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Population predictions"
    if(change.xlab) labx<-plot.opt$xlab
    plot(ppred,wres,pch=plot.opt$pch,col=plot.opt$col,main=plot.opt$main, xlab=labx,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  if(plot.ind) {
    laby<-"Individual weighted residuals"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Individual predictions"
    if(change.xlab) labx<-plot.opt$xlab
    plot(ipred,iwres,pch=plot.opt$pch,col=plot.opt$col,main=plot.opt$main, xlab=labx,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  if(plot.pop & "npde"%in%plot.opt$which.pres) {
    laby<-"NPDE"
    if(change.ylab) laby<-plot.opt$ylab
    labx<-"Population predictions"
    if(change.xlab) labx<-plot.opt$xlab
    plot(ppred,npde,pch=plot.opt$pch,col=plot.opt$col,main=plot.opt$main, xlab=labx,ylab=laby,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    abline(h=0,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=-1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    abline(h=1.96,lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
  }
  }
}

#######################	  	Individual fits		 ########################

saemix.plot.fits<-function(saemixObject,...) {
# Plot of the model fits overlayed on the data
  plot.opt<-saemixObject["prefs"]
  plot.opt$main<-""
  plot.opt$xlab<-paste(saemixObject["data"]["name.X"]," (",saemixObject["data"]["units"]$x,")",sep="")
  plot.opt$ylab<-paste(saemixObject["data"]["name.response"]," (",saemixObject["data"]["units"]$y,")",sep="")
  plot.opt$new<-TRUE
  plot.opt$ilist<-1:saemixObject["data"]["N"]  
  plot.opt$type<-"p"
  plot.opt$level<-c(1)
  plot.opt$ipred.lty<-1
  plot.opt$ppred.lty<-2
  plot.opt<-replace.plot.options(plot.opt,...)
  plot.opt$ilist<-plot.opt$ilist[plot.opt$ilist %in% 1:saemixObject["data"]["N"]]
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)==0) {
    np<-length(plot.opt$ilist)
    if(np>12) np<-12
    n1<-round(sqrt(np))
    n2<-ceiling(np/n1)
    par(mfrow=c(n1,n2),ask=plot.opt$ask)
  } else par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
  }
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  indplot<-(length(grep(1,plot.opt$level))>0)
  popplot<-(length(grep(0,plot.opt$level))>0)
  if(indplot & plot.opt$smooth & length(saemixObject["results"]["map.psi"])==0) {
    cat("Individual parameter estimates should be computed to produce individual plots, conditional means will be used.\n")
  }
  if(indplot & !(plot.opt$smooth) & length(saemixObject["results"]["ipred"])==0) {
    cat("For graphs of predictions, please use predict first.\n") 
    return()
  }
  if(popplot & !(plot.opt$smooth) & length(saemixObject["results"]["ppred"])==0) {
    cat("For graphs of predictions, please use predict first.\n") 
    return()
  }
  logtyp<-""
  if(plot.opt$xlog) logtyp<-paste(logtyp,"x",sep="")
  if(plot.opt$ylog) logtyp<-paste(logtyp,"y",sep="")
  pl.line<-(length(plot.opt$level)>0)
  xind<-saemixObject["data"]["data"][,saemixObject["data"]["name.predictors"], drop=FALSE]
  id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
  if(saemixObject["model"]["error.model"]=="exponential")
    yobs<-saemixObject["data"]["yorig"] else yobs<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  for(i1 in plot.opt$ilist) {
    isuj<-unique(id)[i1]
    if((indplot | popplot) & plot.opt$smooth) {
# If smooth requested and either population and individual predictions
      xvec<-xind[id==isuj, saemixObject["data"]["name.X"]]
      xpred<-seq(min(xvec),max(xvec),length.out=100)
      if(dim(xind)[2]==1) xdep<-matrix(xpred,ncol=1) else {
        x1<-xind[id==isuj,]
# creating an expanded dataframe (warning: will not work with different occasions)
# ECO TODO change this when several occasions
        id1<-unlist(lapply(xpred,function(x,vec) max(which(vec<=x)),vec=xvec))
        xdep<-x1[id1,]
        xdep[,saemixObject["data"]["name.X"]]<-xpred
      }
      idx<-rep(i1,dim(xdep)[1])
      if(indplot) {
# ECO TODO change this when several occasions
        if(length(saemixObject["results"]["map.psi"])>0)
	 ypred<-saemixObject["model"]["model"](saemixObject["results"]["map.psi"][, 2:dim(saemixObject["results"]["map.psi"])[2]],idx,xdep) else {
         psiM<-transphi(saemixObject["results"]["cond.mean.phi"], saemixObject["model"]["transform.par"])
         ypred<-saemixObject["model"]["model"](psiM,idx,xdep)
        }
      }
      if(popplot) {
       psiM<-transphi(saemixObject["results"]["mean.phi"], saemixObject["model"]["transform.par"])
       yppred<-saemixObject["model"]["model"](psiM,idx,xdep)
      }
    } else {
# else, use predictions at each observation time (no smooth)
      xpred<-xind[id==isuj,saemixObject["data"]["name.X"]]
      ypred<-saemixObject["results"]["ipred"][id==isuj]
      yppred<-saemixObject["results"]["ppred"][id==isuj]
    }
    vec<-yobs[id==isuj]
    if(indplot) vec<-c(vec,ypred)
    if(popplot) vec<-c(vec,yppred)
    if(length(plot.opt$ylim)>0) limy<-plot.opt$ylim else {
      if(!plot.opt$ylog) limy<-c(min(vec,na.rm=TRUE),max(vec,na.rm=TRUE)) else limy<-c(min(vec[!is.na(vec) & vec>0]),max(vec[!is.na(vec) & vec>0]))
    }
    main<-paste("Subject",isuj)
    if(change.main) main<-plot.opt$main
    plot(xind[id==isuj,saemixObject["data"]["name.X"]], yobs[id==isuj], xlab=plot.opt$xlab,ylab=plot.opt$ylab,log=logtyp,ylim=limy,type=plot.opt$type, main=main,pch=plot.opt$obs.pch,col=plot.opt$obs.col,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    if(pl.line) {
      if(indplot) lines(xpred,ypred,col=plot.opt$ipred.col, lty=plot.opt$ipred.lty,lwd=plot.opt$ipred.lwd)
      if(popplot) lines(xpred,yppred,col=plot.opt$ppred.col, lty=plot.opt$ppred.lty,lwd=plot.opt$ppred.lwd)
    }
  }
}

#######################	   Advanced GOF plots (VPC, npde) ########################
    
plotnpde<-function(xobs,npde,ypred) {
    nclass<-10
    par(mfrow=c(2,2))
    qqnorm(sort(npde),xlab="Sample quantiles (npde)",ylab="Theoretical Quantiles", cex.lab=1.5,main="Q-Q plot versus N(0,1) for npde")
    qqline(sort(npde))
    #Histogram of npde, with N(0,1) superimposed on the plot
    xh<-hist(npde,nclass=nclass,xlab="npde",main="",cex.lab=1.5)
    xpl<-min(npde)+c(0:100)/100*(max(npde)-min(npde))
    ypl<-dnorm(xpl)
    ypl<-ypl/max(ypl)*max(xh$counts)
    lines(xpl,ypl,lwd=2)
    
    #residuals
    plot(xobs,npde,xlab="X",ylab="npde",cex.lab=1.5)
    abline(h=0,lty=2)
    x1<-qnorm(0.05)
    abline(h=x1,lty=3);abline(h=(-x1),lty=3)
    
    plot(ypred,npde,xlab="Predicted Y",ylab="npde",cex.lab=1.5)
    abline(h=0,lty=2)
    abline(h=x1,lty=3);abline(h=(-x1),lty=3)
}

saemix.plot.npde<-function(saemixObject,...) {
  if(length(saemixObject["results"]["npde"])==0) {
    cat("Please estimate the npde first\n")
    return()
  }
  plotnpde(saemixObject["data"]["data"][,saemixObject["data"]["name.X"]], saemixObject["results"]["npde"],saemixObject["results"]["ypred"])
  y<-testnpde(saemixObject["results"]["npde"])
  return(y)
}

saemix.plot.vpc<-function(saemixObject,npc=FALSE,...) {
  if(length(saemixObject["sim.data"]["nsim"])==0) {
    cat("Please simulate data first, using the simul.saemix function.\n") 
    return()
  }
# Internal function
compute.vpc.pi<-function(ysim,xgrp,idrep,nbin,vpc.pi=0.95) {
  nsim<-length(unique(idrep))
  sim.pi.low<-sim.pi.med<-sim.pi.up<-matrix(0,nrow=nbin,ncol=nsim)
  alpha<-(1-vpc.pi)/2
  i0<-1
  for(irep in unique(idrep)) {
    ysim1<-ysim[idrep==irep]
    l1<-unlist(tapply(ysim1,xgrp,function(vec) quantile(vec,c(alpha,0.5,1-alpha))))
    l1<-matrix(l1,ncol=3,byrow=TRUE)
    sim.pi.low[,i0]<-l1[,1]
    sim.pi.med[,i0]<-l1[,2]
    sim.pi.up[,i0]<-l1[,3]
    i0<-i0+1
  }
  return(list(sim.pi.low=sim.pi.low,sim.pi.med=sim.pi.med,sim.pi.up=sim.pi.up))
}

  plot.opt<-saemixObject["prefs"]
  plot.opt$main<-"Visual Predictive Check"
  plot.opt<-replace.plot.options(plot.opt,...)
  if(plot.opt$new) {
    mfrow<-plot.opt$mfrow
    if(length(mfrow)==0) mfrow<-c(1,1)
    par(mfrow=mfrow,ask=plot.opt$ask)
  }
  logtyp<-""
  if(plot.opt$xlog) logtyp<-paste(logtyp,"x",sep="")
  if(plot.opt$ylog) logtyp<-paste(logtyp,"y",sep="")
  
  if(!is.na(pmatch(plot.opt$vpc.method,"optimal"))) {
    cat("Optimal binning not yet implemented, reverting to equal binning\n")
    plot.opt$vpc.method<-"equal"
  }
  if(!is.na(pmatch(plot.opt$vpc.method,"user")) & is.null(plot.opt$vpc.breaks)) {
    cat("User-defined method specified, but vpc.breaks is empty; reverting to equal binning\n")
    plot.opt$vpc.method<-"equal"
  }
  if(!is.na(pmatch(plot.opt$vpc.method,c("equal","width"))) & is.null(plot.opt$vpc.bin)) {
    plot.opt$vpc.bin<-10
  }

# Binning
  xvec<-saemixObject["data"]["data"][,saemixObject["data"]["name.X"]]
  ydat<-saemixObject["data"]["data"][,saemixObject["data"]["name.response"]]
  ysim<-saemixObject["sim.data"]["datasim"]$ysim
  nbin<-plot.opt$vpc.bin
  alpha<-(1-plot.opt$vpc.interval)/2
# ECO TODO: implement the optimal binning algorithm of Marc
#  if(is.na(plot.opt$vpc.bin)) {
#  } else { # binning by quantiles
#    bnds<-unique(quantile(xvec,seq(0,1,length.out=nbin),type=8))
#    xgrp<-findInterval(xvec,bnds)
    if(!is.na(pmatch(plot.opt$vpc.method,"user"))) {
      bnds<-plot.opt$vpc.breaks
      if(min(bnds)>=min(xvec)) bnds<-c(min(xvec)-1,bnds)
      if(max(bnds)<max(xvec)) bnds<-c(bnds,max(xvec))
    }
    if(!is.na(pmatch(plot.opt$vpc.method,"equal"))) {
      xvec2<-xvec;xvec2[xvec2==min(xvec)]<-min(xvec)-1
      bnds<-unique(quantile(xvec2,(0:nbin)/nbin,type=8))
    }
    if(!is.na(pmatch(plot.opt$vpc.method,"width"))) {
      if(plot.opt$xlog) xvec2<-log(xvec) else xvec2<-xvec
      bnds<-seq(min(xvec2),max(xvec2),length.out=(nbin+1))
      if(plot.opt$xlog) bnds<-exp(bnds)
      bnds[1]<-bnds[1]-1
    }
    if(!is.na(pmatch(plot.opt$vpc.method,c("equal","width","user")))) {
      xgrp<-factor(cut(xvec,bnds,include.lowest=F))
      if(!is.na(pmatch(plot.opt$vpc.method,"equal")) & length(unique(xvec))<=nbin)
        xgrp<-match(xvec,sort(unique(xvec)))
    } else {
      
    }
    nbin<-length(unique(xgrp))
    xpl<-tapply(xvec,xgrp,mean)
    if(!is.na(pmatch(plot.opt$vpc.method,c("equal","width","user")))) {
      tab<-cbind(Interval=names(xpl),Centered.On=format(xpl,digits=2))
      row.names(tab)<-1:dim(tab)[1]
      xnam<-switch(EXPR=plot.opt$vpc.method,equal="binning by quantiles on X", width="equal sized intervals",user="user-defined bins")
      cat("Method used for VPC:",xnam,", dividing into the following intervals\n")
      print(tab,quote=F)
    }
# Observed data
    ypl<-tapply(ydat,xgrp,mean)
    obs.bnd<-cbind(tapply(ydat,xgrp,quantile,alpha),tapply(ydat,xgrp,mean), tapply(ydat,xgrp,quantile,1-alpha))
#  }
  if(plot.opt$vpc.pi) {
    idsim<-saemixObject["sim.data"]["datasim"]$idsim
    idrep<-saemixObject["sim.data"]["datasim"]$irep
    isamp<-sample(1:saemixObject["options"]$nb.sim, saemixObject["options"]$nb.simpred)
    idx<-match(idrep,isamp,nomatch=0)>0
    sbnd<-compute.vpc.pi(ysim[idx],xgrp,idrep[idx],nbin,0.95)
    pi.low<-apply(sbnd$sim.pi.low,1,quantile,c(0.025,0.5,0.975))
    pi.med<-apply(sbnd$sim.pi.med,1,quantile,c(0.025,0.5,0.975))
    pi.up<-apply(sbnd$sim.pi.up,1,quantile,c(0.025,0.5,0.975))
    vec1<-c(pi.low,obs.bnd[,1])
    vec2<-c(obs.bnd[,3],pi.up)
    if(length(grep("y",logtyp))>0) {
      vec1<-vec1[vec1>0]
      vec2<-vec2[vec2>0]
    }
    limy<-c(min(vec1),max(vec2))
    xvec1<-xvec
    if(length(grep("x",logtyp))>0) xvec1<-xvec1[xvec1>0]
    limx<-c(min(xvec1),max(xvec1))

    plot(xpl,ypl,type="n",xlim=limx,ylim=limy,xlab=plot.opt$xlab, ylab=plot.opt$ylab,main=plot.opt$main,log=logtyp,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    polygon(c(xpl,rev(xpl)),c(pi.low[1,],rev(pi.low[3,])), col=plot.opt$col.fillpi,lty=plot.opt$lty.lpi, border=plot.opt$col.lpi)
    polygon(c(xpl,rev(xpl)),c(pi.up[1,],rev(pi.up[3,])), col=plot.opt$col.fillpi,lty=plot.opt$lty.lpi, border=plot.opt$col.lpi)
    polygon(c(xpl,rev(xpl)),c(pi.med[1,],rev(pi.med[3,])), col=plot.opt$col.fillmed,lty=plot.opt$lty.lmed, border=plot.opt$col.lmed)
    lines(xpl,pi.low[2,],lty=plot.opt$lty.lpi, col=plot.opt$col.lpi,lwd=plot.opt$lwd.lpi)
    lines(xpl,pi.med[2,],lty=plot.opt$lty.lmed, col=plot.opt$col.lmed,lwd=plot.opt$lwd.lmed)
    lines(xpl,pi.up[2,],lty=plot.opt$lty.lpi, col=plot.opt$col.lpi,lwd=plot.opt$lwd.lpi)
    lines(xpl,obs.bnd[,2],lty=plot.opt$lty.lobs, col=plot.opt$col.lmed,lwd=plot.opt$lwd.lobs)
    for (icol in c(1,3)) lines(xpl,obs.bnd[,icol],lty=plot.opt$lty.lobs, col=plot.opt$col.lobs,lwd=plot.opt$lwd.lobs)
    if(plot.opt$vpc.obs)
      points(xvec,ydat,pch=plot.opt$pch,col=plot.opt$col.pobs)
  } else {
# Simulated data
    nsim<-length(ysim)/length(ydat)
    id.grp<-rep(xgrp,nsim)
    sim.bnd<-cbind(tapply(ysim,id.grp,quantile,alpha),tapply(ysim,id.grp,mean), tapply(ysim,id.grp,quantile,1-alpha))
    vec1<-c(obs.bnd[,1],sim.bnd[,1])
    vec2<-c(obs.bnd[,3],sim.bnd[,3])
    if(length(grep("y",logtyp))>0) {
      vec1<-vec1[vec1>0]
      vec2<-vec2[vec2>0]
    }
    limy<-c(min(vec1),max(vec2))
    xvec1<-xvec
    if(length(grep("x",logtyp))>0) xvec1<-xvec1[xvec1>0]
    limx<-c(min(xvec1),max(xvec1))
    plot(xpl,ypl,type="n",xlim=limx,ylim=limy,xlab=plot.opt$xlab, ylab=plot.opt$ylab,main=plot.opt$main,log=logtyp,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
    polygon(c(xpl,rev(xpl)),c(sim.bnd[,3],rev(sim.bnd[,1])), col=plot.opt$fillcol,lty=plot.opt$ablinelty, border=plot.opt$ablinecol)
    lines(xpl,sim.bnd[,2],lty=plot.opt$ablinelty, col=plot.opt$ablinecol,lwd=plot.opt$ablinelwd)
    lines(xpl,obs.bnd[,2],lty=plot.opt$lty, col=plot.opt$lcol,lwd=plot.opt$lwd)
    for (icol in c(1,3)) lines(xpl,obs.bnd[,icol],lty=plot.opt$lty, col=plot.opt$lcol,lwd=plot.opt$lwd)
    if(plot.opt$vpc.obs)
      points(xvec,ydat,pch=plot.opt$pch,col=plot.opt$pcol)  
  }
  npc.stat<-c()
  if(npc==TRUE) {
    # ECO TODO: compute NPC - interpolation ? 
  }
  return(npc=npc.stat)
}

#######################	   Distribution of random effects ########################
saemix.plot.correlations<-function(saemixObject,...) {
  plot.opt<-saemixObject["prefs"]
  plot.opt$which.par<-"all"
  plot.opt$main<-"Correlations between random effects"
  plot.opt<-replace.plot.options(plot.opt,...)
  if(plot.opt$which.par[1]=="all") plot.opt$which.par<-saemixObject["model"]["name.modpar"]
  plist<-match(plot.opt$which.par,saemixObject["model"]["name.modpar"])
  plist<-plist[!is.na(match(plist,saemixObject["model"]["indx.omega"]))]
  if(plot.opt$indiv.par=="map" & length(saemixObject["results"]["map.psi"])) {
    indiv.par<-saemixObject["results"]["map.psi"]
  } else {
    if(plot.opt$indiv.par=="map") cat("No MAP estimates, using the conditional means for individual parameters.\n")
    indiv.par<-transphi(saemixObject["results"]["cond.mean.phi"], saemixObject["model"]["transform.par"])
  }
  labs<-saemixObject["model"]["name.modpar"][plist]
  pairs(indiv.par[,plist],labels=labs,panel=panel.smooth,main=plot.opt$main, pch=plot.opt$pch,col=plot.opt$col)
}

saemix.plot.randeff<-function(saemixObject,...) {
  plot.opt<-saemixObject["prefs"]
  plot.opt$which.par<-"all"
  plot.opt$main<-""
  plot.opt<-replace.plot.options(plot.opt,...)
  if(plot.opt$which.par[1]=="all") plot.opt$which.par<-saemixObject["model"]["name.modpar"]
  change.xlab<-FALSE
  if(plot.opt$xlab!=saemixObject["prefs"]$xlab & length(plot.opt$which.par)==1) change.xlab<-TRUE
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)==0) {
    np<-length(plot.opt$which.par)
    n1<-round(sqrt(np))
    n2<-ceiling(np/n1)
    if(n1>5 | n2>5) {
      n1<-3
      n2<-4
#      cat("Changing the plot layout\n")
    }
    par(mfrow=c(n1,n2),ask=plot.opt$ask)
  } else par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
  }
#  if(length(grep("l",plot.opt$line.smooth))>0)
  if(plot.opt$smooth)
    pl.psi<-transpsi(matrix(saemixObject["results"]["fixed.effects"],nrow=1), saemixObject["model"]["transform.par"])
  plist<-match(plot.opt$which.par,saemixObject["model"]["name.modpar"])
  for(ipar in plist) {
    tit<-saemixObject["model"]["name.modpar"][ipar]
    if(change.main) tit<-plot.opt$main    
    labx<-""
    if(change.xlab) labx<-plot.opt$xlab
    boxplot(saemixObject["results"]["cond.mean.phi"][,ipar],xlab=labx,main=tit,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
#    if(length(grep("l",plot.opt$line.smooth))>0)
    if(plot.opt$smooth)
      abline(h=pl.psi[1,ipar],lty=plot.opt$lty,lwd=plot.opt$ablinelwd, col=plot.opt$ablinecol)
  }
}

saemix.plot.distpsi<-function(saemixObject,...) {
# Plots the distribution of the model parameters conditional on covariates 
# plot.opt$cov.value: value of the covariates used
# Adds an histogram of individual parameter estimates if plot.opt$indiv.histo==TRUE
  plot.opt<-saemixObject["prefs"]
  plot.opt$which.par<-"all"
  plot.opt$main<-""
  plot.opt<-replace.plot.options(plot.opt,...)
  if(plot.opt$which.par[1]=="all") plot.opt$which.par<-saemixObject["model"]["name.modpar"]
  change.xlab<-FALSE
  if(plot.opt$xlab!=saemixObject["prefs"]$xlab & length(plot.opt$which.par)==1) change.xlab<-TRUE
  change.ylab<-FALSE
  if(plot.opt$ylab!=saemixObject["prefs"]$ylab & length(plot.opt$which.par)==1) change.ylab<-TRUE
  change.main<-FALSE
  if(plot.opt$main!=saemixObject["prefs"]$main) change.main<-TRUE
  nampar<-saemixObject["model"]["name.modpar"]
  plist<-match(plot.opt$which.par,nampar)
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)==0) {
    np<-length(plist)
    n1<-round(sqrt(np))
    n2<-ceiling(np/n1)
    if(n1>5 | n2>5) {
      n1<-3
      n2<-4
#      cat("Changing the plot layout\n")
    }
    par(mfrow=c(n1,n2),ask=plot.opt$ask)
  } else par(mfrow=plot.opt$mfrow,ask=plot.opt$ask)
  }
  if(plot.opt$indiv.histo) {
    if(plot.opt$indiv.par=="map" & length(saemixObject["results"]["map.psi"])) {
      indiv.par<-saemixObject["results"]["map.psi"]
    } else {
      if(plot.opt$indiv.par=="map") cat("No MAP estimates, using the conditional means for individual parameters.\n")
      indiv.par<-transphi(saemixObject["results"]["cond.mean.phi"], saemixObject["model"]["transform.par"])
    }
  }
  xpl<-c(1:100)/100*4
  xpl<-c(rev(-xpl),0,xpl)
  mpar<-saemixObject["results"]["betas"][saemixObject["model"]["indx.fix"]]
  for(ipar in plist) {
    labx<-nampar[ipar]
    if(change.xlab) labx<-plot.opt$xlab
    tit<-plot.opt$main
    if(colSums(saemixObject["model"]["covariate.model"])[ipar]>0) {
      idcov<-which(saemixObject["model"]["covariate.model"][,ipar]==1)
      for(icov in idcov) {
        xcov<-plot.opt$cov.value[icov]
        nlev<-length(unique(saemixObject["model"]["Mcovariates"][,(icov+1)]))
# covariable binaire
        if(is.na(xcov)) {
        if(nlev==2) 
          xcov<-min(saemixObject["model"]["Mcovariates"][,(icov+1)]) else 
# covariable continue
          xcov<-median(saemixObject["model"]["Mcovariates"][,(icov+1)])
        }
# ECO TODO securiser + dans le cas binaire, faire les 2 distributions si xcov=="all"
# Verifier le code ci-dessous
        ig1<-grep(saemixObject["model"]["name.cov"][idcov], saemixObject["model"]["name.fixed"])
        ig2<-grep(nampar[ipar],saemixObject["model"]["name.fixed"])
        iig<-c(ig1,ig2)
   mpar[ipar]<-mpar[ipar]+xcov*saemixObject["results"]["betas"][iig[duplicated(iig)]]
        if(!change.main) {
        if(tit!="") sep1<-"-" else sep1<-""
	xunit<-saemixObject["data"]["units"]$covariates[icov]
        tit<-paste(tit,sep1,saemixObject["model"]["name.cov"][icov],"=",xcov, ifelse(xunit=="-","",xunit),sep="")
        }
      }
      if(length(idcov)>0 & plot.opt$indiv.histo) cat("Warning: histograms of individual parameter estimates do not make sense since covariates enter the model for parameter",nampar[ipar],"\n")
    }
    xpl1<-mpar[ipar]+xpl*sqrt(diag(saemixObject["results"]["omega"]))[ipar]
    xpl2<-transphi(matrix(xpl1,ncol=1),saemixObject["model"]["transform.par"][ipar])
    if(saemixObject["model"]["transform.par"][ipar]==2) {
      ypl<-pnorm(xpl)*derivphi(matrix(xpl1,ncol=1), saemixObject["model"]["transform.par"][ipar])
    } else
      ypl<-dnorm(xpl)*derivphi(matrix(xpl1,ncol=1), saemixObject["model"]["transform.par"][ipar])
    if(plot.opt$indiv.histo) {
      vec<-c(indiv.par[,(ipar+1)],xpl2)
      limx<-c(min(vec),max(vec))
    } else limx<-c(min(xpl2),max(xpl2))
    if(limx[1]<0) limx[1]<-limx[1]*1.05 else limx[1]<-limx[1]*0.95
    if(limx[2]>0) limx[2]<-limx[2]*1.05 else limx[2]<-limx[2]*0.95
    if(plot.opt$indiv.histo) {
      laby<-"Counts"
      if(change.ylab) laby<-plot.opt$ylab
      h1<-hist(indiv.par[,(ipar+1)],xlim=limx,main=tit,xlab=labx,ylab=laby, col=plot.opt$fillcol)
      ypl<-ypl/max(ypl)*max(h1$counts)
      lines(xpl2,ypl,lty=plot.opt$lty,col=plot.opt$lcol,lwd=plot.opt$lwd)
    } else {
      laby<-"Density"
      if(change.ylab) laby<-plot.opt$ylab
      plot(xpl2,ypl,type="l",xlab=labx,ylab=laby,xlim=limx, main=tit,lty=plot.opt$lty,col=plot.opt$lcol,lwd=plot.opt$lwd,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
      }
  }
}

#######################	   Parameters versus covariates  ##########################

saemix.plot.parcov<-function(saemixObject,...) {
# non-user level function
# parameters versus covariates
  saemix.plot.parcov.aux(saemixObject,partype="p",...)  
}

saemix.plot.randeffcov<-function(saemixObject,...) {
# non-user level function
# random effects versus covariates
  saemix.plot.parcov.aux(saemixObject,partype="r",...)  
}

saemix.plot.parcov.aux<-function(saemixObject,partype="p",...) {
# Plot of parameters (parytype="p") or random effects ("r") versus covariates
  plot.opt<-saemixObject["prefs"]
  plot.opt$which.par<-"all"
  plot.opt$which.cov<-"all"
  plot.opt$main<-""
  plot.opt<-replace.plot.options(plot.opt,...)
  change.xlab<-FALSE
  if(plot.opt$xlab!=saemixObject["prefs"]$xlab) change.xlab<-TRUE
  change.ylab<-FALSE
  if(plot.opt$ylab!=saemixObject["prefs"]$ylab) change.ylab<-TRUE
  nampar<-saemixObject["model"]["name.modpar"]
  namcov<-saemixObject["data"]["name.covariates"]
  if(plot.opt$which.par[1]=="all") plot.opt$which.par<-nampar
  if(plot.opt$which.cov[1]=="all") plot.opt$which.cov<-namcov
  if(!is.integer(plot.opt$which.par)) plist<-match(plot.opt$which.par,nampar)
  plist<-plist[!is.na(plist)]
  if(!is.integer(plot.opt$which.cov)) clist<-match(plot.opt$which.cov,namcov)
  clist<-clist[!is.na(clist)]
  if(length(plist)==0) {
    cat("Cannot find parameter",plot.opt$which.par,", please check parameter names\n")
    return()
  }
  if(length(clist)==0) {
    cat("Cannot find covariates",plot.opt$which.cov,", please check covariate names\n")
    return()
  }
  replot<-FALSE
  mfrow<-plot.opt$mfrow
  if(plot.opt$new) {
    if(length(plot.opt$mfrow)==0) {
      if(length(plist)>1 & length(clist)>1) replot<-TRUE
      if(length(clist)>1) np<-length(clist) else np<-length(plist)   
      n1<-round(sqrt(np))
      n2<-ceiling(np/n1)
      mfrow<-c(n1,n2)
    }
    if(!replot) par(mfrow=mfrow,ask=plot.opt$ask)
  }
# ECO TODO: check that map.eta has a first column=Id
  if(partype=="r") { # random effects versus covariates
  if(tolower(plot.opt$indiv.par)=="map") {
    if(length(saemixObject["results"]["map.eta"])==0) {
      cat("Computing ETA estimates and adding them to fitted object.\n")
      saemixObject<-compute.eta.map(saemixObject)
    }
    param<-saemixObject["results"]["map.eta"]
  } else 
    param<-saemixObject["results"]["cond.mean.phi"]
  } else { # parameters versus covariates
# ECO TODO: check that map.psi has a first column=Id; maybe add one to cond.mean.phi for consistency
    if(tolower(plot.opt$indiv.par)=="map") 
      param<-saemixObject["results"]["map.psi"][, 2:dim(saemixObject["results"]["map.psi"])[2]] else 
      param<-transphi(saemixObject["results"]["cond.mean.phi"], saemixObject["model"]["transform.par"])
  }

# ECO: will not work with IOV  
  id<-saemixObject["data"]["data"][,saemixObject["data"]["name.group"]]
  idlist<-unique(id)
  matcov<-saemixObject["data"]["data"][match(idlist,id), saemixObject["data"]["name.covariates"],drop=FALSE]
  for(ipar in plist) {
    if(replot) par(mfrow=mfrow,ask=plot.opt$ask) # new page for each parameter (only if several covariates & several parameters, & plot.new==TRUE)
    xpar<-param[,ipar]
    laby<-nampar[ipar] 
    if(partype=="r") laby<-paste("ETA(",laby,")",sep="")
    if(change.ylab) laby<-plot.opt$ylab
    for(icov in clist) {
      covar<-matcov[,icov]
      labx<-saemixObject["data"]["name.covariates"][icov]
      if(change.xlab) labx<-plot.opt$xlab
      if(length(unique(covar))<=2) {
        boxplot(xpar~covar,xlab=labx,ylab=laby,main=plot.opt$main,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
        if (length(grep("m",plot.opt$line.smooth))>0) {
        y1<-saemixObject["results"]["fixed.psi"][ipar]
        abline(h=y1,lty=plot.opt$ablinelty,col=plot.opt$ablinecol, lwd=plot.opt$ablinelwd)
        }
       } else {
         plot(covar,xpar,xlab=labx,ylab=laby,main=plot.opt$main,pch=plot.opt$pch, col=plot.opt$col,cex.lab=plot.opt$cex.lab,cex.axis=plot.opt$cex.axis,cex.main=plot.opt$cex.main)
         if (length(grep("l",plot.opt$line.smooth))>0) {
          y1<-lm(xpar~covar)
          abline(y1,lty=plot.opt$ablinelty,col=plot.opt$ablinecol, lwd=plot.opt$ablinelwd)
         }
         if (length(grep("s",plot.opt$line.smooth))>0) {
          lines(lowess(covar[!is.na(covar)],xpar[!is.na(covar)]), lty=plot.opt$ablinelty,col=plot.opt$ablinecol, lwd=plot.opt$ablinelwd)
         }
         if (length(grep("m",plot.opt$line.smooth))>0) {
          y1<-saemixObject["results"]["fixed.psi"][ipar]
          abline(h=y1,lty=plot.opt$ablinelty,col=plot.opt$ablinecol,lwd=plot.opt$lwd)
       }
     }
    }
  }
#  return()
}

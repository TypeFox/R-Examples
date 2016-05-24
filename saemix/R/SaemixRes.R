####################################################################################
####			SaemixRes class - definition				####
####################################################################################

###############################
# Definition with initialise
setClass(
  Class="SaemixRes",
  representation=representation(
    name.fixed="character",	# names of fixed parameters in the model
    name.random="character",	# names of random effects
    name.res="character",	# names of parameters of residual error model
    npar.est="numeric",		# nb of parameters estimated (fixed, random & resid)
    fixed.effects="numeric",	# vector with h(mu) and betas in estimation order
    fixed.psi="numeric",	# h(mu)
    betas="matrix",		# estimated mu
    betaC="numeric",		# estimated fixed effects for covariates
    omega="matrix",		# estimated omega
    respar="numeric",		# estimated residual variability
    fim="matrix",		# Fisher information matrix
    se.fixed="numeric",		# estimated SE for fixed effects
    se.omega="numeric",	# estimated SE for Omega
    se.respar="numeric",	# estimated SE for residual variability
    parpop="matrix",	# population parameters at each iteration
    allpar="matrix",	# all parameters (including covariate effects) at each iteration
    indx.fix="numeric",	# index of mean param estimated (was indx.betaI)
    indx.cov="numeric",	# index of cov param estimated (was indx.betaC)
    indx.omega="numeric",	# index of random param estimated (was i1.omega2)
    indx.res="numeric",	# index of param of residual errors estimated (was indx.res)
    MCOV="matrix",		# 
# Individual parameters
    cond.mean.phi="matrix",	# Cond mean estimates of Phi
    cond.mean.psi="matrix",	# Cond mean estimates of Psi
    cond.var.phi="matrix",
    cond.mean.eta="matrix",	# Cond mean estimates of Eta
    cond.shrinkage="numeric",	# Shrinkage for cond mean estimates of Eta
    mean.phi="matrix",
    map.psi="data.frame",	# MAP estimates of individual parameters
    map.phi="data.frame",	# MAP estimates of phi
    map.eta="matrix",		# ETAs corresponding to the MAP estimates
    map.shrinkage="numeric", 	# shrinkage on MAP estimates
    phi="matrix",
    phi.samp="array",		# nb.chains samples in the individual conditional distributions
    phi.samp.var="array",	# variance of samples
# Statistical criteria
    ll.lin="numeric",		# for each method (linearisation, IS, GQ)
    aic.lin="numeric",		# ll=log-likelihood
    bic.lin="numeric",		# aic= Akaike Information Criterion
    ll.is="numeric",		# bic= Bayesian Information Criterion
    aic.is="numeric",
    bic.is="numeric",
    LL="numeric",		# LL for each iteration in the IS algorithm
    ll.gq="numeric",
    aic.gq="numeric",
    bic.gq="numeric",
# Model predictions and residuals
		predictions="data.frame", # data frame containing all the predictions and residuals below
    ypred="numeric",		# vector of mean population predictions
    ppred="numeric",		# vector of population predictions with MAP
    ipred="numeric",		# vector of individual predictions with MAP
    icpred="numeric",		# vector of individual predictions with conditional estimates
    iwres="numeric",		# vector of individual weighted residuals with MAP
    icwres="numeric",		# vector of individual weighted residuals with conditional estimates
    wres="numeric",		# vector of WRES (population weighted residuals)
    npde="numeric",		# vector of npde
    pd="numeric"		# vector of prediction discrepancies
  ),
  validity=function(object){
#    cat ("--- Checking SaemixRes object ---\n")
    return(TRUE)
  }
)

###############################
# initialize

setMethod(
  f="initialize",
  signature="SaemixRes",
  definition= function(.Object,name.fixed,name.random,name.res,fixed.effects, fixed.psi,betaC,betas,omega,respar,cond.mean.phi,cond.var.phi,mean.phi,phi, phi.samp,parpop, allpar,MCOV){
#    cat ("--- initialising SaemixRes Object --- \n")
    if(missing(name.fixed)) name.fixed<-character(0)
    .Object@name.fixed<-name.fixed
    if(missing(name.random)) name.random<-character(0)
    .Object@name.random<-name.random
    if(missing(name.res)) name.res<-character(0)
    .Object@name.res<-name.res
    if(missing(fixed.effects)) fixed.effects<-numeric(0)
    .Object@fixed.effects<-fixed.effects
    if(missing(fixed.psi)) fixed.psi<-numeric(0)
    .Object@fixed.psi<-fixed.psi
    if(missing(betas)) betas<-matrix(nrow=0,ncol=0)
    .Object@betas<-betas
    if(missing(betaC)) betaC<-numeric(0)
    .Object@betaC<-betaC
    if(missing(omega)) omega<-matrix(nrow=0,ncol=0)
#    if(missing(omega)) omega<-matrix(data=NA,nrow=length(),ncol=length())
    .Object@omega<-omega
    if(missing(respar)) respar<-numeric(0)
    .Object@respar<-respar
    if(missing(cond.mean.phi)) cond.mean.phi<-matrix(nrow=0,ncol=0)
    .Object@cond.mean.phi<-cond.mean.phi
    if(missing(cond.var.phi)) cond.var.phi<-matrix(nrow=0,ncol=0)
    .Object@cond.var.phi<-cond.var.phi
    if(missing(mean.phi)) mean.phi<-matrix(nrow=0,ncol=0)
    .Object@mean.phi<-mean.phi
    if(missing(phi)) phi<-matrix(nrow=0,ncol=0)
    .Object@phi<-phi
    if(missing(phi.samp)) phi.samp<-matrix(nrow=0,ncol=0)
    .Object@phi.samp<-phi.samp
    if(missing(parpop)) parpop<-matrix(nrow=0,ncol=0)
    .Object@parpop<-parpop
    if(missing(allpar)) allpar<-matrix(nrow=0,ncol=0)
    .Object@allpar<-allpar
    if(missing(MCOV)) MCOV<-matrix(nrow=0,ncol=0)
    .Object@MCOV<-MCOV
#    if(missing()) <-
#    .Object@<-
# Object validation
#    validObject(.Object)
    return (.Object )
  }
)

####################################################################################
####			SaemixRes class - accesseur				####
####################################################################################

# Getteur
setMethod(
  f ="[",
  signature = "SaemixRes" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "name.fixed"={return(x@name.fixed)},
    "name.res"={return(x@name.res)},
    "name.random"={return(x@name.random)},
    "npar.est"={return(x@npar.est)},
    "fixed.effects"={return(x@fixed.effects)},
    "fixed.psi"={return(x@fixed.psi)},
    "betas"={return(x@betas)},
    "betaC"={return(x@betaC)},
    "omega"={return(x@omega)},
    "respar"={return(x@respar)},
    "fim"={return(x@fim)},
    "se.fixed"={return(x@se.fixed)},
    "se.omega"={return(x@se.omega)},
    "se.respar"={return(x@se.respar)},
    "parpop"={return(x@parpop)},
    "allpar"={return(x@allpar)},
    "indx.fix"={return(x@indx.fix)},
    "indx.cov"={return(x@indx.cov)},
    "indx.omega"={return(x@indx.omega)},
    "indx.res"={return(x@indx.res)},
    "MCOV"={return(x@MCOV)},
    "cond.mean.phi"={return(x@cond.mean.phi)},
    "cond.mean.psi"={return(x@cond.mean.psi)},
    "cond.var.phi"={return(x@cond.var.phi)},
    "cond.mean.eta"={return(x@cond.mean.eta)},
    "cond.shrinkage"={return(x@cond.shrinkage)},
    "mean.phi"={return(x@mean.phi)},
    "map.psi"={return(x@map.psi)},
    "map.phi"={return(x@map.phi)},
    "map.eta"={return(x@map.eta)},
    "map.shrinkage"={return(x@map.shrinkage)},
    "phi"={return(x@phi)},
    "phi.samp"={return(x@phi.samp)},
    "phi.samp.var"={return(x@phi.samp.var)},
    "ll.lin"={return(x@ll.lin)},
    "aic.lin"={return(x@aic.lin)},
    "bic.lin"={return(x@bic.lin)},
    "ll.is"={return(x@ll.is)},
    "aic.is"={return(x@aic.is)},
    "bic.is"={return(x@bic.is)},
    "LL"={return(x@LL)},
    "ll.gq"={return(x@ll.gq)},
    "aic.gq"={return(x@aic.gq)},
    "bic.gq"={return(x@bic.gq)},
    "predictions"={return(x@predictions)},
    "ypred"={return(x@ypred)},
    "ppred"={return(x@ppred)},
    "ipred"={return(x@ipred)},
    "icpred"={return(x@icpred)},
    "iwres"={return(x@iwres)},
    "icwres"={return(x@icwres)},
    "wres"={return(x@wres)},
    "npde"={return(x@npde)},
    "pd"={return(x@pd)},
    stop("No such attribute\n")
   )
  }
)
#paste("    ",slotNames(saemix.res),"={return(x@",slotNames(saemix.res),")},",sep="")

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixRes" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "name.fixed"={x@name.fixed<-value},
    "name.random"={x@name.random<-value},
    "name.res"={x@name.res<-value},
    "npar.est"={x@npar.est<-value},
    "fixed.effects"={x@fixed.effects<-value},
    "fixed.psi"={x@fixed.psi<-value},
    "betas"={x@betas<-value},
    "betaC"={x@betaC<-value},
    "omega"={x@omega<-value},
    "respar"={x@respar<-value},
    "fim"={x@fim<-value},
    "se.fixed"={x@se.fixed<-value},
    "se.omega"={x@se.omega<-value},
    "se.respar"={x@se.respar<-value},
    "parpop"={x@parpop<-value},
    "allpar"={x@allpar<-value},
    "indx.fix"={x@indx.fix<-value},
    "indx.cov"={x@indx.cov<-value},
    "indx.omega"={x@indx.omega<-value},
    "indx.res"={x@indx.res<-value},
    "MCOV"={x@MCOV<-value},
    "cond.mean.phi"={x@cond.mean.phi<-value},
    "cond.mean.psi"={x@cond.mean.psi<-value},
    "cond.var.phi"={x@cond.var.phi<-value},
    "cond.mean.eta"={x@cond.mean.eta<-value},
    "cond.shrinkage"={x@cond.shrinkage<-value},
    "mean.phi"={x@mean.phi<-value},
    "map.phi"={x@map.phi<-value},
    "map.psi"={x@map.psi<-value},
    "map.eta"={x@map.eta<-value},
    "map.shrinkage"={x@map.shrinkage<-value},
    "phi"={x@phi<-value},
    "phi.samp"={x@phi.samp<-value},
    "phi.samp.var"={x@phi.samp.var<-value},
    "ll.lin"={x@ll.lin<-value},
    "aic.lin"={x@aic.lin<-value},
    "bic.lin"={x@bic.lin<-value},
    "ll.is"={x@ll.is<-value},
    "aic.is"={x@aic.is<-value},
    "bic.is"={x@bic.is<-value},
    "LL"={x@LL<-value},
    "ll.gq"={x@ll.gq<-value},
    "aic.gq"={x@aic.gq<-value},
    "bic.gq"={x@bic.gq<-value},
    "predictions"={x@predictions<-value},
    "ypred"={x@ypred<-value},
    "ppred"={x@ppred<-value},
    "ipred"={x@ipred<-value},
    "icpred"={x@icpred<-value},
    "iwres"={x@iwres<-value},
    "icwres"={x@icwres<-value},
    "wres"={x@wres<-value},
    "npde"={x@npde<-value},
    "pd"={x@pd<-value},
    stop("No such attribute\n")
   )
#   validObject(x)
   return(x)
  }
)

####################################################################################
####			SaemixRes class - method to print/show data		####
####################################################################################

setMethod("print","SaemixRes",
  function(x,digits=2,map=FALSE,...) {
#    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
#    cat("Dataset",x@name.data,"\n")
    if(length(x@betas)==0) {
      cat("No fit performed yet.\n")
      return()
    }
    cat("----------------------------------------------------\n")
    cat("-----------------  Fixed effects  ------------------\n")
    cat("----------------------------------------------------\n")
    if(length(x@se.fixed)==0) {
      tab<-cbind(c(x@name.fixed,x@name.res[x@indx.res]), c(x@fixed.effects,x@respar[x@indx.res]))
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-cbind(c(x@name.fixed,x@name.res[x@indx.res]), c(x@fixed.effects,x@respar[x@indx.res]),c(x@se.fixed,x@se.respar[x@indx.res]))
      tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
      if(length(x@indx.cov)>0) {
      wstat<-as.double(tab[,2])/as.double(tab[,3])
      pval<-rep("-",length(wstat))
      pval[x@indx.cov]<-1-normcdf(abs(wstat[x@indx.cov]))
      tab<-cbind(tab,"p-value"=pval)
      }
      is.not.est<-which(as.double(tab[,3])<=.Machine$double.xmin)
      ncol<-dim(tab)[2]-2
      tab[is.not.est,3:dim(tab)[2]]<-rep("-",ncol)
    }
    if(digits>0) {
      for(i in 2:dim(tab)[2]) {
       xcol<-as.double(as.character(tab[,i]))
       idx<-which(!is.na(xcol))
       tab[idx,i]<-format(xcol[idx],digits=digits)
      }
    }
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("-----------  Variance of random effects  -----------\n")
    cat("----------------------------------------------------\n")
#  cat("   ECO TODO: check if Omega or Omega2 (SD or variances) and can we choose ?\n") => returns omega2, and we can't choose
    if(length(x@se.omega)==0) {
      tab<-cbind(x@name.random,diag(x@omega)[x@indx.omega])
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-cbind(x@name.random,diag(x@omega)[x@indx.omega],x@se.omega[x@indx.omega])
      tab<-cbind(tab,100*as.double(tab[,3])/as.double(tab[,2]))
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
    }
    if(digits>0) {
      for(i in 2:dim(tab)[2]) 
         tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits)
    }
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("------  Correlation matrix of random effects  ------\n")
    cat("----------------------------------------------------\n")
    tab<-cov2cor(x@omega[x@indx.omega,x@indx.omega,drop=FALSE])
    if(digits>0) {
      for(i in 1:dim(tab)[2]) tab[,i]<-format(as.double(as.character(tab[,i])),digits=digits)
    }
    try(colnames(tab)<-rownames(tab)<-x@name.random)
    print(tab,quote=FALSE)
    if(length(x@ll.lin)>0 | length(x@ll.is)>0 | length(x@ll.gq)>0) {
    cat("----------------------------------------------------\n")
    cat("---------------  Statistical criteria  -------------\n")
    cat("----------------------------------------------------\n")
    if(length(x@ll.lin)>0) {
    cat("Likelihood computed by linearisation\n")
    cat("      -2LL=",(-2*x@ll.lin),"\n")
    cat("      AIC =",x@aic.lin,"\n")
    cat("      BIC =",x@bic.lin,"\n")
#  cat("   ECO TODO: verifier si ca renvoie LL ou -2LL\n"): ok renvoie -2LL
    }
    if(length(x@ll.is)>0) {
    cat("\nLikelihood computed by importance sampling\n")
    cat("      -2LL=",(-2*x@ll.is),"\n")
    cat("      AIC =",x@aic.is,"\n")
    cat("      BIC =",x@bic.is,"\n")
    }  
    if(length(x@ll.gq)>0) {
    cat("\nLikelihood computed by Gaussian quadrature\n")
    cat("      -2LL=",(-2*x@ll.gq),"\n")
    cat("      AIC =",x@aic.gq,"\n")
    cat("      BIC =",x@bic.gq,"\n")
    }
    cat("----------------------------------------------------\n")
    }
    if(length(x@map.psi)>0 & map) {
    cat("----------------------------------------------------\n")
    cat("---------------  Individual parameters  ------------\n")
    cat("----------------------------------------------------\n")
      if(dim(x@map.psi)[1]<30)
      print(x@map.psi) else {
      cat("Individual estimates for the first 30 subjects:\n")
      print(x@map.psi[1:30,])
      }
    }
  }
)

setMethod("show","SaemixRes",
  function(object) {
#    cat("Nonlinear mixed-effects model fit by the SAEM algorithm\n")
    cat("Fixed effects\n")
    if(length(object@se.fixed)==0) {
      tab<-cbind(c(object@name.fixed,object@name.res[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]))
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-cbind(c(object@name.fixed,object@name.res[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]), c(object@se.fixed,object@se.respar[object@indx.res]))
      tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
      colnames(tab)<-c("Parameter","Estimate","  SE"," CV(%)")
      if(length(object@indx.cov)>0) {
      wstat<-as.double(tab[,2])/as.double(tab[,3])
      pval<-rep("-",length(wstat))
      pval[object@indx.cov]<-1-normcdf(abs(wstat[object@indx.cov]))
      tab<-cbind(tab,"p-value"=pval)
      }
      is.not.est<-which(as.double(tab[,3])<=.Machine$double.xmin)
      ncol<-dim(tab)[2]-2
      tab[is.not.est,3:dim(tab)[2]]<-rep("-",ncol)
    }
      for(i in 2:dim(tab)[2]) {
       xcol<-as.double(as.character(tab[,i]))
       idx<-which(!is.na(xcol))
       tab[idx,i]<-format(xcol[idx],digits=3)
      }
    rownames(tab)<-rep("",dim(tab)[1])
    print(tab,quote=FALSE)

    cat("\nVariance of random effects\n")
#  cat("   ECO TODO: check if Omega or Omega2 (SD or variances) and can we choose ?\n")
    if(length(object@se.omega)==0) {
      tab<-cbind(object@name.random,diag(object@omega)[object@indx.omega])
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-cbind(object@name.random,diag(object@omega)[object@indx.omega],object@se.omega[object@indx.omega])
      tab<-cbind(tab,100*as.double(tab[,3])/as.double(tab[,2]))
      colnames(tab)<-c("Parameter","Estimate","  SE"," CV(%)")
    }
      for(i in 2:dim(tab)[2]) 
         tab[,i]<-format(as.double(as.character(tab[,i])),digits=3)
    rownames(tab)<-rep("",dim(tab)[1])
    print(tab,quote=FALSE)
    if(length(object@ll.lin)>0 | length(object@ll.is)>0 | length(object@ll.gq)>0) {
    cat("\nStatistical criteria\n")
    }
    mat1<-object@omega
    if(sum(abs(mat1-diag(diag(mat1))))>0) {
    cat("\nCorrelation matrix of random effects\n")
    tab<-cov2cor(object@omega[object@indx.omega,object@indx.omega,drop=FALSE])
    for(i in 1:dim(tab)[2]) 
      tab[,i]<-format(as.double(as.character(tab[,i])),digits=3)
    try(colnames(tab)<-rownames(tab)<-object@name.random)
    print(tab,quote=FALSE)
    }
    if(length(object@ll.lin)>0) {
    cat("Likelihood computed by linearisation\n")
    cat("      -2LL=",(-2*object@ll.lin),"\n")
    cat("       AIC=",object@aic.lin,"\n")
    cat("       BIC=",object@bic.lin,"\n")
#  cat("   ECO TODO: verifier si ca renvoie LL ou -2LL\n"): ok renvoie -2LL
    }
    if(length(object@ll.is)>0) {
    cat("Likelihood computed by importance sampling\n")
    cat("      -2LL=",(-2*object@ll.is),"\n")
    cat("       AIC=",object@aic.is,"\n")
    cat("       BIC=",object@bic.is,"\n")
    }  
    if(length(object@ll.gq)>0) {
    cat("Likelihood computed by Gaussian quadrature\n")
    cat("      -2LL=",(-2*object@ll.gq),"\n")
    cat("       AIC=",object@aic.gq,"\n")
    cat("       BIC=",object@bic.gq,"\n")
    }
#    cat("----------------------------------------------------\n")
  }
)

# Could be print, with only head of data
setMethod("showall","SaemixRes",
  function(object) {
    cat("\n----------------------------------------------------\n")
    cat("-----------------  Fixed effects  ------------------\n")
    cat("----------------------------------------------------\n")
    if(length(object@se.fixed)==0) {
      tab<-cbind(c(object@name.fixed,object@name.res[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]))
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-cbind(c(object@name.fixed,object@name.res[object@indx.res]), c(object@fixed.effects,object@respar[object@indx.res]), c(object@se.fixed,object@se.respar[object@indx.res]))
      tab<-cbind(tab,100*abs(as.double(tab[,3])/as.double(tab[,2])))
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
      if(length(object@indx.cov)>0) {
      wstat<-as.double(tab[,2])/as.double(tab[,3])
      pval<-rep("-",length(wstat))
      pval[object@indx.cov]<-1-normcdf(abs(wstat[object@indx.cov]))
      tab<-cbind(tab,"p-value"=pval)
      }
      is.not.est<-which(as.double(tab[,3])<=.Machine$double.xmin)
      ncol<-dim(tab)[2]-2
      tab[is.not.est,3:dim(tab)[2]]<-rep("-",ncol)
    }
    for(i in 2:dim(tab)[2]) {
       xcol<-as.double(as.character(tab[,i]))
       idx<-which(!is.na(xcol))
       tab[idx,i]<-format(xcol[idx],digits=3)
      }
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("-----------  Variance of random effects  -----------\n")
    cat("----------------------------------------------------\n")
#  cat("   ECO TODO: check if Omega or Omega2 (SD or variances) and can we choose ?\n")
    if(length(object@se.omega)==0) {
      tab<-cbind(object@name.random,diag(object@omega)[object@indx.omega])
      colnames(tab)<-c("Parameter","Estimate")
    } else {
      tab<-cbind(object@name.random,diag(object@omega)[object@indx.omega],object@se.omega[object@indx.omega])
      tab<-cbind(tab,100*as.double(tab[,3])/as.double(tab[,2]))
      colnames(tab)<-c("Parameter","Estimate","SE","CV(%)")
    }
      for(i in 2:dim(tab)[2]) 
         tab[,i]<-format(as.double(as.character(tab[,i])),digits=3)
    print(tab,quote=FALSE)
    cat("----------------------------------------------------\n")
    cat("---------------  Statistical criteria  -------------\n")
    cat("----------------------------------------------------\n")
    if(length(object@ll.lin)>0) {
    cat("Likelihood computed by linearisation\n")
    cat("      -2LL=",(-2*object@ll.lin),"\n")
    cat("      AIC =",object@aic.lin,"\n")
    cat("      BIC =",object@bic.lin,"\n")
#  cat("   ECO TODO: verifier si ca renvoie LL ou -2LL\n"): ok renvoie -2LL
    }
    if(length(object@ll.is)>0) {
    cat("\nLikelihood computed by importance sampling\n")
    cat("      -2LL=",(-2*object@ll.is),"\n")
    cat("      AIC =",object@aic.is,"\n")
    cat("      BIC =",object@bic.is,"\n")
    }  
    if(length(object@ll.gq)>0) {
    cat("\nLikelihood computed by Gaussian quadrature\n")
    cat("      -2LL=",(-2*object@ll.gq),"\n")
    cat("      AIC =",object@aic.gq,"\n")
    cat("      BIC =",object@bic.gq,"\n")
    }
    cat("----------------------------------------------------\n")
  }
)

####################################################################################
####			SaemixRes class - method to plot			####
####################################################################################

####################################################################################

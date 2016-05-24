####################################################################################
####			SaemixModel class - definition				####
####################################################################################

setClass(
  Class="SaemixModel",
  representation=representation(
    model="function", 		# name of model function
    description="character",	# model description
    psi0="matrix",		# CI for parameter estimates
    transform.par="numeric",	# distribution for model parameters
    fixed.estim="numeric",	# 1 for fixed parameters estimated
    error.model="character",	# residual error model
    covariate.model="matrix",	# covariate model
    betaest.model="matrix",	# 1st line=ones, next lines=covariate model
    covariance.model="matrix",	# covariance model
    omega.init="matrix",	# CI for Omega
    error.init="numeric",	# CI for residual error
    nb.parameters="integer",	# nb of parameters in the model
    name.modpar="character",	# name of parameters in the model (columns of psi0)
    name.fixed="character",	# name of fixed parameters
    name.random="character",	# name of random parameters
    name.res="character",	# name of residual parameters (maybe not necessary)
    name.predictors="character",# name of predictors 
    name.X="character",	# name of X 
    name.response="character",	# name of response
    name.cov="character",	# name of covariates
    indx.fix="numeric",		# index of mean param estimated (was indx.betaI)
    indx.cov="numeric",		# index of cov param estimated (was indx.betaC)
    indx.omega="numeric",	# index of random param estimated (was i1.omega2)
    indx.res="numeric",		# index of param of residual errors estimated (was indx.res)
    Mcovariates="data.frame"	# matrix of individual covariates in the model
  ),
  validity=function(object){
#    cat ("--- Checking SaemixModel object ---\n")
    if (dim(object@psi0)[1]==0) {
      cat("[ SaemixModel : Error ] Please provide initial estimates for the fixed effect (a matrix with columns named after the parameters in the model).\n")
      return("Missing psi0")
    }
    isize<-0
    npar<-dim(object@psi0)[2]
    if(npar!=length(object@transform.par)) isize<-1
    if(npar!=length(object@fixed.estim)) isize<-1
    if (npar!=dim(object@covariate.model)[2]) isize<-1
    if (npar!=dim(object@covariance.model)[1]) isize<-1
    if (npar!=dim(object@omega.init)[1]) isize<-1
#    cat("npar=",npar,length(object@transform.par),length(object@fixed.estim), dim(object@covariate.model)[2],dim(object@covariance.model)[1],dim(object@omega.init)[1],"\n")
    if(isize==1) {
      cat("[ SaemixModel : Error ] The number of parameters should be the same in the following elements: psi0 (initial conditions), transform.par, fixed.estim, covariate.model, and the matrices covariance.model and omega.init should be square matrices of size equal to the number of parameters. Please check the input.\n")
      return("Size mismatch")
    }
    if(npar<2) {
      cat("[ SaemixModel : Error ] SAEM needs at least two parameters to work on.\n")
      return("Psi0 has size 1")
    }
		if(sum(object@fixed.estim*mydiag(object@covariance.model))==0) {
			cat("[ SaemixModel : Error ] ")
			if(sum(mydiag(object@covariance.model))==0) cat("At least one parameter with IIV must be included in the model.\n") else cat("At least one parameter with IIV must be estimated and not fixed in the model.\n")
			return("Invalid IIV structure")
		}
    if(is.na(match(object@error.model,c('constant','proportional','combined', 'exponential')))) {
      cat("[ SaemixModel : Error ] Invalid residual error model")
      return("Invalid residual error model")
    }
    return(TRUE)
  }
)

setMethod(
  f="initialize",
  signature="SaemixModel",
  definition=function(.Object,model,description,psi0,transform.par,fixed.estim, error.model,covariate.model,covariance.model,omega.init,error.init,nb.parameters, name.modpar){
#    cat ("--- initialising SaemixModel Object --- \n")
    if(missing(model)) {
#      cat("Error initialising SaemixModel object:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
      return(.Object)
    }
    .Object@model<-model
    if(missing(description)) description<-""
    .Object@description<-description
    if(missing(psi0) || length(psi0)==0) {
      cat("Error initialising SaemixModel object:\n   Please provide initial estimates for the fixed effect (a matrix with columns named after the parameters in the model).\n")
      return(.Object)
    }
    npar<-dim(psi0)[2]
    if(missing(name.modpar) || length(name.modpar)==0) {
      y1<-try(name.modpar<-colnames(psi0))
      if(class(y1)=="try-error") {
        cat("     Can't find parameter names.\n")
        name.modpar<-paste("theta",1:npar)
      }
    }
    if(is.null(colnames(psi0))) {
      y1<-try(colnames(psi0)<-name.modpar)
      if(class(y1)=="try-error") {
        cat("Warning:\n   Problem with names of psi0\n")
        colnames(psi0)<-name.modpar<-paste("theta",1:npar)
      }
    }
    if(dim(psi0)[1]<2 & sum(covariate.model)>0){
    	psi0<-rbind(psi0,rep(0,dim(psi0)[2]))
    }
    if(is.null(rownames(psi0))) {
      rownames(psi0)<-rep("",dim(psi0)[1])
      rownames(psi0)[1]<-"Pop.CondInit"
      if(dim(psi0)[1]>1) rownames(psi0)[2:dim(psi0)[1]]<-"Cov.CondInit"
    }
    .Object@psi0<-psi0    
    .Object@name.modpar<-name.modpar
    if(missing(error.model) || length(error.model)==0) error.model<-"constant"
    .Object@error.model<-error.model
# Checking sizes
    .Object@nb.parameters<-npar
    if(missing(transform.par) || length(transform.par)==0) transform.par<-rep(0,npar)
    .Object@transform.par<-transform.par
    if(missing(fixed.estim) || length(fixed.estim)==0) fixed.estim<-rep(1,npar)
    .Object@fixed.estim<-fixed.estim
    if(missing(covariate.model) || length(covariate.model)==0 || sum(covariate.model)==0) covariate.model<-matrix(nrow=0,ncol=npar)
    if(is.null(colnames(covariate.model))) colnames(covariate.model)<-colnames(psi0)
    .Object@covariate.model<-covariate.model
    if(missing(covariance.model) || length(covariance.model)==0) {
      covariance.model<-diag(nrow=npar,ncol=npar)
    } else {
      if(dim(covariance.model)[1]!=dim(covariance.model)[2]) {
        cat("Error initialising SaemixModel object:\n   The covariance model needs to be a square matrix, please check dimensions.\n")
      return(.Object)
      }
    }
    nomg<-dim(covariance.model)[1]
    if(nomg!=npar) {
      cat("Error initialising SaemixModel object:\n   The covariance model needs to have the same size as the number of parameters.\n")
      return(.Object)
    }
    if(is.null(colnames(covariance.model))) colnames(covariance.model)<-rownames(covariance.model)<-colnames(psi0)
    .Object@covariance.model<-covariance.model
    indx.omega<-which(diag(covariance.model)>0)
    .Object@indx.omega<-indx.omega
    if(!missing(omega.init) && length(omega.init)>0) {
      if(dim(omega.init)[1]!=dim(omega.init)[2]) {
        cat("Warning:   the matrix giving the initial conditions for the covariance model (omega.init) needs to be a square matrix. Changing it to the diagonal matrix\n")
        omega.init<-NULL
      }
    }
    if(missing(omega.init) || length(omega.init)==0) {
      omega.init<-diag(fixed.estim)
      diag.omegi<-rep(1,npar)
      j1<-which(transform.par==0)
      if(length(j1)>0) {
      	diag.omegi[j1]<-sapply(psi0[1,j1]**2,function(x) { x[x<1]<-1; return(x)})
#      for(i in j1) d[i]<-max(psi0[i]^2,1)
      }
      omega.init<-diag(diag.omegi,nrow=npar)
    }
    if(is.null(colnames(omega.init))) colnames(omega.init)<-rownames(omega.init)<-colnames(psi0)
    .Object@omega.init<-omega.init
		if(sum(.Object@fixed.estim*mydiag(.Object@covariance.model))==0) {
			cat("Error initialising SaemixModel object:\n")
#			if(sum(mydiag(.Object@covariance.model))==0) cat("At least one parameter with IIV must be included in the model.\n") else cat("At least one parameter with IIV must be estimated and not fixed in the model.\n")
			return(.Object)
		}

## Residual Error model.
# error models are a + bf described by [a b]
# error models :
#   constant            y = f + a*e
#   proportional        y = f + b*f*e
#   combined            y = f + (a+b*f)*e
#   exponential         y = f*exp(a*e)    ( <=>  log(y) = log(f) + a*e )
    if(missing(error.init) || length(error.init)!=2) {
      error.init<-switch(error.model,
        "constant"=c(1,0),
        "exponential"=c(1,0),
        "proportional"=c(0,1),
        "combined"=c(1,1))
     }
    .Object@error.init<-error.init
    .Object@name.res<-c("a","b")
    if(.Object@error.model=='constant') {
      indx.res<-1
    } else {
      if(.Object@error.model=='proportional') {
        indx.res<-2
      } else {
        if(.Object@error.model=='combined') {
          indx.res<-c(1,2) 
        } else {
          if(.Object@error.model=='exponential') {
           indx.res<-1
           }
        }
      }
    }
    if(length(indx.res<2)) .Object@error.init[-indx.res]<-0
    .Object@indx.res<-indx.res
    .Object@betaest.model<-matrix(c(rep(1,.Object@nb.parameters), c(t(.Object@covariate.model))),ncol=.Object@nb.parameters,byrow=TRUE)
    colnames(.Object@betaest.model)<-colnames(.Object@covariate.model)
    if(!is.null(rownames(.Object@covariate.model))) {
      rownames(.Object@betaest.model)<-c("Fixed",rownames(.Object@covariate.model))
    } else {
      rownames(.Object@betaest.model)<-rep("",dim(.Object@betaest.model)[1])
      rownames(.Object@betaest.model)[1]<-"Fixed"
    }
# Object validation
    validObject(.Object)
    return (.Object)
  }
)

####################################################################################
####			saemixData class - accesseur				####
####################################################################################

# Getteur
setMethod(
  f ="[",
  signature = "SaemixModel" ,
  definition = function (x,i,j,drop ){
  switch (EXPR=i,
    "model"={return(x@model)},
    "description"={return(x@description)},
    "psi0"={return(x@psi0)},
    "transform.par"={return(x@transform.par)},
    "fixed.estim"={return(x@fixed.estim)},
    "error.model"={return(x@error.model)},
    "covariate.model"={return(x@covariate.model)},
    "betaest.model"={return(x@betaest.model)},
    "covariance.model"={return(x@covariance.model)},
    "omega.init"={return(x@omega.init)},
    "error.init"={return(x@error.init)},
    "nb.parameters"={return(x@nb.parameters)},
    "name.modpar"={return(x@name.modpar)},
    "name.fixed"={return(x@name.fixed)},
    "name.random"={return(x@name.random)},
    "name.res"={return(x@name.res)},
    "name.X"={return(x@name.X)},
    "name.response"={return(x@name.response)},
    "name.predictors"={return(x@name.predictors)},
    "name.cov"={return(x@name.cov)},
    "indx.fix"={return(x@indx.fix)},
    "indx.cov"={return(x@indx.cov)},
    "indx.omega"={return(x@indx.omega)},
    "indx.res"={return(x@indx.res)},
    "Mcovariates"={return(x@Mcovariates)},
    stop("No such attribute\n")
   )
  }
)

# Setteur
setReplaceMethod(
  f ="[",
  signature = "SaemixModel" ,
  definition = function (x,i,j,value){
  switch (EXPR=i,
    "model"={x@model<-value},
    "description"={return(x@description)},
    "psi0"={x@psi0<-value},
    "transform.par"={x@transform.par<-value},
    "fixed.estim"={x@fixed.estim<-value},
    "error.model"={x@error.model<-value},
    "covariate.model"={x@covariate.model<-value},
    "betaest.model"={x@betaest.model<-value},
    "covariance.model"={x@covariance.model<-value},
    "omega.init"={x@omega.init<-value},
    "error.init"={x@error.init<-value},
    "nb.parameters"={x@nb.parameters<-value},
    "name.modpar"={x@name.modpar<-value},
    "name.fixed"={x@name.fixed<-value},
    "name.random"={x@name.random<-value},
    "name.res"={x@name.res<-value},
    "name.X"={x@name.X<-value},
    "name.response"={x@name.response<-value},
    "name.predictors"={x@name.predictors<-value},
    "name.cov"={x@name.cov<-value},
    "indx.fix"={x@indx.fix<-value},
    "indx.cov"={x@indx.cov<-value},
    "indx.omega"={x@indx.omega<-value},
    "indx.res"={x@indx.res<-value},
    "Mcovariates"={x@Mcovariates<-value},
    stop("No such attribute\n")
   )
   validObject(x)
   return(x)
  }
)

####################################################################################
####			SaemixModel class - method to print/show data		####
####################################################################################

setMethod("print","SaemixModel",
  function(x,...) {
    cat("Nonlinear mixed-effects model\n")
    distrib<-c("normal","log-normal","probit","logit")
    cat("  Model function")
    if(length(x@description)>0 && nchar(x@description)>0) cat(": ",x@description)
    cat("\n")
    print(x@model)
    cat("  Nb of parameters:",x@nb.parameters,"\n")
    cat("      parameter names: ",x@name.modpar,"\n")
    cat("      distribution:\n")
    tab<-cbind(Parameter=x@name.modpar,Distribution=distrib[x@transform.par+1], Estimated=ifelse(x@fixed.estim==1,"Estimated","Fixed"))
    print(tab,quote=FALSE)
    cat("  Variance-covariance matrix:\n")
    tab<-x@covariance.model
#    try(colnames(tab)<-rownames(tab)<-x@name.modpar)
    print(tab,quote=FALSE)
    st1<-paste(x@name.res,x@error.init,sep="=")
    cat("  Error model:",x@error.model,", initial values:",st1[x@indx.res],"\n")
   if(dim(x@covariate.model)[1]>0) {
      cat("  Covariate model:")
      if(sum(x@covariate.model)==0) cat(" none\n") else {
        cat("\n")
        print(x@covariate.model)
    }
  } else cat("    No covariate in the model.\n")
    cat("    Initial values\n")
    print(x@psi0)
  }
)

setMethod("show","SaemixModel",
  function(object) {
    cat("Nonlinear mixed-effects model\n")
    cat("  Model function")
    if(length(object@description)>0 && nchar(object@description)>0) {
      cat(": ",object@description,"\n")}
    else {
      cat("\n")
      print(object@model)
    }
    fix1<-ifelse(object@fixed.estim==1,""," [fixed]")
    cat("    ",object@nb.parameters,"parameters:", paste(object@name.modpar,fix1,sep=""),"\n")
    cat("     error model:",object@error.model,"\n")
    if(dim(object@covariate.model)[1]>0) {
      cat("     covariate model:\n")
      print(object@covariate.model) 
    } else cat("No covariate\n")
  }
)

setMethod("showall","SaemixModel",
  function(object) {
    cat("Nonlinear mixed-effects model\n")
    distrib<-c("normal","log-normal","probit","logit")
    cat("  Model function")
    if(length(object@description)>0 && nchar(object@description)>0) cat(": ",object@description)
    cat("\n")
    print(object@model)
    cat("  Nb of parameters:",object@nb.parameters,"\n")
    cat("      parameter names: ",object@name.modpar,"\n")
    if(length(object@name.fixed)>0) cat("      fixed parameters: ",object@name.fixed,"\n")
    if(length(object@name.random)>0) cat("      random parameters: ",object@name.random,"\n")
    if(length(object@name.res)>0) cat("      parameters of residual variability: ",object@name.res,"\n")
    if(length(object@name.predictors)>0) cat("      predictors: ",object@name.predictors,"\n")
    if(length(object@name.X)>0) cat("      X predictor: ",object@name.X,"\n")
    if(length(object@name.cov)>0) cat("      covariates: ",object@name.cov,"\n")
    cat("      distribution:\n")
    tab<-cbind(Parameter=object@name.modpar, Distribution=distrib[object@transform.par+1], Estimated=ifelse(object@fixed.estim==1, "Estimated","Fixed"))
    print(tab,quote=FALSE)
    cat("  Variance-covariance matrix:\n")
    tab<-object@covariance.model
    print(tab,quote=FALSE)
    cat("  Initial estimate for variance-covariance matrix:\n")
    print(object@omega.init)
    st1<-paste(object@name.res,object@error.init,sep="=")
    cat("  Error model:",object@error.model,", initial values:",st1[object@indx.res],"\n")
   if(dim(object@covariate.model)[1]>0) {
      cat("  Covariate model:")
      if(sum(object@covariate.model)==0) cat(" none\n") else {
        cat("\n")
        print(object@covariate.model)
    }
  } else cat("  No covariate in the model.\n")
    cat("    Initial values\n")
    print(object@psi0)
    if(length(object@indx.fix)>0) cat("      index for fixed parameters: ", object@indx.fix,"\n")
    if(length(object@indx.cov)>0) cat("      index for covariate parameters: ", object@indx.cov,"\n")
    if(length(object@indx.omega)>0) cat("      index for random parameters: ", object@indx.omega,"\n")
    if(length(object@indx.res)>0) cat("      index for parameters of residual variability: ", object@indx.res,"\n")
    if(length(object@Mcovariates)>0) print(object@Mcovariates)
  }
)


####################################################################################
####				Summary method for SaemixModel			####
####################################################################################

setMethod("summary","SaemixModel",
  function(object, print=TRUE, ...) {
    if(print) {
    	cat("Nonlinear mixed-effects model\n")
    cat("  Model function")
    if(length(object@description)>0 && nchar(object@description)>0) {
      cat(": ",object@description,"\n")}
    else {
      cat("\n")
      print(object@model)
    }
    fix1<-ifelse(object@fixed.estim==1,""," [fixed]")
    cat("    ",object@nb.parameters,"parameters:", paste(object@name.modpar,fix1,sep=""),"\n")
    cat("     error model:",object@error.model,"\n")
    if(dim(object@covariate.model)[1]>0) {
      cat("     covariate model:\n")
      print(object@covariate.model) 
    } else cat("No covariate\n")
    }
    distrib<-c("normal","log-normal","probit","logit")
    tab.par<-data.frame(Parameter=object@name.modpar, Distribution=distrib[object@transform.par+1], Estimated=ifelse(as.numeric(object@betaest.model[1,])==1,"estimated","fixed"), Initial.value=object@psi0[1,])
    tab.res<-data.frame(parameters=object@name.res,Initial.value=object@error.init)   
    res<-list(model=list(model.function=object@model, error.model=object@error.model),parameters=list(fixed=tab.par, residual.error=tab.res),covariance.model=object@covariance.model, covariate.model=object@covariate.model)
    invisible(res)
 }
)

####################################################################################
####			SaemixModel class - method to plot			####
####################################################################################

# Plot simulations from the model

setMethod("plot","SaemixModel",
  function(x,y,range=c(0,1),psi,predictors,...) {
    if(missing(psi)) psi<-x@psi0[1,]
    psi<-matrix(psi,nrow=1)
    npred<-length(x@name.predictors)
    if(npred==0 & missing(predictors)) npred<-1 else {
      if(npred==0 & !missing(predictors)) {
        npred<-1+length(predictors)
      } else {
        if(npred>1 & (missing(predictors) || length(predictors)<(npred-1))) {
        cat("Please provide the value of the predictors other than X\n")
        return()
      }
     }
    }
    npts<-100
    psi<-matrix(rep(psi,npts+1),byrow=T,nrow=(npts+1))
    id<-matrix(rep(1,npts+1),ncol=1)
    xval<-range[1]+(range[2]-range[1])*c(0:100)/100
    if(npred==1) {
      xdep<-matrix(xval,ncol=1)
    } else {
      xdep<-cbind(xval,matrix(rep(predictors[1:(npred-1)],(npts+1)), byrow=T,nrow=(npts+1)))
      if(length(x@name.X)>0) {
        colnames(xdep)<-c(x@name.X,x@name.predictors[x@name.predictors!=x@name.X])
        xdep<-xdep[,match(x@name.predictors,colnames(xdep))]
      } else colnames(xdep)<-paste("Predictor",1:npred)
    }
    ypred<-try(x@model(psi,id,xdep))
    if(!is.numeric(ypred)) {
      cat("Problem when attempting to obtain predictions from the model.\n")
      cat("Usage: plot(x,range=c(0,1),psi,predictors) \n")
      cat("Possible solutions can be:\n")
      cat("   1. provide suitable values for X (option range=c(<lower bound>, <upper bound>))\n")
      cat("   2. provide values for additional predictors (option predictors=c(<value for predictor 1>, <value for predictor 2>, ...)).\n")
      cat("   3. check values for the model parameters (defaults to component psi0[1,] of the model).\n")
      cat("   4. the predictor used the X-axis is assumed to be in the first column; please check your model is written in a compatible way.\n")
    } else {
      if(length(x@name.X)==0 | length(x@name.predictors)==0) cat("Warning: X predictor supposed to be on the first axis\n")
      cat("Plot characteristics:\n")
      if(npred>1) {
        for(j in 1:dim(xdep)[2]) {
    if(length(x@name.X)==0) {
      if(j>1) cat("   predictor:",colnames(xdep)[j],"=",xdep[1,j],"\n")
    } else {
      if(colnames(xdep)[j]!=x@name.X) cat("    predictor:",colnames(xdep)[j],"=",xdep[1,j],"\n")
    }
      }}
      cat("   range for X-axis:",min(xval),"-",max(xval),"\n")
      cat("   parameters used in the simulation:", paste(x@name.modpar,"=",psi[1,],collapse=", "),"\n")
      plot(xval,ypred,type="l",xlab=ifelse(length(x@name.X)==0, "X",x@name.X),ylab=ifelse(length(x@name.response)==0, "Response",x@name.response))
    }
  }
)

####################################################################################
####			SaemixModel class - User-level function			####
####################################################################################

saemixModel<-function(model,psi0,description="",error.model=character(), transform.par=numeric(),fixed.estim=numeric(),covariate.model=matrix(nrow=0,ncol=0), covariance.model=matrix(nrow=0,ncol=0),omega.init=matrix(nrow=0,ncol=0),error.init=numeric(), name.modpar=character()) {
# Creating model from class
  if(missing(model)) {
    cat("Error in saemixModel:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
    return("Creation of SaemixModel failed")  
  }
  xcal<-try(typeof(model))
  if(class(xcal)=="try-error") {
    cat("Error in saemixModel:\n   the model function does not exist.\n")
    return("Creation of SaemixModel failed")  
  }
  if(typeof(model)=="character") {
    if(exists(model)) model<-get(model) else {
      cat("Error in saemixModel:\n   The argument model to saemixModel must be a valid function.\n")
      return("Creation of SaemixModel failed")
    }
  }
  if(!is.function(model)) {
    cat("Error in saemixModel:\n   The argument model to saemixModel must be a valid function.\n")
    return("Creation of SaemixModel failed")
  }
  if(length(formals(model))!=3) {
    cat("Error in saemixModel:\n   The model must be a function, accepting 3 arguments: psi (a vector of parameters), id (a vector of indices) and xidep (a matrix of predictors). Please see the documentation for examples.\n")
    return("Creation of SaemixModel failed")
  }
  if(missing(psi0) || length(psi0)==0) {
    cat("Error in saemixModel:\n   please provide initial estimates psi0 for at least the fixed effects.\n")
    return("Creation of SaemixModel failed")  
  }
  if(is.null(dim(psi0))) {
    psi1<-matrix(psi0,nrow=1)
    if(!is.null(names(psi0))) colnames(psi1)<-names(psi0)
    psi0<-psi1
    cat("Warning: psi0 given as a vector, reshaping it to a matrix.\n")
  }
  if(is.null(colnames(psi0))) {
    cat("Warning: no names given for the parameters in the model, please consider including parameter names.\n")
  }
  xmod<-try(new(Class="SaemixModel",model=model,description=description,psi0=psi0, error.model=error.model, transform.par=transform.par,fixed.estim=fixed.estim, covariate.model=covariate.model,covariance.model=covariance.model, omega.init=omega.init,error.init=error.init,name.modpar=name.modpar))
  if(class(xmod)=="SaemixModel") x1<-try(validObject(xmod),silent=FALSE) else x1<-xmod
  if(class(x1)!="try-error") cat("\n\nThe following SaemixModel object was successfully created:\n\n") else xmod<-"Creation of SaemixModel failed"
  print(xmod)
  return(xmod)
}

####################################################################################

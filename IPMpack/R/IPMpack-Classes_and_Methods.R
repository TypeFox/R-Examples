# TODO: Add comment
# 
# Author: ctg
###############################################################################


### FUNCTIONS FOR PREDICTING GROWTH AND SURVIVAL ###########################
### AS A FUNCTION OF  SIZE, COV, AND GROWTH AND SURV OBJECTS ################

# Register a growth  generic 2
#
# Parameters - (testing update)
#   size - size measurement this time-step
#   sizeNext - size measurement the next time-step (midpoints)
#   cov - the covariate (light, etc) this time-step
#   grow.obj - the growth object 
#
#   Note "growth" takes a SINGLE value for covariate although can take a VECTOR 
#   for size,sizeNext
#
#
# Returns -
#   the growth transition from size to sizeNext
setGeneric("growth",
    function(size, sizeNext, cov, growthObj) standardGeneric("growth"))


# Register a survival  generic
#
#Parameters -
#   size - size measurement this time-step
#   cov - the covariate (light env, etc)
#   surv.obj - the survival object
#
#   Note "surv" takes a SINGLE value for covariate although can take a VECTOR 
#   for size
#
#Returns -
#  The survival probability at size
setGeneric("surv",
    function(size, cov, survObj) standardGeneric("surv"))


# Register a cumulative growth  generic
#
# Parameters -
#   size - size measurement this time-step
#   sizeNext - size measurement the next time-step (midpoints)
#   cov - the covariate (light, etc) this time-step
#   grow.obj - the growth object 
#
#   Note "growth" takes a SINGLE value for covariate although can take a VECTOR 
#   for size,sizeNext
#
#
# Returns -
#   the growth transition from size to sizeNext
setGeneric("growthCum",
    function(size, sizeNext, cov, growthObj) standardGeneric("growthCum"))


### CLASSES UNDERLYING THE GROWTH / SURVIVAL FUNCTIONS ######

## GROWTH OBJECTS ##
# Create a generic growth object containing a lm
setClass("growthObj",
    representation(fit = "lm", sd = "numeric"))

setClass("growthObjPois",
    representation(fit = "glm"))

# Create a generic growth object with normal errors on increment
setClass("growthObjIncr",
    representation(fit = "lm", sd = "numeric"))

# Create a generic growth object with truncated normal errors on increment
setClass("growthObjTruncIncr",
    representation(fit = "list", varcov="matrix"))

# Create a generic growth object with log normal errors on increment
setClass("growthObjLogIncr",
    representation(fit = "lm", sd = "numeric"))

# Create a generic growth object with declining errors 
setClass("growthObjDeclineVar",
    representation(fit = "list"))

# Create a generic growth object with declining errors for increment
setClass("growthObjIncrDeclineVar",
    representation(fit = "list"))

# Create a generic growth object with declining errors for logincrement
setClass("growthObjLogIncrDeclineVar",
    representation(fit = "list"))

# Create a generic growth object containing the Hossfeld parameters 
setClass("growthObjHossfeld",
    representation(paras="numeric",
        sd="numeric", 
        logLik="numeric", 
        hessian="matrix"))

# Create a generic growth object with a poisson model 
setClass("growthObjPois",
    representation(fit = "glm"))

setClass("growthObjNegBin",
    representation(fit = "list"))

## SURVIVAL OBJECTS ##
# Create a generic survival object
setClass("survObj",
    representation(fit = "glm"))

# Create a generic survival object for use where over-dispersion
# modeled, using Diggles approximate correction for the transform
setClass("survObjOverDisp",
    representation(fit = "glm"))



## FECUNDITY OBJECTS ##
# Create a generic fecundity object
setClass("fecObj",
    representation(fitFec = "list",
        fecNames = "character",
        fecConstants = "data.frame",
        offspringSplitter = "data.frame",
        fecByDiscrete = "data.frame",
        vitalRatesPerOffspringType = "data.frame",
        Transform = "character",
        offspringRel = "lm",
        sdOffspringSize = "numeric")
)



# =============================================================================
# =============================================================================
## DISCRETE TRANSITION MATRIX OBJECTS ##
# Create a generic discrete transition matrix object
setClass("discreteTrans",
    representation(discreteTrans = "matrix",
        meanToCont = "matrix",
        sdToCont = "matrix",
        moveToDiscrete = "glm"))

# =============================================================================
# =============================================================================
## INTEGER FECUNDITY OBJECTS ##
# Create a generic fecundity object
setClass("fecObjInteger",
    representation(fitFec = "list",
        fecNames = "character",
        fecConstants = "data.frame",
        offspringSplitter = "data.frame",
        fecByDiscrete = "data.frame",
        vitalRatesPerOffspringType = "data.frame",
        Transform = "character",
        offspringRel = "glm",
        thetaOffspringSize = "numeric",
        distOffspring = "character")
)

# =============================================================================
# =============================================================================
## DISCRETE TRANSITION MATRIX INTEGER OBJECTS ##
# Create a generic discrete transition matrix object
setClass("discreteTransInteger",
    representation(discreteTrans = "matrix",
        meanToCont = "matrix",
        thetaToCont = "matrix",
        moveToDiscrete = "glm",
        distToCont = "character"))





######## DEFINE METHODS #######################################################
# =============================================================================
# =============================================================================
# SURVIVAL METHODS
# Method to obtain probability of survival using
# logistic regression on size with a single covariate
#
#Parameters -
#   size = current size (vector)
#   cov = current discrete covariate (.e.g. light..., single value)
#   survObj = a survival object, containig e.g. a glm fitted to data
#
#Returns -
#  survival probability for given sizes and covariate level
setMethod("surv", 
    c("numeric","data.frame","survObj"),
    function(size,cov,survObj){
      
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              survObj@fit$formula))>0) newd$expsize=exp(size)
      if (length(grep("logsize",
              survObj@fit$formula))>0) newd$logsize=log(size)
      if (length(grep("logsize2",
              survObj@fit$formula))>0) newd$logsize2=(log(size))^2
      
#				print(head(newd))
      
      u <- predict(survObj@fit,newd,type="response")
      return(u);
    })

# =============================================================================
# =============================================================================
# Method to obtain probability of survival using
# logistic regression on size with a single covariate
#  where the logistic regression was modeled with over-dispersion
#  (e.g., using MCMCglmm) - !over-dispersion assumed to be set to 1
#
#Parameters -
#   size = current size (vector)
#   cov = current discrete covariate (.e.g. light..., single value)
#   survObj = a survival object, containig e.g. a glm fitted to data
#
#Returns -
#  survival probability for given sizes and covariate level
setMethod("surv", 
    c("numeric","data.frame","survObjOverDisp"),
    function(size,cov,survObj){
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              survObj@fit$formula))>0) newd$expsize=exp(size)
      if (length(grep("logsize",
              survObj@fit$formula))>0) newd$logsize=log(size)
      if (length(grep("logsize2",
              survObj@fit$formula))>0) newd$logsize2=(log(size))^2
      
      u <- predict(survObj@fit,newd,type="link")
      c2 <- ((16 * sqrt(3))/(15 * pi))^2  #from MCMCglmm course , search for c2
      u <- invLogit(u/sqrt(1+c2)) 
      return(u);
    })

# =============================================================================
# =============================================================================
# GROWTH METHODS
# Method to obtain growth transitions -
# here linear regression to size (with powers) with a single covariate
#
#Parameters -
#   size = size now (vector)
#   sizeNext = size going to  (vector)
#   cov = the covariate (.e.g. light..., single value)
#   growthObj = a growth object
#
#Returns -
#  growth transition probability from size to sizeNext at that covariate level 
setMethod("growth", 
    c("numeric","numeric","data.frame","growthObj"),
    function(size,sizeNext,cov,growthObj){
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      if (length(grep("expsize",
              growthObj@fit$formula))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              growthObj@fit$formula))>0) { newd$logsize <- log(size)}
      
      mux <- predict(growthObj@fit,newd,type="response")
      sigmax <-growthObj@sd
      u <- dnorm(sizeNext,mux,sigmax,log=FALSE)  
      return(u);
    })

# =============================================================================
# =============================================================================
#  growth transition (poisson model) probability from size to sizeNext at that 
#  covariate level 
#	NOTE DO NOT USE THIS WITH AN IPM!!
setMethod("growth", 
    c("numeric","numeric","data.frame","growthObjPois"),
    function(size,sizeNext,cov,growthObj){
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              growthObj@fit$formula))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              growthObj@fit$formula))>0) { newd$logsize <- log(size)}
      
      mux <- predict(growthObj@fit,newd,type="response")
      u <- dpois(sizeNext,mux,log=FALSE)  
      return(u);
    })

# =============================================================================
# =============================================================================
#  growth transition (poisson model) probability from size to sizeNext at that 
#  covariate level 
#	NOTE DO NOT USE THIS WITH AN IPM!!
setMethod("growth", 
    c("numeric","numeric","data.frame","growthObjNegBin"),
    function(size,sizeNext,cov,growthObj){
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              growthObj@fit$formula))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              growthObj@fit$formula))>0) { newd$logsize <- log(size)}
      
      mux <- predict(growthObj@fit[[1]],newd,type="response")
      u <- dnbinom(sizeNext,mu=mux,size=growthObj@fit[[2]],log=FALSE)  
      return(u);
    })

# =============================================================================
# =============================================================================
# growth for predicting next incr 
setMethod("growth", 
    c("numeric","numeric","data.frame","growthObjIncr"),
    function(size,sizeNext,cov,growthObj){
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              growthObj@fit$formula))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              growthObj@fit$formula))>0) newd$logsize=log(size)
      
      mux <- predict(growthObj@fit,newd,type="response")
      
      #print(mux)
      
      sigmax <-growthObj@sd
      u <- dnorm(sizeNext,size+mux,sigmax,log=FALSE)  
      return(u); 
    })


# =============================================================================
# =============================================================================
# growth for predicting next truncated incr 
setMethod("growth", 
    c("numeric","numeric","data.frame","growthObjTruncIncr"),
    function(size,sizeNext,cov,growthObj){
      require(truncnorm)
      
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              growthObj@fit$formula))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
      
      mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
      sigmax <- growthObj@fit$sigma
      u <- dtruncnorm(sizeNext,a=size,b=Inf,mean=size+mux,sd=sigmax)  
      return(u); 
    })


# =============================================================================
# =============================================================================
# growth for predicting next logincr 
setMethod("growth", 
    c("numeric","numeric","data.frame","growthObjLogIncr"),
    function(size,sizeNext,cov,growthObj){
      
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              growthObj@fit$formula))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              growthObj@fit$formula))>0) newd$logsize=log(size)
      
      mux <- predict(growthObj@fit,newd,type="response")
      sigmax <-growthObj@sd
      u <- dlnorm(sizeNext-size,mux,sigmax,log=FALSE)  
      return(u);
    })

# =============================================================================
# =============================================================================
# growth for predicting next logincr with decline var
setMethod("growth", 
    c("numeric", "numeric", "data.frame", "growthObjLogIncrDeclineVar"),
    function(size,sizeNext,cov,growthObj){
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
      
      mux <- .predictMuX(grObj = growthObj, newData = newd, covPred = cov)
      sigmax2 <- growthObj@fit$sigmax2
      var.exp.coef <- growthObj@fit$var.exp.coef
      sigmax2 <- sigmax2 * exp(2 * (var.exp.coef * mux));
      
      u <- dlnorm(sizeNext - size, mux, sqrt(sigmax2), log = FALSE)  
      return(u);
    })

# =============================================================================
# =============================================================================
# growth for predicting next logincr with decline var
setMethod("growthCum", 
    c("numeric", "numeric", "data.frame", "growthObjLogIncrDeclineVar"),
    function(size, sizeNext, cov, growthObj){
      
      newd <- data.frame(cbind(cov, size = size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
      if (length(grep("logsize", names(growthObj@fit$coefficients))) > 0) {
          newd$logsize = log(size)
      }
      mux <- .predictMuX(grObj = growthObj, newData=newd, covPred = cov)
      sigmax2 <- growthObj@fit$sigmax2
      var.exp.coef <- growthObj@fit$var.exp.coef
      sigmax2 <- sigmax2 * exp(2 * (var.exp.coef * mux));
      
      u <- plnorm(sizeNext - size, mux, sqrt(sigmax2), log.p =  FALSE)  
      return(u);
    })


# =============================================================================
# =============================================================================
# Slightly alternative approach, using cumulative probs ##

# growth for predicting next size with a polynomial or log
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
    c("numeric","numeric","data.frame","growthObj"),
    function(size, sizeNext, cov, growthObj){
      
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      if (length(grep("expsize",
              growthObj@fit$formula))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              growthObj@fit$formula))>0) newd$logsize=log(size)
      mux <- predict(growthObj@fit,newd,type="response")
      sigmax <-growthObj@sd
      u <- pnorm(sizeNext,mux,sigmax, log.p =FALSE)  
      return(u);
    })

# =============================================================================
# =============================================================================
# growth for predicting next incr with a polynomial or log
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
    c("numeric","numeric","data.frame","growthObjIncr"),
    function(size,sizeNext,cov,growthObj){
      
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      if (length(grep("expsize",
              growthObj@fit$formula))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              growthObj@fit$formula))>0) newd$logsize=log(size)
      
      mux <- predict(growthObj@fit,newd,type="response")
      sigmax <-growthObj@sd
      u <- pnorm(sizeNext,size+mux,sigmax, log.p =FALSE)  
      return(u); 
    })

# =============================================================================
# =============================================================================
# growth for predicting next truncated incr with cumulative 
setMethod("growthCum", 
    c("numeric","numeric","data.frame","growthObjTruncIncr"),
    function(size,sizeNext,cov,growthObj){
      require(truncnorm)
      
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
      
      mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
      sigmax <- sqrt(growthObj@fit$sigmax)
      
      u <- ptruncnorm(sizeNext,a=size,b=Inf,mean=size+mux,sd=sigmax)  
      return(u); 
    })

# =============================================================================
# =============================================================================
# growth for predicting next logincr with a polynomial or log
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
    c("numeric","numeric","data.frame","growthObjLogIncr"),
    function(size, sizeNext, cov, growthObj){
      
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      if (length(grep("expsize",
              names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
      
      mux <- predict(growthObj@fit,newd,type="response")
      sigmax <-growthObj@sd
      u <- plnorm(sizeNext-size,mux,sigmax,log.p = FALSE)  
      return(u);
    })

# =============================================================================
# =============================================================================
#Simple growth methods, using  declining variance in growth
# using pnorm (i.e. getting cumulative at boundary points and doing difference)
setMethod("growthCum", 
    c("numeric","numeric","data.frame","growthObjDeclineVar"),
    function(size,sizeNext,cov,growthObj){
      
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              names(growthObj@fit$coefficients))) > 0) newd$logsize=log(size)
      
      mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
      sigmax2 <- growthObj@fit$sigmax2
      var.exp.coef<-growthObj@fit$var.exp.coef
      sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
      u <- pnorm(sizeNext,mux,sqrt(sigmax2),log.p =FALSE)  
      return(u);
    })

setMethod("growthCum", 
		c("numeric","numeric","data.frame","growthObjIncrDeclineVar"),
		function(size, sizeNext, cov, growthObj){
			newd <- data.frame(cbind(cov,size=size),
					stringsAsFactors = FALSE)
			newd$size2 <- size^2
			newd$size3 <- size^3
			
			if (length(grep("expsize",
							names(growthObj@fit$coefficients))) > 0) newd$expsize <- exp(size)
			if (length(grep("logsize",
							names(growthObj@fit$coefficients))) > 0) newd$logsize = log(size)
			
			mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
			sigmax2 <- growthObj@fit$sigmax2
			var.exp.coef<-growthObj@fit$var.exp.coef
			sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
			u <- pnorm(sizeNext, size + mux, sqrt(sigmax2), log.p =  FALSE)  
			return(u);
		})


# =============================================================================
# =============================================================================
#Simple growth methods, using  declining variance in growth
setMethod("growth", 
    c("numeric","numeric","data.frame","growthObjDeclineVar"),
    function(size,sizeNext,cov,growthObj){
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              names(growthObj@fit$coefficients)))>0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              names(growthObj@fit$coefficients)))>0) newd$logsize=log(size)
      
      mux <- .predictMuX(grObj=growthObj,newData=newd)
      sigmax2 <- growthObj@fit$sigmax2
      var.exp.coef<-growthObj@fit$var.exp.coef
      sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
      u <- dnorm(sizeNext,mux,sqrt(sigmax2),log=FALSE)  
      return(u);
    })

# =============================================================================
# =============================================================================
#same but with declining variance in growth on incrment
setMethod("growth", 
    c("numeric","numeric","data.frame","growthObjIncrDeclineVar"),
    function(size, sizeNext, cov, growthObj){
      newd <- data.frame(cbind(cov,size=size),
          stringsAsFactors = FALSE)
      newd$size2 <- size^2
      newd$size3 <- size^3
      
      if (length(grep("expsize",
              names(growthObj@fit$coefficients))) > 0) newd$expsize <- exp(size)
      if (length(grep("logsize",
              names(growthObj@fit$coefficients))) > 0) newd$logsize = log(size)
      
      mux <- .predictMuX(grObj=growthObj,newData=newd,covPred=cov)
      sigmax2 <- growthObj@fit$sigmax2
      var.exp.coef<-growthObj@fit$var.exp.coef
      sigmax2<-sigmax2*exp(2*(var.exp.coef*mux));
      u <- dnorm(sizeNext, size + mux, sqrt(sigmax2), log =  FALSE)  
      return(u);
    })

# =============================================================================
# =============================================================================
## Define a new growth method for Hossfeld growth (classic midpoint rule)
setMethod("growth", c("numeric", "numeric", "data.frame", "growthObjHossfeld"), 
    function(size, sizeNext, cov, growthObj) { 
      mux <- size+Hossfeld(size, growthObj@paras) 
      sigmax <- growthObj@sd 
      u <- dnorm(sizeNext, mux, sigmax, log =  FALSE) 
      return(u)
    }) 

# =============================================================================
# =============================================================================
## Define a new growth method for Hossfeld growth (classic midpoint rule)
setMethod("growthCum", c("numeric", "numeric", "data.frame", 
                "growthObjHossfeld"), 
    function(size, sizeNext, cov, growthObj) { 
      mux <- size+Hossfeld(size, growthObj@paras) 
      sigmax <- growthObj@sd 
      u <- pnorm(sizeNext, mux, sigmax, log.p =  FALSE) 
      return(u)
    }) 

### CLASSES AND FUNCTIONS FOR MATRICES (ENV, TMATRIX [compound or not], FMATRIX) 
# =============================================================================
# =============================================================================
#Class for the matrix that holds the env matrix 
setClass("envMatrix",
    representation(nEnvClass = "numeric"), #number of covariate levels
    contains="matrix")

# =============================================================================
# =============================================================================
#Class for the matrix that holds the IPM
setClass("IPMmatrix",
    representation(nDiscrete = "numeric", #number of discrete classes
        nEnvClass = "numeric", #number of covariate levels
        nBigMatrix = "numeric", #the resolution of the IPM
        meshpoints = "numeric",
        env.index = "numeric",
        names.discrete = "character"),
    contains="matrix")

# =============================================================================
# =============================================================================
# Method combining growth and survival for doing outer (not a generic, so that 
# don't have to
# define for all the different classes / growth and surv take care of that)
growSurv <- function(size,sizeNext,cov,growthObj,survObj){
  growth(size, sizeNext, cov, growthObj) * surv(size,cov,survObj)
}

# =============================================================================
# =============================================================================
## Function defining growth using Hossfeld function
#
# Parameters - size - DBH
#            - par - three parameters
#
# Returns - growth increment (not log scale)
Hossfeld <- function(size,par) {
  deltDBH <- (par[2]*par[3]*size^(par[3]-1))/((par[2]+((size^par[3])/par[1]))^2)
  return(deltDBH)
}

# =============================================================================
# =============================================================================
## Function to fit Hossfeld function using optim 
#
# Parameters - par - three parameters
#            - dataf - a data-frame
#
# Returns the SS

wrapHossfeld <- function(par, dataf) { 
  pred <- Hossfeld(dataf$size, par[1:3]) 
  ss <- sum((pred - dataf$incr)^2, na.rm = T)
  return(ss) 
} 

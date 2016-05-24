###################################################
### chunk number 1: 
###################################################

algo.hmm <- function(disProgObj, control = list(range=range, Mtilde=-1, noStates=2, trend=TRUE, noHarmonics=1,covEffectEqual=FALSE, saveHMMs = FALSE, extraMSMargs=list() )){

  # check if the msm package is available
  if (!requireNamespace("msm")) {
      stop("the HMM method requires package ", sQuote("msm"))
  }
  
  # Set the default values if not yet set
  if(is.null(control$Mtilde)){ control$Mtilde <- -1 }
  if(is.null(control$noStates)){ control$noStates <- 2 }
  if(is.null(control$trend)){ control$trend <- TRUE }
  if(is.null(control$noHarmonics)){ control$noHarmonics <- 1 }
  if(is.null(control$covEffectEqual)){ control$covEffectEqual <- FALSE }
  if(is.null(control$saveHMMs)){ control$saveHMMs <- FALSE }
  if(is.null(control$extraMSMargs)){ control$extraMSMargs <- list() }

  #Stop if not enough for estimation
  if(min(control$range) < 2) {
    stop("Error: Too few values as reference values")
  }

  # initialize the necessary vectors
  alarm <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  upperbound <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  control$hmms <- list()


  ##############################################
  #Repeat for each time point to monitor on-line
  ############################################## 
  for (i in 1:length(control$range)) {
    #Function is so slow some sort of perfomance indicator is usually necessary
    cat(paste("i=",i," (out of ",length(control$range),")\n",sep=""))
    #Initialize observations for each round -- can be done sequentally
    first <- ifelse(control$Mtilde== -1, 1, max(control$range[i]-control$Mtilde+1,1))
    t <- first:control$range[i]
    observed <- disProgObj$observed[t]

    #Init data
    counts <- data.frame(observed, t)
    names(counts) <- c("observed","t")
    #Initialize formula
    formulaStr <- ifelse(control$trend, "~ 1 + t ", "~ 1 ")
    #Create formula and add harmonics as covariates -- this could be done recursively?
    for (j in seq_len(control$noHarmonics)) {
      counts[,paste("cos",j,"t",sep="")] <- cos(2*j*pi*(t-1)/disProgObj$freq)
      counts[,paste("sin",j,"t",sep="")] <- sin(2*j*pi*(t-1)/disProgObj$freq)
      formulaStr <- paste(formulaStr,"+ cos",j,"t + sin",j,"t ",sep="")
    }
  
    #Obtain crude inits
    q <- quantile(observed,seq(0,1,length=control$noStates+1))
    lvl <- cut(observed,breaks=q,include.lowest=TRUE)
    crudeMean <- as.numeric(tapply(observed, lvl, mean))

    hcovariates <- list()
    hmodel <- list()
    for (j in seq_len(control$noStates)) {
      hcovariates[[j]] <- as.formula(formulaStr)
      val <- crudeMean[j]
      #Substitution necessary, as hmmPois does lazy evaluation of rate argument
      hmodel[[j]] <- eval(substitute(msm::hmmPois(rate=val),list(val=crudeMean[j])))
    }

    #Any constraints on the parameters of the covariates for the different states
    hconstraint <- list()
    if (control$covEffectEqual) {
      hconstraint <- list(t=rep(1,control$noStates))
      for (j in seq_len(control$noHarmonics)) {
        hconstraint[[paste("sin",j,"t",sep="")]] <- rep(1,control$noStates)
        hconstraint[[paste("cos",j,"t",sep="")]] <- rep(1,control$noStates)
      }
    }

    #Prepare object for msm fitting
    msm.args <- list(formula = observed ~ t,
        data = counts,
        #HMM with "noStates" states having equal initial values
        qmatrix = matrix(1/control$noStates,control$noStates,control$noStates),
        #y|x \sim Po( \mu[t] ) with some initial values
        hmodel = hmodel,
        #Models for \log \mu_t^1 and \log \mu_t^2
        hcovariates = hcovariates,
        #Force the effects of the trend and harmonics to be equal for all states
        hconstraint=hconstraint)
    #Add additional msm arguments
    msm.args <- modifyList(msm.args, control$extraMSMargs)

    # fit the HMM
    hmm <- do.call(what=msm::msm, args=msm.args)
    
    #In case the model fits should be saved.
    if (control$saveHMMs) {
      control$hmms[[i]] <- hmm
    }

    #If most probable state of current time point (i.e. last obs) equals the 
    #highest state then do alarm 
#    print(observed)
#    print(matrix(viterbi.msm(hmm)$fitted,ncol=1))
    alarm[i] <- msm::viterbi.msm(hmm)$fitted[length(t)] == control$noStates

    #Upperbound does not have any meaning -- compute posterior probability!
    upperbound[i] <- 0
  
   
  }

  #Add name and data name to control object.
  control$name <- paste("hmm:", control$trans)
  control$data <- paste(deparse(substitute(disProgObj)))
  #no need for hmm object -- control$hmm  <- hmm

  # return alarm and upperbound vectors
  result <- list(alarm = alarm, upperbound = upperbound, disProgObj=disProgObj,control=control)

  class(result) = "survRes" # for surveillance system result
  return(result)
}




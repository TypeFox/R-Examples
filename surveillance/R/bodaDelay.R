#     ____________________________
#    |\_________________________/|\
#    ||                         || \
#    ||    bodaDelay            ||  \
#    ||                         ||  |
#    ||                         ||  |
#    ||                         ||  |
#    ||                         ||  |
#    ||                         ||  |
#    ||                         ||  /
#    ||_________________________|| /
#    |/_________________________\|/
#       __\_________________/__/|_
#      |_______________________|/ )
#    ________________________    (__
#   /oooo  oooo  oooo  oooo /|   _  )_
#  /ooooooooooooooooooooooo/ /  (_)_(_)
# /ooooooooooooooooooooooo/ /    (o o)
#/C=_____________________/_/    ==\o/==

# Author: M.Salmon

################################################################################
# CONTENTS 
################################################################################
# # MAIN FUNCTION
# Function that manages input and output.
# # FIT GLM FUNCTION
# Function that fits a GLM. 
# # THRESHOLD FUNCTION
# Function that calculates the threshold.
# # DATA GLM FUNCTION
# Function that prepares data for the GLM.
# # FORMULA FUNCTION
# Function that writes the formula for the GLM.
################################################################################
# END OF CONTENTS
################################################################################
################################################################################
# MAIN FUNCTION
################################################################################

bodaDelay <- function(sts, control = list(range = NULL, b = 3, w = 3,
                                          mc.munu=100, mc.y=10,
                                          pastAberrations = TRUE, 
                                          verbose = FALSE,
                                          alpha = 0.01, trend = TRUE,
                                          limit54=c(5,4), 
                                          inferenceMethod=c("asym","INLA"),
                                          quantileMethod=c("MC","MM"),
                                          noPeriods = 1, 
                                          pastWeeksNotIncluded = 26,
                                          delay = TRUE)) {
  
  ######################################################################
  # Use special Date class mechanism to find reference months/weeks/days
  ######################################################################  
  if (is.null( sts@epochAsDate)) {
    epochAsDate <- FALSE
  } else {
    epochAsDate <-    sts@epochAsDate
  }
  
  ######################################################################
  # Fetch observed and population
  ######################################################################  
  # Fetch observed
  observed <- observed(sts)
  freq <- sts@freq
  if (epochAsDate) {
    epochStr <- switch( as.character(freq), "12" = "month","52" =    "week",
                        "365" = "day")
  } else { 
    epochStr <- "none"
  }
  
  # Fetch population (if it exists)
  if (!is.null(population(sts))) {
    population <- population(sts) 
  } else {
    population <- rep(1,length(observed))
  }
  
  ######################################################################
  # Fix missing control options
  ######################################################################
  if (is.null(control[["b",exact=TRUE]])) { control$b = 5 }
  
  if (is.null(control[["w", exact = TRUE]])) { control$w = 3 }
  
  if (is.null(control[["range", exact=TRUE]])) {
    control$range <- (freq*(control$b)+control$w +1):dim(observed)[1] 
  }
  if (is.null(control[["pastAberrations",exact=TRUE]])) {control$pastAberrations=TRUE}
  
  if (is.null(control[["verbose",exact=TRUE]]))    {control$verbose=FALSE}
  
  if (is.null(control[["alpha",exact=TRUE]]))        {control$alpha=0.05}
  
  if (is.null(control[["trend",exact=TRUE]]))        {control$trend=TRUE}
  
  
  # No alarm is sounded
  #    if fewer than cases = 5 reports were received in the past period = 4
  #    weeks. limit54=c(cases,period) is a vector allowing the user to change
  #    these numbers
  
  if (is.null(control[["limit54",exact=TRUE]]))    {control$limit54=c(5,4)}
  
  if (is.null(control[["noPeriods",exact=TRUE]])){control$noPeriods=1}
  
  # Use factors in the model? Depends on noPeriods, no input from the user.
  if (control$noPeriods!=1) {
    control$factorsBool=TRUE
  } else {
    control$factorsBool=FALSE
  }
  
  
  # How many past weeks not to take into account?
  if (is.null(control[["pastWeeksNotIncluded",exact=TRUE]])){ 
    control$pastWeeksNotIncluded=control$w
  }
  
  # Correct for delays?
  if (is.null(control[["delay",exact=TRUE]])) { control$delay = FALSE }
  # Reporting triangle here?
  if (control$delay) {
    if (is.null( sts@control$reportingTriangle$n)) {stop("You have to provide a reporting triangle in control of the sts-object")}
    
    if (!(length(apply(sts@control$reportingTriangle$n,1,sum,na.rm=TRUE))==length(sts@observed)))
    {stop("The reporting triangle number of lines is not the length of the observed slot.")}  
    
    if (!(sum(apply(sts@control$reportingTriangle$n,1,sum,na.rm=TRUE)==sts@observed)==length(sts@observed)))
    {stop("The reporting triangle is wrong: not all cases are in the reporting triangle.")}
  }
  
  # setting for monte carlo integration
  if(is.null(control[["mc.munu",exact=TRUE]])){ control$mc.munu <- 100 }
  if(is.null(control[["mc.y",exact=TRUE]])){ control$mc.y <- 10 }
  ######################################################################
  # Check options
  ######################################################################
  if (!((control$limit54[1] >= 0) && (control$limit54[2] > 0))) {
    stop("The limit54 arguments are out of bounds: cases >= 0 and period > 0.")
  }
  # inference method
  if(is.null(control[["inferenceMethod",exact=TRUE]])){ 
    control$inferenceMethod <- "asym" }
  else {
    control$inferenceMethod <- match.arg(control$inferenceMethod, 
                                         c("asym","INLA"))
  }
  
  if(is.null(control[["quantileMethod",exact=TRUE]])){ 
    control$quantileMethod <- "MC" }
  else {
    control$quantileMethod <- match.arg(control$quantileMethod, 
                                        c("MC","MM"))
  }
  
  #Check if the INLA package is available.
  if (control$inferenceMethod=="INLA"){
    if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("The bodaDelay function requires the INLA package to be installed.\n",
           "  The package is not available on CRAN, but can be easily obtained\n",
           "  from <http://www.r-inla.org/download>.\n",
           "  Alternatively, set inferenceMethod to \"asym\".")
    }
  }
  
  # Define objects
  n <- control$b*(2*control$w+1)
  
  
  
  
  
  # loop over columns of sts
  
  
  #Vector of dates
  if (epochAsDate){
    vectorOfDates <- as.Date(sts@epoch, origin="1970-01-01")
  } else {
    vectorOfDates <- seq_len(length(observed))
  }    
  
  # Parameters
  b <- control$b
  w <- control$w
  noPeriods <- control$noPeriods
  verbose <- control$verbose
  reportingTriangle <- sts@control$reportingTriangle
  timeTrend <- control$trend
  alpha <- control$alpha
  factorsBool <- control$factorsBool
  pastAberrations <- control$pastAberrations
  glmWarnings <- control$glmWarnings
  delay <- control$delay
  k <- control$k
  verbose <- control$verbose
  pastWeeksNotIncluded <- control$pastWeeksNotIncluded
  mc.munu <- control$mc.munu
  mc.y <- control$mc.y
  # Loop over control$range
  for (k in control$range) {
    
    ######################################################################
    # Prepare data for the glm
    ######################################################################
    dayToConsider <- vectorOfDates[k]  		
    diffDates <- diff(vectorOfDates)
    delay <- control$delay
    
    dataGLM <- bodaDelay.data.glm(dayToConsider=dayToConsider, 
                                  b=b, freq=freq, 
                                  epochAsDate=epochAsDate,
                                  epochStr=epochStr,
                                  vectorOfDates=vectorOfDates,w=w,
                                  noPeriods=noPeriods,
                                  observed=observed,population=population,
                                  verbose=verbose,
                                  pastWeeksNotIncluded=pastWeeksNotIncluded,
                                  reportingTriangle=reportingTriangle, 
                                  delay=delay)  
    
    ######################################################################
    # Fit the model 
    ######################################################################      
    argumentsGLM <- list(dataGLM=dataGLM,reportingTriangle=reportingTriangle,
                         timeTrend=timeTrend,alpha=alpha,
                         factorsBool=factorsBool,pastAberrations=pastAberrations,
                         glmWarnings=glmWarnings,
                         verbose=verbose,delay=delay,
                         inferenceMethod=control$inferenceMethod)
    
    model <- do.call(bodaDelay.fitGLM, args=argumentsGLM)
    if(is.null(model)){
      sts@upperbound[k] <- NA
      sts@alarm[k] <- NA
    }
    else{
      ######################################################################
      # Calculate the threshold 
      ###################################################################### 
      quantileMethod <- control$quantileMethod
      argumentsThreshold <- list(model,alpha=alpha,dataGLM=dataGLM,reportingTriangle,
                                 delay=delay,k=k,control=control,mc.munu=mc.munu,mc.y=mc.y,
                                 inferenceMethod=control$inferenceMethod,
                                 quantileMethod=quantileMethod)
      
      threshold <- do.call(bodaDelay.threshold,argumentsThreshold)
      
      ######################################################################
      # Output results if enough cases 
      ###################################################################### 
      sts@upperbound[k] <- threshold
      enoughCases <- (sum(observed[(k-control$limit54[2]+1):k])
                      >=control$limit54[1])
      sts@alarm[k] <- FALSE
      if (is.na(threshold)){sts@alarm[k] <- NA}
      else {
        if (sts@observed[k]>sts@upperbound[k]) {sts@alarm[k] <- TRUE}
      }
      if(!enoughCases){
        sts@upperbound[k] <- NA
        sts@alarm[k] <- NA
      }
    }
 
  } #done looping over all time points
  
  return(sts[control$range,]) 
}
################################################################################
# END OF MAIN FUNCTION
################################################################################

################################################################################
# FIT GLM FUNCTION
################################################################################

bodaDelay.fitGLM <- function(dataGLM,reportingTriangle,alpha,
                             timeTrend,factorsBool,delay,pastAberrations,
                             glmWarnings,verbose,inferenceMethod,...) {
  # Model formula depends on whether to include a time trend or not.
  
  theModel <- formulaGLMDelay(timeBool=timeTrend,factorsBool,delay,outbreak=FALSE)
  
  if(inferenceMethod=="INLA"){
    
    E <- max(0,mean(dataGLM$response, na.rm=TRUE))
    link=1
    model <- INLA::inla(as.formula(theModel),data=dataGLM,
                        family='nbinomial',E=E,
                        control.predictor=list(compute=TRUE,link=link),
                        control.compute=list(cpo=TRUE,config=TRUE),
                        control.inla = list(int.strategy = "grid",dz=1,diff.logdens = 10),
                        control.family = list(hyper = list(theta = list(prior = "normal", param = c(0, 0.001)))))
    
    if (pastAberrations){
      # if we have failures => recompute those manually
      #if (sum(model$cpo$failure,na.rm=TRUE)!=0){
      #   model <- inla.cpo(model)
      #}
      # Calculate the mid p-value
      vpit <- model$cpo$pit
      vcpo <- model$cpo$cpo
      midpvalue <- vpit - 0.5*vcpo
      # Detect the point with a high mid p-value
      
      # outbreakOrNot <- midpvalue
      #outbreakOrNot[midpvalue  <= (1-alpha)] <- 0
      outbreakOrNot <- ifelse(midpvalue  > (1-alpha), 1, 0) 
      outbreakOrNot[is.na(outbreakOrNot)] <- 0# FALSE
      outbreakOrNot[is.na(dataGLM$response)] <- 0#FALSE
      
      # Only recompute the model if it will bring something!
      if (sum(outbreakOrNot)>0){
        dataGLM <- cbind(dataGLM,outbreakOrNot)
        theModel <- formulaGLMDelay(timeBool=timeTrend,factorsBool,delay,outbreak=TRUE)
        
        model <- INLA::inla(as.formula(theModel),data=dataGLM,
                            family='nbinomial',E=E,
                            control.predictor=list(compute=TRUE,link=link),
                            control.compute=list(cpo=FALSE,config=TRUE),
                            control.inla = list(int.strategy = "grid",dz=1,diff.logdens = 10),
                            control.family = list(hyper = list(theta = list(prior = "normal", param = c(0, 0.001)))))
        
        # if we have failures => recompute those manually
        #  if (sum(model$cpo$failure,na.rm=TRUE)!=0){model <- inla.cpo(model)}
        vpit <- model$cpo$pit
        vcpo <- model$cpo$cpo
        midpvalue <- vpit - 0.5*vcpo 
      }
    }
    
  }
  if (inferenceMethod=="asym"){
    
    model <- MASS::glm.nb(as.formula(theModel),data=dataGLM)
    if(!model$converged){
      return(NULL)
    }
  }
  
  return(model)
}
################################################################################
# END OF FIT GLM FUNCTION
################################################################################

################################################################################
# THRESHOLD FUNCTION
################################################################################
bodaDelay.threshold <- function(model, mc.munu,mc.y,alpha,
                                delay,k,control,dataGLM,reportingTriangle,
                                inferenceMethod,quantileMethod...) {
  quantileMethod <- control$quantileMethod
  if (inferenceMethod=="INLA"){
    E <- max(0,mean(dataGLM$response, na.rm=TRUE))
    # Sample from the posterior
    jointSample <- INLA::inla.posterior.sample(mc.munu,model,hyper.user.scale = FALSE)
    
    # take variation in size hyperprior into account by also sampling from it
    theta <- t(sapply(jointSample, function(x) x$hyperpar))
    
    if (delay){
      mu_Tt <- numeric(mc.munu)
      N_Tt <- numeric(mc.munu*mc.y)
      
      # Maximal delay + 1
      Dmax0 <- ncol(as.matrix(reportingTriangle$n))
      # The sum has to be up to min(D,T-t). This is how we find the right indices.
      loopLimit <- min(Dmax0,which(is.na(as.matrix(reportingTriangle$n)[k,]))-1,na.rm=TRUE)
      
      # Find the mu_td and sum
      for (d in 1:loopLimit)
      {
        if(sum(dataGLM$response[dataGLM$delay==d],na.rm=TRUE)!=0){
          mu_Tt <- mu_Tt + exp(t(sapply(jointSample, function(x) x$latent[[nrow(dataGLM)-Dmax0+d]])))
        }
        
      }
      
      # with no delay this is similar to boda.
    } else {
      mu_Tt <- exp(t(sapply(jointSample, function(x) x$latent[[nrow(dataGLM)]])))
      
    }
  }
  
  if (inferenceMethod=="asym"){
    E <- 1
    # Sample from the posterior
    set.seed(1)
    
    
    # take variation in size hyperprior into account by also sampling from it
    theta <- rnorm(n=mc.munu,mean=summary(model)$theta,sd=summary(model)$SE.theta)
    
    if (delay){
      # Maximal delay + 1
      Dmax0 <- ncol(as.matrix(reportingTriangle$n))
      mu_Tt <- numeric(mc.munu)
      newData <- tail(dataGLM,n=Dmax0)
      
      P=predict(model,type="link",se.fit=TRUE,
                newdata=newData)
      
      # The sum has to be up to min(D,T-t). This is how we find the right indices.
      loopLimit <- min(Dmax0,which(is.na(as.matrix(reportingTriangle$n)[k,]))-1,na.rm=TRUE)
      
      # Find the mu_td and sum
      
      for (d in 1:loopLimit)
      {
        if(sum(dataGLM$response[dataGLM$delay==d],na.rm=TRUE)!=0){
          mu_Tt <- mu_Tt + exp(rnorm(n=mc.munu,mean=P$fit[d],sd=P$se.fit[d]))
        }
        
      }
      
      # with no delay this is similar to boda.
    } else {
      
      newData <- tail(dataGLM,n=1)
      
      
      P=try(predict(model,type="link",se.fit=TRUE,
                    newdata=newData),silent=TRUE)
      if (class(P)=="try-error"){P<- NA
      return(NA)}
      set.seed(1)
      mu_Tt <- exp(rnorm(n=mc.munu,mean=P$fit,sd=P$se.fit))
      
    }
    
  }
  
  if(quantileMethod=="MC"){
    N_Tt <- rnbinom(n=mc.y*mc.munu,size=theta,mu=E*mu_Tt)
    # We have to ditch the na values (values for which theta was negative)
    N_Tt <- N_Tt[is.na(N_Tt)==FALSE]
    
    qi <- quantile(N_Tt, probs=(1-alpha), type=3, na.rm=TRUE)
  }
  if(quantileMethod=="MM"){
    mu_Tt <- mu_Tt[mu_Tt>=0&theta>0]
    
    theta <- theta[mu_Tt>=0&theta>0]
    
    minBracket <- qnbinom(p=(1-alpha), 
                          mu=E*min(mu_Tt),
                          size=max(theta))
    
    maxBracket <- qnbinom(p=(1-alpha), 
                          mu=E*max(mu_Tt),
                          size=min(theta))
    
    qi <- qmix(p=(1-alpha), mu=E*mu_Tt, size=theta,
               bracket=c(minBracket, maxBracket))
    

    
  }
  
  
  return(as.numeric(qi))
  
  
}
################################################################################
# END OF THRESHOLD GLM FUNCTION
################################################################################

################################################################################
# DATA GLM FUNCTION
################################################################################
bodaDelay.data.glm <- function(dayToConsider, b, freq, 
                               epochAsDate,epochStr,
                               vectorOfDates,w,noPeriods,
                               observed,population,
                               verbose,pastWeeksNotIncluded,reportingTriangle,delay){
  
  
  # Identify reference time points
  
  # Same date but with one year, two year, etc, lag
  # b+1 because we need to have the current week in the vector
  referenceTimePoints <- algo.farrington.referencetimepoints(dayToConsider,b=b,
                                                             freq=freq,
                                                             epochAsDate=epochAsDate,
                                                             epochStr=epochStr
  )
  
  if (sum((vectorOfDates %in% min(referenceTimePoints)) == rep(FALSE,length(vectorOfDates))) == length(vectorOfDates)){
    warning("Some reference values did not exist (index<1).")
  }
  
  
  # Create the blocks for the noPeriods between windows (including windows)
  # If noPeriods=1 this is a way of identifying windows, actually.
  
  blocks <- blocks(referenceTimePoints,vectorOfDates,epochStr,dayToConsider,
                   b,w,noPeriods,epochAsDate)
  
  # Here add option for not taking the X past weeks into account
  # to avoid adaptation of the model to emerging outbreaks
  blocksID <- blocks
  
  
  # Extract values for the timepoints of interest only
  
  blockIndexes <- which(is.na(blocksID)==FALSE) 
  
  
  # Time
  
  # if epochAsDate make sure wtime has a 1 increment
  if (epochAsDate){
    wtime <- (as.numeric(vectorOfDates[blockIndexes])-
                as.numeric(vectorOfDates[blockIndexes][1]))/as.numeric(diff(vectorOfDates))[1]
  } else {
    wtime <-     as.numeric(vectorOfDates[blockIndexes])
  }
  
  # Factors
  seasgroups <- as.factor(blocks[blockIndexes])
  
  # Observed
  response <- as.numeric(observed[blockIndexes])
  response[length(response)] <- NA
  # Population
  pop <- population[blockIndexes]
  
  if (verbose) { print(response)}
  
  
  # If the delays are not to be taken into account it is like farringtonFlexible
  if (!delay) {
    dataGLM <- data.frame(response=response,wtime=wtime,population=pop,
                          seasgroups=seasgroups,vectorOfDates=vectorOfDates[blockIndexes])
    
    dataGLM$response[(nrow(dataGLM)-pastWeeksNotIncluded):nrow(dataGLM)] <- NA
  }
  # If the delays are to be taken into account we need a bigger dataframe
  else {
    # Delays
    
    delays <- as.factor(0:(dim(reportingTriangle$n)[2]-1))
    
    # Take the subset of the reporting triangle corresponding to the timepoints used for fitting the model
    reportingTriangleGLM <- reportingTriangle$n[rownames(reportingTriangle$n) %in% as.character(vectorOfDates[blockIndexes]),]
    
    # All vectors of data will be this long: each entry will correspond to one t and one d
    lengthGLM <- dim(reportingTriangleGLM)[2]*dim(reportingTriangleGLM)[1]
    
    # Create the vectors for storing data
    responseGLM <- numeric(lengthGLM)
    wtimeGLM <- numeric(lengthGLM)
    seasgroupsGLM <- numeric(lengthGLM)
    popGLM <- numeric(lengthGLM)
    vectorOfDatesGLM <- numeric(lengthGLM)
    delaysGLM <- numeric(lengthGLM) 
    
    # Fill them D by D
    D <- dim(reportingTriangleGLM)[2]
    for (i in (1:dim(reportingTriangleGLM)[1])){
      vectorOfDatesGLM[((i-1)*D+1):(i*D)] <- rep(vectorOfDates[blockIndexes][i],D)
      wtimeGLM[((i-1)*D+1):(i*D)] <- rep(wtime[i],D)
      popGLM[((i-1)*D+1):(i*D)] <- rep(pop[i],D)
      seasgroupsGLM[((i-1)*D+1):(i*D)] <- rep(seasgroups[i],D)
      responseGLM[((i-1)*D+1):(i*D)] <- reportingTriangleGLM[i,]
      delaysGLM[((i-1)*D+1):(i*D)] <- 0:(D-1)
      
    }
    
    responseGLM[((i-1)*D+1):(i*D)] <- rep (NA, D)
    responseGLM[(length(responseGLM)-pastWeeksNotIncluded*D):length(responseGLM)] <- NA
    
    dataGLM <- data.frame(response=responseGLM,wtime=wtimeGLM,population=popGLM,
                          seasgroups=as.factor(seasgroupsGLM),vectorOfDates=as.Date(vectorOfDatesGLM,origin="1970-01-01"),delay=delaysGLM)
    
  }
  
  
  
  
  
  return(as.data.frame(dataGLM))
  
}

################################################################################
# END OF DATA GLM FUNCTION
################################################################################

################################################################################
# FORMULA FUNCTION
################################################################################
# Function for writing the good formula depending on timeTrend,
# and factorsBool

formulaGLMDelay <- function(timeBool=TRUE,factorsBool=FALSE,delay=FALSE,outbreak=FALSE){
  # Description
  # Args:
  #     populationOffset: ---
  # Returns:
  #     Vector of X
  
  # Smallest formula
  formulaString <- "response ~ 1"
  
  # With time trend?
  if (timeBool){
    formulaString <- paste(formulaString,"+wtime",sep ="")}
  
  
  
  # With factors?
  if(factorsBool){
    formulaString <- paste(formulaString,"+as.factor(seasgroups)",sep ="")}
  #   # With delays?
  if(delay){
    formulaString <- paste(formulaString,"+as.factor(delay)",sep ="")}  
  if(outbreak){
    formulaString <- paste(formulaString,"+f(outbreakOrNot,model='linear', prec.linear = 1)",sep ="")}  
  # Return formula as a string
  return(formulaString) 
}
################################################################################
# END OF FORMULA FUNCTION
################################################################################


######################################################################
# CDF of the negbin mixture with different means and sizes
######################################################################
pmix <- function(y, mu, size) {
  PN <- pnbinom(y, mu=mu, size=size)
  lala <- 1/sum(!is.na(PN))*sum(PN,
                                na.rm=TRUE)
  return(lala)
}
######################################################################
# END OF CDF of the negbin mixture with different means and sizes
######################################################################

######################################################################
# Find the root(s) of a 1D function using the bisection method
#
# Params:
# f - the function to minimize or the first derivate of the function to optim
# reltol - relative tolerance epsilon
######################################################################
bisection <- function(f, bracket) {
  ##Boolean for convergence
  convergence <- FALSE
  ##Loop until converged
  while (!convergence) {
    #Half the interval (problem with ints: what uneven number?)
    x <- ceiling(mean(bracket))
    ##Direct hit? -> stop
    if (isTRUE(all.equal(f(x),0))) break
    ##Choose the interval, containing the root
    bracket <- if (f(bracket[1])*f(x) <= 0) c(bracket[1],x) else c(x,bracket[2])
    ##Have we obtained convergence?
    convergence <- (bracket[1]+1) == bracket[2]
  }
  #Return the value of x^{n+1}
  return(ceiling(mean(bracket)))
}

######################################################################
# END OF BISECTION FUNCTION
######################################################################

######################################################################
##Find the p-quantile of the mixture distribution using bisectioning
##
## Parameters:
## p - the q_p quantile is found
## mu - mean vector
## size - size param
## bracket - vector length two, s.t. qmix(bracket[1] < 1-alpha and
## qmix(bracket[2]) > 1-alpha. Exception: if bracket[1]=0
## then qmix(bracket[1] > 1-alpha is ok.
######################################################################
qmix <- function(p, mu, size, bracket=c(0,mu*100)) {
  
  target <- function(y) {
    pmix(y=y,mu=mu,size=size) - p
  }
  
  if (target(bracket[1]) * target(bracket[2]) > 0) {
    
    if ((bracket[1] == 0) & (target(bracket[1]) > 0)) return(0)
    stop("Not a good bracket.")
  }
  bisection(target, bracket=bracket)
}

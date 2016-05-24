### R code from vignette source 'Rnw/algo_farrington.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: algo_farrington.Rnw:25-35
###################################################
anscombe.residuals <- function(m,phi) {
  y <- m$y
  mu <- fitted.values(m)
  
  #Compute raw Anscombe residuals
  a <- 3/2*(y^(2/3) * mu^(-1/6) - mu^(1/2))
  #Compute standardized residuals
  a <- a/sqrt(phi * (1-hatvalues(m)))
  return(a)
}


################################################################################
# WEIGHTS FUNCTION
################################################################################

algo.farrington.assign.weights <- function(s,weightsThreshold=1) {
    #s_i^(-2) for s_i<weightsThreshold and 1 otherwise
    gamma <- length(s)/(sum(    (s^(-2))^(s>weightsThreshold) ))
    omega <- numeric(length(s))
    omega[s>weightsThreshold] <- gamma*(s[s>weightsThreshold]^(-2))
    omega[s<=weightsThreshold] <- gamma
    return(omega)
}

###################################################
### code chunk number 3: algo_farrington.Rnw:136-305
###################################################
algo.farrington.fitGLM <- function(response,wtime,timeTrend=TRUE,reweight=TRUE,...) {
  #Model formula depends on whether to include a time trend or not.
  theModel <- as.formula(ifelse(timeTrend, "response~1+wtime","response~1"))

  #Fit it -- this is slow. An improvement would be to use glm.fit here.
  model <- glm(theModel, family = quasipoisson(link="log"))
    
 #Check convergence - if no convergence we return empty handed.
  if (!model$converged) {
    #Try without time dependence
    if (timeTrend) {
      cat("Warning: No convergence with timeTrend -- trying without.\n")
      #Set model to one without time trend
      theModel <- as.formula("response~1")
      model <- glm(response ~ 1, family = quasipoisson(link="log"))
    } 

    if (!model$converged) {
      cat("Warning: No convergence in this case.\n")
      print(cbind(response,wtime))
      return(NULL)
    }
  }

  #Overdispersion parameter phi
  phi <- max(summary(model)$dispersion,1)
  
  #In case reweighting using Anscome residuals is requested
  if (reweight) {
    s <- anscombe.residuals(model,phi)
    omega <- algo.farrington.assign.weights(s)
    model <- glm(theModel,family=quasipoisson(link="log"),weights=omega)
    #Here, the overdispersion often becomes small, so we use the max
    #to ensure we don't operate with quantities less than 1.
    phi <- max(summary(model)$dispersion,1)
  } # end of refit.
  

  #Add wtime, response and phi to the model
  model$phi <- phi
  model$wtime <- wtime
  model$response <- response
  #Done
  return(model)
}

######################################################################
# The algo.farrington.fitGLM function in a version using glm.fit 
# which is faster than the call using "glm. 
# This saves lots of overhead and increases speed.
#
# Author: Mikko Virtanen (@thl.fi) with minor modifications by Michael Hoehle
# Date:   9 June 2010 
#
# Note: Not all glm results may work on the output. But for the
# necessary ones for the algo.farrington procedure work.
######################################################################

algo.farrington.fitGLM.fast <- function(response,wtime,timeTrend=TRUE,reweight=TRUE, ...) {
  #Create design matrix and formula needed for the terms object
  #Results depends on whether to include a time trend or not.
  if (timeTrend) {
    design<-cbind(intercept=1,wtime=wtime) 
    Formula<-response~wtime 
  } else {
    design<-matrix(1,nrow=length(wtime),dimnames=list(NULL,c("intercept")))
    Formula<-response~1
  }
  
  #Fit it using glm.fit which is faster than calling "glm"
  model <- glm.fit(design,response, family = quasipoisson(link = "log"))
      
   #Check convergence - if no convergence we return empty handed.
   if (!model$converged) {
      #Try without time dependence
     if (timeTrend) {
       cat("Warning: No convergence with timeTrend -- trying without.\n")
       #Drop time from design matrix
       design <- design[,1,drop=FALSE]
       #Refit
       model <- glm.fit(design,response, family = quasipoisson(link = "log"))
       Formula<-response~1
     } 
     #No convergence and no time trend. That's not good.
   }

   #Fix class of output to glm/lm object in order for anscombe.residuals to work
   #Note though: not all glm methods may work for the result
   class(model) <- c("glm","lm")

   #Overdispersion parameter phi
   phi <- max(summary.glm(model)$dispersion,1)
    
   #In case reweighting using Anscome residuals is requested
   if (reweight) {
     s <- anscombe.residuals(model,phi)
     omega <- algo.farrington.assign.weights(s)
     model <- glm.fit(design,response, family = quasipoisson(link = "log"), weights = omega)
     #Here, the overdispersion often becomes small, so we use the max
     #to ensure we don't operate with quantities less than 1.
     phi <- max(summary.glm(model)$dispersion,1)
   } # end of refit.
    
   model$phi <- phi
   model$wtime <- wtime
   model$response <- response
   model$terms <- terms(Formula)
   # cheating a bit, all methods for glm may not work
   class(model)<-c("algo.farrington.glm","glm","lm") # 23/10/2012 (SM):
                                        # added "lm" class to avoid warnings
                                        # from predict.lm about fake object
   #Done
  return(model)
}

######################################################################
# Experimental function to include a population offset in the 
# farrington procedure based on algo.farrington.fitGLM
# Alternative: include populationOffset argument in the two other
# fit functions, but I suspect use of this is not so common
#
# Parameters:
#  takes an additional "population" parameter
######################################################################

algo.farrington.fitGLM.populationOffset <- function(response,wtime,population,timeTrend=TRUE,reweight=TRUE,...) {
  #Model formula depends on whether to include a time trend or not.
  theModel <- as.formula(ifelse(timeTrend, "response~offset(log(population)) + 1 + wtime","response~offset(log(population)) + 1"))

  #Fit it -- this is slow. An improvement would be to use glm.fit here.
  model <- glm(theModel, family = quasipoisson(link="log"))
    
 #Check convergence - if no convergence we return empty handed.
  if (!model$converged) {
    #Try without time dependence
    if (timeTrend) {
     model <- glm(response ~ 1, family = quasipoisson(link="log"))
     cat("Warning: No convergence with timeTrend -- trying without.\n")
    } 

    if (!model$converged) {
      cat("Warning: No convergence in this case.\n")
      print(cbind(response,wtime))
      return(NULL)
    }
  }

  #Overdispersion parameter phi
  phi <- max(summary(model)$dispersion,1)
  
  #In case reweighting using Anscome residuals is requested
  if (reweight) {
    s <- anscombe.residuals(model,phi)
    omega <- algo.farrington.assign.weights(s)
    model <- glm(theModel,family=quasipoisson(link="log"),weights=omega)
    #Here, the overdispersion often becomes small, so we use the max
    #to ensure we don't operate with quantities less than 1.
    phi <- max(summary(model)$dispersion,1)
  } # end of refit.
  

  #Add wtime, response and phi to the model
  model$phi <- phi
  model$wtime <- wtime
  model$response <- response
  model$population <- population
  #Done
  return(model)
}




###################################################
### code chunk number 4: algo_farrington.Rnw:344-370
###################################################

algo.farrington.threshold <- function(pred,phi,alpha=0.01,skewness.transform="none",y) {
  #Fetch mu0 and var(mu0) from the prediction object
  mu0 <- pred$fit
  tau <- phi + (pred$se.fit^2)/mu0
  #Standard deviation of prediction, i.e. sqrt(var(h(Y_0)-h(\mu_0))) 
  switch(skewness.transform,
         "none" = { se <- sqrt(mu0*tau); exponent <- 1},
         "1/2" = { se <- sqrt(1/4*tau); exponent <- 1/2},
         "2/3"  = { se <- sqrt(4/9*mu0^(1/3)*tau); exponent <- 2/3},
         { stop("No proper exponent in algo.farrington.threshold.")})

  #Note that lu can contain NA's if e.g. (-1.47)^(3/2)
  lu <- sort((mu0^exponent + c(-1,1)*qnorm(1-alpha/2)*se)^(1/exponent),na.last=FALSE)

  #Ensure that lower bound is non-negative
  lu[1] <- max(0,lu[1],na.rm=TRUE)

  #Compute quantiles of the predictive distribution based on the
  #normal approximation on the transformed scale
  q <- pnorm( y^(1/exponent) , mean=mu0^exponent, sd=se)
  m <- qnorm(0.5, mean=mu0^exponent, sd=se)^(1/exponent)
  
  #Return lower and upper bounds
  return(c(lu,q=q,m=m))
}


###################################################
### code chunk number 5: algo_farrington.Rnw:412-451
###################################################
######################################################################
# Compute indices of reference value using Date class
#
# Params:
#  t0 - Date object describing the time point
#  b  - Number of years to go back in time
#  w  - Half width of window to include reference values for
#  epochStr - "1 month", "1 week" or "1 day"
#  epochs - Vector containing the epoch value of the sts/disProg object
#
# Details:
#  Using the Date class the reference values are formed as follows:
#   Starting from d0 go i, i in 1,...,b years back in time.
#   
# Returns:
#  a vector of indices in epochs which match
######################################################################
refvalIdxByDate <- function(t0, b, w, epochStr, epochs) {
  refDays <- NULL
  refPoints <- seq( t0, length=b+1, by="-1 year")[-1]

  #Loop over all b-lagged points and append appropriate w-lagged points
  for (j in 1:length(refPoints)) {
      refPointWindow <- c(rev(seq(refPoints[j], length=w+1, by=paste("-",epochStr,sep=""))),
                          seq(refPoints[j], length=w+1, by=epochStr)[-1])
      refDays <- append(refDays,refPointWindow)
  }
  if (epochStr == "1 week") {
      #What weekday is t0 (0=Sunday, 1=Monday, ...)
      epochWeekDay <- as.numeric(format(t0,"%w"))
      #How many days to go forward to obtain the next "epochWeekDay", i.e. (d0 - d) mod 7
      dx.forward <- (epochWeekDay - as.numeric(format(refDays,"%w"))) %% 7
      #How many days to go backward to obtain the next "epochWeekDay", i.e. (d - d0) mod 7
      dx.backward <- (as.numeric(format(refDays,"%w")) - epochWeekDay) %% 7
      #What is shorter - go forward or go backward? 
      #By convention: always go to the closest weekday as t0
      refDays <- refDays + ifelse(dx.forward < dx.backward, dx.forward, -dx.backward)
  }
  if (epochStr == "1 month") {
      #What day of the month is t0 (it is assumed that all epochs have the same value here)
      epochDay <- as.numeric(format(t0,"%d"))
      #By convention: go back in time to closest 1st of month
      refDays <- refDays -  (as.numeric(format(refDays, "%d")) - epochDay)
  }
  #Find the index of these reference values
  wtime <- match(as.numeric(refDays), epochs)
  return(wtime)
} 



###################################################
### code chunk number 6: algo_farrington.Rnw:571-769
###################################################

algo.farrington <- function(disProgObj, control=list(range=NULL, b=3, w=3, reweight=TRUE, verbose=FALSE,alpha=0.01,trend=TRUE,limit54=c(5,4),powertrans="2/3",fitFun=c("algo.farrington.fitGLM.fast","algo.farrington.fitGLM","algo.farrington.fitGLM.populationOffset"))) { 
  #Fetch observed
  observed <- disProgObj$observed
  freq <- disProgObj$freq
  epochStr <- switch( as.character(freq), "12" = "1 month","52" =  "1 week","365" = "1 day")
  #Fetch population (if it exists)
  if (!is.null(disProgObj$populationFrac)) {
    population <- disProgObj$populationFrac } 
  else { 
    population <- rep(1,length(observed))
  }

  ######################################################################
  # Fix missing control options
  ######################################################################
  if (is.null(control$range)) {
    control$range <- (freq*control$b - control$w):length(observed)
  }
  if (is.null(control$b))        {control$b=5}
  if (is.null(control$w))        {control$w=3}
  if (is.null(control$reweight)) {control$reweight=TRUE}
  if (is.null(control$verbose))  {control$verbose=FALSE}
  if (is.null(control$alpha))    {control$alpha=0.05}
  if (is.null(control$trend))    {control$trend=TRUE}
  if (is.null(control$plot))     {control$plot=FALSE}
  if (is.null(control$limit54))  {control$limit54=c(5,4)}
  if (is.null(control$powertrans)){control$powertrans="2/3"}
  if (is.null(control$fitFun))   {
    control$fitFun="algo.farrington.fitGLM.fast"
  } else {
    control$fitFun <- match.arg(control$fitFun, c("algo.farrington.fitGLM.fast","algo.farrington.fitGLM","algo.farrington.fitGLM.populationOffset"))
  }

  #Use special Date class mechanism to find reference months/weeks/days
  if (is.null(disProgObj[["epochAsDate",exact=TRUE]])) { 
    epochAsDate <- FALSE 
  } else { 
    epochAsDate <- disProgObj[["epochAsDate",exact=TRUE]] 
  }
    
  #check options
  if (!((control$limit54[1] >= 0) &  (control$limit54[2] > 0))) {
    stop("The limit54 arguments are out of bounds: cases >= 0 and period > 0.")
  }
  #Check control$range is within bounds.
  if (any((control$range < 1) | (control$range > length(disProgObj$observed)))) {
      stop("Range values are out of bounds (has to be within 1..",length(disProgObj$observed)," for the present data).")
  }
  
  # initialize the necessary vectors
  alarm <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  trend <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  upperbound <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  # predictive distribution
  pd <- matrix(data = 0, nrow = length(control$range), ncol = 2)

  # Define objects
  n <- control$b*(2*control$w+1)
  
  # 2: Fit of the initial model and first estimation of mean and dispersion
  #    parameter
  for (k in control$range) {
    # transform the observed vector in the way
    # that the timepoint to be evaluated is at last position
    #shortObserved <- observed[1:(maxRange - k + 1)]

    if (control$verbose) { cat("k=",k,"\n")}
            
    #Find index of all epochs, which are to be used as reference values
    #i.e. with index k-w,..,k+w 
    #in the years (current year)-1,...,(current year)-b
    if (!epochAsDate) {
      wtimeAll <- NULL
      for (i in control$b:1){
        wtimeAll <- append(wtimeAll,seq(k-freq*i-control$w,k-freq*i+control$w,by=1))
      }
      #Select them as reference values - but only those who exist
      wtime <- wtimeAll[wtimeAll>0]
      if (length(wtimeAll) != length(wtime)) {
          warning("@ range= ",k,": With current b and w then ",length(wtimeAll) - length(wtime),"/",length(wtimeAll), " reference values did not exist (index<1).")
      }
    } else { #Alternative approach using Dates
      t0 <- as.Date(disProgObj$week[k], origin="1970-01-01")
      wtimeAll <- refvalIdxByDate( t0=t0, b=control$b, w=control$w, epochStr=epochStr, epochs=disProgObj$week)
      #Select them as reference values (but only those not being NA!)
      wtime <- wtimeAll[!is.na(wtimeAll)]
      #Throw warning if necessary
      if (length(wtimeAll) != length(wtime)) {
          warning("@ range= ",k,": With current b and w then ",length(wtimeAll) - length(wtime),"/",length(wtimeAll), " reference values did not exist (index<1).")
      }
    }
    
    #Extract values from indices
    response <- observed[wtime]
    pop <- population[wtime]

    if (control$verbose) { print(response)}

    ######################################################################
    #Fit the model with overdispersion -- the initial fit
    ######################################################################
    #New feature: fitFun can now be the fast function for fitting the GLM
    model <- do.call(control$fitFun, args=list(response=response,wtime=wtime,population=pop,timeTrend=control$trend,reweight=control$reweight))

    #Stupid check to pass on NULL values from the algo.farrington.fitGLM proc.
    if (is.null(model)) return(model)

    ######################################################################
    #Time trend
    #
    #Check whether to include time trend, to do this we need to check whether
    #1) wtime is signifcant at the 95lvl
    #2) the predicted value is not larger than any observed value
    #3) the historical data span at least 3 years.
    doTrend <- control$trend
#Bug discovered by Julia Kammerer and Sabrina Heckl: Only investigate trend if it actually was part of the GLM
    #if (control$trend) {
    if ("wtime" %in% names(coef(model))){
      #is the p-value for the trend significant (0.05) level
      p <- summary.glm(model)$coefficients["wtime",4]
      significant <- (p < 0.05)
      #prediction for time k
      mu0Hat <- predict.glm(model,data.frame(wtime=c(k),population=population[k]),type="response")
      #have to use at least three years of data to allow for a trend
      atLeastThreeYears <- (control$b>=3)
      #no horrible predictions
      noExtrapolation <- mu0Hat <= max(response)
     
      #All 3 criteria have to be met in order to include the trend. Otherwise
      #it is removed. Only necessary to check this if a trend is requested.
      if (!(atLeastThreeYears && significant && noExtrapolation)) {
        doTrend <- FALSE
        model <- do.call(control$fitFun, args=list(response=response,wtime=wtime,population=pop,timeTrend=FALSE,reweight=control$reweight))
      }
    } else {
      doTrend <- FALSE
    }
    #done with time trend
    ######################################################################
    
    ######################################################################
    # Calculate prediction & confidence interval                         #
    ######################################################################
    #Predict value - note that the se is the mean CI
    #and not the prediction error of a single observation
    pred <- predict.glm(model,data.frame(wtime=c(k),population=population[k]),dispersion=model$phi,
                        type="response",se.fit=TRUE)
    #Calculate lower and upper threshold
    lu <- algo.farrington.threshold(pred,model$phi,skewness.transform=control$powertrans,alpha=control$alpha, observed[k])

    ######################################################################
    # If requested show a plot of the fit.
    ######################################################################
    if (control$plot) {
      #Compute all predictions
      data <- data.frame(wtime=seq(min(wtime),k,length=1000))
      preds <- predict(model,data,type="response",dispersion=model$phi)

      #Show a plot of the model fit.
      plot(c(wtime, k), c(response,observed[k]),ylim=range(c(observed[data$wtime],lu)),,xlab="time",ylab="No. infected",main=paste("Prediction at time t=",k," with b=",control$b,",w=",control$w,sep=""),pch=c(rep(1,length(wtime)),16))
      #Add the prediction
      lines(data$wtime,preds,col=1,pch=2)

      #Add the thresholds to the plot
      lines(rep(k,2),lu[1:2],col=3,lty=2)
    }


    ######################################################################
    #Postprocessing steps
    ######################################################################

    #Compute exceedance score unless less than 5 reports during last 4 weeks.
    #Changed in version 0.9-7 - current week is included now
    enoughCases <- (sum(observed[(k-control$limit54[2]+1):k])>=control$limit54[1])

    #18 May 2006: Bug/unexpected feature found by Y. Le Strat. 
    #the okHistory variable meant to protect against zero count problems,
    #but instead it resulted in exceedance score == 0 for low counts. 
    #Now removed to be concordant with the Farrington 1996 paper.
    X <- ifelse(enoughCases,(observed[k] - pred$fit) / (lu[2] - pred$fit),0)

    #Do we have an alarm -- i.e. is observation beyond CI??
    #upperbound only relevant if we can have an alarm (enoughCases)
    trend[k-min(control$range)+1] <- doTrend
    alarm[k-min(control$range)+1] <- (X>1)
    upperbound[k-min(control$range)+1] <- ifelse(enoughCases,lu[2],0)
    #Compute bounds of the predictive
    pd[k-min(control$range)+1,] <- lu[c(3,4)]

  }#done looping over all time points

  #Add name and data name to control object.
  control$name <- paste("farrington(",control$w,",",0,",",control$b,")",sep="")
  control$data <- paste(deparse(substitute(disProgObj)))
  #Add information about predictive distribution
  control$pd   <- pd

  # return alarm and upperbound vectors 
  result <- list(alarm = alarm, upperbound = upperbound, trend=trend, 
                 disProgObj=disProgObj, control=control) 
  class(result) <- "survRes" 

  #Done
  return(result) 
}




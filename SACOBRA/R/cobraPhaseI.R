#
#
#Samineh Bagheri, Patrick Koch, Wolfgang Konen
#Cologne Univeristy of Applied Science
#
#27.April.2014
#cobraPhaseI.R
#phase I
#Find a feasible point
#
#
#

#'  Find a feasible solution.
#'
#'  Find a feasible solution using the COBRA optimizer phase I by searching new infill points 
#'  with the help of RBF surrogate models
#'
#'  @param cobra an object of class COBRA, this is a (long) list containing all settings
#'      from \code{\link{cobraInit}}
#'
#'  @return \code{cobra}, an object of class COBRA
#'  
#' @seealso   \code{\link{cobraPhaseII}}, \code{\link{cobraInit}}
#' @author Wolfgang Konen, Samineh Bagheri, Patrick Koch, Cologne Univeristy of Applied Sciences
#' @export
#' 
cobraPhaseI <- function(cobra){
  
  print("no Feasible point available")
  print("PHASE I Started")
  phase<-"PHASE I"
  

  ##################################Initializing cobra parameters####################################
  #xbestIndex<-which(cobra$Fres==cobra$fbest)    # /WK/ we have cobra$ibest
  fn=cobra$fn
  dimension=cobra$dimension
  nConstraints <- cobra$nConstraints
  EPS <- rep(cobra$epsilonInit,cobra$nConstraints) # Initializing Margin for all constraints
  n<-nrow(cobra$A)
  iteration<-cobra$iteration
  feasibleSolutionExists<-(0 %in% cobra$numViol)
  #predY initialization should move to init code for phase I and phase II
  predY = rep(NA,cobra$initDesPoints) # structure to store surrogate optimization results
  cobra$predC = matrix(nrow=cobra$initDesPoints,ncol=cobra$nConstraints) # matrix to store predicted constraint values
  constraintPrediction = NULL # actual constraint value prediction
  optimizerConvergence = rep(1,cobra$initDesPoints) # vector to store optimizer convergence
  cobra$optimizationTime <- rep(0,cobra$initDesPoints)
  feas <- sapply(1:nrow(cobra$Gres), FUN = function(i) !any(cobra$Gres[i,]>0)) # feasibility of initial design
  #feval initialization should move to init code for phase I and phase II
  feval <- rep(NA,cobra$initDesPoints) # structure to store function evaluations on surrogate
  fbestArray <- rep(cobra$fbest,cobra$initDesPoints)
  penaF <- cobra$penaF    # /WK/
  sigmaD <- cobra$sigmaD;  # /WK/
  cobra$Ffeas <- cobra$fbest      # /WK/ needed?
  cobra$Fall <- min(cobra$Fres)   # /WK/ needed?
  
  
  ########################################################################################################
  # STEP4:                                                                                               #
  #  While a feasible Point is not found do the loop                                                     #
  ########################################################################################################
  
  while(!feasibleSolutionExists && n<cobra$feval){   # /WK/ bug fix
    gc()
    iteration<-iteration+1
    
    ################################################
    # STEP4.1: Update RBF for m constraints       #
    ################################################
    constraintSurrogates <- list()  
    cobra$A <- as.matrix(cobra$A)
    cobra$Gres <- as.matrix(cobra$Gres)
    constraintSurrogates <- trainCubicRBF(cobra$A,cobra$Gres,squares=cobra$squaresC)
    fitnessSurrogate <- trainCubicRBF(cobra$A,cobra$Fres,squares=cobra$squaresF)
    
    ################################################
    # STEP4.2:Determine distance requirement(XI) #
    ################################################
    
    gama<-cobra$XI[(nrow(cobra$A) %% length(cobra$XI))+1]     # /WK/ bug fix
    ro<-gama*cobra$l
    
    ################################################
    # STEP4.3: select next point                   #
    ################################################
    #optimization sub problem
    
    # surrogate penalty function for unconstraint optimization methods
    subProb<-function(x){ 
      distance <- distLine(x,cobra$A)
      subC<-pmax((ro-distance),0)
      penalty1<-sum(subC)
      #cat(">>> Entering interpRBF\n")
      constraintPrediction <<- interpRBF(x,constraintSurrogates)+EPS^2       # why to the power of two? #
      maxViolation <- sapply(1:length(constraintPrediction), FUN=function(i)max(0,constraintPrediction[i]))
      #penalty2 <- sum(maxViolation)
      # Constraint handling: Coit et al. (1996)
      NFT <- rep(0.1,length(maxViolation))
      kappa <- 2
      
      penalty2 <- (length(maxViolation)*sum(maxViolation)*n)*sum(sapply(1:length(maxViolation),FUN=function(i)(maxViolation[i]/NFT[i])^kappa))
      # y<-interpRBF(x, fitnessSurrogate) + (penalty1*sigmaD[1] + penalty2)*penaF[1]
      #y<- sum(max(interpRBF(x,constraintSurrogates),0)^2) + (penalty1*sigmaD[1] + penalty2)*penaF[1]
      f<-interpRBF(x, fitnessSurrogate)
      y <- f  + penalty2 #+ penalty1*sigmaD[1]
      #cat("<<< leaving interpRBF\n")
      return(y)
    }
    
    # surrogate evaluation of 'f' for constraint optimization methods
    subProb2 <- function(x){
      #y<-predict.RBFinter(constraintSurrogates,matrix(x,ncol=dimension))
      y<- sum(max(interpRBF(x,constraintSurrogates),0)^2)
      return(y)
    }
    
    
    
    # surrogate evaluation of '\vec{g}' for constraint optimization methods
    gCOBRA <- function(x) {
      h <- c()
      distance <- distLine(x,cobra$A)
      subC<-pmax((ro-distance),0)
      h[1] <- sum(subC)
      #cat(">>> Entering interpRBF\n")
      constraintPrediction <<- interpRBF(x,constraintSurrogates)+EPS^2
      h <- (-1.0)*c(h[1], constraintPrediction) # TODO -1* ... is required for COBYLA constraints, maybe also for other optimizers? 
      #cat("<<< leaving interpRBF\n")
      return(h)
    }
    
    
    ptm <- proc.time()
    subMin <- list()
    cat(cobra$seqOptimizer, " optimization on surrogate...\n")
    
    switch(cobra$seqOptimizer,
           COBYLA={ subMin<-nloptr::cobyla(cobra$xbest,fn=subProb2,lower=cobra$lower,upper=cobra$upper,hin=gCOBRA,control=list(maxeval=cobra$seqFeval)); subMin$feval=subMin$iter },
           ISRES ={subMin<- isres2(cobra$xbest,fn=subProb2,lower=cobra$lower,upper=cobra$upper,hin=gCOBRA, maxeval=cobra$seqFeval); subMin$feval=subMin$iter},
           HJKB = { subMin<-dfoptim::hjkb(cobra$xbest,fn=subProb,lower=cobra$lower,upper=cobra$upper,control=list(maxfeval=cobra$seqFeval)) },
           NMKB = { subMin<-nmkb2(cobra$xbest,fn=subProb,lower=cobra$lower,upper=cobra$upper, control=list(maxfeval=cobra$seqFeval,tol=cobra$seqTol)) },
           #ACTIVECMA  = { subMin <- ActiveOnePlusOneCMAES(subProb, length(cobra$xbest), cobra$xbest, opts=list(esname="ActiveCMAES", lb=cobra$lower, ub=cobra$upper, maxFunEvals=cobra$seqFeval)); 
           #               subMin$convergence <- 1;
           #             }
    )
    cobra$optimizationTime <- c(cobra$optimizationTime, (proc.time()-ptm)[3])
    cat("Optimization time for", subMin$feval, " iterations:", cobra$optimizationTime[length(cobra$optimizationTime)], "seconds\n")  
    cat("Predicted infill value:",subMin$value,"\n")
    if (penaF[1]*penaF[2] < penaF[3])
      penaF[1] <- penaF[1]*penaF[2]
    
    CHECKDIST=T               
    if (CHECKDIST) {
      # check whether the solution returned fulfills the distance requirement
      # 
      ed = ro - distLine(subMin$par,constraintSurrogates$xp)
      violatedDist = which(ed>0)
      if (length(violatedDist)>0) {
        #
        # If distance requirement is not fulfilled, increase sigmaD[1] (the empirical penalty factor 
        # in fitFuncPenalRBF or subProb). This influences currently only NMKB (not COBYLA),  
        # because sigmaD is only used in fitFuncPenalRBF or subProb.
        #
        if(sigmaD[1]*sigmaD[2]<sigmaD[3]) {
          sigmaD[1] <- sigmaD[1]*sigmaD[2]
          cat("***   Increasing sigmaD to: ",sigmaD[1],"at iteration",n,"  ***\n")          
        }
      }
    }

    
    ################################################
    # STEP4.4: Evaluate real functions             #
    ################################################
    cobra$fe<-cobra$fe+1
    predY <- c(predY,subMin$value)
    feval <- c(feval,subMin$feval)
    optimizerConvergence <- c(optimizerConvergence,subMin$convergence)
    cobra$predC <- rbind(cobra$predC,constraintPrediction)
    xNew<-subMin$par
    #cat("X:",xNew,"\n")
    xNew <- unlist(sapply(1:cobra$dimension, FUN=function(i)max(cobra$lower[i],xNew[i])))
    xNew <- unlist(sapply(1:cobra$dimension, FUN=function(i)min(cobra$upper[i],xNew[i])))
    
    ##rescaling
    if(cobra$rescale){      # /WK/ bug fix: added missing cobra$rescale clause
      xNewTemp<-sapply(1:dimension , function(i){  
        y<-scales::rescale(xNew[i],to=c(cobra$originalLower[i],cobra$originalUpper[i]),from=c(cobra$lower[i],cobra$upper[i]))
        return(y)
      })
    } else {
      xNewTemp <- xNew
    }
    xNewEval<-fn(xNewTemp)
    #cat("F/Max(C):",xNewEval[1],"/",max(xNewEval[2:length(xNewEval)]),"\n")
    newNumViol<-length(which((unlist(xNewEval)[-1])>0)) # Calculating number of constraint Violations for new point
    
    
    
    if( newNumViol < 1 ) feas = c(feas, TRUE)
    else feas = c(feas, FALSE)
    newMaxViol<-max(0,max((unlist(xNewEval)[-1])) )            # Calculating maximum violation
    
    
    
    ################################################
    # STEP4.5: Update Information                  #
    ################################################
    cobra$A<-rbind(cobra$A,xNew)
    cobra$Fres <- c(cobra$Fres, xNewEval[1])
    newConstraintLine <- xNewEval[2:(cobra$nConstraints+1)]#sapply(2:(m+1), FUN=function(i)(c(Gres[i],xNewEval[i])))
    cobra$Gres = rbind(cobra$Gres,newConstraintLine)
    cobra$numViol<-c(cobra$numViol,newNumViol)
    cobra$maxViol<-c(cobra$maxViol,newMaxViol)
    cobra$phase<-c(cobra$phase,phase)
    xNewIndex<-length(cobra$numViol)
    s<-sprintf("%s.[%d]: %f %f | %f | %f" , phase ,nrow(cobra$A), cobra$A[xNewIndex,1] ,cobra$A[xNewIndex,2] , xNewEval[1] , newMaxViol)
    print(s)
    
    
    ################################################
    # STEP4.6: Update best Point                   #
    ################################################
    if( (cobra$numViol[xNewIndex] < cobra$numViol[cobra$ibest]) || 
          ((cobra$numViol[xNewIndex]==cobra$numViol[cobra$ibest]) 
           && (cobra$maxViol[xNewIndex]< cobra$maxViol[cobra$ibest])) ){
      
      cobra$xbest<-xNew
      cobra$fbest<-cobra$Fres[xNewIndex]
      cobra$ibest<-xNewIndex
    }
    
    
    
    cobra$fbestArray<-c(cobra$fbestArray,cobra$fbest)
    cobra$xbestArray<-rbind(cobra$xbestArray,cobra$xbest)
    feasibleIndices <- which(sapply(1:nrow(cobra$Gres),FUN=function(i)(all(cobra$Gres[i,]<0))))
    # xbestIndex<-which.min(cobra$Fres[feasibleIndices])                      # finding index of the best point so far
    
    n<-nrow(cobra$A)
    
    
    # result data frame
    df <- data.frame(cobra$Fres, 
                     predY, 
                     feas, 
                     feasPred=rep(NA,length(feas)),   # /WK/ bug fix
                     cobra$numViol,
                     cobra$maxViol,
                     feval, 
                     cobra$fbestArray,
                     #cobra$A,  
                     #cobra$Gres, 
                     #predC, 
                     cobra$seqOptimizer,
                     cobra$optimizationTime,
                     optimizerConvergence, 
                     row.names=NULL)
    colnames(df) = c("y", 
                     "predY", 
                     "feasible",
                     "feasPred",                # /WK/ bug fix
                     "nViolations",
                     "maxViolation",
                     "FEval",
                     "Best",
                     #sprintf("x%i",1:cobra$dimension),
                     #sprintf("c%i",1:cobra$nConstraints), 
                     #sprintf("predC%i",1:cobra$nConstraints),
                     "optimizer",
                     "optimizationTime",
                     "conv")
    df <- cbind(iter=1:nrow(df),df)
    df <- cbind(df,seed=cobra$cobraSeed)
    cobra$df <- df
    cobra$phase1DesignPoints <- nrow(df)
    
    if (cobra$saveSurrogates) {
      cobra$constraintSurrogates = constraintSurrogates;
      cobra$fitnessSurrogate = fitnessSurrogate;
    }
    
    if (cobra$saveIntermediate) {
      # save intermediate results
      # cobraResult = list(cobra=cobra, df=df, constraintSurrogates=constraintSurrogates, fn=fn) 
      cobraResult = cobra
      if (is.na(file.info("results")$isdir)) dir.create("results")    # if directory "results" does not exist, create it
      save(cobraResult, file=sprintf("results/cobra-%s-%s-%i.RData",cobra$fName,cobra$seqOptimizer,cobra$cobraSeed))
    }
    

    
    feasibleSolutionExists<-(0 %in% cobra$numViol)
  }
  
  print("END of PHASE I")
  
  return(cobra)
  
  
}
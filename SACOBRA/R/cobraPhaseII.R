#
#Samineh Bagheri, Patrick Koch, Wolfgang Konen
#Cologne Univeristy of Applied Sciences
#
#April,2014 - May,2015
#cobraPhaseII.R
#
#'  Improve the feasible solution by searching new infill points
#'
#'  Improve the feasible solution using the COBRA optimizer phase II
#'  by searching new infill points with the help of RBF surrogate models. 
#'  May be even called if no feasible solution is found yet, then phase II will try to find
#'  feasible solutions.
#'
#'  @param cobra an object of class COBRA, this is a (long) list containing all settings
#'      from \code{\link{cobraInit}}
#'
#'  @return \code{cobra}, an object of class COBRA from \code{\link{cobraInit}}, 
#'    enhanced here by the following elements (among others):
#'      \item{\code{fn}}{ function returning an (m+1)-vector \code{c(objective,g1,...,gm)}. This
#'            function may be a rescaled and plog-transformed version of the original \code{fn} 
#'            passed into \code{\link{cobraInit}}. The original \code{fn} is in 
#'            \code{cobra$originalFn}. }
#'      \item{\code{df}}{  data frame with summary of the optimization run (see below)}
#'      \item{\code{df2}}{  data frame with additional summary information (see below)}
#'      \item{\code{A}}{ (feval x dim)-matrix containing all evaluated points 
#'            in input space. If rescale==TRUE, all points are in \strong{rescaled} input space. }
#'      \item{\code{Fres}}{ a vector of the objective values of all evaluated points }
#'      \item{\code{Gres}}{ a matrix of the constraint values of all evaluated points }
#'      \item{\code{fbest}}{ the best feasible objective value found }
#'      \item{\code{xbest}}{ the point in input space yielding the best feasible objective value }
#'      \item{\code{ibest}}{ the corresponding iteration number (row of cobra$df, of cobra$A)}
#'      \item{\code{PLOG}}{ If TRUE, then the objective surrogate model is trained on the 
#'            \code{\link{plog}}-transformed objective function. }
#'      
#'   Note that \code{cobra$Fres}, \code{cobra$fbest}, \code{cobra$fbestArray} and similar contain 
#'   always the objective values of the orignial function \code{cobra$fn[1]}. (The surrogate models 
#'   may be trained on a \code{\link{plog}}-transformed version of this function.)
#'   
#'   The data frame \code{cobra$df} contains one row per iteration with columns 
#'   \describe{
#'      \item{iter}{  }
#'      \item{y}{   true objective value \code{Fres} }
#'      \item{predY}{  surrogate objective value. Note: The surrogate may be trained on  
#'            plog-transformed training data, but \code{predY} is transformed back to the original 
#'            objective range. NA for the initial design points.}
#'      \item{predSolu}{  surrogate objective value at best-known solution \code{cobra$solu}, if given. 
#'            If \code{cobra$solu} is NULL, take the current point instead. Note: The surrogate may be trained on  
#'            plog-transformed training data, but \code{predSolu} is transformed back to the original 
#'            objective range. NA for the initial design points.}
#'      \item{feasible}{  }
#'      \item{feasPred}{  }
#'      \item{nViolations}{  }
#'      \item{maxViolation}{  }
#'      \item{FEval}{  number of function evaluations in sequential optimizer. NA if it was a repair step }
#'      \item{Best}{  ever-best feasible objective value \code{fbest}. As long as there is 
#'            no feasible point, take among those with minimum number of violated constraints the
#'            one with minimum Fres. }
#'      \item{optimizer}{ e.g. "COBYLA"  }
#'      \item{optimizationTime}{  in sec}
#'      \item{conv}{  }
#'      \item{seed}{  }
#'   }
#'
#'   The data frame \code{cobra$df2} contains one row per phase-II-iteration with columns 
#'   \describe{
#'      \item{iter}{  }
#'      \item{predY}{  surrogate objective value. Note: The surrogate may be trained on  
#'            plog-transformed training data, but \code{predY} is transformed back to the original 
#'            objective range. NA for the initial design points.}
#'      \item{predVal}{   surrogate objective value + penalty }
#'      \item{predSolu}{   surrogate objective value at true solution (see \code{cobra$df$predSolu}) }
#'      \item{predSoluPenal}{   surrogate objective value + penalty at true solution (only diagnostics)}
#'      \item{sigmaD}{  }
#'      \item{penaF}{  }
#'      \item{XI}{  the DRC element used in the current iteration }
#'      \item{EPS}{  }
#'   }
#'
#' @seealso   \code{\link{cobraPhaseI}}, \code{\link{cobraInit}}
#' @author Wolfgang Konen, Samineh Bagheri, Patrick Koch, Cologne Univeristy of Applied Sciences
#' @export

##########################################################################################################
# Some rules about the COBRA-II-code:
#
# - cobra$df contains one row for each iteration, including initial points
#   cobra$df2 contains one row for each phase-II iteration only
# - cobra$PLOG is set by adFit (SACOBRA.R) and adFit is called at the start of each iteration
#   (see trainSurrogates here in phase II)
# - cobra$solu is always in original input space. But when it is used (for diagnostics) in 
#   updateSaveCobra, then a local copy \code{solu} is made, and - if cobra$rescale==TRUE - 
#   \code{solu} is transformed to the rescaled space.
# - the rescaled space has the bounds [rep(cobra$newlower,d),rep(cobra$newupper,d)], (usually -1 and 1) 
# - cobra$xbest,cobra$fbest,cobra$ibest refer always to the same infill point (same iteration).
# - What is cobra$fbest before a feasible infill is found? - See \code{updateSaveCobra}:
#   If the new infill point is feasible, take its fn[1]-value (of course!). If it is not feasible,
#   leave it at the setting from cobraInitial.R#400: from all points of the initial design
#   with minimum number of violated constraints, take the one with smallest Fres.
# 
##########################################################################################################

cobraPhaseII <- function(cobra){
  gc()
  
  verboseprint(cobra$verbose, important=FALSE,"There is at least one feasible point in the population Or PhaseI is skipped")
  verboseprint(2, important=TRUE,"PHASE II Started")
  phase<-"PHASE II"
  
  testit::assert("cobraPhaseII: cobra$fbest is NULL!",!is.null(cobra$fbest))
  testit::assert("cobraPhaseII: cobra$ibest is NULL!",!is.null(cobra$ibest))
  
  ########################################################################################################
  # STEP5: 
  # Initializing the parameters and
  #  Initialize the margin  and counters                                                                 #
  ########################################################################################################
  #xbest = cobra$xbest
  #fbest = cobra$fbest
  #xbestArray <- cobra$xbestArray
  #fbestArray <- cobra$fbestArray
  newErr1<-0
  newErr2<-0
  err1<-c()
  err2<-c()
  fn=cobra$fn
  dimension=cobra$dimension
  #Gres=cobra$Gres
  #Fres=cobra$Fres
  nConstraints <- cobra$nConstraints
  #numViol=cobra$numViol
  #maxViol=cobra$maxViol
  #A <- cobra$A
  CHECKDIST=T       
  Tfeas <- cobra$Tfeas
  Tinfeas <- cobra$Tinfeas
  Cfeas<-0                # Starting Counters
  Cinfeas<-0
  EPS <- rep(cobra$epsilonInit,cobra$nConstraints) # Initializing margin for all constraints
  #eps <- cobra$epsilonInit        # Initializing margin for fitness function
  n <- nrow(cobra$A)
  nRepair<- 0
  if(n==cobra$initDesPoints){
    predY = rep(NA,cobra$initDesPoints) # structure to store surrogate optimization results
    predVal = rep(NA,cobra$initDesPoints) 
    if(cobra$nConstraints!=0){
      cobra$predC = matrix(nrow=cobra$initDesPoints,ncol=cobra$nConstraints) # matrix to store predicted constraint values
      colnames(cobra$predC) <- paste("C",1:cobra$nConstraints,sep="") 
      feas <- sapply(1:nrow(cobra$Gres), FUN = function(i) !any(cobra$Gres[i,]>0)) # feasibility of initial design
      
    }
    
    #browser()
    feasPred <- rep(NA,cobra$initDesPoints)
    optimizerConvergence = rep(1,cobra$initDesPoints) # vector to store optimizer convergence
    cobra$optimizationTime <- rep(0,cobra$initDesPoints)
    feval <- rep(NA,cobra$initDesPoints) # structure to store function evaluations on surrogate
    fbestArray <- rep(cobra$fbest,cobra$initDesPoints)
  }else{
    predY<-cobra$df$predY
    predVal<-cobra$df$predVal
    #cobra$predC<-cobra$predC
    feas<-cobra$df$feasible
    feasPred<-cobra$df$feasPred
    optimizerConvergence <- cobra$df$conv
    cobra$optimizationTime <- cobra$optimizationTime
    feval <- cobra$df$FEval
    fbestArray <- cobra$fbestArray
  }
  constraintSurrogates <- NULL # have them as variables on the global level of cobraPhaseII  
  fitnessSurrogate     <- NULL # such that inner trainSurrogates() can access them with "<<-"
  fitnessSurrogate1    <- NULL # built model according to the original value of the fitness values
  fitnessSurrogate2    <- NULL # built model according to the plog transformed of the fitness values
  
  constraintPrediction = NULL # actual constraint value prediction
  penaF <- cobra$penaF    
  sigmaD <- cobra$sigmaD;  
  cobra$important<-FALSE
  nCobyla <- 0;               # only for ISRESCOBY
  nIsres <- 0;                # only for ISRESCOBY
  

  if (n >= cobra$feval) warning("n is after Phase I equal or larger than cobra$feval")
  
  # -----------------------------------------------------------------------------------------------
  # ---------------- helper functions cobraPhaseII ------------------------------------------------
  # -----------------------------------------------------------------------------------------------
  
  # surrogate penalty function for unconstraint optimization methods
  subProb<-function(x){ 
    #cat(">>> Entering subProb\n")
    distance <- distLine(x,cobra$A)
    subC<-pmax((ro-distance),0)
    penalty1<-sum(subC)
    
    constraintPrediction <<- interpRBF(x,constraintSurrogates)+EPS^2
    if (cobra$trueFuncForSurrogates) constraintPrediction <<- fn(x)[-1]+EPS^2
    maxViolation <- sapply(1:length(constraintPrediction), FUN=function(i)max(0,constraintPrediction[i]))
    penalty2 <- sum(maxViolation)
    
    
    switch(cobra$constraintHandling,
           
           # J.A. Joines and C.R. Houck:  On the use of non-stationary penalty functions to solve 
           # nonlinear constrained optimization problems with GA's. In Proceedings of the First IEEE 
           # Conference on Evolutionary Computation, p. 579-584, (1994)
           JOINESHOUCK={
             C=0.5
             alpha=2
             beta=2
             y<-interpRBF(x, fitnessSurrogate) + (C * cobra$feval)^alpha * sum(maxViolation)^beta 
             #browser()
           },
           
           # A.E. Smith and D.M. Tate: Genetic optimization using a penalty function. In Proceedings of the 
           # Fifth International Conference on Genetic Algorithms p. 499-505, (1993)
           SMITHTATE={
             #TODO
             lambda = 1
             beta1 = 3
             beta2 = 2
             kappa = 1
             d = sum(maxViolation)^kappa
             y <- interpRBF(x, fitnessSurrogate) + penalty
           },
           
           # D.W. Coit, A.E. Smith and D.M. Tate: Adaptive penalty methods for genetic optimization of 
           # constrained combinatorial problems. In INFORMS Journal on Computing 8(2), (1996)
           COIT={
             fFeas <- cobra$fbest
             fAll <- min(cobra$Fres)
             NFT <- rep(0.05, cobra$nConstraints)
             kappa <- 1
             y <- interpRBF(x, fitnessSurrogate) + (fFeas - fAll) * sum((maxViolation/NFT)^kappa)
           },
           
           # T. Baeck and S. Khuri: An evolutionary heuristic for the maximum independent set problem 
           # In Proceedings of the First IEEE Conference on Evolutionary Computation, p. 531-535, (1994)
           BAECKKHURI={
             #TODO
             K=10^9
             p = cobra$nConstraints
             s <- length(cobra$numViol)
             penalty <- K - s*(K/nConstraints)
             y <- interpRBF(x, fitnessSurrogate) + penalty
             #browser()
           },
           
           DEFAULT={
             y<-interpRBF(x, fitnessSurrogate) + (penalty1*sigmaD[1] + penalty2)*penaF[1]
           }
    )
    
    #PKDEBUG=TRUE|FALSE
    #if(PKDEBUG){# Joines and Houck 1994
    #  C=5
    #  alpha=2
    #  beta=2
    #  y<-interpRBF(x, fitnessSurrogate) + (C * cobra$feval)^alpha * sum(maxViolation)^beta 
    #}else{y<-interpRBF(x, fitnessSurrogate) + (penalty1*sigmaD[1] + penalty2)*penaF[1]}
    
    if(any(is.nan(x))){
      warning("subProb: x value is NaN, returning Inf")
      return(Inf)
    }
    #cat(">>SubProb: ", interpRBF(x, fitnessSurrogate), " / ", (C * cobra$feval)^alpha * sum(maxViolation)^beta , " ||", (C * cobra$feval)^alpha, " / ", sum(maxViolation)^beta ,"\n")
    if (cobra$trueFuncForSurrogates) y<-fn(x)[1] + (penalty1*sigmaD[1] + penalty2)*penaF[1]
    #cat("<<< leaving subProb\n")
    #browser()
    return(y)
  }
  
  # surrogate evaluation of 'f' for constraint optimization methods
  subProb2 <- function(x){
    if(any(is.nan(x))){
      warning("subProb2: x value is NaN, returning Inf")
      return(Inf)
    }
    
    if (cobra$trueFuncForSurrogates) {
      y<-fn(x)[1]
    } else {
      y<-predict.RBFinter(fitnessSurrogate,matrix(x,ncol=dimension))
    }
    
    return(y)
  }
  
  fitFuncPenalRBF <- function(x) { 
    distRequirement <- function(x,fitnessSurrogate,ro) {
      ed = ro - distLine(x,fitnessSurrogate$xp)      # distLine: euclidean distances between x and each xp
      violatedDist = which(ed>0)
      return(sum(ed[violatedDist]))    
    }
    
    if(any(is.nan(x))){
      warning("fitFuncPenalRBF: x value is NaN, returning Inf")
      return(Inf)
    }
    
    y = interpRBF(x, fitnessSurrogate)
    if (cobra$trueFuncForSurrogates) y<-fn(x)[1]
    constraintPrediction <<- cstr <- interpRBF(x,constraintSurrogates) +EPS^2
    if (cobra$trueFuncForSurrogates) constraintPrediction <<- cstr <- fn(x)[-1]+EPS^2
    violatedConstraints = which(cstr>0)
    penalty = sum(cstr[violatedConstraints])
    
    DISTREQUIREMENT=T
    if (DISTREQUIREMENT) {
      penalty = penalty + distRequirement(x,fitnessSurrogate,ro)*sigmaD[1]
    }
    
    #print(c(distRequirement(x,reg.model,gammaD)*343*penaF[1],x[10:12]))
    return(y + penalty*penaF[1])
  }
  
  # surrogate evaluation of '\vec{g}' for constraint optimization methods
  gCOBRA <- function(x) {
    h <- c()
    distance <- distLine(x,cobra$A)
    subC<-pmax((ro-distance),0)
    h[1] <- sum(subC)*cobra$drFactor
    #cat(">>> Entering interpRBF\n")
    constraintPrediction <<- interpRBF(x,constraintSurrogates)+EPS^2
    if (cobra$trueFuncForSurrogates) constraintPrediction <<-  fn(x)[-1]+EPS^2
    #/WK/ Bug fix: the above line for cobra$trueFuncForSurrogates was missing before
    
    DBG=FALSE
    if (DBG & h[1]>0) {
      cat("gCOBRA: ",h,max(constraintPrediction),"\n")
      if (h<770) browser()
    }
    
    if(any(is.nan(x))){
      warning("gCOBRA: x value is NaN, returning Inf")
      return(c(Inf,rep(Inf,cobra$nConstraints)))
    }
    
    h <- (-1.0)*c(h[1], constraintPrediction) # TODO -1* ... is required for COBYLA constraints, maybe also for other optimizers? 
    #cat("<<< leaving interpRBF\n")
    return(h)
  }
  
  
  # check whether the solution returned in subMin fulfills the distance requirement.
  # If not, increase sigmaD[1]. Return sigmaD.
  # 
  checkDistanceReq <- function(subMin,constraintSurrogates,sigmaD,CHECKDIST) {
    if (CHECKDIST) {
      ed = ro - distLine(subMin$par,constraintSurrogates$xp)
      violatedDist = which(ed>0)
      if (length(violatedDist)>0) {
        #
        # If distance requirement is not fulfilled, increase sigmaD[1] (the empirical penalty factor 
        # in fitFuncPenalRBF or subProb). This influences currently only NMKB and HJKB (not COBYLA),  
        # because sigmaD[1] is only used in fitFuncPenalRBF or subProb.
        #
        if(sigmaD[1]*sigmaD[2]<sigmaD[3]) {
          sigmaD[1] <- sigmaD[1]*sigmaD[2]
          verboseprint(cobra$verbose, important=FALSE,paste("***   Increasing sigmaD to: ",sigmaD[1],"at iteration",n,"  ***"))          
        }
      }
    }
    return(sigmaD)
  }
   
  # update cobra information (A, Fres, Gres)
  # update counters Cfeas, Cinfeas on the global level of cobraPhaseII
  #
  updateInfoAndCounters <- function(cobra,xNew,xNewEval,newNumViol,newMaxViol,phase)
  {
    cobra$A<-rbind(cobra$A,xNew)
    cobra$Fres <- c(cobra$Fres, xNewEval[1])
    cobra$Gres = rbind(cobra$Gres,xNewEval[-1])
    cobra$numViol<-c(cobra$numViol,newNumViol)
    cobra$maxViol<-c(cobra$maxViol,newMaxViol)
    cobra$phase<-c(cobra$phase,phase)
    if(nrow(cobra$A) %% cobra$verboseIter == 0){#important to print
      #browser()
      cobra$important=TRUE
    }else{
      cobra$important=FALSE
    }
    xNewIndex<-length(cobra$numViol)
    
    
    verboseprint(cobra$verbose, important=FALSE,(sprintf("%s.[%d]: %f %f | %f | %f" 
                  , phase ,nrow(cobra$A), cobra$A[xNewIndex,1] ,cobra$A[xNewIndex,2] , xNewEval[1] , newMaxViol)))
    #SB: I added another print line which prints the best found solution after several iterations, if we do not have any interest for this the following lines can be commented
    realXbest<-sapply(1:length(cobra$xbest) , function(i){scales::rescale(cobra$xbest[i],from=c(cobra$newlower,cobra$newupper),to=c(cobra$lower[i],cobra$upper[i]))})

    verboseprint(cobra$verbose, important=cobra$important,(sprintf("%s.[%d]: %f %f | %f | %f" 
                                                                   , "Best Result" ,nrow(cobra$A), realXbest[1] ,realXbest[2] , cobra$fbest[1] , cobra$maxViol[cobra$ibest])))   
    
    if(cobra$numViol[xNewIndex]==0){ #check if the new point is feasible
      Cfeas <<- Cfeas+1
      Cinfeas <<- 0
    }else{
      Cinfeas <<- Cinfeas+1
      Cfeas <<- 0
    }
    
    return(cobra)
  }
  
  # adjust margins (EPS)
  # may change EPS, Cfeas, and Cinfeas on the global level of cobraPhaseII
  #
  adjustMargins <- function(Cfeas,Tfeas,Cinfeas,Tinfeas,EPS,epsMax) {
    if(Cfeas >= Tfeas){
      EPS <<- EPS/2

      verboseprint(cobra$verbose, important=FALSE,sprintf("reducing epsilon to %f",EPS[1]))
      Cfeas <<- 0
    }
    
    if(Cinfeas>=Tinfeas){
      EPS <<- pmin(2*EPS,epsMax) 
      verboseprint(cobra$verbose, important=FALSE,sprintf("increasing epsilon to %f",EPS[1]))
      Cinfeas <<- 0
    }
  }

  # update the information in list 'cobra', including data frames df and df2, 
  # and - if cobra$saveIntermediate==TRUE - save it in subdir 'results/'
  #
  updateSaveCobra <- function(cobra,xNew,feas,feasPred,feval,optimizerConvergence,
                              predY,predVal,subMin,sigmaD,penaF,gama,EPS)
  {
    if (cobra$DEBUG_XI) {
      df_fxStart <- c(cobra$df$fxStart,cobra$fn(cobra$xStart)[1])
      df_fxbest <-  c(cobra$df$fxbest,cobra$fn(cobra$xbest)[1])
      df_RS  <-  c(cobra$df$RS,!all(cobra$xbest==cobra$xStart))
      df_RS2 <-  c(cobra$df$RS2,cobra$DEBUG_RS)
      if (is.null(cobra$df)) {
        df_XI <- c(rep(NA,cobra$initDesPoints),gama)
      } else {
        df_XI <- c(cobra$df$XI,gama)
      }
    }

    xNewIndex<-length(cobra$numViol)
    
    #testit::assert(cobra$fbest<=min(cobra$fbestArray[which(cobra$numViol==0)]))
    
    ### /WK/ Bug fix: the 3 commented lines below were mind-buggingly complex and not 
    ### correct. In some cases, cobra$fbestArray (a.k.a df$Best) would *increase* in value!!
    ### Changed it to a simpler logic with the help of cobra$ibest, which is - if set - 
    ### the index to the so-far-best feasible solution. If there is no feasible solution yet 
    ### it is NULL.
    ###
    #    index<-which(cobra$Fres==cobra$fbest)
    #    if(cobra$numViol[tail(index,1)]==0){       #If the so-far-best is a feasible point
    #      if((cobra$numViol[xNewIndex]==0 )&&(cobra$Fres[xNewIndex] < cobra$fbest )){  #If xNew is feasible and even better
    #SB: this condition is always true therefore replaced with Samineh's suggestion   
    #if (!is.null(cobra$ibest)){   # if we have an everBestFeasible at index ibest and fitness value fbest
      #testit::assert(cobra$fn(cobra$A[cobra$ibest,])[1]==cobra$fbest)
      #browser()
    #Samineh's suggestion:
    if(cobra$numViol[cobra$ibest]==0){    # if the so-far best is feasible...
      testit::assert(fn(cobra$A[cobra$ibest,])[1]==cobra$fbest)
      testit::assert(cobra$Fres[cobra$ibest]==cobra$fbest)
      if((cobra$numViol[xNewIndex]==0 )&&(cobra$Fres[xNewIndex] < cobra$fbest )){  #If xNew is feasible and even better
        cobra$xbest<-xNew
        cobra$fbest<-cobra$Fres[xNewIndex]
        cobra$ibest<-xNewIndex    
        #print(cobra$ibest) #do we need this?
        cobra$progressCount<-0
      }
    }else{                                # if we have not yet an everBestFeasible ...
      if(cobra$numViol[xNewIndex]==0){
        cobra$xbest<-xNew                 # ... take xNew, if it is feasible
        cobra$fbest<-cobra$Fres[xNewIndex]
        cobra$ibest<-xNewIndex
        cobra$progressCount<-0
      }
      
    }
    # If we do not have a feasible point AND xNew is not feasible, then leave the triple
    # (xbest,fbest,ibest) at the setting of cobraInitial.R, line 400: From all points 
    # with minimum number of violated constraints, take the one with smallest Fres.
    
    cobra$fbestArray<-c(cobra$fbestArray,cobra$fbest)
    cobra$xbestArray<-rbind(cobra$xbestArray,cobra$xbest)
    feasibleIndices <- which(sapply(1:nrow(cobra$Gres),FUN=function(i)(all(cobra$Gres[i,]<0))))
    xbestIndex<-which.min(cobra$Fres[feasibleIndices])                      # finding index of the best point so far
    
    # only diagnostics, needed for cobra$df2 /WK/
    solu <- cobra$solu; 
    if (is.null(solu)) {
      solu=subMin$par;
    } else {
      if (cobra$rescale) 
        if (is.matrix(solu)) {
          solu <- t(sapply(1:nrow(solu),function(i){ forwardRescale(solu[i,],cobra)}))
        } else {
          solu <- forwardRescale(solu,cobra);
          ### OLD: forwardRescale(solu,cobra$originalL,cobra$originalU);
        }
    }
    # now solu is always in *rescaled* input space
    
    predSoluFunc <- function(x)getPredY(x,fitnessSurrogate,cobra);
    if (is.matrix(solu)) {      # in case of multiple global optima in solu:
      predSolu <- sapply(1:nrow(solu),function(i){ predSoluFunc(solu[i,])}) ;
      predSoluPenal <- sapply(1:nrow(solu),function(i){ fitFuncPenalRBF(solu[i,])}) ;
    } else {
      predSolu <- predSoluFunc(solu);      
      predSoluPenal <- fitFuncPenalRBF(solu);      
    }
    predSolu <- min(predSolu)   # Why min? - In case of multiple global optima: predSolu is the 
                                # value of fitFuncPenalRBF at the best solution solu
    predSoluPenal <- min(predSoluPenal) 
    if (is.null(cobra$df)) {
      df_predSolu <- c(rep(NA,cobra$initDesPoints),predSolu)
    } else {
      df_predSolu <- c(cobra$df$predSolu,predSolu)
    
    }

    # result data frame
    testit::assert("predY",length(cobra$Fres)==length(predY))
    testit::assert("feas",length(cobra$Fres)==length(feas))
    testit::assert("feasPred",length(cobra$Fres)==length(feasPred))
    testit::assert("cobra$numViol",length(cobra$Fres)==length(cobra$numViol))
    testit::assert("cobra$maxViol",length(cobra$Fres)==length(cobra$maxViol))
    testit::assert("cobra$fbestArray",length(cobra$Fres)==length(cobra$fbestArray))
    testit::assert("cobra$optimizationTime",length(cobra$Fres)==length(cobra$optimizationTime))
    testit::assert("optimizerConvergence",length(cobra$Fres)==length(optimizerConvergence))
    df <- data.frame(y=cobra$Fres, 
                     predY=predY,           # surrogate fitness
                     #predSolu=interpRBF(solu,fitnessSurrogate),  # OLD and WRONG
                     predSolu=df_predSolu,
                     feasible=feas, 
                     feasPred=feasPred,
                     nViolations=cobra$numViol,
                     maxViolation=cobra$maxViol,
                     FEval=feval, 
                     Best=cobra$fbestArray,
                     optimizer=rep(cobra$seqOptimizer,length(cobra$Fres)),
                     optimizationTime=cobra$optimizationTime,
                     conv=optimizerConvergence,
                     row.names=NULL
    )
    if (cobra$DEBUG_XI) {
      optimum <- cobra$fn(solu)[1];
      df$fxStart=df_fxStart # objective function at xStart
      df$fxbest=df_fxbest   # objective function at xbest
      df$exbest=df_fxbest - optimum   # error (objective function - optimum) at xbest
      df$RS=df_RS           # TRUE, if it is an iteration with random start point
      df$RS2=df_RS2         # the same
      df$iter2=1:nrow(df)
      df$errFy=df$y - optimum # the error of the optimizer result in every iteration
      df$XI=df_XI
      #if(tail(df_RS,1)==TRUE) browser()
      #browser()
      testit::assert(df$RS==df_RS2)
      if (any(df$fxbest[!df$RS]!=df$fxStart[!df$RS])) {
        browser()
        df$fxbest[!df$RS]=df$fxStart[!df$RS]  # symptomatic fix for the next assert
      }
      testit::assert(df$fxbest[!df$RS]==df$fxStart[!df$RS])
    }
    df <- cbind(iter=1:nrow(df),df)
    df <- cbind(df,seed=cobra$cobraSeed)
    cobra$df <- df
    cobra$df2 <- rbind(cobra$df2,data.frame(
      iter=tail(df$iter,1),
      predY=tail(predY,1),           # surrogate fitness at current point xNew  
      predVal=tail(predVal,1),       # surrogate fitness + penalty at xNew
      predSolu=predSolu,             # surrogate fitness at solu (only diagnostics). 
      predSoluPenal=predSoluPenal,   # surrogate fitness + penalty at solu (only diagnostics). 
      sigmaD=sigmaD[1],
      penaF=penaF[1],
      XI=gama,
      fBest=tail(df$Best,1),
      EPS=EPS[1]
      #,fBest2=fitFuncPenalRBF(xNew)    # not the same as predVal, since penaF or sigmaD might have changed (!)
      #,fSolu=fitFuncPenalRBF(min(solu))# not the same as predSolu for the same reason  
      #,feas=feas,
      #,data.frame(xNew,row.names=NULL)
    ))
    
    # consistency check for data frames df and df2:
    ninit = cobra$initDesPoints # length(which(is.na(cobra$df$predY)))
    msg <- "updateSaveCobra: wrong nrow for df and df2";
    if (is.null(cobra$phase1DesignPoints)) {
      testit::assert(msg,nrow(cobra$df)==nrow(cobra$df2)+cobra$initDesPoints)
    } else {
      testit::assert(msg,nrow(cobra$df)==nrow(cobra$df2)+cobra$phase1DesignPoints)      
    }
    
    if (cobra$saveSurrogates) {
      cobra$constraintSurrogates = constraintSurrogates;
      cobra$fitnessSurrogate = fitnessSurrogate;
      #browser()
    }
    
    if (cobra$saveIntermediate) {
      # save intermediate results
      # cobraResult = list(cobra=cobra, df=df, constraintSurrogates=constraintSurrogates, fn=fn) 
      cobraResult = cobra
      if (is.na(file.info("results")$isdir)) dir.create("results")    # if directory "results" does not exist, create it
      save(cobraResult, file=sprintf("results/cobra-%s-%s-%i.RData",cobra$fName,cobra$seqOptimizer,cobra$cobraSeed))
    }
    
    return(cobra)
  } # updateSaveCobra()
  
  trainSurrogates <- function(cobra) {
    verboseprint(cobra$verbose,important=FALSE,paste(">> Training" ,cobra$RBFmodel,"surrogates","..."))
    ptm <- proc.time()
    cobra$A <- as.matrix(cobra$A)
    cobra$Gres <- as.matrix(cobra$Gres)
    cobra$Fres <- as.vector(cobra$Fres)
    #WK: added an option to build the surrogates only from the lowest 0.9-quantile of Fres
    #    but SKIP_HIGH is now deprecated, since adaptive plog (see below) is better
    if (cobra$SKIP_HIGH) {
      quFres <- quantile(cobra$Fres,0.9)
      ind=which(cobra$Fres<=quFres)
    } else {
      ind=1:length(cobra$Fres)  # all
    }
    A=cobra$A[ind,]
     
       
    #added option of adaptive plog
    if(cobra$sac$aFF){
     # print("adjusting fitness function")
      cobra<-adFit(cobra,ind)
      Fres<-cobra$SurrogateInput
        }else{
      Fres=cobra$Fres[ind]
    }
    
    if(cobra$DOSAC >0){
      if(cobra$PLOG[length(cobra$PLOG)] && printP){
        verboseprint(cobra$verbose, important=TRUE,"PLOG transformation is done")
      }   
    } 
    printP<<-FALSE
    
   # Fres=cobra$Fres[ind]
    #Fres[-ind]=quFres    # test: clip the highest quantile (does not work yet)
    #adjust fitness function with adaptive plog
   # Fres <- switch(as.character(cobra$sac$aFF),
    #       "TRUE" =  adFit(cobra,ind),
    #       "FALSE"=  cobra$Fres[ind]
    #       )
    # --- OLD ---
    #     if(cobra$sac$aFF){
    #       print("adjusting fitness function")
    #       Fres<-adFit(cobra,ind)
    #         }else{
    #       Fres=cobra$Fres[ind]
    #     }
    Gres=cobra$Gres[ind,]  
   
    #SB: added a switch to select type of the RBF model
    sw=switch(cobra$RBFmodel,
              "cubic" =     {constraintSurrogates <<- trainCubicRBF(A,Gres,squares=cobra$squaresC,rho=cobra$RBFrho)
                             fitnessSurrogate <<- trainCubicRBF(A,Fres,squares=cobra$squaresF,rho=cobra$RBFrho)},
              "Gaussian"=   {constraintSurrogates <<- trainGaussRBF(A,Gres,cobra$RBFwidth,squares=cobra$squaresC,RULE=cobra$RULE,rho=cobra$RBFrho);
                             fitnessSurrogate <<- trainGaussRBF(A,Fres,cobra$RBFwidth,squares=cobra$squaresF,RULE=cobra$RULE,rho=cobra$RBFrho)},
              "InvalidRBFmodel"
    )
   
   if (cobra$DEBUG_RBF) {
     if (cobra$dimension!=2) stop("cobra$DEBUG_RBF currently only for d=2")
     
     N=100
     X=outer(rep(1,N),seq(cobra$originalL[1],cobra$originalU[1],length=N))  
     Y=t(X)   
     fac <- (cobra$originalU[1]-cobra$originalL[1])/
            (cobra$newupper-cobra$newlower)
     newWin <- (cobra$initDesPoints==nrow(A))
     ZS <- drawSurrogate3d(fitnessSurrogate,cobra$fName,X,Y,newWindow=newWin,
                           lower=cobra$originalL,upper=cobra$originalU)
     rgl::points3d(cbind(fac*A,Fres),color="white",size=10)
     if (is.null(cobra$TrueZ)) {
       fitFunc <- function(x) {
         z=cobra$originalfn(x)
         return(z[1])
       }
       Z = X*0
       for (i in 1:N)
         for (j in 1:N)
           Z[i,j] = fitFunc(c(X[i,j],Y[i,j]))      
       cobra$TrueZ <- Z
     }
     approxErr = sqrt(sum((cobra$TrueZ-ZS)^2))/(N*N)
     newpnt <- inverseRescale(A[nrow(A),],cobra)
     cat(sprintf("N=%d, new pnt  at (%7.2f,%7.2f), surrogate range = %7.2f, approx err=%7.2f\n",
                 fitnessSurrogate$npts,newpnt[1],newpnt[2],max(ZS)-min(ZS),approxErr))
     wm <- which.min(ZS)
     minpnt <- c(X[wm],Y[wm])
     cat(sprintf("N=%d, surr min at (%7.2f,%7.2f)\n",
                 fitnessSurrogate$npts,minpnt[1],minpnt[2]))
     DO_SNAPSHOTS=F
     if (DO_SNAPSHOTS) {
       if (fitnessSurrogate$npts %% 2 == 0) {
         if (fitnessSurrogate$npts==6) browser()
         if (!file.exists("images.d")) dir.create("images.d")
         rgl::rgl.snapshot(sprintf("images.d/%s-03b-%03d.png",cobra$fName,
                                   fitnessSurrogate$npts))       
       }
     }
   }
   
   #SB: added the possibilty to measure p-effect after every 10 iterations
   if((cobra$sac$onlinePLOG && nrow(cobra$A)%%cobra$sac$onlineFreqPLOG==0 )|| nrow(cobra$A)==cobra$initDesPoints){
    
       Fres1<-cobra$Fres             #without plog
       Fres2<-sapply(cobra$Fres,plog)#with plog
     
       #two models are built after each onlineFreqPLOG iterations:
       #                                            fitnessSurrogate1-> 
     sw=switch(cobra$RBFmodel,
               "cubic" =     {
                              fitnessSurrogate1 <<- trainCubicRBF(A,Fres1,squares=cobra$squaresF,rho=cobra$RBFrho)
                              fitnessSurrogate2 <<- trainCubicRBF(A,Fres2,squares=cobra$squaresF,rho=cobra$RBFrho)},
               "Gaussian"=   {
                              fitnessSurrogate1 <<- trainGaussRBF(A,Fres1,cobra$RBFwidth,squares=cobra$squaresF,RULE=cobra$RULE,rho=cobra$RBFrho)
                              fitnessSurrogate2 <<- trainGaussRBF(A,Fres2,cobra$RBFwidth,squares=cobra$squaresF,RULE=cobra$RULE,rho=cobra$RBFrho)},
               "InvalidRBFmodel"
     ) 
     
   }
   
    testit::assert(sprintf("Wrong value %s for cobra$RBFmodel",sw),sw!="InvalidRBFmodel")
    # --- OLD ---
    #   constraintSurrogates <- trainCubicRBF(cobra$A,cobra$Gres,squares=cobra$squaresC)
    #   fitnessSurrogate <- trainCubicRBF(cobra$A,cobra$Fres,squares=cobra$squaresF)
    DO_ASSERT=F
    if (DO_ASSERT) {
      #might need adjust due to rescale /WK/  
      conFunc <- {function(x)fn(x)[-1];}
      Gres = t(sapply(1:nrow(cobra$A),function(i){conFunc(cobra$A[i,])}))
      testit::assert("Gres-assertion failed",all(Gres==cobra$Gres))
      testit::assert("cobra$A-assertion failed",all(cobra$A==constraintSurrogates$xp))
      Gres = t(sapply(1:nrow(cobra$A),function(i){interpRBF(cobra$A[i,],constraintSurrogates)}))
      for (i in 1:ncol(cobra$Gres)) {
        z = (Gres[,i]-cobra$Gres[,i]) / (max(Gres[,i])-min(Gres[,i]))
        if (max(abs((z))>1e-5)) {
          verboseprint(cobra$verbose, important=FALSE,paste("interpRBF(..,constraintSurrogates)-assertion failed for constraint",i,""))
          print(max(abs(z))) #do we need this?
          browser()
        }        
      }
      cat("All assertions passed\n")
    }
   verboseprint(cobra$verbose, important=FALSE,paste(" finished (",(proc.time() - ptm)[3],"sec )"))
    return(cobra);
    
  } # trainSurrogates()
  
  getPredY <- function(xNew,fitnessSurrogate,cobra) {
    predy <- interpRBF(xNew, fitnessSurrogate)      
    if (cobra$PLOG[length(cobra$PLOG)]) {
      predy <- plogReverse(predy,tail(cobra$pShift,1))
    }
    return (predy)
  }
  getPredY1 <- function(xNew,fitnessSurrogate,cobra) {
    predy <- interpRBF(xNew, fitnessSurrogate)      
    #if (cobra$PLOG[length(cobra$PLOG)]) {
    #  predy <- plogReverse(predy,tail(cobra$pShift,1))
    #}
    return (predy)
  }

  isresCobyla <- function(xStart,fn=subProb2,hin=gCOBRA, cobra) {
    maxeval=cobra$seqFeval; subMin$feval=subMin$iter
    subMin1 <- isres2(xStart,fn=subProb2,lower=cobra$lower,upper=cobra$upper,hin=gCOBRA, maxeval=cobra$seqFeval); 
    subMin2 <- nloptr::cobyla(xStart,fn=subProb2,lower=cobra$lower,upper=cobra$upper,hin=gCOBRA, control=list(maxeval=cobra$seqFeval,xtol_rel=cobra$seqTol)); 
    if (subMin1$value < subMin2$value) {
      nIsres <<- nIsres+1;
      return (subMin1);
    } else {
      nCobyla <<- nCobyla+1;
      return (subMin2);
    }
  }
  
  # -----------------------------------------------------------------------------------------------
  # ---------------- end helper functions cobraPhaseII --------------------------------------------
  # -----------------------------------------------------------------------------------------------
  
  ########################################################################################################
  # STEP6:                                                                                               #
  #  Improve the feasible point                                                                          #
  ########################################################################################################
  printP<-TRUE
  
  while(n < cobra$feval){

    
    ##########################################################
    # STEP6.1: UPDATE RBF MODEL for fitness and constratnit  #
    ##########################################################
    cobra <- trainSurrogates(cobra);  # side effect: constraintSurrogates, fitnessSurrogate
    
    # browser()
    #xp=paste("[",,"]",sep="")
    ##set new surrogate fitness function
    #  this<-paste("xp=",,"squ=",,"coef=",,)
    # evaluate(matlab, this)
    #   if(isOpen(matlab)){
    #    setFunction(matlab, " \
    #     function f=mySFFun(x) \
    #     z = ones(1,size(xp,1)).'*x - xp; \
    #     "); 
    #   }

    ##########################################################
    # STEP6.2: Determine Distance requirement                #
    ##########################################################
    ##gama<-cobra$XI[(cobra$iteration %% length(cobra$XI))+1]
    ##gama<-cobra$XI[(nrow(cobra$A)   %% length(cobra$XI))+1]           # /WK/ bug fix
    gama<-cobra$XI[((nrow(cobra$A)-nRepair) %% length(cobra$XI))+1]     # /WK/ 2nd bug fix
    
    ### only debug: for searching the reason why long DRC is worse than short DRC ###
    #if (n>40) {
    #  # force a switch to short DRC after iteration 40
    #  XI2 = c(0.001,0.0)
    #  gama<-XI2[((nrow(cobra$A)-nRepair) %% length(XI2))+1]  
    #}
    ### --- ###
    
    ro<-gama*cobra$l
    
    
    ##########################################################
    # STEP6.3: Select Iterate                                #
    ##########################################################
    ptm <- proc.time()
    subMin <- list()
    # print_level = 2  # { 0 | 2 } print no or more diagnositc information --> only avail for nloptr, not for wrapper cobyla 
    verboseprint(cobra$verbose, important=FALSE,paste(cobra$seqOptimizer, " optimization on surrogate ..."))
    
    #### /WK/ - only debug - ####
    #     if (nrow(cobra$A)>495) {
    #       save(list=ls(all=TRUE),file="WK-cobyla-freeze.Rdata")
    #       cat("Snapshot before COBYLA saved to WK-cobyla-freeze.Rdata\n")
    #       source("cobyla-WK.R")   # another cobyla-wrapper with print_level=2 and maxtime=120 (seconds)
    #       cat("Sourced a new cobyla with print_level=2 and maxtime=120 (seconds)\n")
    #      
    #     }
    #### /WK/ - only debug - ####
    
    #### /SB/-Random Start algorithm
    if(cobra$sac$RS){
      cobra<-RandomStart(cobra)
      xStart<-cobra$xStart
      #if(any(cobra$xbest!=xStart)) cobra$progressCount<-0
    }else{
      xStart<-cobra$xbest
    }
    
    switch(cobra$seqOptimizer,
           
           #RANDOMLHS={subMin<-RS(fun=subProb,lb=cobra$lower, ub=cobra$upper, n=cobra$seqFeval)},
           COBYLA={ subMin<-nloptr::cobyla(xStart,fn=subProb2,lower=cobra$lower,upper=cobra$upper,hin=gCOBRA,control=list(maxeval=cobra$seqFeval,xtol_rel=cobra$seqTol)); subMin$feval=subMin$iter },
           ISRES ={ subMin<-isres2(xStart,fn=subProb2,lower=cobra$lower,upper=cobra$upper,hin=gCOBRA, maxeval=cobra$seqFeval); subMin$feval=subMin$iter},
           HJKB = { subMin<-dfoptim::hjkb(xStart,fn=subProb,lower=cobra$lower,upper=cobra$upper,control=list(maxfeval=cobra$seqFeval)) },
           NMKB = { subMin<-nmkb2(xStart,fn=subProb,lower=cobra$lower,upper=cobra$upper, control=list(maxfeval=cobra$seqFeval,tol=cobra$seqTol)) },
           #ACTIVECMA  = { subMin <- ActiveOnePlusOneCMAES(xStart, subProb, length(cobra$xbest), opts=list(esname="ActiveCMAES", lb=cobra$lower, ub=cobra$upper, 
           #                maxFunEvals=cobra$seqFeval, 
           #                mu=cobra$seqMu, 
           #                lambda=cobra$seqLambda, 
           #                sigma=cobra$seqStepSize)); 
           #                subMin$convergence <- 1;
           #},
           #RANDOMSEARCH = { subMin <- randomSearch(cobra$xbest, fn=subProb2, lower=cobra$lower,upper=cobra$upper,control=list(maxfeval=cobra$seqFeval, sd=0.05))},
           ISRESCOBY = { subMin<-isresCobyla(xStart,fn=subProb2,hin=gCOBRA, cobra); subMin$feval=subMin$iter},
           DEOPTIM={subMin<-DEoptim::DEoptim(fn=subProb,lower=cobra$lower,upper=cobra$upper,DEoptim::DEoptim.control(itermax=300,trace=Inf)); subMin$feval=subMin$optim$nfeval ; subMin$par=subMin$optim$bestmem;subMin$convergence=NA}
    )
    cobra$optimizationTime <- c(cobra$optimizationTime, (proc.time()-ptm)[3])
   # browser()
    verboseprint(cobra$verbose, important=FALSE,paste(" finished (",subMin$feval,"iterations,",tail(cobra$optimizationTime,1),"sec )"))
    #cat("Optimization time for", subMin$feval, " iterations:", tail(cobra$optimizationTime,1), "seconds\n")  
    #cat("Predicted infill value:",subMin$value,"\n")
    
    if (cobra$DEBUG_XI) {
      # If cobra$DEBUG_XI==TRUE, then print xStart in every iteration 
      # and add later some extra debug info to cobra$df (columns fxStart, fxbest, errY, XI, RS)
      cat("** xStart =",xStart," **\n")
      cobra$xStart = xStart
      cobra$DEBUG_RS = (!all(xStart==cobra$xbest)) 
    }
    
    #if (n==(cobra$feval-1)) browser()   # /WK/ only for analysis G06
    
    
    if (penaF[1]*penaF[2] < penaF[3])
      penaF[1] <- penaF[1]*penaF[2]
    
    
    sigmaD <- checkDistanceReq(subMin,constraintSurrogates,sigmaD,CHECKDIST)    


    ##########################################################
    # STEP6.4: Evaluate real functions                       #
    ##########################################################
    cobra$progressCount<-cobra$progressCount+1
    cobra$fe<-cobra$fe+1
   #browser()
   #cobra$radi<-c(cobra$radi,cobra$radi[length(cobra$radi)])
   #browser()
   
    xNew<-subMin$par
    #cat("X:",xNew,"\n")
    xNew <- pmax(xNew,cobra$lower)  
    xNew <- pmin(xNew,cobra$upper)  
    newPredY <- getPredY(xNew,fitnessSurrogate,cobra) 
    if (cobra$trueFuncForSurrogates) newPredY<-fn(xNew)[1]   
    predY <- c(predY,newPredY)      # bug fix: now predY is the fitness surrogate value /WK/
    predVal <- c(predVal,subMin$value)   # fitness + penalty (in case of NMKB et al.) /WK/
    feval <- c(feval,subMin$feval)
    optimizerConvergence <- c(optimizerConvergence,subMin$convergence)
   
    #cobra$predC <- rbind(cobra$predC,constraintPrediction)
    cobra$predC <- rbind(cobra$predC,interpRBF(xNew,constraintSurrogates))
    #browser()

 
    xNewTemp <- xNew
    xNewEval<-fn(xNewTemp)
    #cat("F/Max(C):",xNewEval[1],"/",max(xNewEval[2:length(xNewEval)]),"\n")
   #SB: calculating the pEffect
   if((nrow(cobra$A)%%cobra$sac$onlineFreqPLOG==0 && cobra$sac$onlinePLOG) || nrow(cobra$A)==cobra$initDesPoints){
    # print(paste("n:",nrow(cobra$A)))
     
     newPredY1 <- getPredY1(xNew,fitnessSurrogate1,cobra)
     newPredY2 <- getPredY1(xNew,fitnessSurrogate2,cobra)
     newErr1<-abs(newPredY1-xNewEval[1])
     newErr2<-abs(plogReverse(newPredY2)-xNewEval[1])
     #newErr2<-abs(newPredY2-xNewEval[1])
     
     err1<-c(err1,newErr1)
     err2<-c(err2,newErr2)
     
     errRatio<-err1/err2
     
     if(is.infinite(newErr2)){
       errRatio[length(errRatio)]<-0 
     }else if(is.infinite(newErr1)){
       errRatio[length(errRatio)]<-Inf
     }

     cobra$pEffect<-log10(quantile(errRatio,na.rm = TRUE)[4])   
   }
    
    newNumViol<-length(which((unlist(xNewEval)[-1]) > cobra$conTol)) # Calculating number of constraint Violations for new point #0 change to conTol
    feas = c(feas, newNumViol < 1 )
    newNumPred<-length(which(cobra$predC[nrow(cobra$predC),] > cobra$conTol)) # the same on constraint surrogates
    feasPred = c(feasPred, newNumPred < 1 ) 
                         
    if((max(0,max((unlist(xNewEval)[-1])) )) > cobra$conTol){ # Calculating maximum violation
      newMaxViol<-max(0,max((unlist(xNewEval)[-1])) )  
    }else{
      newMaxViol<-0
    }
    
     
    
    ##########################################################
    # STEP6.5 & 6.6: Update Information and Counters         #
    ##########################################################
    cobra <- updateInfoAndCounters(cobra,xNew,xNewEval,newNumViol,newMaxViol,phase)
    
    ##########################################################
    # STEP6.7: Adjust Margins                                #
    ##########################################################
    adjustMargins(Cfeas,Tfeas,Cinfeas,Tinfeas,EPS,cobra$epsilonMax)
    
    n<-nrow(cobra$A)
    

    ##########################################################
    # STEP6.8: Update and save cobra                         #
    ##########################################################
    cobra <- updateSaveCobra(cobra,xNew,feas,feasPred,feval,optimizerConvergence,
                             predY,predVal,subMin,sigmaD,penaF,gama,EPS)
    
    
    ##########################################################
    # STEP6.9: Repair Infeasible                             #
    ##########################################################
    
    if(cobra$repairInfeas==TRUE) 
    {
      xNewHasEverBestFitness = (cobra$Fres[length(cobra$numViol)] < cobra$fbest + cobra$ri$marFres) 
      if (!cobra$ri$repairOnlyFresBetter) xNewHasEverBestFitness = TRUE
      # /WK/ if repairOnlyFresBetter=T, repair only iterates with fitness < so-far-best-fitness  
      
      if( (newNumViol>=1) && 
          (xNewHasEverBestFitness) &&   
          (newMaxViol < cobra$ri$repairMargin) &&
          (n < cobra$feval) )            # /WK/ bug fix: no repair, if the last iteration would issue a repair
      {
        if(n %% cobra$verboseIter ==0){#important to print
          cobra$important=TRUE
        }else{
          cobra$important=FALSE
        }
        
        # Build surrogate anew, based on current cobra$A, cobra$Gres
        # This is important for accurate constraint surrogates models near current infeasible point
        cobra <- trainSurrogates(cobra);  # side effect: constraintSurrogates, fitnessSurrogate
        
        repairInfeasible <- repairInfeasRI2
        if (cobra$ri$RIMODE==3) repairInfeasible <- repairChootinan
        #---this is the normal repairInfeasible call (checkit=FALSE)---
        xNewRepaired<-repairInfeasible(subMin$par,xNewEval[-1], constraintSurrogates,cobra)
        #---this is the repairInfeasible call with checkit=TRUE (debug, extra printout) ---
        #xNewRepaired<-repairInfeasible(subMin$par,xNewEval[-1], constraintSurrogates,cobra,TRUE,function(x){fn(x)[-1]})
        
        nRepair <- nRepair +1 
        if(all(xNew==xNewRepaired)){
          print("cannot repair infeasible solution")
        }
        else {
          xNew<-xNewRepaired
          
          ##########################################################
          # STEP R6.4: Evaluate real functions                     #
          ##########################################################
          cobra$fe<-cobra$fe+1
          #browser()
          cobra$radi<-c(cobra$radi,cobra$radi[length(cobra$radi)])
         # browser()
          xNew <- pmax(xNew,cobra$lower)   
          xNew <- pmin(xNew,cobra$upper)  
          newPredY <- getPredY(xNew,fitnessSurrogate,cobra)
          if (cobra$trueFuncForSurrogates) newPredY<-fn(xNew)[1]   
          predY <- c(predY,newPredY)      # bug fix: now predY is the fitness surrogate value /WK/
          predVal <- c(predVal, fitFuncPenalRBF(xNew))   # fitness + penalty (in case of NMKB et al.) /WK/
          feval <- c(feval,NA)
          optimizerConvergence <- c(optimizerConvergence,NA)
          cobra$optimizationTime <- c(cobra$optimizationTime, NA)
          #cobra$predC <- rbind(cobra$predC,constraintPrediction)
          cobra$predC <- rbind(cobra$predC,interpRBF(xNew,constraintSurrogates))
          
         
          
          xNewTemp <- xNew
          xNewEval<-fn(xNewTemp)
          newNumViol<-length(which((unlist(xNewEval)[-1]) > cobra$conTol)) # Calculating number of constraint Violations for new point #0 change to conTol
          feas = c(feas, newNumViol < 1 ) 
          newNumPred<-length(which(cobra$predC[nrow(cobra$predC),] > cobra$conTol)) # the same on constraint surrogates
          feasPred = c(feasPred, newNumPred < 1 ) 
          
          if((max(0,max((unlist(xNewEval)[-1])) )) > cobra$conTol){ # Calculating maximum violation
            newMaxViol<-max(0,max((unlist(xNewEval)[-1])) )  
          }else{
            newMaxViol<-0
          }
          
          ##########################################################
          # STEP R6.5 & 6.6: Update Information and Counters         #
          ##########################################################
          cobra <- updateInfoAndCounters(cobra,xNew,xNewEval,newNumViol,newMaxViol,phase)
          
          ##########################################################
          # STEP R6.7: Adjust Margins                                #
          ##########################################################
          adjustMargins(Cfeas,Tfeas,Cinfeas,Tinfeas,EPS,cobra$epsilonMax)
          
          n<-nrow(cobra$A)
          
          ##########################################################
          # STEP R6.8: Update and save cobra                         #
          ##########################################################
          cobra <- updateSaveCobra(cobra,xNew,feas,feasPred,feval,optimizerConvergence,
                                   predY,predVal,subMin,sigmaD,penaF,gama,EPS)
        } # else of 'all(xNew==xNewRepaired)'
      } 
    } # if(cobra$repairInfeas==TRUE)
    #######################################################################
    
    if(cobra$TrustRegion){ #TRUST REGION: is not wroking yet (in debugging phase)
      cobra<-trustRegion(cobra)
      xNewRefined<-cobra$refinedX

      if(cobra$TRDONE)if(any(cobra$xbest!=xNewRefined)){ xNew<-xNewRefined
      
      
      ##########################################################
      # STEP TRA6.4: Evaluate real functions                     #
      ##########################################################
      cobra$fe<-cobra$fe+1
      xNew <- pmax(xNew,cobra$lower)   
      xNew <- pmin(xNew,cobra$upper)  
      newPredY <- getPredY(xNew,fitnessSurrogate,cobra)
      if (cobra$trueFuncForSurrogates) newPredY<-fn(xNew)[1]   
      predY <- c(predY,newPredY)      # bug fix: now predY is the fitness surrogate value /WK/
      predVal <- c(predVal, fitFuncPenalRBF(xNew))   # fitness + penalty (in case of NMKB et al.) /WK/
      feval <- c(feval,NA)
      optimizerConvergence <- c(optimizerConvergence,NA)
      cobra$optimizationTime <- c(cobra$optimizationTime, NA)
      cobra$radi<-c(cobra$radi,cobra$radi[length(cobra$radi)])
      #cobra$predC <- rbind(cobra$predC,constraintPrediction)
      cobra$predC <- rbind(cobra$predC,interpRBF(xNew,constraintSurrogates))
      
      
      xNewTemp <- xNew
      
      xNewEval<-fn(xNewTemp)
      
      newNumViol<-length(which((unlist(xNewEval)[-1]) > cobra$conTol)) # Calculating number of constraint Violations for new point #0 change to conTol
      feas = c(feas, newNumViol < 1 ) 
      newNumPred<-length(which(cobra$predC[nrow(cobra$predC),] > cobra$conTol)) # the same on constraint surrogates
      feasPred = c(feasPred, newNumPred < 1 ) 
      
      if((max(0,max((unlist(xNewEval)[-1])) )) > cobra$conTol){ # Calculating maximum violation
        newMaxViol<-max(0,max((unlist(xNewEval)[-1])) )  
      }else{
        newMaxViol<-0
      }
      
      ##########################################################
      # STEP TRA6.5 & 6.6: Update Information and Counters         #
      ##########################################################
      cobra <- updateInfoAndCounters(cobra,xNew,xNewEval,newNumViol,newMaxViol,phase)
      
      ##########################################################
      # STEP TRA6.7: Adjust Margins                                #
      ##########################################################
      adjustMargins(Cfeas,Tfeas,Cinfeas,Tinfeas,EPS,cobra$epsilonMax)
      
      n<-nrow(cobra$A)
      ##########################################################
      # STEP TRA6.8: Update and save cobra                         #
      ##########################################################
      cobra <- updateSaveCobra(cobra,xNew,feas,feasPred,feval,optimizerConvergence,
                               predY,predVal,subMin,sigmaD,penaF,gama,EPS)
      }
     #cobra<-adaptRadi(cobra)
    }#TRUST REGION: is not wroking yet (in debugging phase)
    
  } # while(n)
  
  # only info for ISRESCOBY optimizer
  cobra$nCobyla <- nCobyla;
  cobra$nIsres <- nIsres;
  
  return(cobra)
}

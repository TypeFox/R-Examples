#
#Samineh Bagheri, Patrick Koch, Wolfgang Konen
#Cologne Univeristy of Applied Sciences
#
#April,2014 - May,2015
#cobraInitial.R
#


######################################################################################
# cobraInit
#
#'   Initial phase for COBRA optimizer
#'
#'   Constraint-based optimization initialization
#'  
#Detail: 
#'   If \code{epsilonInit} or \code{epsilonMax} are NULL on input, then \code{cobra$epsilonInit}
#'   and \code{cobra$epsilonMax},  resp., are set to \code{0.005*l} where \code{l} is the smallest 
#'   side of the search box.\cr
#'   Note that the parameters \code{penaF}, \code{sigmaD}, \code{constraintHandling} are only 
#'   relevant for penalty-based internal optimizers NMKB or HJKB.
#'   
#'   @param xStart        a vector containing the starting point for the optimization problem
#'   @param fn            objective function that is to be minimized, should return a vector of the 
#'                        objective function value and the constraint values
#'   @param fName         file name .Rdata where the results of \code{\link{cobraPhaseII}} are saved
#'   @param lower         lower bound of search space, same dimension as \code{xStart}
#'   @param upper         upper bound of search space, same dimension as \code{xStart}
#'   @param nConstraints  number of constraints
#'   @param feval         maximum number of function evaluations
#'   @param initDesign    ["RANDOM"] one out of ["RANDOM","LHS","BIASED","OPTIMIZED","OPTBIASED"]
#'   @param initDesPoints number of initial points, must be smaller than feval
#'   @param initDesOptP   [NULL] only for initDesign=="OPTBIASED": number of points for the "OPT"
#'                        phase. If NULL, take initDesPoints.
#'   @param initBias      [0.005] bias for normal distribution in "OPTBIASED" and "BIASED"
#'   @param seqOptimizer  ["COBYLA"] string defining the optimization method for COBRA phases I 
#'                        and II, one out of ["COBYLA","ISRES","HJKB","NMKB","ISRESCOBY"]
#'   @param seqFeval      maximum number of function evaluations on the surrogate model
#'   @param seqTol        [1e-6] Convergence tolerance, see param \code{tol} in \code{\link[dfoptim]{nmkb}}
#    @param seqMu         (deprecated, for ACTIVECMA only)
#    @param seqLambda     (deprecated, for ACTIVECMA only)
#    @param seqStepSize   (deprecated, for ACTIVECMA only) initial (global) step size for optimization on the surrogate model 
#'   @param penaF         [c(3,1.7,3e5)] parameters for dynamic penalty factor (fct subProb in cobraPhaseII): \code{c(start,augment,max)}
#'   @param sigmaD        [c(3,2.0,100)] parameters for dynamic distance factor (fct subProb in cobraPhaseII): \code{c(start,augment,max)}
#'   @param squaresF      [FALSE] set to TRUE for \code{fitnessSurrogate <- trainCubicRBF(..., squares=T)}
#'   @param squaresC      [FALSE] set to TRUE for \code{constraintSurrogates <- trainCubicRBF(..., squares=T)}
#'   @param XI            magic parameters for the distance requirement (DR)
#'   @param drFactor      [1.0] factor multiplied to the DR-constraint (gCOBRA)
#'   @param epsilonInit   initial constant added to each constraint to maintain a certain margin to boundary
#'   @param epsilonMax    maximum for constant added to each constraint 
#'   @param cobraSeed     seed for random number generator
#'   @param conTol        [0.0] constraint violation tolerance
#'   @param repairInfeas  [FALSE] if TRUE, try to repair infeasible solutions
#'   @param repairMargin  -- deprecated --
#'   @param ri            [\code{\link{defaultRI}()}] list with other parameters for 
#'                        \code{\link{repairInfeasRI2}}
#'   @param saveSurrogates [FALSE] if TRUE, then cobraPhaseII returns the last surrogate models in
#'                        cobra$fitnessSurrogate and cobra$constraintSurrogates
#'   @param saveIntermediate [FALSE] if TRUE, then cobraPhaseI + II save intermediate results
#'                        in dir 'results/' (create it, if necessary)
#'   @param RBFmodel      ["cubic"] a string assigning type of the RBF model, "cubic" or "Gaussian"
#'   @param RBFwidth      [-1] only relevant for Gaussian RBF model, see \code{\link{trainGaussRBF}}   
#'   @param RBFrho        [0.0] experimental: 0: interpolating, >0, approximating (spline-like) Gaussian RBFs
#'   @param GaussRule     ["One"] only relevant for Gaussian RBF model, see \code{\link{trainGaussRBF}}                                     
#'   @param trueFuncForSurrogates  [FALSE] if TRUE, use the true (constraint & fitness) functions
#'                        instead of surrogates (only for debug analysis)
#'   @param skipPhaseI    [FALSE] if TRUE, then skip \code{\link{cobraPhaseI}}
#'   @param rescale       [TRUE] if TRUE change \code{[lower,upper]} to hypercube \code{[newlower,newupper]^d}
#'   @param newlower      [-1] lower bound of each rescaled input space dimension, if rescale==TRUE
#'   @param newupper      [+1] upper bound of each rescaled input space dimension, if rescale==TRUE
#'   @param solu          [NULL] the best-known solution (only for diagnostics). This is normally a 
#'                        vector of length d. If there are multiple solutions, it is a matrix with d
#'                        columns (each row is a solution). If NULL, then the current best point
#'                        will be used in \code{\link{cobraPhaseII}}. 
#'                        \code{solu} is given in original input space.
#'   @param TrustRegion   [FALSE] if TRUE, an embedded  trust region algorithm \code{\link{trustRegion}} is performed. 
#'   @param TRlist        [\code{\link{defaultTR}()}] a list of parameters, needed only 
#'                        in case \code{TrustRegion==TRUE}.  
#'   @param sac           [\code{\link{defaultSAC}(DOSAC)}] list with other parameters for SACOBRA  
#'   @param DOSAC         [0|1|2] if >0, any elements of \code{sac} not set by the user are set to \code{defaultSAC(DOSAC)}.
#'                        0: COBRA-R settings, 1: SACOBRA settings, 2: SACOBRA settings with fewer parameters
#'                        and more online adujustement (aFF and aCF are done parameter free).
#'   @param constraintHandling ["DEFAULT"] (other choices: "JOINESHOUCK", "SMITHTATE", "COIT", "BAECKKHURI";
#'                        experimental, only for penalty-based internal optimizers NMKB or HJKB, 
#'                        see the code in function \code{subProb} in \code{\link{cobraPhaseII}})          
#'   @param DEBUG_XI      [FALSE] if TRUE, then print in cobraPhaseII extra debug information: 
#'                        xStart in every iteration to console and add some extra debug 
#'                        columns to cobra$df
#'   @param DEBUG_RBF     [FALSE] visualize RBF (only for dimension==2)                      
#'   @param SKIP_HIGH     [FALSE] (deprecated) if TRUE, then build the surrogate models in cobraPhaseII by 
#'                        skipping the data points having Fres (objective function) in the highest decile.
#'   @param verbose       [1] one out of [0|1|2], how much output to print
#'   @param verboseIter   [10]                   
#'                        
#'   @return \code{cobra}, an object of class COBRA, this is a (long) list containing most
#'   of the argument settings (see above) and in addition (among others)
#'      \item{\code{A}}{ (feval x dim)-matrix containing the initial design points in input .  
#'            space. If rescale==TRUE, all points are in  \strong{rescaled} input space. }
#'      \item{\code{Fres}}{ a vector of the objective values of the initial design points }
#'      \item{\code{Gres}}{ a matrix of the constraint values of the initial design points }
#'      \item{\code{Tfeas}}{ the threshhold parameter for the number of consecutive iterations 
#'            that yield feasible solutions before margin epsilon is reduced }
#'      \item{\code{Tinfeas}}{ the threshhold parameter for the number of consecutive iterations 
#'            that yield infeasible solutions before margin epsilon is increased }
#'      \item{\code{numViol}}{ number of constraint violations }
#'      \item{\code{maxViol}}{ maximum constraint violation }
#'      \item{\code{refinedX}}{A vector of all refined solutions generated by trust region algorithm (see \code{trustRegion})}
#'   
#'   Note that \code{cobra$Fres}, \code{cobra$fbest}, \code{cobra$fbestArray} and similar contain 
#'   always the objective values of the orignial function \code{cobra$fn[1]}. (The surrogate models 
#'   may be trained on a \code{\link{plog}}-transformed version of this function.)
#' 
#' @seealso   \code{\link{startCobra}}, \code{\link{cobraPhaseI}}, \code{\link{cobraPhaseII}}
#' @author Wolfgang Konen, Samineh Bagheri, Patrick Koch, Cologne Univeristy of Applied Sciences
#' @export
#' 
#' 
#' 
######################################################################################
cobraInit <- function(xStart, fn, fName, lower, upper, nConstraints, feval, 
                      initDesign="RANDOM", 
                      initDesPoints=2*length(xStart)+1, initDesOptP=NULL, initBias=0.005,
                      seqOptimizer="COBYLA", seqFeval=1000, seqTol=1e-6, 
                      # seqMu=2, seqLambda=10, seqStepSize=0.05,   # deprecated, for ACTIVECMA only
                      penaF=c(3.0, 1.7, 3e5), squaresF=TRUE, squaresC=TRUE, conTol=0.0,
                      constraintHandling="DEFAULT",
                      sigmaD=c(3.0,2.0,100), 
                      repairInfeas=FALSE, repairMargin=NULL,ri=defaultRI(),
                      DOSAC=1, sac=defaultSAC(DOSAC), 
                      epsilonInit=NULL, epsilonMax=NULL, solu=NULL,
                      saveIntermediate=FALSE, 
                      saveSurrogates=FALSE, RBFmodel="cubic", RBFwidth=-1,GaussRule="One",
                      RBFrho=0.0,
                      skipPhaseI=TRUE,trueFuncForSurrogates=FALSE, drFactor=1.0, XI=DRCL,
                      # teta=c(0.1,0.05,0.01,0.005,0.001,0.0005), 
                      #teta=DRCL,
                      rescale=TRUE,newlower=-1,newupper=1,  
                      DEBUG_XI=FALSE, SKIP_HIGH=FALSE, DEBUG_RBF=FALSE,
                      TrustRegion=FALSE,TRlist=defaultTR(),
                      verbose=1,verboseIter=10,cobraSeed ){
  ##
  ## /WK/ TODO: RBFwidth=-1 would lead to strange behaviour, if used with RBFmodel="Gaussian"
  ## 
  if (!is.null(repairMargin)) {
    stop("Parameter repairMargin no longer supported! Please use parameter ri=defaultRI(repairMargin).")
  }
  originalfn<- fn
  originalL <- lower
  originalU <- upper
  phase<-"init"
  
  dimension<-length(xStart)         # number of parameters
  #browser()
  if(rescale){
    lb<-rep(newlower,dimension)
    up<-rep(newupper,dimension)
    xStart<-sapply(1:dimension , function(i){scales::rescale(xStart[i],to=c(lb[i],up[i]),from=c(originalL[i],originalU[i]))
    })
    fn<-rescaleWrapper(fn,originalL,originalU,dimension,newlower,newupper)
    lower<-lb
    upper<-up
  }
  #browser()
  l <- min(upper - lower) # length of smallest side of search space
  if (is.null(epsilonInit)) epsilonInit<- 0.005*l
  if (is.null(epsilonMax))  epsilonMax <- 2*0.005*l
  if (is.null(initDesOptP)) initDesOptP <- initDesPoints
  
  testit::assert("", initDesPoints<feval)
  testit::assert("", length(fn(xStart))==(nConstraints+1))
  
  set.seed(cobraSeed)
  dimension<-length(xStart)         # number of parameters
  iteration<-0  
  # We just have one starting point
  xStart <- xStart
  
  I<-matrix()
  Gres <- matrix()
  Fres <- c()
  randomInitialPoints <- 10 # additional random points, required for optimized design to avoid crashing of SVD in the next phases

  
 
  # helper for switch(initDesign): evaluate random solutions in I on fn
  randomResultsFactory <- function(I,fn,dimension) {
    sapply(1:nrow(I), function(i){
      verboseprint(verbose=0,important=FALSE,sprintf("iteration %03d: ",i))
      x<-I[i,]    
      res<-fn(x)
      return(res)
    })
  }
  
  # helper for switch(initDesign): clip all entries in I[,j] to be between lower[j] and upper[j] 
  clipLowerUpper <- function(I,dimension) {
    I <- t(unlist(sapply(1:nrow(I), FUN=function(i)pmax(I[i,],lower))))
    I <- t(unlist(sapply(1:nrow(I), FUN=function(i)pmin(I[i,],upper))))
    #     I1 = I*0
    #     for(i in 1:nrow(I)){
    #       I1[i,] <- unlist(sapply(1:dimension, FUN=function(j)max(lower[j],I[i,j]))) # lower-->lower[j]
    #       I1[i,] <- unlist(sapply(1:dimension, FUN=function(j)min(upper[j],I1[i,j]))) # upper-->upper[j]
    #     }
    #     testit::assert("check I",all(I1==I))
    return(I)    
  }
  
  cat("\n")
  sw=switch(initDesign,
         "RANDOM" = {
           set.seed(cobraSeed)
        
           I <- t(as.matrix(sapply(1:(initDesPoints-1), FUN=function(i)stats::runif(dimension,lower,upper))))    
                     
           I <- rbind(I,xStart)
           I <- clipLowerUpper(I,dimension)
           
           randomResults<-randomResultsFactory(I,fn,dimension)
           # update structures for random solutions
           Gres <- matrix(t(randomResults[2:nrow(randomResults),]),ncol=nConstraints,nrow=initDesPoints) #Changed to be adapted to case of having 1 constraint
           Fres <- as.vector(randomResults[1,])                                                          #should be chnaged for other initialization as well
           names(Fres) <- NULL
           colnames(I)=NULL
           rownames(I)=NULL
           "ok"
         },
         
         "LHS" = {      # latin hypercube sampling + xStart
           #require(lhs)
           set.seed(cobraSeed)
           I <- lhs::randomLHS(initDesPoints-1,dimension)    # LHS with values uniformly in [0,1]
           for (k in 1:ncol(I)) {
             I[,k] <- lower[k] + I[,k]*(upper[k]-lower[k])
           }
           I <- rbind(I,xStart)
           I <- clipLowerUpper(I,dimension)
           
           randomResults<-randomResultsFactory(I,fn,dimension)
           # update structures for random solutions
           #SB: In order to adapt the code ato address unconstraint problems
           if(nConstraints==0){
             Gres<-NULL
             Fres <- as.vector(randomResults)                                                          
             
           }else{
             Gres <- matrix(t(randomResults[2:nrow(randomResults),]),ncol=nConstraints,nrow=initDesPoints) #Changed to be adapted to case of having 1 constraint  
             Fres <- as.vector(randomResults[1,])                                                          
           }
           names(Fres) <- NULL
           colnames(I)=NULL
           rownames(I)=NULL
           "ok"
         },
         
         "BIASED" = {  # biased initial population based on starting point
           I <- t(as.matrix(sapply(1:(initDesPoints-1), FUN=function(i)rnorm(dimension,as.vector(xStart),initBias))))
           names(xStart)=NULL
           I <- rbind(I,xStart)
           I <- clipLowerUpper(I,dimension)
           
           randomResults<-randomResultsFactory(I,fn,dimension)
           # update structures for random solutions
           Gres <- matrix(t(randomResults[-1,]),ncol=nConstraints,nrow=ncol(randomResults))
           Fres <- as.vector(randomResults[1,])
           names(Fres) <- NULL
           "ok"
         },
         
         "BIASEDINF" = {  # infeasible initial population based on starting point
           I <- t(as.matrix(sapply(1:(initDesPoints-1), FUN=function(i)rnorm(dimension,as.vector(xStart),initBias))))
           names(xStart)=NULL
           I <- rbind(I,xStart)
           I <- clipLowerUpper(I,dimension)
           
           randomResults<-randomResultsFactory(I,fn,dimension)
           # update structures for random solutions
           Gres <- matrix(t(randomResults[-1,]),ncol=nConstraints,nrow=ncol(randomResults))
           Fres <- as.vector(randomResults[1,])
           names(Fres) <- NULL
           #browser()
           infeasibleSolution = sapply(1:initDesPoints, FUN=function(i)any(Gres[i,]>0))
           ind = which(infeasibleSolution)
           I <- I[ind,]
           Gres <- as.matrix(Gres[ind,])
           Fres <- Fres[ind]
           nInfeasible <- nrow(I)
           
           while(nInfeasible < initDesPoints){
             x <- runif(dimension,0,1)
             I <- rbind(I,x)
             nInfeasible <- nInfeasible + 1
             #evaluate random solution on fn
             res<-fn(x) 
             Gres <- rbind(Gres,res[2:length(res)])
             Fres <- c(Fres,res[1])
           }
           Gres <- as.matrix(Gres)
           colnames(I)=NULL
           rownames(I)=NULL
           names(Fres) <- NULL
           "ok"
         },
         
         "OPTIMIZED" =  { # optimized (HJKB) initial points    
           # start Hooke & Jeeves initial search
           initialHookeJeeves <- initialHjkb(xStart,fn=fn,lower=lower,upper=upper,control=list(maxfeval=initDesPoints) )
           I <- as.matrix(initialHookeJeeves$historyData$xArchive)
           duplicatedIndices <- which(duplicated(I))
           if(length(duplicatedIndices)!=0){
             I <- I[-duplicatedIndices,]
           }
           
           I <- clipLowerUpper(I,dimension)
           colnames(I)=NULL
           rownames(I)=NULL
           Gres <- as.matrix(initialHookeJeeves$historyData$constraintArchive)
           if(length(duplicatedIndices)!=0){
           Gres <- Gres[-duplicatedIndices,]
           }
           Fres <- initialHookeJeeves$historyData$yArchive
           if(length(duplicatedIndices)!=0){
             Fres <- Fres[-duplicatedIndices]
           }

           names(Fres) = NULL
           #TODO probably need to reset initial design size, because Hooke&Jeeves exceeds budget:
           initDesPoints <- length(Fres)
           "ok"           
         },
            
         "OPTCOBYLA"=, "OPTBIASED"= {  # COBYLA-optimized initial points + BIASED design
              # start COBYLA initial search
              #
              #resetSoluArchive()    # this was necessary for fnArchive_OLD.R
              #setArchiveFunc(fn)
              fnArchiveF <- fnArchiveFactory(fn);
              initialCobyla <- nloptr::cobyla(xStart,fn=function(x){fnArchiveF(x)[1]},lower=lower,upper=upper
                                      , hin=function(x){-fn(x)[-1]} 
                                      , control=list(maxeval=initDesOptP,xtol_rel = 1e-9))
              I <- (environment(fnArchiveF))$getSoluArchive()
              I <- clipLowerUpper(I,dimension)
              
              
              randomResults<-randomResultsFactory(I,fn,dimension)
              # update structures for random solutions
              Gres <- matrix(t(randomResults[-1,]),ncol=nConstraints,nrow=ncol(randomResults))
              Fres <- as.vector(randomResults[1,])
              
              if (initDesign=="OPTBIASED") {
                # find the best feasible point (if any, else the best infeasible point) ...
                feas<-apply(Gres,1,function(x)all(x<=0))
                if (any(feas==TRUE)) Fres[!feas] <- Inf
                xStart <- I[which.min(Fres),]
                #
                # ... and construct initDesPoints 'BIASED' design points around it
                I <- t(as.matrix(sapply(1:(initDesPoints-1), FUN=function(i)rnorm(dimension,as.vector(xStart),initBias))))
                names(xStart)=NULL
                I <- rbind(I,xStart)
                I <- clipLowerUpper(I,dimension)
                randomResults<-randomResultsFactory(I,fn,dimension)
                # update structures for random solutions
                Gres <- matrix(t(randomResults[-1,]),ncol=nConstraints,nrow=ncol(randomResults))
                Fres <- as.vector(randomResults[1,])
              }

              names(Fres) <- NULL
              #testit::assert("",initDesPoints==length(Fres))
              #TODO probably need to reset initial design size, because COBYLA varies budget:
              initDesPoints <- length(Fres)
              "ok"           
            },
            
         "InvalidInitDesign"
  ) # switch
  cat("\n")
  testit::assert(sprintf("Wrong value %s for initDesign",initDesign),sw!="InvalidInitDesign")
  
  # parameters for cycling distance requirement
  #teta<-teta   #distance requirement phaseI
  XI<-XI#c(0.01, 0.001, 0.0005)      #distance requirement phaseII
  
  Tfeas<-floor(2*sqrt(dimension) )  # The threshhold parameter for the number of consecutive iterations that yield feasible solution before the margin is reduced
  Tinfeas<-floor(2*sqrt(dimension)) # The threshold parameter for the number of consecutive iterations that yield infeasible solutions before the margin is increased
  Nmax<-feval
  
  
  fe<-nrow(I) #number of function evaluations
  fbest<-c() # best function value found
  xbest<-c() # best point found # a numeric vector with the length of dimension
  fbestArray<-c()
  xbestArray<-matrix()
    
  ########################################################################################################
  # Update structures                                                                              #
  ########################################################################################################
  # browser()
  #testit::assert("Gres is not a matrix!",is.matrix(Gres) || nConstraints!=0)
  numViol<-sapply(1:initDesPoints , function(i){ # number of initial Constraint violations
    return(sum(Gres[i,]>0))
  })
  
  maxViol<-sapply(1:initDesPoints , function(i){
    y<-max(0,Gres[i,])
    return(y)
  })
  
  A<-I  # A contains all evaluated points
  n<-nrow(A)

  # determining best feasible objective value (fbest) so far and best point (xbest)
  if(0 %in% numViol){
    fbest<-min(Fres[which(numViol==0)])
    xbest <- I[which.min(fbest),]
    ibest <- which(Fres==fbest)
  }else{
    # if there is no feasibe point yet, take one from the points with min number of violated constraints:
    minNumIndex<-which(numViol==min(numViol))
    # /WK/ the folllowing lines can select pretty big fbest values which lead to a very bad pShift:
    #index<-minNumIndex[which(maxViol[minNumIndex]==min(maxViol[minNumIndex]))]
    #fbest <- Fres[index[1]]
    #xbest <- A[index[1],]
    #ibest <- index[1]
    # /WK/ alternatve: select among the min-num-viol point the one with smallest Fres
    FresMin <- Fres[minNumIndex]
    ind <- which.min(FresMin)
    index <- minNumIndex[ind]
    
    fbest <- Fres[index[1]]
    xbest <- A[index[1],]
    ibest <- index[1]
  }
  fbestArray<-rep(fbest,initDesPoints)
  xbestArray<-xbest
  
  for(i in c(2:n)){
    xbestArray<-rbind(xbestArray,xbest)
  }
  cat("*** Starting run with seed",cobraSeed,"***\n")
  print("Initialization is done")
  
  cobra<-list(fn=fn,
            xStart=xbest,
            fName=fName,
            dimension=dimension, 
            nConstraints=nConstraints, 
            lower=lower,
            upper=upper,
            newlower=newlower,
            newupper=newupper,
            originalL=originalL,
            originalU=originalU,
            originalfn=originalfn,
            rescale=rescale,
            feval=feval, 
            A=A, 
            fbestArray=fbestArray, 
            xbestArray=xbestArray, 
            xbest=xbest, 
            fbest=fbest,
            ibest=ibest,
            Fres=Fres, 
            Gres=Gres, 
            numViol=numViol,
            maxViol=maxViol,
            epsilonInit=epsilonInit, 
            epsilonMax=epsilonMax, 
            XI=XI, 
            drFactor=drFactor,
            Tinfeas=Tinfeas,
            Tfeas=Tfeas,
            iteration=iteration, 
            initDesPoints=initDesPoints,
            seqOptimizer=seqOptimizer,
            seqFeval=seqFeval,
            seqTol=seqTol,
            #seqMu=seqMu,           # deprecated, only for ACTIVECMA
            #seqLambda=seqLambda,
            #seqStepSize=seqStepSize,
            penaF=penaF,
            sigmaD=sigmaD,
            squaresF=squaresF,
            squaresC=squaresC,
            cobraSeed=cobraSeed,
            conTol=conTol,
            constraintHandling=constraintHandling,
            l=l, 
            repairInfeas=repairInfeas,
            #repairMargin=repairMargin,   # now in cobra$ri
            ri=ri,
            fe=fe,
            saveIntermediate=saveIntermediate,
            saveSurrogates=saveSurrogates,
            RBFmodel=RBFmodel,
            RBFwidth=RBFwidth,
            RBFrho=RBFrho,
            RULE=GaussRule,
            skipPhaseI=skipPhaseI,
            trueFuncForSurrogates=trueFuncForSurrogates,
            solu=solu,
            TrustRegion=TrustRegion,
            TRlist=TRlist,
            radi=rep(TRlist$radiInit,initDesPoints),
            sac=sac,
            TRDONE=rep(FALSE,initDesPoints),
            TRind=c(),
            refinedX=c(),
            fCount=0,
            sCount=0,
            DOSAC=DOSAC,
            PLOG=FALSE,
            pShift=0,
            pEffect=NA,
            progressCount=0,
            DEBUG_XI=DEBUG_XI,
            DEBUG_RBF=DEBUG_RBF,
            SKIP_HIGH=SKIP_HIGH,
            verbose=verbose,
            verboseIter=verboseIter,
            phase=rep(phase,initDesPoints)
            )
 # browser()
  ## /WK/ if cobra$ri
  cobra$ri <- setOpts(cobra$ri,defaultRI())
  ## /SB/: extension to SACOBRA
  if(DOSAC>0) {     
    verboseprint(cobra$verbose, important=FALSE,"Parameter and function adjustment phase")
    cobra$sac<-setOpts(cobra$sac,defaultSAC(DOSAC))
    cobra$pEffect<-cobra$sac$pEffectInit
    cobra$TRlist<-setOpts(cobra$TRlist,defaultTR())
   # if (any(XI!=DRCL) & cobra$sac$aDRC==TRUE) {
   #   warning("XI is different from default (DRCL), but sac$aDRC==TRUE, so XI will be set by automatic DRC adjustment!")
   # }
    if(cobra$sac$aDRC){
      if(length(XI)!=length(DRCL)){
        warning("XI is different from default (DRCL), but sac$aDRC==TRUE, so XI will be set by automatic DRC adjustment!")  
      }else if(any(XI!=DRCL)){
        warning("XI is different from default (DRCL), but sac$aDRC==TRUE, so XI will be set by automatic DRC adjustment!")
      }
      
      verboseprint(cobra$verbose,important=FALSE,"adjusting DRC")
      DRC<-adDRC(max(cobra$Fres),min(cobra$Fres))
      cobra$XI<-DRC
    }
    
    # --- adFit is now called in *each* iteration of cobraPhaseII (adaptive plog) ---
    #
    # if(res$sac$aFF){
    #   print("adjusting fitness function")
    #   res<-adFit(res)
    #   res$fn<-newfn 
    # }
    
    if(cobra$sac$aCF && nConstraints!=0 ){
      verboseprint(cobra$verbose,important=FALSE,"adjusting constraint functions")
      cobra<-adCon(cobra)
      #cobra$fn<-fn  
    }
  } # (DOSAC)
  
  class(cobra) <- c("COBRA","list")
  
  return(cobra)
}

verboseprint<-function(verbose, important, message){
  if(verbose!=0){ # if verbose == 0 do not print anything
                  # if verbose == 1 print only important messages 
                  # if verbose == 2 print everything
    if((verbose==2) || ((verbose==1) && (important==TRUE)) ){
      print(message)
    }
  }
}
verbosecat<-function(verbose, message, important=FALSE){
  if(verbose!=0){ # if verbose== 0 do not cat anything
                  # if verbose ==1 cat only important messages 
                  # if verbose ==2 cat everything
    if((verbose==2) || ((verbose==1) && (important==TRUE)) ){
      cat(message)
    }
  }
}


#' 
#' 
#' 

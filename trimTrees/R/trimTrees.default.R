# ----------------------------------------------------------------
# R function: hitRate
# ----------------------------------------------------------------
# This function calculates the average hit rate for a group of 
# forecasts. The function has as its arguments the forcasts PIT values
# and the interval of interest. 
# ----------------------------------------------------------------

"hitRate" <- function(matrixPIT, interval=c(0.25,0.75)){
  if(any(is.na(matrixPIT))) stop("NAs not permitted in matrixPIT.")
  matrixPIT <- as.matrix(matrixPIT)
  numberOfForecasters <- ncol(matrixPIT)
  numberOfForecasts <- nrow(matrixPIT)
  HR <- rep(0,numberOfForecasters)
  for(i in 1:numberOfForecasters){
    HR[i] <- sum((matrixPIT[,i]>=interval[1]) & (matrixPIT[,i]<=interval[2]))/numberOfForecasts
  }
  return(HR)
}

# ----------------------------------------------------------------
# R function: trimTrees
# ----------------------------------------------------------------

"trimTrees" <-
  function(xtrain,
           ytrain, 
           xtest,
           ytest=NULL,
           ntree=500,
           mtry=if (!is.null(ytrain) && !is.factor(ytrain))
             max(floor(ncol(xtrain)/3), 1) else floor(sqrt(ncol(xtrain))),
           nodesize=if (!is.null(ytrain) && !is.factor(ytrain)) 5 else 1,
           trim=0.0, 
           trimIsExterior=TRUE,
           uQuantiles = seq(0.05,0.95,0.05), 
           methodIsCDF=TRUE) {
    
    ## Some checks  
    if(is.null(nrow(xtest)) || is.null(ncol(xtest)))
      stop("xtest contains no data.")
    
    if(is.null(ytest)){
      ytest <- rep(0,nrow(xtest))
      isytestNull <- TRUE
    } else isytestNull <- FALSE
    
    if(! class(ytest) %in% c("numeric","integer", "factor") )
      stop("ytest must be numeric or factor.")
       
    testdat <- !is.null(xtest)
    if (testdat) {
      if (ncol(xtrain) != ncol(xtest))
        stop("xtrain and xtest must have same number of columns.")
    }
    
    if( nrow(xtest) != length(ytest) )
      stop("xtest and ytest must contain the same number of samples.")
    
    if (any(is.na(xtest))) stop("NA not permitted in xtest.")
    if (any(is.na(ytest))) stop("NA not permitted in ytest.")
    
    if(! class(trim) %in% c("numeric","integer"))
      stop("trim must be numeric.")

    if(trim <0 | trim > 0.5)
      stop("trim must be greater than or equal to 0 and less than or equal to 0.5.")

    if(! class(uQuantiles) %in% c("numeric","integer") )
      stop("uQuantiles must be numeric")
    
    if(sum(uQuantiles <=0) + sum(uQuantiles>= 1) != 0)
      stop("uQuantiles must be greater than 0 and less than 1.")
    
    if(is.unsorted(uQuantiles, strictly=TRUE))
      stop("uQuantiles must be a strictly increasing sequence.")
    
          ntrain <- nrow(xtrain)
          ntest  <- nrow(xtest)
          # Run cinbag (i.e., a modified version of random forest) to get the inbag counts.
          rf <- cinbag(xtrain, ytrain, ntree=ntree, nodesize=nodesize, mtry=mtry, keep.inbag=TRUE)
          inbagCount <- rf$inbagCount
          # Find the terminal node numbers for the training and testing rows.
          # You can use predict on the object rf.  cinbag outputs a "random forest" object.
          termNodetrain <- attr(predict(rf,xtrain,nodes=TRUE),"nodes")
          termNodeNewX <- attr(predict(rf,xtest,nodes=TRUE),"nodes")
          
          # Initialize inputs.
          trainingDataSort <- sort(as.double(ytrain), method = "shell", index.return = TRUE)
          ytrainSortedIndex <- trainingDataSort$ix
          inbagCountSorted <- inbagCount[ytrainSortedIndex,]
          termNodetrainSorted <- termNodetrain[ytrainSortedIndex,]
          ytrainSorted <- trainingDataSort$x
          forestSupport <- unique(ytrainSorted)
          nSupport <- length(forestSupport)
          nQuantiles <- length(uQuantiles)
          indextest <- rep(0,ntree)

          # Initialize tree outputs.
          treeValues <- matrix(0,ntrain,ntree)
          treeCounts <- matrix(0,nSupport,ntree)
          treeCumCounts <- matrix(0,nSupport + 1,ntree)
          treeCDFs <- matrix(0,nSupport + 1,ntree)
          treePMFs <- matrix(0,nSupport, ntree)
          treeMeans <- matrix(0,ntest,ntree)
          treeVars <- matrix(0,ntest,ntree)
          treePITs <- matrix(0,ntest,ntree)
          treeQuantiles <- matrix(0,ntree,nQuantiles)
          treeFirstPMFValues <- matrix(0,ntest,ntree)
              
          # Initialize ensemble outputs.
          bracketingRate <- rep(0,ntest)
          bracketingRateAllPairs <- matrix(0,ntree,ntree)
          trimmedEnsembleCDFs <- matrix(0,ntest,nSupport + 1)
          trimmedEnsemblePMFs <- matrix(0,ntest,nSupport)
          trimmedEnsembleMeans <- rep(0,ntest)
          trimmedEnsembleVars <- rep(0,ntest)
          trimmedEnsemblePITs <- rep(0,ntest) 
          trimmedEnsembleQuantiles <- rep(0,nQuantiles)
          trimmedEnsembleComponentScores <- matrix(0,nQuantiles,2)
          trimmedEnsembleScores <- matrix(0,ntest,4)
          untrimmedEnsembleCDFs <- matrix(0,ntest,nSupport + 1)
          untrimmedEnsemblePMFs <- matrix(0,ntest,nSupport)
          untrimmedEnsembleMeans <- rep(0,ntest)
          untrimmedEnsembleVars <- rep(0,ntest)
          untrimmedEnsemblePITs <- rep(0,ntest) 
          untrimmedEnsembleQuantiles <- rep(0,nQuantiles)
          untrimmedEnsembleComponentScores <- matrix(0,nQuantiles,2)
          untrimmedEnsembleScores <- matrix(0,ntest,4)  
          
          tol <- 5*.Machine$double.eps
          
          if(methodIsCDF){
            ttout <- .C("trimTreesCDF",
                        
                        # random forest inputs.
                        as.integer(inbagCountSorted),
                        as.integer(termNodetrainSorted),
                        as.integer(ntree), 
                        as.double(ytrainSorted),
                        as.integer(ntrain),
                        forestSupport=as.double(forestSupport),
                        as.integer(nSupport),           
                        as.integer(termNodeNewX), 
                        as.double(ytest),
                        as.integer(ntest),
                        
                        # user inputs.
                        as.double(trim),
                        as.logical(trimIsExterior),
                        as.double(uQuantiles),
                        as.integer(nQuantiles),
                        
                        # tree outputs.
                        treeValues=as.double(treeValues),
                        treeCounts=as.double(treeCounts),
                        treeCumCounts=as.double(treeCumCounts),
                        treeCDFs=as.double(treeCDFs),
                        treePMFs=as.double(treePMFs),
                        treeMeans=as.double(treeMeans),
                        treeVars=as.double(treeVars),
                        treePITs=as.double(treePITs),
                        treeQuantiles=as.double(treeQuantiles),
                        treeFirstPMFValues=as.double(treeFirstPMFValues),      
                        
                        # ensemble outputs.
                        bracketingRate=as.double(bracketingRate),     
                        bracketingRateAllPairs=as.double(bracketingRateAllPairs),
                 
                        trimmedEnsembleCDFs=as.double(trimmedEnsembleCDFs),
                        trimmedEnsemblePMFs=as.double(trimmedEnsemblePMFs),
                        trimmedEnsembleMeans=as.double(trimmedEnsembleMeans),
                        trimmedEnsembleVars=as.double(trimmedEnsembleVars),
                        trimmedEnsemblePITs=as.double(trimmedEnsemblePITs),
                        trimmedEnsembleQuantiles=as.double(trimmedEnsembleQuantiles),
                        trimmedEnsembleComponentScores=as.double(trimmedEnsembleComponentScores),
                        trimmedEnsembleScores=as.double(trimmedEnsembleScores),
                        
                        untrimmedEnsembleCDFs=as.double(untrimmedEnsembleCDFs),
                        untrimmedEnsemblePMFs=as.double(untrimmedEnsemblePMFs),
                        untrimmedEnsembleMeans=as.double(untrimmedEnsembleMeans),
                        untrimmedEnsembleVars=as.double(untrimmedEnsembleVars),
                        untrimmedEnsemblePITs=as.double(untrimmedEnsemblePITs),
                        untrimmedEnsembleQuantiles=as.double(untrimmedEnsembleQuantiles),
                        untrimmedEnsembleComponentScores=as.double(untrimmedEnsembleComponentScores),
                        untrimmedEnsembleScores=as.double(untrimmedEnsembleScores),
                        
                        tol=as.double(tol),
                        PACKAGE="trimTrees"
            )
          
          }
          else{
            ttout <- .C("trimTreesMA",
                        
                        # random forest inputs.
                        as.integer(inbagCountSorted),
                        as.integer(termNodetrainSorted),
                        as.integer(ntree), 
                        as.double(ytrainSorted),
                        as.integer(ntrain),
                        forestSupport=as.double(forestSupport),
                        as.integer(nSupport),           
                        as.integer(termNodeNewX), 
                        as.double(ytest),
                        as.integer(ntest),
                        
                        # user inputs.
                        as.double(trim),
                        as.logical(trimIsExterior),
                        as.double(uQuantiles),
                        as.integer(nQuantiles),
                        
                        # tree outputs.
                        treeValues=as.double(treeValues),
                        treeCounts=as.double(treeCounts),
                        treeCumCounts=as.double(treeCumCounts),
                        treeCDFs=as.double(treeCDFs),
                        treePMFs=as.double(treePMFs),
                        treeMeans=as.double(treeMeans),
                        treeVars=as.double(treeVars),
                        treePITs=as.double(treePITs),
                        treeQuantiles=as.double(treeQuantiles),
                        treeFirstPMFValues=as.double(treeFirstPMFValues),      
                        
                        # ensemble outputs.
                        bracketingRate=as.double(bracketingRate),     
                        bracketingRateAllPairs=as.double(bracketingRateAllPairs),
                        
                        trimmedEnsembleCDFs=as.double(trimmedEnsembleCDFs),
                        trimmedEnsemblePMFs=as.double(trimmedEnsemblePMFs),
                        trimmedEnsembleMeans=as.double(trimmedEnsembleMeans),
                        trimmedEnsembleVars=as.double(trimmedEnsembleVars),
                        trimmedEnsemblePITs=as.double(trimmedEnsemblePITs),
                        trimmedEnsembleQuantiles=as.double(trimmedEnsembleQuantiles),
                        trimmedEnsembleComponentScores=as.double(trimmedEnsembleComponentScores),
                        trimmedEnsembleScores=as.double(trimmedEnsembleScores),
                        
                        untrimmedEnsembleCDFs=as.double(untrimmedEnsembleCDFs),
                        untrimmedEnsemblePMFs=as.double(untrimmedEnsemblePMFs),
                        untrimmedEnsembleMeans=as.double(untrimmedEnsembleMeans),
                        untrimmedEnsembleVars=as.double(untrimmedEnsembleVars),
                        untrimmedEnsemblePITs=as.double(untrimmedEnsemblePITs),
                        untrimmedEnsembleQuantiles=as.double(untrimmedEnsembleQuantiles),
                        untrimmedEnsembleComponentScores=as.double(untrimmedEnsembleComponentScores),
                        untrimmedEnsembleScores=as.double(untrimmedEnsembleScores),
                        
                        tol=as.double(tol),
                        PACKAGE="trimTrees"
            )
       }
        if(isytestNull){
          ttout$treePITs <- rep(NA,ntest*ntree)
          ttout$bracketingRate <- rep(NA, ntest)
          ttout$bracketingRateAllPairs <- rep(NA, ntree*ntree)
          ttout$trimmedEnsemblePITs <- rep(NA, ntest)
          ttout$trimmedEnsembleComponentScores <- rep(NA,nQuantiles*2)
          ttout$trimmedEnsembleScores <- rep(NA,ntest*4)
          ttout$untrimmedEnsemblePITs <- rep(NA, ntest)
          ttout$untrimmedEnsembleComponentScores <-rep(NA,nQuantiles*2)
          ttout$untrimmedEnsembleScores <- rep(NA, ntest*4)
        }
          out <- list(forestSupport=forestSupport,
                      treeValues = matrix(ttout$treeValues,ntrain,ntree),
                      treeCounts = matrix(ttout$treeCounts,nSupport,ntree),
                      treeCumCounts = matrix(ttout$treeCumCounts,nSupport + 1,ntree),
                      treeCDFs = matrix(ttout$treeCDFs,nSupport + 1,ntree),
                      treePMFs = matrix(ttout$treePMFs,nSupport, ntree),
                      treeMeans = matrix(ttout$treeMeans,ntest,ntree),
                      treeVars = matrix(ttout$treeVars,ntest,ntree),
                      treePITs = matrix(ttout$treePITs,ntest,ntree),
                      treeQuantiles = matrix(ttout$treeQuantiles,ntree,nQuantiles),
                      treeFirstPMFValues = matrix(ttout$treeFirstPMFValues,ntest,ntree),
                      bracketingRate = ttout$bracketingRate,
                      bracketingRateAllPairs = matrix(ttout$bracketingRateAllPairs,ntree,ntree),
                      trimmedEnsembleCDFs = matrix(ttout$trimmedEnsembleCDFs,ntest,nSupport + 1),
                      trimmedEnsemblePMFs = matrix(ttout$trimmedEnsemblePMFs,ntest,nSupport),
                      trimmedEnsembleMeans = ttout$trimmedEnsembleMeans,
                      trimmedEnsembleVars = ttout$trimmedEnsembleVars,
                      trimmedEnsemblePITs = ttout$trimmedEnsemblePITs,
                      trimmedEnsembleQuantiles = ttout$trimmedEnsembleQuantiles,
                      trimmedEnsembleComponentScores = matrix(ttout$trimmedEnsembleComponentScores,nQuantiles,2),
                      trimmedEnsembleScores = matrix(ttout$trimmedEnsembleScores,ntest,4, dimnames=list(seq(1,ntest,1), c("LinQuanS", "LogQuanS", "RPS", "TMS"))),
                      untrimmedEnsembleCDFs = matrix(ttout$untrimmedEnsembleCDFs,ntest,nSupport + 1),
                      untrimmedEnsemblePMFs = matrix(ttout$untrimmedEnsemblePMFs,ntest,nSupport),
                      untrimmedEnsembleMeans = ttout$untrimmedEnsembleMeans,
                      untrimmedEnsembleVars = ttout$untrimmedEnsembleVars,
                      untrimmedEnsemblePITs = ttout$untrimmedEnsemblePITs,
                      untrimmedEnsembleQuantiles = ttout$untrimmedEnsembleQuantiles,
                      untrimmedEnsembleComponentScores = matrix(ttout$untrimmedEnsembleComponentScores,nQuantiles,2),
                      untrimmedEnsembleScores = matrix(ttout$untrimmedEnsembleScores,ntest,4,dimnames=list(seq(1,ntest,1), c("LinQuanS", "LogQuanS", "RPS", "TMS")))       
                      )     
          class(out) <- "trimTree"
          return(out)
  }        
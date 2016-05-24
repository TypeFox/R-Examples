
#################################################################
# Manipulation of the "ComparisonResults" objects
#################################################################


# =====================================================
# Function that joins  experimental results objects.
# The joining is carried out by some specified dimension,
# the most common being joining experiments carried out in
# different data sets (dimension 4), or experiments with
# different learners (dimension 3) on the same data sets.
# =====================================================
# Luis Torgo, Aug 2013
# =====================================================
# Example runs:
# > bestScores(mergeEstimationRes(subset(earth,stats='e1',vars=1:3),
#                   subset(nnet,stats='e1',vars=4:6),by=3))
# > bestScores(mergeEstimationRes(nnet,earth,rf,rpartXse,svm,by=3))
#
mergeEstimationRes <- function(...,by='tasks') {

  s <- list(...)
  if (length(s) < 2) return(s[[1]])
  
  if ((! by %in% c('metrics','workflows','tasks')) &&
      (! by %in% 1:3))
    stop('mergeEstimationRes:: invalid value on "by" argument!')
  for(i in 2:length(s)) 
    if (!identical(s[[i]][[1]][[1]]@estTask@method,s[[1]][[1]][[1]]@estTask@method))
      stop('mergeEstimationRes:: trying to join performance estimation objects with different estimation settings!')

  if (!is.numeric(by))
    by <- match(by,c('metrics','workflows','tasks'))

  if (by == 1) {
      sameTasks <- sapply(s[2:length(s)],function(x) identical(names(s[[1]]),names(x)))
      if (!all(sameTasks)) stop("mergeEstimationRes:: to join by metrics all objects need to address the same tasks!")
      sameWFs <- sapply(s[2:length(s)],function(x) identical(names(s[[1]][[1]]),names(x[[1]])))
      if (!all(sameWFs)) stop("mergeEstimationRes:: to join by metrics all objects need to use the same workflows!")
      for(e in s[2:length(s)]) 
          for(t in 1:length(e))
              for(w in 1:length(e[[t]])) {
                  s[[1]][[t]][[w]]@iterationsScores <- cbind(s[[1]][[t]][[w]]@iterationsScores,
                                                                   e[[t]][[w]]@iterationsScores)
                  s[[1]][[t]][[w]]@estTask@metrics <- c(s[[1]][[t]][[w]]@estTask@metrics,e[[t]][[w]]@estTask@metrics)
              }

  } else if (by == 2) {
      sameTasks <- sapply(s[2:length(s)],function(x) identical(names(s[[1]]),names(x)))
      if (!all(sameTasks)) stop("mergeEstimationRes:: to join by workflows all objects need to address the same tasks!")
      sameStats <- sapply(s[2:length(s)],function(x) identical(colnames(s[[1]][[1]][[1]]@iterationsScores),colnames(x[[1]][[1]]@iterationsScores)))
      if (!all(sameStats)) stop("mergeEstimationRes:: to join by workflows all objects need to estimate the same metrics!")
      for(e in s[2:length(s)]) 
          for(t in 1:length(e))
              s[[1]][[t]] <- c(s[[1]][[t]],e[[t]])
  } else if (by == 3) {
      sameWFs <- sapply(s[2:length(s)],function(x) identical(names(s[[1]][[1]]),names(x[[1]])))
      if (!all(sameWFs)) stop("mergeEstimationRes:: to join by tasks all objects need to use the same workflows!")
      sameStats <- sapply(s[2:length(s)],function(x) identical(colnames(s[[1]][[1]][[1]]@iterationsScores),colnames(x[[1]][[1]]@iterationsScores)))
      if (!all(sameStats)) stop("mergeEstimationRes:: to join by tasks all objects need to estimate the same metrics!")
      for(e in s[2:length(s)]) 
          s[[1]] <- c(s[[1]],e)
      
  }
  return(ComparisonResults(s[[1]]))
}

# =====================================================
# Small auxiliary functions to obtain information from 
# ComparisonResults objects.
# =====================================================
# Luis Torgo, Aug 2013
# =====================================================
taskNames      <- function(o) names(o)
workflowNames  <- function(o) names(o[[1]])
metricNames    <- function(o) o[[1]][[1]]@estTask@metrics



  
## ======================================================================
## Obtaining some information from the estimation experiments
## ======================================================================
## Luis Torgo, Jan 2009, 2014
## ======================================================================


## ----------------------------------------------------------
## The scores on all iterations
## ----------------------------------------------------------
getScores <- function(results,workflow,task) {
  if (!inherits(results,"ComparisonResults")) stop(results," is not of class 'ComparisonResults''.\n")

  results[[task]][[workflow]]@iterationsScores
}


## ----------------------------------------------------------
## The scores on all iterations
## ----------------------------------------------------------
getIterationsInfo <- function(obj,workflow=1,task=1,rep,fold,it) {
    if (missing(rep) && missing(fold) && missing(it)) return(obj[[task]][[workflow]]@iterationsInfo)
    if ((missing(rep) || missing(fold)) && missing(it))
        stop("getIterationsInfo:: to get the results of a particular iteration you need to supply both 'rep' and 'fold', or simply 'it'")
    if (!missing(it)) {
        if (it > nrow(obj[[task]][[workflow]]@iterationsScores)) stop(paste("getIterationsInfo:: only",nrow(obj[[task]][[workflow]]@iterationsScores),"iterations available.\n"))
        obj[[task]][[workflow]]@iterationsInfo[[it]]
    } else {
        if (rep >  obj[[task]][[workflow]]@estTask@method@nReps || fold > obj[[task]][[workflow]]@estTask@method@nFolds) stop(paste("getIterationsInfo:: only",obj[[task]][[workflow]]@estTask@method@nReps,"repetitions and",obj[[task]][[workflow]]@estTask@method@nFolds,"folds available.\n"))
        obj[[task]][[workflow]]@iterationsInfo[[(rep-1)*obj[[task]][[workflow]]@estTask@method@nFolds+fold]]
    }
}

## ----------------------------------------------------------
## The scores on all iterations
## ----------------------------------------------------------
getIterationsPreds <- function(obj,workflow=1,task=1,rep,fold,it,predComp="preds") {
    if (missing(rep) && missing(fold) && missing(it))
        return(sapply(obj[[task]][[workflow]]@iterationsInfo, function(x) x[[predComp]]))
    if ((missing(rep) || missing(fold)) && missing(it))
        stop("getPredictionsInfo:: to get the results of a particular iteration you need to supply both 'rep' and 'fold', or simply 'it'")
    if (!missing(it)) {
        if (it > nrow(obj[[task]][[workflow]]@iterationsScores)) stop(paste("getIterationInfo:: only",nrow(obj[[task]][[workflow]]@iterationsScores),"iterations available.\n"))
        obj[[task]][[workflow]]@iterationsInfo[[it]][[predComp]]
    } else {
        if (rep >  obj[[task]][[workflow]]@estTask@method@nReps || fold > obj[[task]][[workflow]]@estTask@method@nFolds) stop(paste("getIterationInfo:: only",obj[[task]][[workflow]]@estTask@method@nReps,"repetitions and",obj[[task]][[workflow]]@estTask@method@nFolds,"folds available.\n"))
        obj[[task]][[workflow]]@iterationsInfo[[(rep-1)*obj[[task]][[workflow]]@estTask@method@nFolds+fold]][[predComp]]
    }
}



# ======================================================================
# Get some summary statistics of all evaluation metrics of the performance
# of a workflow on a certain task.
# =====================================================
# Luis Torgo, Jan 2009, 2014
# =====================================================
estimationSummary <- function(results,workflow,task) {
  if (!inherits(results,"ComparisonResults")) stop(results," is not of class 'ComparisonResults''.\n")
  .scores2summary(results[[task]][[workflow]])
}
  


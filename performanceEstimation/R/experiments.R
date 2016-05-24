################################################################# 
## THIS FILE CONTAINS FUNCTIONS THAT ARE RELATED TO RUNNING     #
## EXPERIMENTS WITH MODELLING TOOLS                             #
#################################################################
## Author : Luis Torgo (ltorgo@dcc.fc.up.pt)     Date: Aug 2013 #
## License: GPL (>= 2)                                          #
#################################################################





###################################################################
## FUNCTIONS FOR CARRYING OUT PERFORMANCE ESTIMATION EXPERIMENTS ## 
###################################################################


## ==============================================================
##
## ==============================================================
## Luis Torgo, Aug 2013
## ==============================================================
##
performanceEstimation <- function(tasks,workflows,estTask,...) {

  if (!is(tasks,'list')) tasks <- list(tasks)
  if (!is(workflows,'list')) workflows <- list(workflows)
  
  if (is.null(names(workflows)))
    names(workflows) <- paste('var',1:length(workflows),sep='.')
  
  nTasks <- length(tasks)
  taskNames <- sapply(tasks,function(x) x@taskName)
  nWFs <- length(workflows)
  wfNames <- sapply(workflows,function(x) x@name)
  
  allRes <- vector("list",nTasks)
  names(allRes) <- taskNames
  
  cat('\n\n##### PERFORMANCE ESTIMATION USING ',
      switch(class(estTask@method),
             CV='CROSS VALIDATION',
             Holdout='HOLD OUT',
             MonteCarlo='MONTE CARLO',
             Bootstrap='BOOTSTRAP',
             LOOCV='LOOCV',
             ),
      ' #####')
  
  for(d in 1:length(tasks)) {

    cat('\n\n** PREDICTIVE TASK ::',tasks[[d]]@taskName)

    taskRes <- vector("list",nWFs)
    names(taskRes) <- wfNames
    ##rr <- NULL
    for (s in 1:length(workflows)) {

      cat('\n\n++ MODEL/WORKFLOW ::',workflows[[s]]@name,"\n")
      
      taskRes[[s]] <- do.call(
               switch(class(estTask@method),
                      CV='cvEstimates',
                      Holdout='hldEstimates',
                      Bootstrap='bootEstimates',
                      MonteCarlo='mcEstimates',
                      LOOCV='loocvEstimates'
                      ),

               c(list(workflows[[s]],
                    tasks[[d]],
                    estTask),...)
                         )
    }
    
    allRes[[d]] <- taskRes

  }
  
  ComparisonResults(allRes)
}



#################################################################
## Cross Validation Experiments
#################################################################



## =====================================================
## Function that performs a cross validation experiment
## of a system on a given data set.
## The function is completely generic. The generality comes
## from the fact that the function that the user provides
## as the system to evaluate, needs in effect to be a
## user-defined function that takes care of the learning,
## testing and calculation of the statistics that the user
## wants to estimate through cross validation. 
## =====================================================
## Luis Torgo, Jan 2009
## =====================================================
## Example runs:
## x <- cvEstimates(Workflow('cv.rpartXse',se=2),
##                      PredTask(medv~.,Boston),
##                      CvTask(1,10,1234))
##
cvEstimates <- function(wf,task,estTask,cluster) {
    ## the function to use on the loop (either sequential or parallel)
    `%fun%` <- if (!missing(cluster)) foreach::`%dopar%` else foreach::`%do%`

    ## registering (and eventually creating) the parallel backend
    if (!missing(cluster)) {
        if (is(cluster,"cluster")) doParallel::registerDoParallel(cluster)
        else {
            cl <- parallel::makeCluster(parallel::detectCores()%/%2)
            doParallel::registerDoParallel(cl)
            on.exit(parallel::stopCluster(cl))
        }
        cat(sprintf('cvEstimates: Running in parallel with %d worker(s)\n', foreach::getDoParWorkers()))
    }
            
        
    show(estTask)

    ## Did the user supplied the data splits for all folds and repetitions?
    userSplit <- !is.null(estTask@method@dataSplits)
    
    n <- nrow(eval(task@dataSource))
    if (!userSplit) n.each.part <- n %/% estTask@method@nFolds
    
    nits <- estTask@method@nFolds*estTask@method@nReps
    itsInfo <- vector("list",nits)

    ## Stratified sampling stuff
    if (!userSplit && estTask@method@strat) { 
        respVals <- responseValues(task@formula,eval(task@dataSource))
        regrProb <- is.numeric(respVals)
        if (regrProb) {  # regression problem
            ## the bucket to which each case belongs  
            b <- cut(respVals,10)  # this 10 should be parametrizable
        } else {
            b <- respVals
        }
        ## how many on each bucket
        bc <- table(b)
        ## how many should be on each test partition
        bct <- bc %/% estTask@method@nFolds

    }

    permutation <- 1:n
    it <- NULL  # dummy assignment due to Note on cran-check
    if (missing(cluster)) cat("Iteration :")
    itsInfo <- foreach::foreach(it=1:nits,
                                .packages=.loadedPackages(),
                                .export=if (class(task@dataSource) == "data.frame") NULL else c(as.character(task@dataSource))
                                ) %fun%  {
        
        if (missing(cluster)) cat(" ",it)
        nfold <- (it - 1) %% estTask@method@nFolds + 1
        nrep <- (it - 1) %/% estTask@method@nFolds + 1
        
        if (!userSplit) {
            set.seed(estTask@method@seed*nrep)
            permutation <- sample(n)
        }
        
        if (!userSplit) {
            if (estTask@method@strat) {
                out.fold <- c()
                for(x in seq(along=levels(b))) 
                    if (bct[x]) out.fold <- c(out.fold,which(b[permutation] == levels(b)[x])[((nfold-1)*bct[x]+1):((nfold-1)*bct[x]+bct[x])])
            } else {
                out.fold <- ((nfold-1)*n.each.part+1):(nfold*n.each.part)
            }
        } else out.fold <- outFold(estTask@method@dataSplits,it)
        
        it.res <- runWorkflow(wf,
                              task@formula,
                                        #perm.data[-out.fold,],
                              eval(task@dataSource)[permutation[-out.fold],],
                                        #perm.data[out.fold,])
                              eval(task@dataSource)[permutation[out.fold],])

        c(it.res,list(train=permutation[-out.fold]))
        
    }
    
    ## randomize the number generator to avoid undesired
    ## problems caused by inner set.seed()'s
    set.seed(prod(as.integer(unlist(strsplit(strsplit(date()," ")[[1]][4],":")))))
    
    ## Calculate the metrics estimation
    scores <- .scoresIts(task,estTask,itsInfo)
    
    EstimationResults(task,wf,estTask,scores,itsInfo)
}

#################################################################
# Hold Out Experiments
#################################################################



# =====================================================
# Function that performs a hold out experiment
# of a system on a given data set.
# The function is completely generic. The generality comes
# from the fact that the function that the user provides
# as the system to evaluate, needs in effect to be a
# user-defined function that takes care of the learning,
# testing and calculation of the statistics that the user
# wants to estimate through hold out. A few example
# functions are provided (cv.rpartXse, cv.lm, cv.nnet)
# =====================================================
# Luis Torgo, Feb 2010
# =====================================================
# Example runs:
# x <- hldEstimates(learner('cv.rpartXse',list(se=2)),
#              dataset(medv~.,Boston),
#              hldTask(4,0.25,1234))
#
hldEstimates <- function(wf,task,estTask,cluster) {
    ## the function to use on the loop (either sequential or parallel)
    `%fun%` <- if (!missing(cluster)) foreach::`%dopar%` else foreach::`%do%`

    ## registering (and eventually creating) the parallel backend
    if (!missing(cluster)) {
        if (is(cluster,"cluster")) doParallel::registerDoParallel(cluster)
        else {
            cl <- parallel::makeCluster(parallel::detectCores()%/%2)
            doParallel::registerDoParallel(cl)
            on.exit(parallel::stopCluster(cl))
        }
        cat(sprintf('hldEstimates: Running in parallel with %d worker(s)\n', foreach::getDoParWorkers()))
    }

    show(estTask)

    ## Did the user supplied the data splits for all folds and repetitions?
    userSplit <- !is.null(estTask@method@dataSplits)

    n <- nrow(eval(task@dataSource))
    if (!userSplit) n.test <- as.integer(n * estTask@method@hldSz)

    itsInfo <- vector("list",estTask@method@nReps)

    if (!userSplit & estTask@method@strat) {  # stratified sampling
        respVals <- responseValues(task@formula,eval(task@dataSource))
        regrProb <- is.numeric(respVals)
        if (regrProb) {  # regression problem
            ## the bucket to which each case belongs  
            b <- cut(respVals,10)  # this 10 should be parameterizable
        } else {
            b <- respVals
        }
        ## how many on each bucket
        bc <- table(b)
        ## how many should be on each test partition
        bct <- as.integer(bc * estTask@method@hldSz)
    }

    permutation <- 1:n

    r <- NULL  # dummy assignment due to Note on cran-check
    if (missing(cluster)) cat("Iteration :")
    itsInfo <- foreach::foreach(r=1:estTask@method@nReps,
                                .packages=.loadedPackages(),
                                .export=if (class(task@dataSource) == "data.frame") NULL else c(as.character(task@dataSource))
                                ) %fun%  {


        if (missing(cluster)) cat(' ',r)

        if (!userSplit) {
            set.seed(estTask@method@seed*r)
            permutation <- sample(n)
        } 


        if (!userSplit) {
            if (estTask@method@strat) {
                out.fold <- c()
                for(x in seq(along=levels(b))) 
                    if (bct[x]) out.fold <- c(out.fold,which(b[permutation] == levels(b)[x])[1:bct[x]])
            } else {
                out.fold <- 1:n.test
            }
        } else out.fold <- outFold(estTask@method@dataSplits,r)

        it.res <- runWorkflow(wf,
                              task@formula,
                              eval(task@dataSource)[permutation[-out.fold],],
                              eval(task@dataSource)[permutation[out.fold],])

        c(it.res,list(train=permutation[-out.fold]))
          
    }
    cat('\n')
  
    ## randomize the number generator to avoid undesired
    ## problems caused by inner set.seed()'s
    set.seed(prod(as.integer(unlist(strsplit(strsplit(date()," ")[[1]][4],":")))))

    ## Calculate the metrics estimation
    scores <- .scoresIts(task,estTask,itsInfo)
    
    EstimationResults(task,wf,estTask,scores,itsInfo)
}





#################################################################
# Leave One Out Cross Validation (LOOCV) Experiments
#################################################################



# =====================================================
# Function that performs a LOOCV experiment
# of a system on a given data set.
# The function is completely generic. The generality comes
# from the fact that the function that the user provides
# as the system to evaluate, needs in effect to be a
# user-defined function that takes care of the learning,
# testing and calculation of the statistics that the user
# wants to estimate through hold out. 
# =====================================================
# Luis Torgo, Mar 2010
# =====================================================
# Example runs:
# x <- loocvEstimates(learner('cv.rpartXse',list(se=2)),
#            dataset(medv~.,Boston))
#
loocvEstimates <- function(wf,task,estTask,verbose=FALSE,cluster) {
    ## the function to use on the loop (either sequential or parallel)
    `%fun%` <- if (!missing(cluster)) foreach::`%dopar%` else foreach::`%do%`

    ## registering (and eventually creating) the parallel backend
    if (!missing(cluster)) {
        if (is(cluster,"cluster")) doParallel::registerDoParallel(cluster)
        else {
            cl <- parallel::makeCluster(parallel::detectCores()%/%2)
            doParallel::registerDoParallel(cl)
            on.exit(parallel::stopCluster(cl))
        }
        cat(sprintf('loocvEstimates: Running in parallel with %d worker(s)\n', foreach::getDoParWorkers()))
    }

    show(estTask)

    ## Did the user supplied the data splits for all folds and repetitions?
    userSplit <- !is.null(estTask@method@dataSplits)

    n <- nrow(eval(task@dataSource))

    itsInfo <- vector("list",n)

    r <- NULL  # dummy assignment due to Note on cran-check
    if (verbose && missing(cluster)) cat("Iteration :")
    itsInfo <- foreach::foreach(r=1:n,
                                .packages=.loadedPackages(),
                                .export=if (class(task@dataSource) == "data.frame") NULL else c(as.character(task@dataSource))
                                ) %fun%  {

        if (verbose && missing(cluster)) cat('*')

        if (!userSplit) {
            set.seed(estTask@method@seed*r)
            out.fold <- r
        } else out.fold <- outFold(estTask@method@dataSplits,r)

        it.res <- runWorkflow(wf,
                              task@formula,
                              eval(task@dataSource)[-out.fold,],
                              eval(task@dataSource)[out.fold,])
        
        c(it.res,list(train=(1:n)[-out.fold]))
    }
    if (verbose && missing(cluster)) cat('\n')
    
    ## randomize the number generator to avoid undesired
    ## problems caused by inner set.seed()'s
    set.seed(prod(as.integer(unlist(strsplit(strsplit(date()," ")[[1]][4],":")))))
    
    ## Calculate the metrics estimation
    scores <- .scoresIts(task,estTask,itsInfo)
    
    EstimationResults(task,wf,estTask,scores,itsInfo)

}




#################################################################
# Bootstrap Experiments
#################################################################



# =====================================================
# Function that performs a bootstrap experiment
# of a system on a given data set.
# The function is completely generic. The generality comes
# from the fact that the function that the user provides
# as the system to evaluate, needs in effect to be a
# user-defined function that takes care of the learning,
# testing and calculation of the statistics that the user
# wants to estimate through cross validation. 
# =====================================================
# Luis Torgo, Apr 2010
# =====================================================
# Example runs:
# x <- bootEstimates('cv.rpartXse',list(se=2)),
#                      dataset(medv~.,Boston),
#                      bootTask(1234,10))
#
bootEstimates <- function(wf,task,estTask,cluster) {
    ## the function to use on the loop (either sequential or parallel)
    `%fun%` <- if (!missing(cluster)) foreach::`%dopar%` else foreach::`%do%`

    ## registering (and eventually creating) the parallel backend
    if (!missing(cluster)) {
        if (is(cluster,"cluster")) doParallel::registerDoParallel(cluster)
        else {
            cl <- parallel::makeCluster(parallel::detectCores()%/%2)
            doParallel::registerDoParallel(cl)
            on.exit(parallel::stopCluster(cl))
        }
        cat(sprintf('bootEstimates: Running in parallel with %d worker(s)\n', foreach::getDoParWorkers()))
    }

    show(estTask)

    if (estTask@method@type == '.632')
        resub <- runWorkflow(wf,task@formula,eval(task@dataSource),eval(task@dataSource))
    
    ## Did the user supplied the data splits for all folds and repetitions?
    userSplit <- !is.null(estTask@method@dataSplits)

    n <- nrow(eval(task@dataSource))
    
    itsInfo <- vector("list",estTask@method@nReps)

    r <- NULL  # dummy assignment due to Note on cran-check
    if (missing(cluster)) cat("Iteration :")
    itsInfo <- foreach::foreach(r=1:estTask@method@nReps,
                                .packages=.loadedPackages(),
                                .export=if (class(task@dataSource) == "data.frame") NULL else c(as.character(task@dataSource))
                                ) %fun%  {


        if (missing(cluster)) cat(' ',r)

        if (!userSplit) {
            set.seed(estTask@method@seed*r)
            idx <- sample(n,n,replace=T)
            it.res <- runWorkflow(wf,
                                  task@formula,
                                  eval(task@dataSource)[idx,],
                                  eval(task@dataSource)[-idx,])
            c(it.res,list(train=idx))
        } else {
            it.res <- runWorkflow(wf,
                                  task@formula,
                                  eval(task@dataSource)[outFold(estTask@method@dataSplits,r,"train"),],
                                  eval(task@dataSource)[outFold(estTask@method@dataSplits,r),])
            c(it.res,list(train=outFold(estTask@method@dataSplits,r,"train")))

        }
      
    }
    cat('\n')

    ## randomize the number generator to avoid undesired
    ## problems caused by inner set.seed()'s
    set.seed(prod(as.integer(unlist(strsplit(strsplit(date()," ")[[1]][4],":")))))
  
    ## Calculate the metrics estimation
    if (estTask@method@type == ".632") {  # this method is different from all others
        trReq <- estTask@trainReq || any(estTask@metrics %in% c("nmse","nmae","theil")) || (is.regression(task) && is.null(estTask@metrics))
        
        nIts <- length(itsInfo)
        
        standEval <- if (estTask@evaluator == "" ) TRUE else FALSE
        if (standEval) 
            evalFunc <- if (is.classification(task)) "classificationMetrics" else "regressionMetrics"
        else
            evalFunc <- estTask@evaluator

        wts <- intersect(estTask@metrics,c("trTime","tsTime","totTime"))
        predMs <- setdiff(estTask@metrics,wts)
        stats <- if (!is.null(predMs)) list(stats=predMs) else NULL

        ## getting the resubstitution scores
        trR <- if (trReq) list(train.y=eval(task@dataSource)[1:n,task@target]) else NULL
        fstArgs <- if (standEval) list(trues=resub$trues,preds=resub$preds) else resub
        resubScores <- do.call(evalFunc,
                               c(fstArgs,
                                 stats,
                                 trR,
                                 estTask@evaluator.pars))

        ## The structure holding all scores
        if (is.null(estTask@metrics)) {
            ncols <- length(resubScores)
            namcols <- names(resubScores)
        } else {
            ncols <- length(estTask@metrics)
            namcols <- estTask@metrics
        }
        scores <- matrix(NA,nrow=nIts,ncol=ncols,
                         dimnames=list(1:nIts,namcols))

        for(i in 1:nIts) {
            trR <- if (trReq) list(train.y=eval(task@dataSource)[itsInfo[[i]]$train,task@target]) else NULL
            fstArgs <- if (standEval) list(trues=itsInfo[[i]]$trues,preds=itsInfo[[i]]$preds) else itsInfo[[i]]
            ss <- 0.632*do.call(evalFunc,
                                c(fstArgs,
                                  stats,
                                  trR,
                                  estTask@evaluator.pars)
                                ) + 0.368*resubScores
            
            if (is.null(predMs))  scores[i,] <- ss else scores[i,predMs] <- ss

            if (length(wts)) {
                allts <- as.numeric(itsInfo[[i]]$times)
                scores[i,wts] <- c(trTime=allts[1],tsTime=allts[2],
                                   totTime=allts[1]+allts[2])[wts]
            }
        }
        
    } else scores <- .scoresIts(task,estTask,itsInfo)
    
    
    EstimationResults(task,wf,estTask,scores,itsInfo)

}



#################################################################
## Monte Carlo Experiments
#################################################################



## =====================================================
## Function that performs a Monte Carlo experiment of a 
## system on a given data set.
## The function is completely generic. The generality comes
## from the fact that the function that the user provides
## as the system to evaluate, needs in effect to be a
## user-defined function that takes care of the learning,
## testing and calculation of the statistics that the user
## wants to estimate through this experiment. A few example
## functions are provided.
## =====================================================
## Luis Torgo, Aug 2009
## =====================================================

mcEstimates <- function(wf, task, estTask, verbose=TRUE, cluster) {
    ## the function to use on the loop (either sequential or parallel)
    `%fun%` <- if (!missing(cluster)) foreach::`%dopar%` else foreach::`%do%`

    ## registering (and eventually creating) the parallel backend
    if (!missing(cluster)) {
        if (is(cluster,"cluster")) doParallel::registerDoParallel(cluster)
        else {
            cl <- parallel::makeCluster(parallel::detectCores()%/%2)
            doParallel::registerDoParallel(cl)
            on.exit(parallel::stopCluster(cl))
        }
        cat(sprintf('mcEstimates: Running in parallel with %d worker(s)\n', foreach::getDoParWorkers()))
    }


    show(estTask)

    ## Did the user supplied the data splits for all  repetitions?
    userSplit <- !is.null(estTask@method@dataSplits)
    
    itsInfo <- vector("list",estTask@method@nReps)

    n <- NROW(eval(task@dataSource))

    if (!userSplit) {
        train.size <- if (estTask@method@szTrain < 1) as.integer(n*estTask@method@szTrain) else estTask@method@szTrain
        test.size <- if (estTask@method@szTest < 1) as.integer(n*estTask@method@szTest) else estTask@method@szTest
        if (n-test.size+1 <= train.size+1) stop('mcEstimates:: Invalid train/test sizes.',call.=FALSE)
    } else {
        train.size <- length(estTask@method@dataSplits[[1]]$train)
        test.size <- length(estTask@method@dataSplits[[1]]$test)
    }
  
    set.seed(estTask@method@seed)

    if (!userSplit) {
        selection.range <- (train.size+1):(n-test.size+1)
        starting.points <- sort(sample(selection.range,estTask@method@nReps))
    } else {
        starting.points <- sapply(estTask@method@dataSplits, function(d) d$test[1])
    }

    it <- NULL  # dummy assignment due to Note on cran-check
    itsInfo <- foreach::foreach(it=seq(along=starting.points),
                                .packages=.loadedPackages(),
                                .export=if (class(task@dataSource) == "data.frame") NULL else c(as.character(task@dataSource))
                                ) %fun%  {



        start <- starting.points[it]
        if (missing(cluster)) cat('Repetition ',it,'\n\t start test = ',
                                  start,'; test size = ',test.size,'\n')


        if (!userSplit) {
            rep.res <- runWorkflow(wf,
                                   task@formula,
                                   eval(task@dataSource)[(start-train.size):(start-1),],
                                   eval(task@dataSource)[start:(start+test.size-1),])
        } else {
            rep.res <- runWorkflow(wf,
                                   task@formula,
                                   eval(task@dataSource)[estTask@method@dataSplits[[it]]$train,],
                                   eval(task@dataSource)[estTask@method@dataSplits[[it]]$test,])

        }

        c(rep.res,list(train=(start-train.size):(start-1)))

    }
    cat('\n')

    ## randomize the number generator to avoid undesired
    ## problems caused by inner set.seed()'s
    set.seed(prod(as.integer(unlist(strsplit(strsplit(date()," ")[[1]][4],":")))))
    
    ## Calculate the metrics estimation
    scores <- .scoresIts(task,estTask,itsInfo)
    
    EstimationResults(task,wf,estTask,scores,itsInfo)
}




# =====================================================
# Small utility functions 
# =====================================================


is.regression <- function(task) task@type == 'regr'

is.classification <- function(task) task@type == 'class'

responseValues <- function(formula,data,na=NULL) model.response(model.frame(formula,data,na.action=na))


## ----------------------------------------
## Internal functions that are not exported

outFold <- function(ds,it,what="test") if (is.list(ds[[1]])) ds[[it]][[what]] else ds[[it]]

.scores2summary <- function(obj)
    apply(obj@iterationsScores,2,function(x)
          c(avg=mean(x,na.rm=TRUE),std=sd(x,na.rm=TRUE),
            med=median(x,na.rm=TRUE),iqr=IQR(x,na.rm=TRUE),
            min=min(x,na.rm=TRUE),max=max(x,na.rm=TRUE),
            invalid=sum(is.na(x)))
          )


.scores2long <- function(itRes) {
    d <- data.frame(rep=1:nrow(itRes),itRes)
    s <- reshape(d,direction='long',varying=list(2:(ncol(itRes)+1)),idvar='rep',v.names='score')
    colnames(s)[2] <- 'stat'
    s[,2] <- factor(s[,2],labels=colnames(d)[2:(ncol(itRes)+1)])
    s
}


.statScores <- function(compRes,stat=1) {
    r <- list()
    for(t in compRes) {
        ws <- NULL
        for(w in t)
            ws <- cbind(ws,t(.scores2summary(w)[stat,,drop=FALSE]))
        colnames(ws) <- names(t)
        r <- c(r,list(ws))
    }
    names(r) <- names(compRes)
    r
}

## Though simpler and more elegant this one fails due to over-simplification of
## sapply when we have only one metric (and it did not worked with simplify=FALSE
## on sapply)
## .statScores.old <- function(compRes,stat=1) lapply(compRes,function(t) sapply(t,function(w) .scores2summary(w)[stat,,drop=FALSE]))


## calculates the scores of all iterations of an estimation exp
.scoresIts <- function(task,estTask,its) {
    trReq <- estTask@trainReq || any(estTask@metrics %in% c("nmse","nmae","theil")) || (is.regression(task) && is.null(estTask@metrics))

    nIts <- length(its)

    standEval <- if (estTask@evaluator == "" ) TRUE else FALSE
    if (standEval) 
        evalFunc <- if (is.classification(task)) "classificationMetrics" else "regressionMetrics"
    else
        evalFunc <- estTask@evaluator

    scores <- NULL
    wts <- intersect(estTask@metrics,c("trTime","tsTime","totTime"))
    predMs <- setdiff(estTask@metrics,wts)
    stats <- if (!is.null(predMs)) list(stats=predMs) else NULL
    
    for(i in 1:nIts) {
        trR <- if (trReq) list(train.y=eval(task@dataSource)[its[[i]]$train,task@target]) else NULL
        fstArgs <- if (standEval) list(trues=its[[i]]$trues,preds=its[[i]]$preds) else its[[i]]
        ss <- do.call(evalFunc,
                      c(fstArgs,
                        stats,
                        trR,
                        estTask@evaluator.pars))
        if (is.null(scores)) {
            if (is.null(estTask@metrics)) {
                ncols <- length(ss)
                namcols <- names(ss)
            } else {
                ncols <- length(estTask@metrics)
                namcols <- estTask@metrics
            }
            scores <- matrix(NA,nrow=nIts,ncol=ncols,
                             dimnames=list(1:nIts,namcols))
        }
        if (is.null(predMs))  scores[i,] <- ss else scores[i,predMs] <- ss
        
        if (length(wts)) {
            allts <- as.numeric(its[[i]]$times)
            scores[i,wts] <- c(trTime=allts[1],tsTime=allts[2],totTime=allts[1]+allts[2])[wts]
        }
    }
    scores
}


.loadedPackages <- function(bases=c("datasets","grDevices","stats","utils","base","graphics","methods")) setdiff(sapply(strsplit(search()[grep("package",search())],":"),function(x) x[2]),bases)

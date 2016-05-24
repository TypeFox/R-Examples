################################################################# 
# THIS FILE CONTAINS FUNCTIONS THAT ARE RELATED TO RUNNING      #
# EXPERIMENTS WITH MODELLING TOOLS                              #
# IT IS A PART OF THE PACKAGE DMwR                              #
#################################################################
# Author : Luis Torgo (ltorgo@inescporto.pt)     Date: Jan 2009 #
# License: GPL (>= 2)                                           #
#################################################################




#################################################################
# GENERIC EXPERIMENTAL COMPARISONS
#################################################################

# =====================================================
# This function is intended for carrying out comparisons
# of a set of learning systems over a set of data sets.
# The type of experimental methodology that is used depends
# on the class of the experimental settings object (3rd arg.):
# if it is of class "cvSettings" then a Cross Validation
# experimental comparison is carried out; if it is of class
# "mcSettings" then a Monte Carlo comparison is done.
#
# The user provides the data sets (a list), the learning
# system (another list) and the experimental comparison
# settings in an object of the adequate class (see above).
#
# The result is a complex data structure with all
# information on the experiment.
#
# There are a few extra utility functions to work over
# the resulting data structure, e.g.
# bestScores(res)
# getVariant("cv.nnet-v2",res)
# =====================================================
# Luis Torgo, Jan-Aug 2009
# =====================================================
# Example calls:
# z <- experimentalComparison(
#       # data sets
#       list(dataset(medv~.,Boston,'b1'),
#         dataset(medv ~ lstat + rm + dis, Boston,'b2')
#         ),
#       # learners
#       list(sld5.se1=learner('slidingWindowTest',
#                 pars=list(learner=learner("MC.rpartXse",pars=list(se=1)),
#                   relearn.step=5,
#                   eval.func='erros')
#                 ),
#         sld10.se0=learner('slidingWindowTest',
#                 pars=list(learner=learner("MC.rpartXse",pars=list(se=0)),
#                   relearn.step=10,
#                   eval.func='erros')
#                 )
#         ),
#       # experiment settings (in this case a Monte Carlo experiment)
#       mcSettings(5,100,100,1234)
#                   )
#
# y <- experimentalComparison(
#       # data sets
#       c(dataset(medv~.,Boston,'b1'),
#         dataset(medv ~ lstat + rm + dis, Boston,'b2')
#         ),
#       # learners
#       c(variants('cv.rpartXse',se=c(0,0.5,1))),
#       # experiment settings (in this case a 3x10-fold Cross Validation)
#       cvSettings(3,10,1234)
#                   )
#
experimentalComparison <- function(datasets,systems,setts,...) {
  ##require(abind,quietly=T)

  if (!is(datasets,'list')) datasets <- list(datasets)
  if (!is(systems,'list')) systems <- list(systems)
  
  if (is.null(names(systems)))
    names(systems) <- paste('var',1:length(systems),sep='.')
  
  results <- compExp(systems,
                     lapply(datasets,function(x) as(x,'task')),
                     setts,array())

  r <- NULL
  
  cat('\n\n##### ',
      switch(class(setts),
             cvSettings='CROSS VALIDATION',
             hldSettings='HOLD OUT',
             mcSettings='MONTE CARLO',
             bootSettings='BOOTSTRAP',
             loocvSettings='LOOCV',
             ),
      ' EXPERIMENTAL COMPARISON #####')
  
  for(d in 1:length(datasets)) {

    cat('\n\n** DATASET ::',datasets[[d]]@name)

    rr <- NULL
    for (s in 1:length(systems)) {

      cat('\n\n++ LEARNER ::',systems[[s]]@func,
          ' variant -> ',names(systems)[s],'\n')

      var.res <- do.call(
               switch(class(setts),
                      cvSettings='crossValidation',
                      hldSettings='holdOut',
                      bootSettings='bootstrap',
                      mcSettings='monteCarlo',
                      loocvSettings='loocv'
                      ),
#               if(is(setts,'cvSettings')) 'crossValidation' else 'monteCarlo',
               c(list(systems[[s]],
                    datasets[[d]],
                    setts),...)
                         )
      
      rr <- abind(rr,var.res@foldResults,along=3)
      
    }

    r <- abind(r,rr,along=4)

  }
  results@foldResults <- r
  dimnames(results@foldResults)[3:4] <-
    list(names(systems),unlist(lapply(datasets,function(x) x@name)))

  results
}



#################################################################
# Cross Validation Experiments
#################################################################



# =====================================================
# Function that performs a cross validation experiment
# of a system on a given data set.
# The function is completely generic. The generality comes
# from the fact that the function that the user provides
# as the system to evaluate, needs in effect to be a
# user-defined function that takes care of the learning,
# testing and calculation of the statistics that the user
# wants to estimate through cross validation. 
# =====================================================
# Luis Torgo, Jan 2009
# =====================================================
# Example runs:
# x <- crossValidation(learner('cv.rpartXse',list(se=2)),
#                      dataset(medv~.,Boston),
#                      cvSettings(1,10,1234))
#
crossValidation <- function(sys,ds,sets,itsInfo=F) {

  show(sets)

  n <- nrow(ds@data)
  n.each.part <- n %/% sets@cvFolds

  itsI <- results <- NULL

  if (sets@strat) {  # stratified sampling
    respVals <- resp(ds@formula,ds@data)
    regrProb <- is.numeric(respVals)
    if (regrProb) {  # regression problem
      # the bucket to which each case belongs  
      b <- cut(respVals,10)  # this 10 should be parametrizable
    } else {
      b <- respVals
    }
    # how many on each bucket
    bc <- table(b)
    # how many should be on each test partition
    bct <- bc %/% sets@cvFolds
    # still missing (due to rounding effects of the previous statement)
    #rem <- n.test-sum(bct)
    #ib <- 1
    #nb <- length(bct)
    #while (rem) 
    #  if (bct[ib] < bc[ib]) {
    #    bct[ib] <- bct[ib]+1
    #    rem <- rem-1
    #    ib <- ib %% nb + 1
    #  }

  }
  
  for(r in 1:sets@cvReps) {
    cat('Repetition ',r,'\nFold:')

    set.seed(sets@cvSeed*r)
    permutation <- sample(n)
    perm.data <- ds@data[permutation,]

    for(i in seq(sets@cvFolds)) {
      cat(' ',i)
      
      if (sets@strat) {
          out.fold <- c()
          for(x in seq(along=levels(b))) 
              if (bct[x]) out.fold <- c(out.fold,which(b == levels(b)[x])[((i-1)*bct[x]+1):((i-1)*bct[x]+bct[x])])
      } else {
          out.fold <- ((i-1)*n.each.part+1):(i*n.each.part)
      }

      it.res <- runLearner(sys,
                           ds@formula,
                           perm.data[-out.fold,],
                           perm.data[out.fold,])

      
     if (itsInfo && !is.null(tmp <- attr(it.res,'itInfo'))) itsI <- c(itsI,tmp)
     results <- rbind(results,it.res)
      
    }
    cat('\n')
  }
  rownames(results) <- 1:nrow(results)
  colnames(results) <- names(it.res)
  
  # randomize the number generator to avoid undesired
  # problems caused by inner set.seed()'s
  set.seed(prod(as.integer(unlist(strsplit(strsplit(date()," ")[[1]][4],":")))))

  if (itsInfo) return(structure(cvRun(sys,as(ds,'task'),sets,results),itsInfo=itsI))
  else         return(cvRun(sys,as(ds,'task'),sets,results))

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
# x <- holdOut(learner('cv.rpartXse',list(se=2)),
#              dataset(medv~.,Boston),
#              hldSettings(4,0.25,1234))
#
holdOut <- function(sys,ds,sets,itsInfo=F) {

  show(sets)

  n <- nrow(ds@data)
  n.test <- as.integer(n * sets@hldSz)

  itsI <- results <- NULL
  
  if (sets@strat) {  # stratified sampling
    respVals <- resp(ds@formula,ds@data)
    regrProb <- is.numeric(respVals)
    if (regrProb) {  # regression problem
      # the bucket to which each case belongs  
      b <- cut(respVals,10)  # this 10 should be parameterizable
    } else {
      b <- respVals
    }
    # how many on each bucket
    bc <- table(b)
    # how many should be on each test partition
    bct <- as.integer(bc * sets@hldSz)
    # still missing (due to rounding effects of the previous statement)
    #rem <- n.test-sum(bct)
    #ib <- 1
    #nb <- length(bct)
    #while (rem) 
    #  if (bct[ib] < bc[ib]) {
    #    bct[ib] <- bct[ib]+1
    #    rem <- rem-1
    #    ib <- ib %% nb + 1
    #  }
  }
  
  for(r in 1:sets@hldReps) {
    cat('Repetition ',r)

    set.seed(sets@hldSeed*r)
    permutation <- sample(n)
    perm.data <- ds@data[permutation,]

    if (sets@strat) {
        out.fold <- c()
        for(x in seq(along=levels(b))) 
            if (bct[x]) out.fold <- c(out.fold,which(b == levels(b)[x])[1:bct[x]])
    } else {
        out.fold <- 1:n.test
    }

    it.res <- runLearner(sys,
                         ds@formula,
                         perm.data[-out.fold,],
                         perm.data[out.fold,])
    
    if (itsInfo && !is.null(tmp <- attr(it.res,'itInfo'))) itsI <- c(itsI,tmp)
    results <- rbind(results,it.res)
      
    cat('\n')
  }
  rownames(results) <- 1:nrow(results)
  colnames(results) <- names(it.res)
  
  # randomize the number generator to avoid undesired
  # problems caused by inner set.seed()'s
  set.seed(prod(as.integer(unlist(strsplit(strsplit(date()," ")[[1]][4],":")))))

  if (itsInfo) return(structure(hldRun(sys,as(ds,'task'),sets,results),itsInfo=itsI))
  else         return(hldRun(sys,as(ds,'task'),sets,results))
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
# x <- loocv(learner('cv.rpartXse',list(se=2)),
#            dataset(medv~.,Boston))
#
loocv <- function(sys,ds,sets,itsInfo=F,verbose=F) {

  show(sets)

  n <- nrow(ds@data)

  itsI <- results <- NULL

  if (verbose) cat('Iteration: ')
  for(r in 1:n) {
    if (verbose) cat('*')

    set.seed(sets@loocvSeed*r)

    it.res <- runLearner(sys,
                         ds@formula,
                         ds@data[-r,],
                         ds@data[r,])
    
    if (itsInfo && !is.null(tmp <- attr(it.res,'itInfo'))) itsI <- c(itsI,tmp)
    results <- rbind(results,it.res)
      
  }
  if (verbose) cat('\n')
  rownames(results) <- 1:nrow(results)
  colnames(results) <- names(it.res)
  

  if (itsInfo) return(structure(loocvRun(sys,as(ds,'task'),sets,results),itsInfo=itsI))
  else         return(loocvRun(sys,as(ds,'task'),sets,results))
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
# x <- bootstrap('cv.rpartXse',list(se=2)),
#                      dataset(medv~.,Boston),
#                      bootSettings(1234,10))
#
bootstrap <- function(sys,ds,sets,itsInfo=F,verbose=T) {

  show(sets)

  n <- nrow(ds@data)

  itsI <- results <- NULL

  if (verbose) cat('Repetition: ')
  for(r in 1:sets@bootReps) {
    if (verbose) cat(' ',r)

    set.seed(sets@bootSeed*r)
    idx <- sample(n,n,replace=T)
    
    it.res <- runLearner(sys,
                         ds@formula,
                         ds@data[idx,],
                         ds@data[-idx,])

    if (itsInfo && !is.null(tmp <- attr(it.res,'itInfo'))) itsI <- c(itsI,tmp)
    results <- rbind(results,it.res)
      
  }
  if (verbose) cat('\n')

  rownames(results) <- 1:nrow(results)
  colnames(results) <- names(it.res)
  
  # randomize the number generator to avoid undesired
  # problems caused by inner set.seed()'s
  set.seed(prod(as.integer(unlist(strsplit(strsplit(date()," ")[[1]][4],":")))))
  
  if (itsInfo) return(structure(bootRun(sys,as(ds,'task'),sets,results),itsInfo=itsI))
  else         return(bootRun(sys,as(ds,'task'),sets,results))

}



#################################################################
# Monte Carlo Experiments
#################################################################



# =====================================================
# Function that performs a Monte Carlo experiment of a 
# system on a given data set.
# The function is completely generic. The generality comes
# from the fact that the function that the user provides
# as the system to evaluate, needs in effect to be a
# user-defined function that takes care of the learning,
# testing and calculation of the statistics that the user
# wants to estimate through this experiment. A few example
# functions are provided.
# =====================================================
# Luis Torgo, Aug 2009
# =====================================================

monteCarlo <- function(learner,
                       data.set,
                       mcSet,
                       itsInfo=F,verbose=T) {

  show(mcSet)
  
  itsI <- results <- NULL

  n <- NROW(data.set@data)

  train.size <- if (mcSet@mcTrain < 1) as.integer(n*mcSet@mcTrain) else mcSet@mcTrain
  test.size <- if (mcSet@mcTest < 1) as.integer(n*mcSet@mcTest) else mcSet@mcTest
  if (n-test.size+1 <= train.size+1) stop('monteCarlo:: Invalid train/test sizes.')
  
  set.seed(mcSet@mcSeed)
  selection.range <- (train.size+1):(n-test.size+1)
  starting.points <- sort(sample(selection.range,mcSet@mcReps))


  # main loop over all repetitions
  for(it in seq(along=starting.points)) {
    start <- starting.points[it]

    if (verbose)  cat('Repetition ',it,'\n\t start test = ',
                      start,'; test size = ',test.size,'\n')

    itDS <- dataset(data.set@formula,
                    data.set@data[(start-train.size):(start+test.size-1),])
    
    rep.res <- runLearner(learner,
                          data.set@formula,
                          data.set@data[(start-train.size):(start-1),],
                          data.set@data[start:(start+test.size-1),])

    if (itsInfo && !is.null(tmp <- attr(rep.res,'itInfo'))) itsI <- c(itsI,tmp)
    results <- rbind(results,rep.res)

  }
  if (verbose) cat('\n')
  rownames(results) <- 1:nrow(results)

  # randomize the number generator to avoid undesired
  # problems caused by inner set.seed()'s
  set.seed(prod(as.integer(unlist(strsplit(strsplit(date()," ")[[1]][4],":")))))

  if (itsInfo) return(structure(mcRun(learner,as(data.set,'task'),mcSet,results),itsInfo=itsI))
  else         return(mcRun(learner,as(data.set,'task'),mcSet,results))


}





#################################################################
# Manipulation of the "compExps" objects
#################################################################


# =====================================================
# Function that joins  experimental results objects.
# The joining is carried out by some specified dimension,
# the most common being joining experiments carried out in
# different data sets (dimension 4), or experiments with
# different learners (dimension 3) on the same data sets.
# =====================================================
# Luis Torgo, Aug 2009
# =====================================================
# Example runs:
# > bestScores(join(subset(earth,stats='e1',vars=1:3),
#                   subset(nnet,stats='e1',vars=4:6),by=3))
# > bestScores(join(nnet,earth,rf,rpartXse,svm,by=3))
#
join <- function(...,by='datasets') {

  s <- list(...)
  if (length(s) < 2) return(s[1])
  
  if ((! by %in% c('iterations','statistics','variants','datasets')) &&
      (! by %in% 1:4))
    stop('join:: invalid value on "by" argument!')
  for(i in 2:length(s)) 
    if (!identical(s[[i]]@settings,s[[1]]@settings))
      stop('join:: trying to join experimental comparisons with different settings!')

  ##require(abind)
  if (!is.numeric(by))
    by <- match(by,c('iterations','statistics','variants','datasets'))

  r <- s[[1]]
  if (by == 3) for(i in 2:length(s)) r@learners <- c(r@learners,s[[i]]@learners)
  if (by == 4) for(i in 2:length(s)) r@datasets <- c(r@datasets,s[[i]]@datasets)
  for(i in 2:length(s))
    r@foldResults <- abind(r@foldResults,s[[i]]@foldResults,along=by)
  r
}

# =====================================================
# Small auxiliary functions to obtain information from 
# compExp objects.
# =====================================================
# Luis Torgo, Mar 2011
# =====================================================
dsNames      <- function(res) dimnames(res@foldResults)[[4]]
learnerNames <- function(res) dimnames(res@foldResults)[[3]]
statNames    <- function(res) dimnames(res@foldResults)[[2]]


# =====================================================
# This function receives a data structure resulting from
# the experimentalComparison() function and provides list with as
# many components as there are data sets and for each data
# set, a data frame with the best system and respective
# score for each statistic.
# By default the best is obtained as the minimum average
# score over all folds, but through a second parameter
# the user may indicate that some statistics are to be maximized
# and not minimized.
# =====================================================
# Luis Torgo, Jan-Aug 2009
# =====================================================
# Example calls:
# 1) a situation with only minimizations
# > bestScores(results)
#
# 2) a situation where the 2nd statistic is to be maximized
# > bestScores(results,c(F,T,F,F))
#
bestScores <- function(compRes,maxs=rep(F,dim(compRes@foldResults)[2])) {
  if (!inherits(compRes,'compExp')) stop(compRes,' needs to be of class "compExp".\n')
  if (length(maxs) == 1) maxs <- rep(maxs,dim(compRes@foldResults)[2])
  else if (length(maxs) != dim(compRes@foldResults)[2]) stop('"maxs" needs to have the same size as the number of evaluation statistics.\n')

  avgsDS <- apply(compRes@foldResults,4,
                  function(d) {
                    m <- matrix(apply(d,3, function(x) apply(x,2,mean,na.rm=T)),
                                nrow=dim(compRes@foldResults)[2])
                    dimnames(m)[2] <- dimnames(compRes@foldResults)[3]
                    data.frame(m)
                  }
                  )
  
  nms <- lapply(avgsDS,
                function(p)
                   sapply(1:nrow(p),
                          function(l,tab,f)
                             if (f[l]) which.max(tab[l,]) else which.min(tab[l,]),p,maxs))
                
  nms <- lapply(nms,function(x) names(x))

  scs <- lapply(avgsDS,
                function(p)
                   sapply(1:nrow(p),
                          function(l,tab,f)
                             if (f[l]) max(tab[l,],na.rm=T) else min(tab[l,],na.rm=T),p,maxs))

  
  
  res <- lapply(1:length(nms),
                function(x,p,s) data.frame(system=p[[x]],score=s[[x]],row.names=dimnames(compRes@foldResults)[[2]],stringsAsFactors=F),
                nms,scs)
  names(res) <- dimnames(compRes@foldResults)[[4]]
  res
    
}


# =====================================================
# This function receives a data structure resulting from
# the experimentalComparison() function and provides a list with as
# many components as there are data sets and for each data
# set, another list with as many components as there are
# evaluation statistics. For each component of this list a top X
# is given the the names and scores of the best X systems on
# that data set / statistic.
# =====================================================
# Luis Torgo, Nov 2009
# =====================================================
# Example calls:
#
rankSystems <- function(compRes,top=5,maxs=rep(F,dim(compRes@foldResults)[2])) {
  if (!inherits(compRes,'compExp')) stop(compRes,' needs to be of class "compExp".\n')
  if (length(maxs) == 1) maxs <- rep(maxs,dim(compRes@foldResults)[2])
  else if (length(maxs) != dim(compRes@foldResults)[2]) stop('"maxs" needs to have the same size as the number of evaluation statistics.\n')
  dss <- list()
  for(d in 1:dim(compRes@foldResults)[4]) {
    dss[[d]] <- list()
    for(s in 1:dim(compRes@foldResults)[2]) {
      m <- apply(compRes@foldResults[,s,,d],2,mean,na.rm=T)
      t <- m[order(m,decreasing=maxs[s])[1:top]]
      dss[[d]][[s]] <- data.frame(system=names(t),score=t,row.names=1:top,stringsAsFactors=F)
      names(dss[[d]])[s] <- dimnames(compRes@foldResults)[[2]][s]
    }

  }
  names(dss) <- dimnames(compRes@foldResults)[[4]]
  dss
}
      
# =====================================================
# This function receives a data structure resulting from
# the experimentalComparison() function and provides a list with as
# many components as there are data sets and for each data
# set, the scores of all systems on a single statistic. The scores are
# by default the mean but any other numeric summarization function
# can be provided in the argument summary.
# =====================================================
# Luis Torgo, Nov 2009
# =====================================================
# Example calls:
#
# > statScores(svmR,'NTrades')
# > statScores(svmR,'NTrades','max')
# 
statScores <- function(compRes,stat,summary='mean') {
  if (!inherits(compRes,'compExp'))
    stop(compRes,' needs to be of class "compExp".\n')
  if (length(stat)>1)
    stop('"stat" should contain a single statistic.\n')
  
  dss <- list()
  for(d in 1:dim(compRes@foldResults)[4]) 
    dss[[d]]  <- apply(compRes@foldResults[,stat,,d],2,summary,na.rm=T)

  names(dss) <- dimnames(compRes@foldResults)[[4]]
  dss
}



# ======================================================================
# Construction of comparative analysis tables based on the results of 
# comparative experiment obtained with experimentalComparison() (i.e. based
# on a compExp object).
# The first argument is the compExp object that resulted from the 
# experiments. Then we have the system against which all remaning are
# compared to (defaults to first in the structure). Finally we can
# provide a vector of the names of the statistics we which to get a
# table (defaults to all).
# =====================================================
# Luis Torgo, Jan-Aug 2009
# =====================================================
compAnalysis <- function(comp,
                       against=dimnames(comp@foldResults)[[3]][1],
                       stats=dimnames(comp@foldResults)[[2]],
                       datasets=dimnames(comp@foldResults)[[4]],
                       show=T) {
  if (!inherits(comp,'compExp')) stop(comp,' is not of class compExp.\n')

  d <- dim(comp@foldResults)
  dn <- dimnames(comp@foldResults)
  
  n.ds <- length(datasets)
  
  n.methods <- d[3]
  if (n.methods < 2)
    stop('More than one method is necessary for significance comparisons.\n')
  methods <- dn[[3]]

  ag <- which(methods==against)
  theOne <- methods[ag]
  theRest <- methods[-ag]
  
  sig.res <- c('  ','+ ','++','- ','--')
  
  res <- list()
  for(s in stats) {
    for(d in seq(along=datasets)) {
      la <- ls <- list()

      la[[theOne]] <- getSummaryResults(comp,theOne,datasets[d])[1,s] # 1= avg
      ls[[theOne]] <- getSummaryResults(comp,theOne,datasets[d])[2,s] # 2= std

      for(m in theRest) {
        
        la[[m]] <- getSummaryResults(comp,m,datasets[d])[1,s]
        ls[[m]] <- getSummaryResults(comp,m,datasets[d])[2,s]


        t <- try(w <- wilcox.test(comp@foldResults[,s,m,d],
                                  comp@foldResults[,s,theOne,d],
                                  paired=T)
                 )
        if (inherits(t,"try-error")) i <- NA
        else if (w$p.value  > 0.05)  i <- '  '
        else if (w$p.value > 0.001)  i <- if (la[[m]] > la[[theOne]]) '+ ' else '- '
        else i <- if (la[[m]] > la[[theOne]]) '++' else '--'
        
        la[[paste('sig.',m,sep='')]] <- factor(i,levels=sig.res)
        ls[[paste('sig.',m,sep='')]] <- factor('  ',levels=sig.res)
      }
      
      if (exists('temp')) {temp[2*d-1,] <- la; temp[2*d,] <- ls}
      else {temp <- data.frame(la); temp[2,] <- ls}
    }
    row.names(temp)[seq(1,by=2,len=n.ds)] <- datasets
    res[[s]] <- temp
    rm('temp')
  }
  if (show) {
    cat('\n== Statistical Significance Analysis of Comparison Results ==\n')
    cat('\nBaseline Learner::\t',against,' (Learn.1)\n')
    for(s in stats) {
      cat('\n** Evaluation Metric::\t',s,'\n')
      for(d in seq(along=datasets)) showTab(datasets[d],res[[s]][(d-1)*2+1:2,])
    }
    cat('\nLegends:\nLearners -> ')
#    cat('\t|')
    for(s in seq(along=c(theOne,theRest))) cat(paste('Learn',s,sep='.'),'=',c(theOne,theRest)[s],'; ')
    cat("\nSignif. Codes -> 0 '++' or '--' 0.001 '+' or '-' 0.05 ' ' 1\n")

    return(if (length(stats)==1) invisible(res[[1]]) else invisible(res))
  } else  if (length(stats)==1) res[[1]] else res
}

showTab <- function(d,x,dig=3) {
  rownames(x) <- c('AVG','STD')
  colnames(x)[1] <- c('Learn.1')
  for(i in 2*(1:(ncol(x)%/%2))) colnames(x)[c(i,i+1)] <- c(paste('Learn',i/2+1,sep='.'),paste('sig',i/2+1,sep='.'))
  cat('\n- Dataset:',d,'\n')
  print(x)
}
  
# ======================================================================
# Get the fold results of a certain variant (learning system) over a
# certain dataset. Both can be specified by name or number.
# You get a matrix with as many rows as there are folds and as many
# columns as there are evaluation statistics.
# =====================================================
# Luis Torgo, Jan 2009
# =====================================================
getFoldsResults <- function(results,learner,dataSet) {
  if (!inherits(results,'compExp')) stop(results,' is not of class compExp.\n')

  results@foldResults[,,learner,dataSet]
}


# ======================================================================
# Get some summary statistics of the performance of a variant on a certain
# data set, for all evaluation statistics.
# =====================================================
# Luis Torgo, Jan 2009
# =====================================================
getSummaryResults <- function(results,learner,dataSet) {
  if (!inherits(results,'compExp')) stop(results,' is not of class compExp.\n')

  apply(results@foldResults[,,learner,dataSet,drop=F],2,function(x)
        c(avg=mean(x,na.rm=T),std=sd(x,na.rm=T),
          min=if (all(is.na(x))) NA else min(x,na.rm=T),
          max=if (all(is.na(x))) NA else max(x,na.rm=T),
          invalid=sum(is.na(x)))
        )
}
  



#################################################################
# Learners and Variants of Learners
#################################################################


# =====================================================
# This function generates a named list of objects of class
# learner.
# It is used for easily generating a list of variants of
# a learning system that is usable by the experimentalComparison()
# function.
# If you give only the learning system name it will generate
# a learner wit the default parameters of that system.
# The names of the components are generated automatically
# and have the form <sys>-v<x> where <x> is an increasing
# integer. In the case of using defaults the name is
# <sys>-defaults
# =====================================================
# Luis Torgo, Jan 2009
# =====================================================
# Example call:
# ex1 <- variants('cv.rpartXse',se=c(0,0.5,1))
# ex2 <- variants('nnet')
#
variants <- function(sys,varsRootName=sys,as.is=NULL,...) {

    vars <- list(...)
    if (!length(vars)) {
        vars <- c(learner(sys,list()))
        names(vars)[1] <- paste(varsRootName,'.v1',sep='')
        return(vars)
    }
  
    ## the parameters not involved in variants generation
    ## their names:
    allnovar <- c(as.is,names(vars)[which(sapply(vars,length)==1)])
    ## their positions in vars:
    toExcl <- which(names(vars) %in% allnovar)
    ## their number:
    nExcl <- length(toExcl)
    
    ## checking how many variants per parameter and generate the grid
    nvarsEach <- rep(1,length(vars))
    varying <- if (nExcl) (1:length(vars))[-toExcl] else 1:length(vars)
    if (length(varying)) nvarsEach[varying] <- lapply(vars[varying],length)
    idxsEach <- lapply(nvarsEach,function(x) 1:x)
    theVars <- expand.grid(idxsEach)
    
    ## now go for generating the different variants
    vs <- list()
    for(i in 1:nrow(theVars)) {
        ## start
        varPars <- list()
        for(k in 1:ncol(theVars)) {
            if (nExcl & (k %in% toExcl)) varPars <- c(varPars,vars[k])
            else varPars <- c(varPars,vars[[k]][theVars[i,k]])
        }
        names(varPars) <- names(vars)
        vs <- c(vs,learner(sys,varPars))
    }
    names(vs) <- paste(varsRootName,'.v',1:nrow(theVars),sep='')
    vs
}



# =====================================================
# This function obtains the parameter settings associated
# to a certain variant name in the context of the variants
# of an experimental comparison
# =====================================================
# Luis Torgo, Jan 2009
# =====================================================
# Example Call:
# > getVariant('cv.nnet-v6',cvResults)
#
# Note: The result of this can then be "run" as follows,
# > runLearner(getVariant('cv.nnet-v6',cvResults),
#              medv~.,Boston[1:100,],Boston[-(1:100),])
#
getVariant <- function(var,ExpsData) 
  ExpsData@learners[[which(names(ExpsData@learners) == var)]]



# =====================================================
# Function that can be used to call a learning system
# whose information is stored in an object of class learner.
# =====================================================
# Luis Torgo, Fev 2009
# =====================================================
# Example run:
# l  <- learner('nnet',pars=list(size=4,linout=T))
# runLearner(l,medv ~ ., Boston)
#
runLearner <- function(l,...) {
  if (!inherits(l,'learner')) stop(l,' is not of class "learner".')
  do.call(l@func,c(list(...),l@pars))
}






#################################################################
# Sliding and Growing Windows Approaches to Learn+Test a Model
#################################################################

slidingWindowTest <- function(learner,
                               form,train,test,
                               relearn.step=1,verbose=T) {

  data <- rbind(train,test)
  n <- NROW(data)
  train.size <- NROW(train)
  sts <- seq(train.size+1,n,by=relearn.step)

  preds <- vector()
  for(s in sts) {

    if (verbose) cat('*')

    ps <- runLearner(learner,
                     form=form,
                     train=data[(s-train.size):(s-1),],
                     test=data[s:min((s+relearn.step-1),n),]
                     )

    preds <- c(preds,ps)
  }

  
  if (verbose) cat('\n')
  if (is.factor(resp(form,train))) return(factor(preds,levels=1:3,labels=levels(resp(form,train)))) else return(preds)
}


growingWindowTest <- function(learner,
                               form,train,test,
                               relearn.step=1,verbose=T) {

  data <- rbind(train,test)
  n <- NROW(data)
  train.size <- NROW(train)
  sts <- seq(train.size+1,n,by=relearn.step)

  preds <- vector()
  for(s in sts) {

    if (verbose) cat('*')

    ps <- runLearner(learner,
                     form=form,
                     train=data[1:(s-1),],
                     test=data[s:min((s+relearn.step-1),n),]
                     )

    preds <- c(preds,ps)
  }

  
  if (verbose) cat('\n')
  if (is.factor(resp(form,train))) return(factor(preds,levels=1:3,labels=levels(resp(form,train)))) else return(preds)

}




#################################################################
## Functions calculating evaluation metrics
#################################################################


# =====================================================================
# Function to calculate some standard regression evaluation statistics
# ---------------------------------------------------------------------
# L. Torgo (2009)
#
# Examples:
# s <- regr.eval(tr,ps,train.y=data[,'Y'])
# s <- regr.eval(tr,ps,stats=c('mse','mae'))
#
regr.eval <- function(trues,preds,
                      stats=if (is.null(train.y)) c('mae','mse','rmse','mape') else c('mae','mse','rmse','mape','nmse','nmae'),
                      train.y=NULL)
{
  allSs <- c('mae','mse','rmse','mape','nmse','nmae')
  if (any(c('nmse','nmad') %in% stats) && is.null(train.y))
    stop('regr.eval:: train.y parameter not specified.',call.=F)
  if (!all(stats %in% allSs))
    stop("regr.eval:: don't know how to calculate -> ",call.=F,
         paste(stats[which(!(stats %in% allSs))],collapse=','))
  N <- length(trues)
  sae <- sum(abs(trues-preds))
  sse <- sum((trues-preds)^2)
  r <- c(mae=sae/N,mse=sse/N,rmse=sqrt(sse/N),mape=sum(abs((trues-preds)/trues))/N)
  if (!is.null(train.y)) r <- c(r,c(nmse=sse/sum((trues-mean(train.y))^2),nmae=sae/sum(abs(trues-mean(train.y)))))
  return(r[stats])
}

# =====================================================================
# Function to calculate some standard classification evaluation statistics
# ---------------------------------------------------------------------
# L. Torgo (2012)
#
# Examples:
# s <- class.eval(tr,ps)
# s <- class.eval(tr,ps,benMtrx=matrix(c(2,-13,-4,5),2,2))
#
class.eval <- function(trues,preds,
                       stats=if (is.null(benMtrx)) c('err') else c('err','totU'),
                       benMtrx=NULL,
                       allCls=levels(factor(trues)))

  {
    preds <- factor(preds,levels=allCls)
    trues <- factor(trues,levels=allCls)
    allSs <- c('acc','err','totU')
    if (any(c('totU') %in% stats) && is.null(benMtrx))
      stop('class.eval:: benMtrx parameter not specified.',call.=F)
    if (!all(stats %in% allSs))
      stop("class.eval:: don't know how to calculate -> ",call.=F,
           paste(stats[which(!(stats %in% allSs))],collapse=','))
    N <- length(trues)
    cm <- as.matrix(table(trues,preds))
    a <- sum(diag(cm))/N
    r <- c(acc=a,err=1-a)
    if (!is.null(benMtrx))
      if (!all(dim(cm)==dim(benMtrx)))
        stop("class.eval:: dimensions of confusion and benefits metrices do not match",call.=F)
      else r <- c(r,totU=sum(cm*benMtrx))
    
    return(r[stats])
  }   


# =====================================================================
# Function to calculate some standard  evaluation statistics for time series
# problems
# ---------------------------------------------------------------------
# L. Torgo (2013)
#
# Examples:
# s <- ts.eval(tr,ps,train.y=data[,'Y'])
# s <- ts.eval(tr,ps,stats=c('mse','mae'))
#
ts.eval <- function(trues,preds,
                    stats=if (is.null(train.y)) c('mae','mse','rmse','mape') else c('mae','mse','rmse','mape','nmse','nmae','theil'),
                    train.y=NULL)
{
  print(stats)
  cat(length(trues),'\t',length(preds),'\n')
  r <- if (!is.null(train.y))  c(regr.eval(trues,preds,setdiff(stats,'theil'),train.y),theil=sum((trues-preds)^2)/sum((c(train.y[length(train.y)],trues[-length(trues)])-preds)^2)) else regr.eval(trues,preds,setdiff(stats,'theil'),train.y)
  return(r[stats])
}





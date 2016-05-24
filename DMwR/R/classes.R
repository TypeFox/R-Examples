################################################################# 
# THIS FILE CONTAINS THE CLASSES AND RESPECTIVE METHODS OF THE  #
# PACKAGE DMwR                                                  #
#################################################################
# Author : Luis Torgo (ltorgo@dcc.fc.up.pt)     Date: Jan 2009 #
# License: GPL (>= 2)                                           #
#################################################################

# Defined classes :
#   learner, task, dataset,
#   cvSettings, cvRun,
#   mcSettings, mcRun,
#   hldSettings, hldRun
#   loocvSettings, loocvRun
#   bootSettings, bootRun
#   expSettings, compExp,
#   tradeRecord



# ==============================================================
# CLASS: learner
# ==============================================================
# Luis Torgo, Jan 2009
# ==============================================================


# --------------------------------------------------------------
# class def
setClass("learner", representation(func="character",pars="list"))


# --------------------------------------------------------------
# constructor function
learner <- function(func,pars=list()) {
  if (missing(func)) stop("\nYou need to provide a function name.\n")
  new("learner",func=func,pars=pars)
}



# --------------------------------------------------------------
# Methods:


# show
setMethod("show","learner",
          function(object) {
            cat('\nLearner:: ',deparse(object@func),'\n\nParameter values\n')
            for(n in names(object@pars))
              cat('\t',n,' = ',deparse(object@pars[[n]]),'\n')
            cat('\n\n')
          }
          )




# ==============================================================
# CLASS: task
# ==============================================================
# Luis Torgo, Jan 2009
# ==============================================================


# --------------------------------------------------------------
# class def
setClass("task", representation(name="character",formula="formula"))


# --------------------------------------------------------------
# constructor function
task <- function(formula,data,name=NULL) {
  if (missing(formula) || missing(data))
    stop('\nYou need to provide a formula and a data frame.\n')
  if (inherits(try(model.frame(formula,data),T),"try-error"))
    stop('\nInvalid formula for the given data frame.\n')
  if (is.null(name)) {
    m <- match.call()
    name <- m$data
  }
  new("task",name=name,formula=formula)
}


# --------------------------------------------------------------
# Methods:

# show
setMethod("show","task",
          function(object) {
            cat('\nTask Name :: ',object@name)
            cat('\nFormula   :: ')
            print(object@formula)
            cat('\n')
          }
          )






# ==============================================================
# CLASS: dataset
# ==============================================================
# Luis Torgo, Jan 2009
# ==============================================================


# --------------------------------------------------------------
# class def
setClass("dataset", representation(data="data.frame"),contains='task')


# --------------------------------------------------------------
# constructor function
dataset <- function(formula,data,name=NULL) {
  if (missing(formula) || missing(data))
    stop('\nYou need to provide a formula and a data frame.\n')
  if (inherits(try(model.frame(formula,data),T),"try-error"))
    stop('\nInvalid formula for the given data frame.\n')

  if (is.null(name)) {
    m <- match.call()
    name <- deparse(m$data)
  }
  new("dataset",name=name,
      formula=formula,data=as.data.frame(model.frame(formula,data)))
}




# --------------------------------------------------------------
# Methods:

# show

setMethod("show","dataset",
          function(object) {
            cat('\nTask Name :: ',object@name)
            cat('\nFormula   :: ')
            print(object@formula)
            cat('Task Data ::\n\n')
            str(object@data,give.attr=F)
            cat('\n')
          }
          )





# ==============================================================
# CLASS: cvSettings
# ==============================================================
# Luis Torgo, Jan 2009
# ==============================================================


# --------------------------------------------------------------
# class def
setClass("cvSettings",
         representation(cvReps='numeric',  # nr. of repetitions
                        cvFolds='numeric', # nr. of folds of each rep.
                        cvSeed='numeric',  # seed of the random generator
                        strat='logical')   # is the sampling stratified?
         )


# --------------------------------------------------------------
# constructor function
cvSettings <- function(r=1,f=10,s=1234,st=F) 
  new("cvSettings",cvReps=r,cvFolds=f,cvSeed=s,strat=st)



# --------------------------------------------------------------
# Methods:


# show
setMethod("show","cvSettings",
          function(object) {
           cat(ifelse(object@strat,'\n Stratified ','\n'),
               object@cvReps,'x',object@cvFolds,
               '- Fold Cross Validation run with seed = ',
               object@cvSeed,'\n')
         })



# ==============================================================
# CLASS: cvRun
# ==============================================================
# Luis Torgo, Jan 2009
# ==============================================================


# --------------------------------------------------------------
# class def
setClass("cvRun",
         representation(learner="learner",
                        dataset="task",
                        settings="cvSettings",
                        foldResults="matrix")
         )


# --------------------------------------------------------------
# constructor function
cvRun <- function(l,t,s,r) {
  o <- new("cvRun")
  o@learner <- l
  o@dataset <- t
  o@settings <- s
  o@foldResults <- r
  o
}


# --------------------------------------------------------------
# Methods:

summary.cvRun <- function(object,...) {
  cat('\n== Summary of a Cross Validation Experiment ==\n')
  print(object@settings)
  cat('\n* Data set :: ',object@dataset@name)
  cat('\n* Learner  :: ',object@learner@func,' with parameters ')
  for(x in names(object@learner@pars))
    cat(x,' = ',
        object@learner@pars[[x]],' ')
  cat('\n\n* Summary of Experiment Results:\n\n')
  apply(object@foldResults,2,function(x)
        c(avg=mean(x,na.rm=T),std=sd(x,na.rm=T),
          min=min(x,na.rm=T),max=max(x,na.rm=T),
          invalid=sum(is.na(x)))
        )
}

plot.cvRun <- function(x,y,...) {
  sum <- apply(x@foldResults,2,function(x)
               c(avg=mean(x,na.rm=T),std=sd(x,na.rm=T),
                 min=min(x,na.rm=T),max=max(x,na.rm=T),
                 invalid=sum(is.na(x)))
               )
  nstats <- ncol(sum)
  layout(matrix(c(1:nstats,nstats,nstats+1:nstats,2*nstats),ncol=2),
         widths=c(3,1),heights=c(rep(5,nstats),1))
  
  par(mar=c(0,4.1,0,0))
  for(s in 1:(nstats-1)) {
    plot(x@foldResults[,s],type='b',xlab='',ylab=colnames(x@foldResults)[s],
         main='',xaxt='n')
    t <- axTicks(2)
    abline(h=t,lty=3,col='gray')
    abline(h=sum[1,s],lty=2)
  }
  par(mar=c(4.1,4.1,0,0))
  plot(x@foldResults[,nstats],type='b',xlab='Folds',
       ylab=colnames(x@foldResults)[nstats], main='')
  t <- axTicks(2)
  abline(h=t,lty=3,col='gray')
  abline(h=sum[1,nstats],lty=2)
  
  par(mar=c(0,0,0,0))
  for(s in 1:(nstats-1)) {
    boxplot(x@foldResults[,s],type='b',xlab='',ylab='',
            main='',xaxt='n',yaxt='n')
    abline(h=sum[1,s],lty=2)
    t <- axTicks(2)
    abline(h=t,lty=3,col='gray')
    
  }
  
  par(mar=c(4.1,0,0,0))
  boxplot(x@foldResults[,nstats],type='b',xlab='',
          ylab='', main='',xaxt='n',yaxt='n')
  abline(h=sum[1,nstats],lty=2)
  t <- axTicks(2)
  abline(h=t,lty=3,col='gray')
  mtext(paste('DATASET:',x@dataset@name),1,line=1,cex=0.8)
  mtext(paste('LEARNER:',x@learner@func),1,line=2,cex=0.8)
  
}




# ==============================================================
# CLASS: hldSettings
# ==============================================================
# Luis Torgo, Feb 2010
# ==============================================================


# --------------------------------------------------------------
# class def
setClass("hldSettings",
         representation(hldReps='numeric', # number of repetitions
                        hldSz='numeric',   # the size (0..1) of the holdout
                        hldSeed='numeric', # the random number seed
                        strat='logical')   # is the sampling stratified?
         )


# --------------------------------------------------------------
# constructor function
hldSettings <- function(r=1,sz=0.3,s=1234,str=F) 
  new("hldSettings",hldReps=r,hldSz=sz,hldSeed=s,strat=str)



# --------------------------------------------------------------
# Methods:


# show
setMethod("show","hldSettings",
          function(object) {
           cat(ifelse(object@strat,'\n Stratified ','\n'),
               object@hldReps,'x',
               100*(1-object@hldSz),'%/',100*object@hldSz,
               '% Holdout run with seed = ',
               object@hldSeed,'\n')
         })


# ==============================================================
# CLASS: hldRun
# ==============================================================
# Luis Torgo, Feb 2010
# ==============================================================


# --------------------------------------------------------------
# class def
#
# Note: foldResults is a matrix with as many columns as there are
# evaluation metrics being tested, and as many rows as there are
# iterations of the experiment (in this case the number of repetitions
# of the hold out experiment.
setClass("hldRun",
         representation(learner="learner",      # the learner
                        dataset="task",         # the data set
                        settings="hldSettings", # the settings of the exp.
                        foldResults="matrix")   # the results
         )


# --------------------------------------------------------------
# constructor function
hldRun <- function(l,t,s,r) {
  o <- new("hldRun")
  o@learner <- l
  o@dataset <- t
  o@settings <- s
  o@foldResults <- r
  o
}


# --------------------------------------------------------------
# Methods:

summary.hldRun <- function(object,...) {
  cat('\n== Summary of a Hold Out Experiment ==\n')
  print(object@settings)
  cat('\n* Data set :: ',object@dataset@name)
  cat('\n* Learner  :: ',object@learner@func,' with parameters:')
  for(x in names(object@learner@pars)) {
    k <- object@learner@pars[[x]]
    k <- paste(k,collapse=' ')
    k <- if (nchar(k) > 5) paste(substr(k,1,5),' ...') else k
    cat('\n\t',x,' = ',k,' ')
  }
  cat('\n\n* Summary of Experiment Results:\n\n')
  apply(object@foldResults,2,function(x)
        c(avg=mean(x,na.rm=T),std=sd(x,na.rm=T),
          min=min(x,na.rm=T),max=max(x,na.rm=T),
          invalid=sum(is.na(x)))
        )     
}

plot.hldRun <- function(x,y,...) {
  sum <- apply(x@foldResults,2,function(x)
               c(avg=mean(x,na.rm=T),std=sd(x,na.rm=T),
                 min=min(x,na.rm=T),max=max(x,na.rm=T),
                 invalid=sum(is.na(x)))
               )
  nstats <- ncol(sum)
  layout(matrix(c(1:nstats,nstats,nstats+1:nstats,2*nstats),ncol=2),
         widths=c(3,1),heights=c(rep(5,nstats),1))
  
  par(mar=c(0,4.1,0,0))
  for(s in 1:(nstats-1)) {
    plot(x@foldResults[,s],type='b',xlab='',ylab=colnames(x@foldResults)[s],
         main='',xaxt='n')
    t <- axTicks(2)
    abline(h=t,lty=3,col='gray')
    abline(h=sum[1,s],lty=2)
  }
  par(mar=c(4.1,4.1,0,0))
  plot(x@foldResults[,nstats],type='b',xlab='Repetitions',
       ylab=colnames(x@foldResults)[nstats], main='')
  t <- axTicks(2)
  abline(h=t,lty=3,col='gray')
  abline(h=sum[1,nstats],lty=2)
  
  par(mar=c(0,0,0,0))
  for(s in 1:(nstats-1)) {
    boxplot(x@foldResults[,s],type='b',xlab='',ylab='',
            main='',xaxt='n',yaxt='n')
    abline(h=sum[1,s],lty=2)
    t <- axTicks(2)
    abline(h=t,lty=3,col='gray')
    
  }
  
  par(mar=c(4.1,0,0,0))
  boxplot(x@foldResults[,nstats],type='b',xlab='',
          ylab='', main='',xaxt='n',yaxt='n')
  abline(h=sum[1,nstats],lty=2)
  t <- axTicks(2)
  abline(h=t,lty=3,col='gray')
  mtext(paste('DATASET:',x@dataset@name),1,line=1,cex=0.8)
  mtext(paste('LEARNER:',x@learner@func),1,line=2,cex=0.8)
  
}







# ==============================================================
# CLASS: loocvSettings
# ==============================================================
# Luis Torgo, Mar 2010
# ==============================================================


# --------------------------------------------------------------
# class def
setClass("loocvSettings",
         representation(loocvSeed='numeric',  # seed of the random generator
                        verbose='logical') # function used to evalute preds
         )


# --------------------------------------------------------------
# constructor function
loocvSettings <- function(seed=1234,verbose=F)
  new("loocvSettings",loocvSeed=seed,verbose=verbose)



# --------------------------------------------------------------
# Methods:


# show
setMethod("show","loocvSettings",
          function(object) {
           cat('\n LOOCV experiment with verbose = ',
               ifelse(object@verbose,'TRUE','FALSE'),' and seed =',
               object@loocvSeed,'\n')
         })


# ==============================================================
# CLASS: loocvRun
# ==============================================================
# Luis Torgo, Mar 2010
# ==============================================================


# --------------------------------------------------------------
# class def
#
# Note: casePreds is a matrix with as many columns as there are
# evaluation metrics being tested, and as many rows as there are
# test cases (i.e. the number of rows of the dataset)
setClass("loocvRun",
         representation(learner="learner",      # the learner
                        dataset="task",         # the data set
                        settings="loocvSettings", # the settings of the exp.
                        foldResults="matrix")   # the results
         )


# --------------------------------------------------------------
# constructor function
loocvRun <- function(l,t,s,r) {
  o <- new("loocvRun")
  o@learner <- l
  o@dataset <- t
  o@settings <- s
  o@foldResults <- r
  o
}


# --------------------------------------------------------------
# Methods:

summary.loocvRun <- function(object,...) {
  cat('\n== Summary of a Leave One Out Cross Validation  Experiment ==\n')
  print(object@settings)
  cat('\n* Data set :: ',object@dataset@name)
  cat('\n* Learner  :: ',object@learner@func,' with parameters:')
  for(x in names(object@learner@pars)) {
    k <- object@learner@pars[[x]]
    k <- paste(k,collapse=' ')
    k <- if (nchar(k) > 5) paste(substr(k,1,5),' ...') else k
    cat('\n\t',x,' = ',k,' ')
  }
  cat('\n\n* Summary of Experiment Results:\n\n')
  apply(object@foldResults,2,function(x)
        c(avg=mean(x,na.rm=T),std=sd(x,na.rm=T),
          min=min(x,na.rm=T),max=max(x,na.rm=T),
          invalid=sum(is.na(x)))
        ) 
}






# ==============================================================
# CLASS: bootSettings
# ==============================================================
# Luis Torgo, Apr 2010
# ==============================================================


# --------------------------------------------------------------
# class def
setClass("bootSettings",
         representation(bootSeed='numeric',  # seed of the random generator
                        bootReps='numeric') # number of repetitions
         )


# --------------------------------------------------------------
# constructor function
bootSettings <- function(seed=1234,nreps=50)
  new("bootSettings",bootSeed=seed,bootReps=nreps)



# --------------------------------------------------------------
# Methods:


# show
setMethod("show","bootSettings",
          function(object) {
           cat('\n Bootstrap experiment settings\n\t Seed = ',
               object@bootSeed,'\n\t Nr. repetitions = ',object@bootReps,'\n')
         })


# ==============================================================
# CLASS: bootRun
# ==============================================================
# Luis Torgo, Apr 2010
# ==============================================================


# --------------------------------------------------------------
# class def
#
setClass("bootRun",
         representation(learner="learner",      # the learner
                        dataset="task",         # the data set
                        settings="bootSettings", # the settings of the exp.
                        foldResults="matrix")   # the results
         )


# --------------------------------------------------------------
# constructor function
bootRun <- function(l,t,s,r) {
  o <- new("bootRun")
  o@learner <- l
  o@dataset <- t
  o@settings <- s
  o@foldResults <- r
  o
}


# --------------------------------------------------------------
# Methods:

summary.bootRun <- function(object,...) {
  cat('\n== Summary of a Bootstrap Experiment ==\n')
  print(object@settings)
  cat('\n* Data set :: ',object@dataset@name)
  cat('\n* Learner  :: ',object@learner@func,' with parameters:')
  for(x in names(object@learner@pars)) {
    k <- object@learner@pars[[x]]
    k <- paste(k,collapse=' ')
    k <- if (nchar(k) > 5) paste(substr(k,1,5),' ...') else k
    cat('\n\t',x,' = ',k,' ')
  }
  cat('\n\n* Summary of Experiment Results:\n\n')
  apply(object@foldResults,2,function(x)
        c(avg=mean(x,na.rm=T),std=sd(x,na.rm=T),
          min=min(x,na.rm=T),max=max(x,na.rm=T),
          invalid=sum(is.na(x)))
        )
}






# ==============================================================
# CLASS: mcSettings
# ==============================================================
# Luis Torgo, Aug 2009
# ==============================================================


# --------------------------------------------------------------
# class def
setClass("mcSettings",
         representation(mcReps='numeric',
                        mcTrain='numeric',mcTest='numeric',
                        mcSeed='numeric')
         )


# --------------------------------------------------------------
# constructor function
mcSettings <- function(r=10,tr=0.25,ts=0.25,s=1234)
  new("mcSettings",mcReps=r,mcTrain=tr,mcTest=ts,mcSeed=s)



# --------------------------------------------------------------
# Methods:


# show
setMethod("show","mcSettings",
          function(object) {
           cat('\n',object@mcReps,
               ' repetitions Monte Carlo Simulation using:',
               '\n\t seed = ', object@mcSeed,
               '\n\t train size = ',object@mcTrain,
               ifelse(object@mcTrain<1,'x NROW(DataSet)',' cases'),
               '\n\t test size = ',object@mcTest,
               ifelse(object@mcTest<1,'x NROW(DataSet)',' cases'),
               '\n'
               )
         })


          
# ==============================================================
# CLASS: mcRun
# ==============================================================
# Luis Torgo, Aug 2009
# ==============================================================
# Note: I'm sure this class and its methods would be generalizable
# together with the CV counter-parts to a single more generic class
# and respective methods. However, this is being made incrementally
# and I'm just felling too lazy and too pressed by deadlines to give
# this task some time... probably in a future version of the package...


# --------------------------------------------------------------
# class def
setClass("mcRun",
         representation(learner="learner",
                        dataset="task",
                        settings="mcSettings",
                        foldResults="matrix")
         )


# --------------------------------------------------------------
# constructor function
mcRun <- function(l,t,s,r) {
  o <- new("mcRun")
  o@learner <- l
  o@dataset <- t
  o@settings <- s
  o@foldResults <- r
  o
}


# --------------------------------------------------------------
# Methods:

summary.mcRun <- function(object,...) {
  cat('\n== Summary of a Monte Carlo Simulation Experiment ==\n')
  print(object@settings)
  cat('\n* Data set :: ',object@dataset@name)
  cat('\n* Learner  :: ',deparse(object@learner@func),' with parameters \n')
  for(x in names(object@learner@pars))
    cat('\t',x,' = ',
        deparse(object@learner@pars[[x]]),'\n')
  cat('\n\n* Summary of Experiment Results:\n\n')
  apply(object@foldResults,2,function(x)
        c(avg=mean(x,na.rm=T),std=sd(x,na.rm=T),
          min=min(x,na.rm=T),max=max(x,na.rm=T),
          invalid=sum(is.na(x)))
        )
}


plot.mcRun <- function(x,y,...) {
  sum <- apply(x@foldResults,2,function(x)
               c(avg=mean(x,na.rm=T),std=sd(x,na.rm=T),
                 min=min(x,na.rm=T),max=max(x,na.rm=T),
                 invalid=sum(is.na(x)))
               )
  nstats <- ncol(sum)
  layout(matrix(c(1:nstats,nstats,nstats+1:nstats,2*nstats),ncol=2),
         widths=c(3,1),heights=c(rep(5,nstats),1))
  
  par(mar=c(0,4.1,0,0))
  for(s in 1:(nstats-1)) {
    plot(x@foldResults[,s],type='b',xlab='',ylab=colnames(x@foldResults)[s],
         main='',xaxt='n')
    t <- axTicks(2)
    abline(h=t,lty=3,col='gray')
    abline(h=sum[1,s],lty=2)
  }
  par(mar=c(4.1,4.1,0,0))
  plot(x@foldResults[,nstats],type='b',xlab='Folds',
       ylab=colnames(x@foldResults)[nstats], main='')
  t <- axTicks(2)
  abline(h=t,lty=3,col='gray')
  abline(h=sum[1,nstats],lty=2)
  
  par(mar=c(0,0,0,0))
  for(s in 1:(nstats-1)) {
    boxplot(x@foldResults[,s],type='b',xlab='',ylab='',
            main='',xaxt='n',yaxt='n')
    abline(h=sum[1,s],lty=2)
    t <- axTicks(2)
    abline(h=t,lty=3,col='gray')
    
  }
  
  par(mar=c(4.1,0,0,0))
  boxplot(x@foldResults[,nstats],type='b',xlab='',
          ylab='', main='',xaxt='n',yaxt='n')
  abline(h=sum[1,nstats],lty=2)
  t <- axTicks(2)
  abline(h=t,lty=3,col='gray')
  mtext(paste('DATASET:',x@dataset@name),1,line=1,cex=0.8)
  mtext(paste('LEARNER:',x@learner@func),1,line=2,cex=0.8)
  
}




# ==============================================================
# CLASS UNION: expSettings
# ==============================================================
# Luis Torgo, Aug 2009
# ==============================================================
setClassUnion("expSettings", c("cvSettings", "mcSettings", "hldSettings","loocvSettings","bootSettings"))




# ==============================================================
# CLASS: compExp
# ==============================================================
# Luis Torgo, Jan 2009
# ==============================================================


# --------------------------------------------------------------
# class def
setClass("compExp",
         representation(learners="list",
                        datasets="list",
                        settings="expSettings",
                        foldResults="array")
         )

# --------------------------------------------------------------
# constructor function
compExp <- function(l,t,s,r) {
  o <- new("compExp")
  o@learners <- l
  o@datasets <- t
  o@settings <- s
  o@foldResults <- r
  o
}



# --------------------------------------------------------------
# Methods:

plot.compExp <- function(x,y,stats=dimnames(x@foldResults)[[2]],...) {
  ##require(lattice)
  ##require(grid)
  
  # Function that transforms a 4-dimension array into a data frame
  # it's similar to reshape() goals but I was unable to use the latter...
  flattenRes <- function(foldRes) {
    dim <- prod(dim(foldRes)[c(1,3,4)])
    m <- matrix(NA,nrow=dim,ncol=dim(foldRes)[2])
    m[1:dim(foldRes)[1],] <- foldRes[,,1,1]
    for(d in 1:dim(foldRes)[4])
      for(v in 1:dim(foldRes)[3])
        m[(d-1)*dim(foldRes)[1]*dim(foldRes)[3]+(v-1)*dim(foldRes)[1]+1:dim(foldRes)[1],] <- foldRes[,,v,d]
    
    colnames(m) <- dimnames(foldRes)[[2]]
    d <- data.frame(m,
                    fold=rep(1:dim(foldRes)[1],prod(dim(foldRes)[3:4])),
                    var=rep(dimnames(foldRes)[[3]],each=dim(foldRes)[1],dim(foldRes)[4]),
                    ds=rep(dimnames(foldRes)[[4]],each=prod(dim(foldRes)[c(1,3)]))
                    )
    
  }
  z <- flattenRes(x@foldResults)
  nstats <- length(stats)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nstats, 1)))
  for(s in seq(along=stats)) {
    form <- as.formula(paste('var ~',stats[s],'| ds'))
    gr <- bwplot(form,data=z,
                 panel=function(x,y) {
                   panel.grid(h=0,v=-1)
                   panel.bwplot(x,y)
                 },... )
    
    pushViewport(viewport(layout.pos.col=1,layout.pos.row=s))
    print(gr,newpage=F)
    upViewport()
  }
  
}

summary.compExp <- function(object,...) {
  cat('\n== Summary of a ',
      switch(class(object@settings),
             cvSettings='Cross Validation',
             hldSettings='Hold Out',
             bootSettings='Bootstrap',
             mcSettings='Monte Carlo'
             ),
      ' Experiment ==\n')
  print(object@settings)
  cat('\n* Data sets :: ',
      paste(sapply(object@datasets,function(x) x@name),collapse=', '))
  cat('\n* Learners  :: ',paste(names(object@learners),collapse=', '))
  
  cat('\n\n* Summary of Experiment Results:\n\n')
  ld <- list()
  
  for(d in 1:dim(object@foldResults)[4]) {
    lv <- list()
    cat("\n-> Datataset: ",dimnames(object@foldResults)[[4]][d],'\n')
    for(v in 1:dim(object@foldResults)[3]) {
      cat("\n\t*Learner:",dimnames(object@foldResults)[[3]][v],"\n")
      tab <- apply(object@foldResults[,,v,d,drop=F],2,function(x)
                   c(avg=mean(x,na.rm=T),std=sd(x,na.rm=T),
                     min=if (all(is.na(x))) NA else min(x,na.rm=T),
                     max=if (all(is.na(x))) NA else max(x,na.rm=T),
                     invalid=sum(is.na(x)))
                   )
      print(tab)
      lv <- c(lv,list(tab))
    }
    cat('\n')
    names(lv) <- dimnames(object@foldResults)[[3]]
    ld <- c(ld,list(lv))
  }
  names(ld) <- dimnames(object@foldResults)[[4]]
  invisible(ld)
  
}


# show
setMethod("show","compExp",
          function(object) {
            cat('\n== ',
                switch(class(object@settings),
                       cvSettings='Cross Validation',
                       hldSettings='Hold Out',
                       bootSettings='Bootstrap',
                       mcSettings='Monte Carlo'
                       ),
                ' Experiment ==\n')
            print(object@settings)
            cat(length(object@learners),' learning systems\n')
            cat('tested on ',length(object@datasets),' data sets\n')
          })


# subset
#
# =====================================================
# Method that selects a subset of the experimental 
# results in a object.
# The subsetting criteria can be one of the four dimensions
# of the foldResults array, i.e. the iterations, the statistcs,
# the learner variants, and the data sets, respectively.
# Subsetting expressions can be provided as numbers or as
# dimension names.
# =====================================================
# Luis Torgo, Aug 2009
# =====================================================
# Example runs:
# > plot(subset(nnet,stats='e1',vars=1:4))
#
setMethod("subset",
          signature(x='compExp'),
          function(x,
                   its=1:dim(x@foldResults)[1],
                   stats=1:dim(x@foldResults)[2],
                   vars=1:dim(x@foldResults)[3],
                   dss=1:dim(x@foldResults)[4])
          {
            rr <- x
            if (!identical(vars,1:dim(x@foldResults)[3])) {
              if (is.character(vars) && length(vars) == 1)
                vars <- grep(vars,names(rr@learners))
              rr@learners <- rr@learners[vars]
            }
            if (!identical(dss,1:dim(x@foldResults)[4])) {
              if (is.character(dss) && length(dss) == 1)
                dss <- grep(dss,dimnames(rr@foldResults)[[4]])
              rr@datasets <- rr@datasets[dss]
            }
            if (is.character(stats) && length(stats) == 1)
              stats <- grep(stats,dimnames(rr@foldResults)[[2]])
            rr@foldResults <- rr@foldResults[its,stats,vars,dss,drop=F]
            rr
  
          }
          )



# ==============================================================
# CLASS: tradeRecord
# ==============================================================
# Luis Torgo, Nov 2009
# ==============================================================


# --------------------------------------------------------------
# class def
setClass("tradeRecord",
         representation(trading="zoo",
                        positions="matrix",
                        trans.cost="numeric",
                        init.cap="numeric",
                        policy.func="character",
                        policy.pars="list")
         )


# --------------------------------------------------------------
# constructor function
tradeRecord <- function(t,p,tc,c,pf,pp) {
  o <- new("tradeRecord")
  o@trading <- t
  o@positions <- p
  o@trans.cost <- tc
  o@init.cap <- c
  o@policy.func <- pf
  o@policy.pars <- pp
  o
}



# --------------------------------------------------------------
# Methods:

plot.tradeRecord <- function(x,y,verbose=T,...) {
  
  market <- cbind(y,coredata(x@trading)[,c('Equity','N.Stocks')])
  candleChart(market,
              TA=c(.addEq(),.addSt()),
              ...)
  if (verbose)
    cat('Rentability = ',100*(coredata(market[nrow(market),'Equity'])/
                              coredata(market[1,'Equity'])-1),'%\n')
}


summary.tradeRecord <- function(object,...) {
  cat('\n== Summary of a Trading Simulation with ',nrow(object@trading),' days ==\n')
  cat('\nTrading policy function : ',object@policy.func,'\n')
  cat('Policy function parameters:\n')
  for(x in names(object@policy.pars))
    cat('\t',x,' = ',deparse(object@policy.pars[[x]]),'\n')
  cat('\n')
  cat('Transaction costs : ',object@trans.cost,'\n')
  cat('Initial Equity    : ',round(object@init.cap,1),'\n')
  cat('Final Equity      : ',round(object@trading[nrow(object@trading),'Equity'],1),'  Return : ',
      round(100*(object@trading[nrow(object@trading),'Equity']/object@init.cap - 1),2),'%\n')
  cat('Number of trading positions: ',NROW(object@positions),'\n')
  cat('\nUse function "tradingEvaluation()" for further stats on this simulation.\n\n')
}


# show
setMethod("show","tradeRecord",
          function(object) {
            cat('\nObject of class tradeRecord with slots:\n\n')
            cat('\t trading: <xts object with a numeric ',dim(object@trading)[1],'x',dim(object@trading)[2],' matrix>\n')
            cat('\t positions: <numeric ',dim(object@positions)[1],'x',dim(object@positions)[2],' matrix>\n')
            cat('\t init.cap : ',object@init.cap,'\n')
            cat('\t trans.cost : ',object@trans.cost,'\n')
            cat('\t policy.func : ',object@policy.func,'\n')
            cat('\t policy.pars : <list with ',length(object@policy.pars),' elements>\n\n')
          })




# This function plots the trading record
# It requires as a second argument the market quotes during the test period
#.Eq <- function(p) p[,'Equity']
#.St <- function(p) p[,'N.Stocks']
#.addEq <- newTA(FUN = .Eq, col = 'red', legend = "Equity")
#.addSt <- newTA(FUN = .St, col = 'green', legend = "N.Stocks")

#plotTradeRecord <- function(o,market,verbose=F,...) {
#  market <- cbind(market,coredata(o@trading)[,c('Equity','N.Stocks')])
#  candleChart(market,TA=c(.addEq(),.addSt()),...)
#  if (verbose) cat('Rentability = ',100*(coredata(market[nrow(market),'Equity'])/coredata(mar#ket[1,'Equity'])-1),'%\n')
#  
#}



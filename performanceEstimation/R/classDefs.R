################################################################# 
## THIS FILE CONTAINS THE CLASSES AND RESPECTIVE CONSTRUCTORS   #
## THAT ARE DEFINED IN THIS PACKAGE                             #
#################################################################
## Author : Luis Torgo (ltorgo@dcc.fc.up.pt)     Date: Nov 2013 #
## License: GPL (>= 2)                                          #
#################################################################


setClassUnion("DataSource",c("call","name","data.frame"))
setClassUnion("OptList",c("list","NULL"))
setClassUnion("OptMatrix",c("matrix","NULL"))
setClassUnion("OptString",c("character","NULL"))
setClassUnion("NameOrCall",c("call","name"))

## ==============================================================
## CLASS: PredTask
##
## Class for storing information concerning a prediction task
## ==============================================================


## --------------------------------------------------------------
## Definition
##
setClass("PredTask",
         slots=c(formula="formula",
                 dataSource="DataSource",
                 taskName="character",
                 type="character",
                 target="character")
         )


## --------------------------------------------------------------
## constructor
##
PredTask <- function(form,data,taskName=NULL,type=NULL,copy=FALSE) {
  if (missing(form) || missing(data))
    stop('\nYou need to provide a formula and a data frame.\n',call.=FALSE)

  tgt <- deparse(form[[2]])

  if (!copy) {
      data <- substitute(data)
      eval.env <- 2
      if (is.null(taskName)) taskName <- paste(deparse(data),tgt,sep=".")
  } else {
      eval.env <- 1
      if (is.null(taskName)) taskName <-  paste(deparse(substitute(data)),tgt,sep=".")
  }

  if (inherits(try(mf <- model.frame(form,eval(data,envir=eval.env),na.action=NULL),TRUE),"try-error"))
    stop('\nInvalid formula for the given data frame.\n',call.=FALSE)

  
  if (is.null(type)) {
      taskType <- if (is.factor(eval(data,envir=eval.env)[,tgt])) "class" else "regr"
  } else {
      if (!(type %in% c("class","regr","ts")))
          stop(paste("PredTask::",type,"tasks not implemented."),call.=FALSE)
      taskType <- type
  }

#  if (taskType == "ts" && !is.numeric(eval(data,envir=eval.env)[[tgt]]))
#      stop("PredTask:: time series tasks should have numeric target.",call.=FALSE)
  
  new("PredTask", formula=form, dataSource=data, taskName=taskName, type=taskType, target=tgt)
}



## ==============================================================
## CLASS: Workflow
##
## Class for storing a workflow to solve a predictive task
## ==============================================================


## --------------------------------------------------------------
## Definition
##
setClass("Workflow",
         slots=c(name="character",
                 func="character",
                 pars="list",
                 deps="OptList"))



## --------------------------------------------------------------
## Constructor
##

Workflow <- function(wf, ..., wfID, deps=NULL) {

    wf.pars <- list(...)
    
    ## if no ID was provided then it is one of the standard workflows
    if (missing(wf) || wf == "standardWF" || wf == "timeseriesWF") {
        if ("type" %in% names(wf.pars)) {
            n <- paste(wf.pars[["learner"]],wf.pars[["type"]],sep='.')
            f <- "timeseriesWF"
        } else {
            n <- wf.pars[["learner"]]
            f <- "standardWF"
        }
    ## using a user-defined workflow
    } else {
        n <- f <- wf
    }
    if (!missing(wfID)) n <- wfID
    
    new("Workflow", name=n, func=f, pars=wf.pars, deps=deps)

}



## ==============================================================
## CLASS: EstCommon
##
## A class containing the common estimation settings
## ==============================================================
setClass("EstCommon",
         slots=c(seed='numeric',         # seed of the random generator
                 dataSplits='OptList')   # user supplied data splits
         )



## ==============================================================
## CLASS: CV
##
## A class containing the settings of a cross validation experiment
## ==============================================================


## --------------------------------------------------------------
## Definition
##
setClass("CV",
         slots=c(nReps='numeric',      # nr. of repetitions
                 nFolds='numeric',     # nr. of folds of each rep.
                 strat='logical'),     # is the sampling stratified?
         contains="EstCommon"
         )


## --------------------------------------------------------------
## constructor
##
CV <- function(nReps=1,nFolds=10,
                       seed=1234,strat=FALSE,
                       dataSplits=NULL) {
    new("CV",
        nReps=nReps,
        nFolds=if (is.null(dataSplits)) nFolds else length(dataSplits)/nReps,
        seed=seed,strat=strat,dataSplits=dataSplits)
}


## ==============================================================
## CLASS: Holdout
##
## A class containing the settings of a holdout experiment
## ==============================================================


## --------------------------------------------------------------
## Definition
##
setClass("Holdout",
         slots=c(nReps='numeric', # number of repetitions
                 hldSz='numeric', # the size (0..1) of the holdout
                 strat='logical'),# is the sampling stratified?
         contains="EstCommon"
         )


## --------------------------------------------------------------
## Constructor
##
Holdout <- function(nReps=1,hldSz=0.3,
                        seed=1234,strat=FALSE,
                        dataSplits=NULL) {
    new("Holdout",
        nReps=if (is.null(dataSplits)) nReps else length(dataSplits),
        hldSz=if (is.null(dataSplits)) hldSz else length(dataSplits[[1]]),
        seed=seed,strat=strat,dataSplits=dataSplits)
}


## ==============================================================
## CLASS: LOOCV
##
## A class containing the settings of a leave one out cross validation
## experiment
## ==============================================================


## --------------------------------------------------------------
## Definition
##
setClass("LOOCV",
         contains="EstCommon"
         )


## --------------------------------------------------------------
## constructor
LOOCV <- function(seed=1234,dataSplits=NULL)
  new("LOOCV",
      seed=seed,
      dataSplits=dataSplits)



## ==============================================================
## CLASS: Bootstrap
##
## A class containing the settings of a boostrap experiment
## ==============================================================


## --------------------------------------------------------------
## Definition
##
setClass("Bootstrap",
         slots=c(type='character', # type of boostrap ("e0" or ".632")
                 nReps='numeric'), # number of repetitions
         contains="EstCommon"
         )


## --------------------------------------------------------------
## constructor
##
Bootstrap <- function(type='e0',nReps=200,seed=1234,dataSplits=NULL) {
     new("Bootstrap",
         type=type,
         nReps=if (is.null(dataSplits)) nReps else length(dataSplits),
         seed=seed,dataSplits=dataSplits)
}



## ==============================================================
## CLASS: MonteCarlo
##
## A class containing the settings of a monte carlo experiment
## ==============================================================


## --------------------------------------------------------------
## Definition
##
setClass("MonteCarlo",
         slots=c(nReps='numeric',
                 szTrain='numeric',
                 szTest='numeric'),
         contains="EstCommon"
         )


## --------------------------------------------------------------
## constructor
##
MonteCarlo <- function(nReps=10,szTrain=0.25,szTest=0.25,
                       seed=1234,dataSplits=NULL)
    new("MonteCarlo",
        nReps= if (is.null(dataSplits)) nReps else length(dataSplits),
        szTrain=szTrain,szTest=szTest,
        seed=seed,dataSplits=dataSplits)



## ==============================================================
## CLASS UNION: EstimationMethod
##
## A class encapsulating all types of Estimation Experiments
## ==============================================================

setClassUnion("EstimationMethod",
              c("CV", "MonteCarlo", "Holdout",
                "LOOCV","Bootstrap"))





## ==============================================================
## CLASS: EstimationTask
##
## A class containing the information on a estimation task, i.e.
## the metrics (and eventually the function to calculate them) and
## the estimation methodology to use
## ==============================================================
setClass("EstimationTask",
         slots=c(metrics='OptString',        # the metrics to be estimated
                 method="EstimationMethod",  # the estimation method to use
                 evaluator='character',      # function used to calculate the metrics
                 evaluator.pars='OptList',   # pars to this function
                 trainReq='logical'          # is the training data required?
         )
         )

## --------------------------------------------------------------
## constructor
##
EstimationTask <- function(metrics=NULL,method=CV(),
                           evaluator="",evaluator.pars=NULL,
                           trainReq=FALSE) {
    new("EstimationTask",
        metrics=metrics,method=method,
        evaluator=evaluator,evaluator.pars=evaluator.pars,
        trainReq=trainReq)
}


    
## ==============================================================
## CLASS: EstimationResults
##
## A class containing the results of a single experiment, i.e. the
## the estimation results of applying a single workflow to a single
## predictive task
## ==============================================================


## --------------------------------------------------------------
## Constructor
##
setClass("EstimationResults",
         slots=c(task             = "PredTask",
                 workflow         = "Workflow",
                 estTask          = "EstimationTask",
                 iterationsScores = "matrix",   # nIts x nStats
                 iterationsInfo   = "list"      # list of nIts lists with comps preds, info and train
                )
         )


## --------------------------------------------------------------
## constructor
##
EstimationResults <- function(t,w,et,sc,e) {
  o                  <- new("EstimationResults")
  o@task             <- t
  o@workflow         <- w
  o@estTask          <- et
  o@iterationsScores <- sc
  ## classification tasks, code back predictions to class labels
  if (et@evaluator=="" & is.factor(model.response(model.frame(t@formula,eval(t@dataSource))))) {
      for (i in 1:length(e))
          e[[i]]$preds <- factor(e[[i]]$preds,levels=levels(responseValues(t@formula,eval(t@dataSource))))
  }
  o@iterationsInfo  <- e
  o
}



## ==============================================================
## CLASS: ComparisonResults
##
## A class containing the results of a full experimental comparison,
## i.e. the  estimation results of applying several workflows to 
## several predictive tasks
## ==============================================================


## --------------------------------------------------------------
## Definition
##
setClass("ComparisonResults",contains="list")

## --------------------------------------------------------------
## constructor
##
ComparisonResults <- function(t) new("ComparisonResults",t)



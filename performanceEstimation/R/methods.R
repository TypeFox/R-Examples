################################################################# 
## THIS FILE CONTAINS THE METHODS OF THE CLASSES THAT WERE      #
## DEFINED IN THIS PACKAGE                                      #
#################################################################
## Author : Luis Torgo (ltorgo@dcc.fc.up.pt)     Date: Nov 2013 #
## License: GPL (>= 2)                                          #
#################################################################




################################################################
## PredTask methods
################################################################


setMethod("show",
          "PredTask",
          function(object) {
            cat('Prediction Task Object:\n')
            cat('\tTask Name         ::',object@taskName,"\n")
            cat('\tTask Type         ::',
                if (object@type == "class") "classification" else "regression","\n")
            cat('\tTarget Feature    ::',object@target,"\n")
            cat('\tFormula           :: ')
            print(object@formula)
            cat('\tTask Data Source  :: ')
            if (is.data.frame(object@dataSource)) cat('internal ',paste(dim(object@dataSource),collapse="x"),"data frame.")
            else cat(deparse(object@dataSource))
            cat("\n")
          }
          )




################################################################
## Workflow methods
################################################################

setMethod("show",
          "Workflow",
          function(object) {
            cat('Workflow Object:\n\tWorkflow ID       :: ',object@name,
                                '\n\tWorkflow Function :: ',object@func)
            if (length(object@pars)) {
                cat('\n\t     Parameter values:\n')
                for(n in names(object@pars)) {
                    if (!is.null(object@pars[[n]])) {
##                    cat('\t\t',n,' = ',deparse(object@pars[[n]]),'\n')
                        k <- object@pars[[n]]
                        k <- if (is.list(k)) paste(paste(names(k),k,sep='='),collapse=" ") else paste(k,collapse=' ')
                        k <- if (nchar(k) > 20) paste(substr(k,1,20),' ...') else k
                        cat('\t\t',n,' -> ',k,'\n')
                    }
                }
            }
            cat('\n')
          }
          )


setMethod("summary",
          "Workflow",
          function(object) {
              print(object)
              cat("\nTo apply the workflow on some predictive task use the function 'runWorkflow()'.",
                  "\nExample: 'runWorkflow(myWFobject,Y ~ .,trainData,testData)'\n")
          }
          )





################################################################
## CV Methods:
################################################################
setMethod("show",
          "CV",
          function(object) {
            userSplit <- !is.null(object@dataSplits)
            cat(ifelse(!userSplit & object@strat,'Stratified ',''),
                object@nReps,'x',object@nFolds,
                '- Fold Cross Validation\n')
            if (!userSplit)
              cat('\t Run with seed = ', object@seed,'\n')
            else
              cat('\t User-supplied data splits\n')
          })







################################################################
## Holdout methods:
################################################################

setMethod("show",
          "Holdout",
          function(object) {
            userSplit <- !is.null(object@dataSplits)
            cat(ifelse(!userSplit & object@strat,'Stratified ',''),
                object@nReps,'x')
            
            if (object@hldSz < 1) cat(100*(1-object@hldSz),'%/',100*object@hldSz,'% Holdout\n')
            else cat(object@hldSz," cases Holdout\n")
            
            if (!userSplit)
              cat('\t Run with seed = ',object@seed,'\n')
            else
              cat('\t User-supplied data splits\n')
          })




################################################################
## LOOCV methods:
################################################################


setMethod("show","LOOCV",
          function(object) {
              userSplit <- !is.null(object@dataSplits)
              cat('LOOCV experiment\n')
              if (!userSplit)
                  cat('\t Run with seed = ',object@seed,'\n')
              else
                  cat('\t Run user-supplied data splits\n')
         })





################################################################
## Bootstrap methods:
################################################################


setMethod("show",
          "Bootstrap",
          function(object) {
            userSplit <- !is.null(object@dataSplits)
            cat(object@nReps,' repetitions of ',ifelse(object@type=='e0','e0','.632'),
                ' Bootstrap experiment\n')
            if (!userSplit)
              cat('\t Run with seed = ', object@seed,'\n')
            else
              cat('\t User-supplied data splits\n')
         })



################################################################
## MonteCarlo methods:
################################################################


setMethod("show",
          "MonteCarlo",
          function(object) {
            userSplit <- !is.null(object@dataSplits)
            cat(object@nReps,
               ' repetitions Monte Carlo Simulation')
            if (userSplit) {
              cat(' using user-supplied data splits\n')
            } else {
              cat(' using:',
                  '\n\t seed = ', object@seed,
                  '\n\t train size = ',object@szTrain,
                  ifelse(object@szTrain<1,'x NROW(DataSet)',' cases'),
                  '\n\t test size = ',object@szTest,
                  ifelse(object@szTest<1,'x NROW(DataSet)',' cases'),
                  '\n'
                  )
            }
         })




################################################################
## EstimationTask methods:
################################################################


setMethod("show",
          "EstimationTask",
          function(object) {
              if (is.null(object@metrics)) cat("Task for estimating all metrics of the selected evaluation function using\n")
              else cat("Task for estimating ",paste(object@metrics,collapse=",")," using\n")
              print(object@method)
          }
         )





################################################################
## EstimationResults methods:
################################################################

setMethod("show",
          "EstimationResults",
          function(object) {
            print(object@estTask)
            cat("\nTask    ID :: ",object@task@taskName,"\nWorflow ID :: ",object@workflow@name,"\n")
            cat("\nOverview of the Scores Estimates:\n")
            print(.scores2summary(object)[1:2,,drop=FALSE])
            cat("\n")
            })



setMethod("summary",
          "EstimationResults",
          function(object) {
              cat('\n*** Summary of a ',
                  switch(class(object@estTask@method),
                         CV='Cross Validation',
                         Holdout='Hold Out',
                         MonteCarlo='Monte Carlo',
                         Bootstrap='Bootstrap',
                         LOOCV='Loocv',
                         ),
                  ' Estimation Experiment ***\n\n')

              print(object@estTask)
              cat('\n* Predictive Task ID :: ',object@task@taskName)
              cat('\n\tTask Type         ::',
                  if (object@task@type == "class") "classification" else "regression","\n")
              cat('\tTarget Feature    ::',object@task@target,"\n")
              cat('\tFormula           :: ')
              print(object@task@formula)
              cat('\tTask Data Source  ::',deparse(object@task@dataSource),"\n")
              cat('\n* Workflow        ID :: ',object@workflow@name,
                  '\n\tWorkflow Function :: ',object@workflow@func)
              if (length(object@workflow@pars)) {
                  cat('\n\t     Parameter values:\n')
                  for(n in names(object@workflow@pars)) {
                      if (!is.null(object@workflow@pars[[n]])) {
                          k <- object@workflow@pars[[n]]
                          k <- if (is.list(k)) paste(paste(names(k),k,sep='='),collapse=" ") else paste(k,collapse=' ')
                          k <- if (nchar(k) > 20) paste(substr(k,1,20),' ...') else k
                          cat('\t\t',n,' -> ',k,'\n')
                      }
                  }
              }
              cat('\n\n* Summary of Score Estimation Results:\n\n')
              print(.scores2summary(object))
          })


if (!isGeneric("plot"))  setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setMethod("plot",
          "EstimationResults",
          function(x,y,...) {
              
              nstats <- ncol(x@iterationsScores)
              
              tit <- paste(x@workflow@name,
                  switch(class(x@estTask@method),
                         CV='Cross Validation',
                         Holdout='Hold Out',
                         LOOCV='Leave One Out',
                         Bootstrap='Bootstrap',
                         MonteCarlo='Monte Carlo'
                         ),
                           "estimation on",x@task@taskName,sep=" "
                  )
              if (nstats == 1) {
                  plt <- ggplot2::qplot(1:nrow(x@iterationsScores),
                               x@iterationsScores[,1],
                               main=tit,
                               xlab='Estimation Iterations',
                               ylab=colnames(x@iterationsScores)[1]) +
                         ggplot2::geom_smooth(method='loess',size=1) +
                         ggplot2::geom_line(stat="hline",yintercept="mean",color="red") +
                         ggplot2::scale_x_discrete()
                  ##print(plt)
              } else {
                  dt <- .scores2long(x@iterationsScores)
                  plt <- ggplot2::ggplot(dt,ggplot2::aes_string(x="rep",y="score")) + 
                      ggplot2::ggtitle(tit) +
                      ggplot2::ylab('Metrics Scores') + ggplot2::xlab('Estimation Iterations')+
                      ggplot2::geom_smooth(ggplot2::aes_string(group="stat"),method='loess',size=1) +
                      ggplot2::geom_line(stat="hline",yintercept="mean",color="red") +
                      ggplot2::scale_x_discrete() +ggplot2::theme(axis.text.x=ggplot2::element_text(angle=270,size=10,vjust=0.5,hjust=0))+
                      ggplot2::facet_grid( stat ~ .,scales = "free_y")
                  ##print(plt)
              }
              plt
              
          }
          )



################################################################
# Comparison Results methods
################################################################

setMethod("plot",
          "ComparisonResults",
          function(x,y,...) {
              
              allRes <- NULL
              taskNames <- names(x)
              for(t in 1:length(x)) {
                  task <- taskNames[t]
                  sysNames <- names(x[[t]])
                  for(s in 1:length(x[[t]])) {
                      d <- .scores2long(x[[t]][[s]]@iterationsScores)
                      d <- cbind(d,sys=sysNames[s],task=taskNames[t])
                      allRes <- rbind(allRes,d)
                  }
              }

              tlt <- paste(switch(class(x[[1]][[1]]@estTask@method),
                                  CV='Cross Validation',
                                  Holdout='Hold Out',
                                  LOOCV='Leave One Out',
                                  Bootstrap='Bootstrap',
                                  MonteCarlo='Monte Carlo'
                                  ),"Performance Estimation Results")
              plt <- ggplot2::ggplot(allRes,ggplot2::aes_string(y="score",x="sys")) +
                     ggplot2::geom_boxplot(ggplot2::aes_string(group="sys")) +
                     ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .15,height=0),color="red",size=3,alpha=0.25) + 
                     ggplot2::ggtitle(tlt) +
                     ggplot2::ylab("Distribution of Statistics Scores") +
                     ggplot2::xlab("Alternative Workflows") +
                     ggplot2::facet_grid(stat ~ task,scales="free_y") +
                     ggplot2::theme(axis.text.x=ggplot2::element_text(angle=270,size=10,vjust=0.5))
              #print(plt)
              plt
                      
  })



setMethod("summary",
          "ComparisonResults",
          function(object) {
              cat('\n== Summary of a ',
                  switch(class(object[[1]][[1]]@estTask@method),
                         CV='Cross Validation',
                         Holdout='Hold Out',
                         LOOCV='Leave One Out',
                         Bootstrap='Bootstrap',
                         MonteCarlo='Monte Carlo'
                         ),
                  'Performance Estimation Experiment ==\n\n')
              print(object[[1]][[1]]@estTask)
              cat('\n* Predictive Tasks :: ',
                  paste(names(object),collapse=', '))
              cat('\n* Workflows  :: ',paste(names(object[[1]]),collapse=', '),"\n")
              
              ##cat('\n\n* Summary of Experiment Results:\n')
              ld <- list()
              for(d in 1:length(object)) {
                  lv <- list()
                  cat("\n-> Task: ",names(object)[d])
                  for(v in 1:length(object[[d]])) {
                      cat("\n  *Workflow:",names(object[[d]])[v],"\n")
                      ss <- .scores2summary(object[[d]][[v]])
                      print(ss)
                      lv <- c(lv,list(ss))
                  }
                  ##cat('\n')
                  names(lv) <- names(object[[d]])
                  ld <- c(ld,list(lv))
              }
              names(ld) <- names(object)
              invisible(ld)
              
          })



setMethod("show",
          "ComparisonResults",
          function(object) {
            cat('\n== ',
                switch(class(object[[1]][[1]]@estTask@method),
                       CV='Cross Validation',
                       Holdout='Hold Out',
                       Bootstrap='Bootstrap',
                       LOOCV='Leave One Out',
                       MonteCarlo='Monte Carlo'
                       ),
                'Performance Estimation Experiment ==\n\n')
            print(object[[1]][[1]]@estTask)
            cat("\n",length(object[[1]]),' workflows applied to ',
                length(object),' predictive tasks\n')
          })



# =====================================================
# Method that selects a subset of the experimental 
# results in a object.
# The subsetting criteria can be one of the four dimensions
# of the foldResults array, i.e. the iterations, the statistcs,
# the workflow variants, and the data sets, respectively.
# Subsetting expressions can be provided as numbers or as
# dimension names.
# =====================================================
# Luis Torgo, Aug 2009
# =====================================================
# Example runs:
# > plot(subset(nnet,stats='e1',vars=1:4))
#
setMethod("subset",
          signature(x='ComparisonResults'),
          function(x,
                   tasks=1:length(x),
                   workflows=1:length(x[[1]]),
                   metrics=1:dim(x[[1]][[1]]@iterationsScores)[2],
                   partial=TRUE) {
              mf <- if (partial) "grep" else "match"
              rr <- x
              if (!identical(workflows,1:length(x[[1]]))) {
                  if (is.character(workflows))
                      workflows <- unlist(lapply(workflows,function(w) do.call(mf,list(w,names(rr[[1]])))))
                  rr <- lapply(rr,function(t) t[workflows])
              }
              if (!identical(tasks,1:length(x))) {
                  if (is.character(tasks))
                      tasks <- unlist(lapply(tasks,function(t) do.call(mf,list(t,names(rr)))))
                  rr <- rr[tasks]
              }
              if (is.character(metrics)) 
                  metrics <- unlist(lapply(metrics,function(m) do.call(mf,list(m,colnames(x[[1]][[1]]@iterationsScores)))))
              rr <- lapply(rr,function(t) lapply(t,function(s) {sn <- s; sn@iterationsScores <- s@iterationsScores[,metrics,drop=F] ; sn@estTask@metrics <- metricNames(rr)[metrics] ; sn}))
              
              ComparisonResults(rr)
          }
          )






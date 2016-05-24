#####################################################################################
## Author: Daniel Sabanes Bove [sabanesd *a*t* roche *.* com],
##         Wai Yin Yeung [ w *.* yeung1 *a*t* lancaster *.* ac *.* uk]
## Project: Object-oriented implementation of CRM designs
##
## Time-stamp: <[Simulations-methods.R] by DSB Fre 16/01/2015 13:41>
##
## Description:
## Methods for handling the simulations output.
##
## History:
## 19/02/2014   file creation
## 30/07/2014   added in methods for pseudo models simulations
###################################################################################


##' @include Simulations-class.R
##' @include helpers.R
{}


##' Plot simulations
##'
##' Summarize the simulations with plots
##'
##' This plot method can be applied to \code{\linkS4class{GeneralSimulations}}
##' objects in order to summarize them graphically. Possible \code{type}s of
##' plots at the moment are: \describe{ \item{trajectory}{Summary of the
##' trajectory of the simulated trials} \item{dosesTried}{Average proportions of
##' the doses tested in patients} } You can specify one or both of these in the
##' \code{type} argument.
##'
##' @param x the \code{\linkS4class{GeneralSimulations}} object we want
##' to plot from
##' @param y missing
##' @param type the type of plots you want to obtain.
##' @param \dots not used
##' @return A single \code{\link[ggplot2]{ggplot}} object if a single plot is
##' asked for, otherwise a \code{\link{gridExtra}{gTree}} object.
##'
##' @importFrom ggplot2 ggplot geom_step geom_bar aes xlab ylab
##' scale_linetype_manual
##' @importFrom gridExtra arrangeGrob
##' 
##' @example examples/Simulations-method-plotSIMsingle.R
##' @export
##' @keywords methods
setMethod("plot",
          signature=
          signature(x="GeneralSimulations",
                    y="missing"),
          def=
          function(x,
                   y,
                   type=
                   c("trajectory",
                     "dosesTried"),
                   ...){

              ## which plots should be produced?
              type <- match.arg(type,
                                several.ok=TRUE)
              stopifnot(length(type) > 0L)

              ## start the plot list
              plotList <- list()
              plotIndex <- 0L

              ## summary of the trajectories
              if("trajectory" %in% type)
              {
                  ## get a matrix of the simulated dose trajectories,
                  ## where the rows correspond to the simulations and
                  ## the columns to the patient index:

                  simDoses <- lapply(x@data,
                                     slot,
                                     "x")

                  maxPatients <- max(sapply(simDoses, length))

                  simDosesMat <- matrix(data=NA,
                                        nrow=length(simDoses),
                                        ncol=maxPatients)

                  for(i in seq_along(simDoses))
                  {
                      simDosesMat[i, seq_along(simDoses[[i]])] <-
                          simDoses[[i]]
                  }


                  ## extract statistics
                  stats <- c("Minimum",
                             "Lower Quartile",
                             "Median",
                             "Upper Quartile",
                             "Maximum")
                  traj.df <-
                      data.frame(patient=
                                 rep(seq_len(maxPatients), each=5L),
                                 Statistic=
                                 factor(rep(stats,
                                            maxPatients),
                                        levels=stats),
                                 traj=
                                 c(apply(simDosesMat, 2L, quantile,
                                         na.rm=TRUE)))

                  ## linetypes for the plot
                  lt <- c("Median" = 1,
                          "Lower Quartile" = 2,
                          "Upper Quartile" = 2,
                          "Minimum" = 4,
                          "Maximum"= 4)

                  ## save the plot
                  plotList[[plotIndex <- plotIndex + 1L]] <-
                      ggplot() +
                          geom_step(aes(x=patient,
                                        y=traj,
                                        group=Statistic,
                                        linetype=Statistic),
                                    size=1.2, colour="blue", data=traj.df) +
                                        ## scale_linetype_manual(values=lt) +
                                            xlab("Patient") + ylab("Dose Level")
              }

              ## average distribution of the doses tried
              if("dosesTried" %in% type)
              {
                  ## get the doses tried
                  simDoses <- lapply(x@data,
                                     slot,
                                     "x")

                  ## get the dose distributions by trial
                  doseDistributions <-
                      sapply(simDoses,
                             function(s){
                                 prop.table(table(factor(s,
                                                         levels=
                                                         x@data[[1]]@doseGrid)))
                             })

                  ## derive the average dose distribution across trial
                  ## simulations
                  averageDoseDist <- rowMeans(doseDistributions)

                  ## get in data frame shape
                  dat <- data.frame(dose=as.numeric(names(averageDoseDist)),
                                    perc=averageDoseDist * 100)

                  ## produce and save the plot
                  plotList[[plotIndex <- plotIndex + 1L]] <-
                      ggplot() +
                          geom_bar(data=as.data.frame(dat),
                                   aes(x=dose, y=perc),
                                   stat="identity",
                                   position="identity",
                                   width=1) +
                                       xlab("Dose level") +
                                           ylab("Average proportion [%]")
              }
              ## then finally plot everything

              ## if there is only one plot
              if(identical(length(plotList),
                           1L))
              {
                  ## just return it
                  return(plotList[[1L]])
              } else {
                  ## otherwise arrange them
                  ret <- do.call(gridExtra::arrangeGrob,
                                 plotList)
                  return(ret)
              }
          })


##' Plot dual-endpoint simulations
##'
##' This plot method can be applied to \code{\linkS4class{DualSimulations}}
##' objects in order to summarize them graphically. In addition to the standard
##' plot types, there is
##' \describe{
##' \item{sigma2W}{Plot a boxplot of the final biomarker variance estimates in
##' the simulated trials}
##' \item{rho}{Plot a boxplot of the final correlation estimates in
##' the simulated trials}
##' }
##'
##' @param x the \code{\linkS4class{DualSimulations}} object we want
##' to plot from
##' @param y missing
##' @param type the type of plots you want to obtain.
##' @param \dots not used
##' @return A single \code{\link[ggplot2]{ggplot}} object if a single plot is
##' asked for, otherwise a \code{\link{gridExtra}{gTree}} object.
##'
##' @importFrom ggplot2 qplot coord_flip scale_x_discrete
##' @importFrom gridExtra arrangeGrob
##' 
##' @example examples/Simulations-method-plot-DualSimulations.R
##' @export
##' @keywords methods
setMethod("plot",
          signature=
          signature(x="DualSimulations",
                    y="missing"),
          def=
          function(x,
                   y,
                   type=
                   c("trajectory",
                     "dosesTried",
                     "sigma2W",
                     "rho"),
                   ...){

              ## start the plot list
              plotList <- list()
              plotIndex <- 0L

              ## which plots should be produced?
              type <- match.arg(type,
                                several.ok=TRUE)
              stopifnot(length(type) > 0L)

              ## substract the specific plot types for
              ## dual-endpoint simulation results
              typeReduced <- setdiff(type,
                                     c("sigma2W", "rho"))

              ## are there more plots from general?
              moreFromGeneral <- (length(typeReduced) > 0)

              ## if so, then produce these plots
              if(moreFromGeneral)
              {
                  genPlot <- callNextMethod(x=x, y=y, type=typeReduced)
              }

              ## now to the specific dual-endpoint plots:

              ## biomarker variance estimates boxplot
              if("sigma2W" %in% type)
              {
                  ## save the plot
                  plotList[[plotIndex <- plotIndex + 1L]] <-
                      qplot(factor(0), y=y, data=data.frame(y=x@sigma2West), geom="boxplot",
                            xlab="", ylab="Biomarker variance estimates") +
                                coord_flip() + scale_x_discrete(breaks=NULL)
              }

              ## correlation estimates boxplot
              if("rho" %in% type)
              {
                  ## save the plot
                  plotList[[plotIndex <- plotIndex + 1L]] <-
                      qplot(factor(0), y=y, data=data.frame(y=x@rhoEst), geom="boxplot",
                            xlab="", ylab="Correlation estimates") +
                                coord_flip() + scale_x_discrete(breaks=NULL)
              }

              ## then finally plot everything
              if(identical(length(plotList),
                           0L))
              {
                  return(genPlot)
              } else if(identical(length(plotList),
                                  1L))
              {
                  ret <- plotList[[1L]]
              } else {
                  ret <- do.call(gridExtra::arrangeGrob,
                                 plotList)
              }

              if(moreFromGeneral)
              {
                  ret <- gridExtra::arrangeGrob(genPlot, ret)
              }

              return(ret)
          })



##' Summarize the simulations, relative to a given truth
##'
##' @param object the \code{\linkS4class{GeneralSimulations}} object we want to
##' summarize
##' @param truth a function which takes as input a dose (vector) and returns the
##' true probability (vector) for toxicity
##' @param target the target toxicity interval (default: 20-35\%) used for the
##' computations
##' @param \dots Additional arguments can be supplied here for \code{truth}
##' @return an object of class \code{\linkS4class{GeneralSimulationsSummary}}
##'
##' @export
##' @keywords methods
setMethod("summary",
          signature=
          signature(object="GeneralSimulations"),
          def=
          function(object,
                   truth,
                   target=c(0.2, 0.35),
                   ...){

              ## extract dose grid
              doseGrid <- object@data[[1]]@doseGrid

              ## evaluate true toxicity at doseGrid
              trueTox <- truth(doseGrid, ...)

              ## find dose interval corresponding to target tox interval
              targetDoseInterval <-
                  sapply(target,
                         function(t){
                             ## we have to be careful because it might be
                             ## that in the range of the dose grid, no
                             ## doses can be found that match the target
                             ## interval boundaries!
                             ## In that case we want to return NA
                             r <- try(uniroot(f=function(x){truth(x, ...) - t},
                                              interval=
                                              range(doseGrid))$root,
                                      silent=TRUE)
                             if(inherits(r, "try-error"))
                             {
                                 return(NA_real_)
                             } else {
                                 return(r)
                             }
                         })

              ## what are the levels above target interval?
              xAboveTarget <- which(trueTox > target[2])

              ## proportion of DLTs in a trial:
              propDLTs <- sapply(object@data,
                                 function(d){
                                     mean(d@y)
                                 })

              ## mean toxicity risk
              meanToxRisk <- sapply(object@data,
                                    function(d){
                                        mean(trueTox[d@xLevel])
                                    })

              ## doses selected for MTD
              doseSelected <- object@doses

              ## replace NA by 0
              doseSelected[is.na(doseSelected)] <- 0

              ## dose most often selected as MTD
              doseMostSelected <-
                  as.numeric(names(which.max(table(doseSelected))))
              xMostSelected <-
                  match(doseMostSelected,
                        table=doseGrid)

              ## observed toxicity rate at dose most often selected
              ## Note: this does not seem very useful!
              ## Reason: In case of a fine grid, few patients if any
              ## will have been treated at this dose.
              tmp <-
                  sapply(object@data,
                         function(d){
                             whichAtThisDose <- which(d@x == doseMostSelected)
                             nAtThisDose <- length(whichAtThisDose)
                             nDLTatThisDose <- sum(d@y[whichAtThisDose])
                             return(c(nAtThisDose=nAtThisDose,
                                      nDLTatThisDose=nDLTatThisDose))
                         })

              obsToxRateAtDoseMostSelected <-
                  mean(tmp["nDLTatThisDose",]) / mean(tmp["nAtThisDose",])

              ## number of patients overall
              nObs <- sapply(object@data,
                             slot,
                             "nObs")

              ## number of patients treated above target tox interval
              nAboveTarget <- sapply(object@data,
                                     function(d){
                                         sum(d@xLevel %in% xAboveTarget)
                                     })

              ## Proportion of trials selecting target MTD
              toxAtDoses <- truth(doseSelected, ...)
              propAtTarget <- mean((toxAtDoses > target[1]) &
                                   (toxAtDoses < target[2]))

              ## give back an object of class GeneralSimulationsSummary,
              ## for which we then define a print / plot method
              ret <-
                  .GeneralSimulationsSummary(
                      target=target,
                      targetDoseInterval=targetDoseInterval,
                      nsim=length(object@data),
                      propDLTs=propDLTs,
                      meanToxRisk=meanToxRisk,
                      doseSelected=doseSelected,
                      doseMostSelected=doseMostSelected,
                      obsToxRateAtDoseMostSelected=obsToxRateAtDoseMostSelected,
                      nObs=nObs,
                      nAboveTarget=nAboveTarget,
                      toxAtDosesSelected=toxAtDoses,
                      propAtTarget=propAtTarget,
                      doseGrid=doseGrid)

              return(ret)
          })


##' Summarize the model-based design simulations, relative to a given truth
##'
##' @param object the \code{\linkS4class{Simulations}} object we want to
##' summarize
##' @param truth a function which takes as input a dose (vector) and returns the
##' true probability (vector) for toxicity
##' @param target the target toxicity interval (default: 20-35\%) used for the
##' computations
##' @param \dots Additional arguments can be supplied here for \code{truth}
##' @return an object of class \code{\linkS4class{SimulationsSummary}}
##'
##' @example examples/Simulations-method-summary.R
##' @export
##' @keywords methods
setMethod("summary",
          signature=
          signature(object="Simulations"),
          def=
          function(object,
                   truth,
                   target=c(0.2, 0.35),
                   ...){

              ## call the parent method
              start <- callNextMethod(object=object,
                                      truth=truth,
                                      target=target,
                                      ...)

              doseGrid <- object@data[[1]]@doseGrid

              ## dose level most often selected as MTD
              xMostSelected <-
                  match(start@doseMostSelected,
                        table=doseGrid)

              ## fitted toxicity rate at dose most often selected
              fitAtDoseMostSelected <-
                  sapply(object@fit,
                         function(f){
                             f$middle[xMostSelected]
                         })

              ## mean fitted toxicity (average, lower and upper quantiles)
              ## at each dose level
              ## (this is required for plotting)
              meanFitMatrix <- sapply(object@fit,
                                      "[[",
                                      "middle")
              meanFit <- list(truth=
                              truth(doseGrid, ...),
                              average=
                              rowMeans(meanFitMatrix),
                              lower=
                              apply(meanFitMatrix,
                                    1L,
                                    quantile,
                                    0.025),
                              upper=
                              apply(meanFitMatrix,
                                          1L,
                                    quantile,
                                    0.975))

              ## give back an object of class SimulationsSummary,
              ## for which we then define a print / plot method
              ret <- .SimulationsSummary(
                  start,
                  fitAtDoseMostSelected=fitAtDoseMostSelected,
                  meanFit=meanFit)

              return(ret)
          })


##' Summarize the dual-endpoint design simulations, relative to given true
##' dose-toxicity and dose-biomarker curves
##'
##' @param object the \code{\linkS4class{DualSimulations}} object we want to
##' summarize
##' @param trueTox a function which takes as input a dose (vector) and returns the
##' true probability (vector) for toxicity.
##' @param trueBiomarker a function which takes as input a dose (vector) and
##' returns the true biomarker level (vector).
##' @param target the target toxicity interval (default: 20-35\%) used for the
##' computations
##' @param \dots Additional arguments can be supplied here for \code{trueTox}
##' and \code{trueBiomarker}
##' @return an object of class \code{\linkS4class{DualSimulationsSummary}}
##'
##' @example examples/Simulations-method-summary-DualSimulations.R
##' @export
##' @keywords methods
setMethod("summary",
          signature=
          signature(object="DualSimulations"),
          def=
          function(object,
                   trueTox,
                   trueBiomarker,
                   target=c(0.2, 0.35),
                   ...){

              ## call the parent method
              start <- callNextMethod(object=object,
                                      truth=trueTox,
                                      target=target,
                                      ...)

              doseGrid <- object@data[[1]]@doseGrid

              ## dose level most often selected as MTD
              xMostSelected <-
                  match(start@doseMostSelected,
                        table=doseGrid)

              ## fitted biomarker level at dose most often selected
              biomarkerFitAtDoseMostSelected <-
                  sapply(object@fitBiomarker,
                         function(f){
                             f$middleBiomarker[xMostSelected]
                         })

              ## mean fitted biomarker curve (average, lower and upper quantiles)
              ## at each dose level
              ## (this is required for plotting)
              meanBiomarkerFitMatrix <- sapply(object@fitBiomarker,
                                               "[[",
                                               "middleBiomarker")
              meanBiomarkerFit <- list(truth=
                                       trueBiomarker(doseGrid, ...),
                                       average=
                                       rowMeans(meanBiomarkerFitMatrix),
                                       lower=
                                       apply(meanBiomarkerFitMatrix,
                                             1L,
                                             quantile,
                                             0.025),
                                       upper=
                                       apply(meanBiomarkerFitMatrix,
                                             1L,
                                             quantile,
                                             0.975))

              ## give back an object of class DualSimulationsSummary,
              ## for which we then define a print / plot method
              ret <- .DualSimulationsSummary(
                  start,
                  biomarkerFitAtDoseMostSelected=biomarkerFitAtDoseMostSelected,
                  meanBiomarkerFit=meanBiomarkerFit)

              return(ret)
          })

##' A Reference Class to represent sequentially updated reporting objects.
##' @name Report
##' @field object The object from which to report
##' @field df the data frame to which columns are sequentially added
##' @field dfNames the names to which strings are sequentially added
Report <-
    setRefClass("Report",
                fields =
                list(object = "ANY",
                     df = "data.frame",
                     dfNames = "character"),
                methods = list(
                dfSave =
                function(res, name) {
                    df <<- cbind(df, res)
                    dfNames <<- c(dfNames, name)
                    return(res)
                },
                report =
                function(slotName,
                         description,
                         percent=TRUE,
                         digits=0,
                         quantiles=c(0.1, 0.9)) {
                    vals <- slot(object, name=slotName)
                    if(percent)
                    {
                        unit <- " %"
                        vals <- vals * 100
                    } else {
                        unit <- ""
                    }

                    res <- paste(round(mean(vals), digits),
                                 unit,
                                 " (",
                                 paste(round(quantile(vals,
                                                      quantiles,
                                                      na.rm=TRUE),
                                             digits),
                                       unit,
                                       collapse=", ",
                                       sep=""),
                                 ")",
                                 sep="")

                    ## print result to the buffer
                    cat(description, ":",
                        "mean",
                        dfSave(res, slotName),
                        "\n")
                }))


##' Show the summary of the simulations
##'
##' @param object the \code{\linkS4class{GeneralSimulationsSummary}} object we want
##' to print
##' @return invisibly returns a data frame of the results with one row and
##' appropriate column names
##'
##' @export
##' @keywords methods
setMethod("show",
          signature=
          signature(object="GeneralSimulationsSummary"),
          def=
          function(object){

              r <- Report$new(object=object,
                              df=
                              as.data.frame(matrix(nrow=1,
                                                   ncol=0)),
                              dfNames=character())

              cat("Summary of",
                  r$dfSave(object@nsim, "nsim"),
                  "simulations\n\n")

              cat("Target toxicity interval was",
                  r$dfSave(paste(round(object@target * 100),
                               collapse=", "),
                         "targetToxInterval"),
                  "%\n")
              cat("Target dose interval corresponding to this was",
                  r$dfSave(paste(round(object@targetDoseInterval, 1),
                               collapse=", "),
                         "targetDoseInterval"),
                  "\n")
              cat("Intervals are corresponding to",
                  "10 and 90 % quantiles\n\n")

              r$report("nObs",
                       "Number of patients overall",
                       percent=FALSE)
              r$report("nAboveTarget",
                       "Number of patients treated above target tox interval",
                       percent=FALSE)
              r$report("propDLTs",
                       "Proportions of DLTs in the trials")
              r$report("meanToxRisk",
                       "Mean toxicity risks for the patients")
              r$report("doseSelected",
                       "Doses selected as MTD",
                       percent=FALSE, digits=1)
              r$report("toxAtDosesSelected",
                       "True toxicity at doses selected")
              cat("Proportion of trials selecting target MTD:",
                  r$dfSave(object@propAtTarget * 100,
                         "percentAtTarget"),
                  "%\n")
              cat("Dose most often selected as MTD:",
                  r$dfSave(object@doseMostSelected,
                         "doseMostSelected"),
                  "\n")
              cat("Observed toxicity rate at dose most often selected:",
                  r$dfSave(round(object@obsToxRateAtDoseMostSelected * 100),
                         "obsToxRateAtDoseMostSelected"),
                  "%\n")

              ## finally assign names to the df
              ## and return it invisibly
              names(r$df) <- r$dfNames
              invisible(r$df)
          })


##' Show the summary of the simulations
##'
##' @param object the \code{\linkS4class{SimulationsSummary}} object we want
##' to print
##' @return invisibly returns a data frame of the results with one row and
##' appropriate column names
##'
##' @example examples/Simulations-method-show-SimulationsSummary.R
##' @export
##' @keywords methods
setMethod("show",
          signature=
          signature(object="SimulationsSummary"),
          def=
          function(object){

              ## call the parent method
              df <- callNextMethod(object)
              dfNames <- names(df)

              ## start report object
              r <- Report$new(object=object,
                              df=df,
                              dfNames=dfNames)

              ## add one reporting line
              r$report("fitAtDoseMostSelected",
                       "Fitted toxicity rate at dose most often selected")

              ## and return the updated information
              names(r$df) <- r$dfNames
              invisible(r$df)
          })

##' Show the summary of the dual-endpoint simulations
##'
##' @param object the \code{\linkS4class{DualSimulationsSummary}} object we want
##' to print
##' @return invisibly returns a data frame of the results with one row and
##' appropriate column names
##'
##' @example examples/Simulations-method-show-DualSimulationsSummary.R
##' @export
##' @keywords methods
setMethod("show",
          signature=
          signature(object="DualSimulationsSummary"),
          def=
          function(object){

              ## call the parent method
              df <- callNextMethod(object)
              dfNames <- names(df)

              ## start report object
              r <- Report$new(object=object,
                              df=df,
                              dfNames=dfNames)

              ## add one reporting line
              r$report("biomarkerFitAtDoseMostSelected",
                       "Fitted biomarker level at dose most often selected",
                       percent=FALSE,
                       digits=1)

              ## and return the updated information
              names(r$df) <- r$dfNames
              invisible(r$df)
          })


##' Plot summaries of the general simulations
##'
##' Graphical display of the general simulation summary
##'
##' This plot method can be applied to
##' \code{\linkS4class{GeneralSimulationsSummary}} objects in order to
##' summarize them graphically. Possible \code{type}s of plots at the moment
##' are:
##'
##' \describe{
##' \item{nObs}{Distribution of the number of patients in the simulated trials}
##' \item{doseSelected}{Distribution of the final selected doses in the trials.
##' Note that this can include zero entries, meaning that the trial was stopped
##' because all doses in the dose grid appeared too toxic.}
##' \item{propDLTs}{Distribution of the proportion of patients with DLTs in the
##' trials}
##' \item{nAboveTarget}{Distribution of the number of patients treated at doses
##' which are above the target toxicity interval (as specified by the
##' \code{truth} and \code{target} arguments to
##' \code{\link{summary,GeneralSimulations-method}})}
##' }
##' You can specify any subset of these in the \code{type} argument.
##'
##' @param x the \code{\linkS4class{GeneralSimulationsSummary}} object we want
##' to plot from
##' @param y missing
##' @param type the types of plots you want to obtain.
##' @param \dots not used
##' @return A single \code{\link[ggplot2]{ggplot}} object if a single plot is
##' asked for, otherwise a \code{\link{gridExtra}{gTree}} object.
##'
##' @importFrom ggplot2 geom_histogram ggplot aes xlab ylab geom_line
##' scale_linetype_manual scale_colour_manual
##' @importFrom gridExtra arrangeGrob
##' @export
##' @keywords methods
setMethod("plot",
          signature=
          signature(x="GeneralSimulationsSummary",
                    y="missing"),
          def=
          function(x,
                   y,
                   type=
                   c("nObs",
                     "doseSelected",
                     "propDLTs",
                     "nAboveTarget"),
                   ...){

              ## convenience function to make histograms
              myHist <- function(x, description)
              {
                  dat <- data.frame(x=x)
                  ggplot() +
                      geom_histogram(aes(x=x, y=100*..density..),
                                     data=dat, binwidth=1, origin=-0.5) +
                                         xlab(description)+
                                             ylab("Percent")
              }


              ## which plots should be produced?
              type <- match.arg(type,
                                several.ok=TRUE)
              stopifnot(length(type) > 0L)

              ## start the plot list
              plotList <- list()
              plotIndex <- 0L

              ## distribution of overall sample size
              if("nObs" %in% type)
              {
                  plotList[[plotIndex <- plotIndex + 1L]] <-
                      myHist(x=x@nObs,
                             description="Number of patients in total")
              }

              ## distribution of final MTD estimate
              if("doseSelected" %in% type)
              {
                  plotList[[plotIndex <- plotIndex + 1L]] <-
                      myHist(x=x@doseSelected,
                             description="MTD estimate")
              }

              ## distribution of proportion of DLTs
              if("propDLTs" %in% type)
              {
                  plotList[[plotIndex <- plotIndex + 1L]] <-
                      myHist(x=x@propDLTs * 100,
                             description="Proportion of DLTs [%]")
              }

              ## distribution of number of patients treated at too much tox
              if("nAboveTarget" %in% type)
              {
                  plotList[[plotIndex <- plotIndex + 1L]] <-
                      myHist(x=x@nAboveTarget,
                             description="Number of patients above target")
              }

              ## first combine these small plots
              if(length(plotList))
              {
                  ret <-
                      ## if there is only one plot
                      if(identical(length(plotList),
                                   1L))
                      {
                          ## just use that
                          plotList[[1L]]
                      } else {
                          ## multiple plots in this case
                          do.call(gridExtra::arrangeGrob,
                                  plotList)
                      }
              }

              ## then return
              ret
          })


##' Plot summaries of the model-based design simulations
##'
##' Graphical display of the simulation summary
##'
##' This plot method can be applied to \code{\linkS4class{SimulationsSummary}}
##' objects in order to summarize them graphically. Possible \code{type} of
##' plots at the moment are those listed in
##' \code{\link{plot,GeneralSimulationsSummary,missing-method}} plus:
##' \describe{
##' \item{meanFit}{Plot showing the average fitted dose-toxicity curve across
##' the trials, together with 95\% credible intervals, and comparison with the
##' assumed truth (as specified by the \code{truth} argument to
##' \code{\link{summary,Simulations-method}})}
##' }
##' You can specify any subset of these in the \code{type} argument.
##'
##' @param x the \code{\linkS4class{SimulationsSummary}} object we want
##' to plot from
##' @param y missing
##' @param type the types of plots you want to obtain.
##' @param \dots not used
##' @return A single \code{\link[ggplot2]{ggplot}} object if a single plot is
##' asked for, otherwise a \code{\link{gridExtra}{gTree}} object.
##'
##' @importFrom ggplot2 geom_histogram ggplot aes xlab ylab geom_line
##' scale_linetype_manual scale_colour_manual
##' @importFrom gridExtra arrangeGrob
##' 
##' @example examples/Simulations-method-plot-SimulationsSummary.R
##' @export
##' @keywords methods
setMethod("plot",
          signature=
          signature(x="SimulationsSummary",
                    y="missing"),
          def=
          function(x,
                   y,
                   type=
                   c("nObs",
                     "doseSelected",
                     "propDLTs",
                     "nAboveTarget",
                     "meanFit"),
                   ...){

              ## which plots should be produced?
              type <- match.arg(type,
                                several.ok=TRUE)
              stopifnot(length(type) > 0L)

              ## substract the specific plot types for model-based
              ## designs
              typeReduced <- setdiff(type,
                                     "meanFit")

              ## are there more plots from general?
              moreFromGeneral <- (length(typeReduced) > 0)

              ## if so, then produce these plots
              if(moreFromGeneral)
              {
                  ret <- callNextMethod(x=x, y=y, type=typeReduced)
              }

              ## is the meanFit plot requested?
              if("meanFit" %in% type)
              {
                  ## which types of lines do we have?
                  linetype <- c("True toxicity",
                                "Average estimated toxicity",
                                "95% interval for estimated toxicity")

                  ## create the data frame, with
                  ## true tox, average estimated tox, and 95% (lower, upper)
                  ## estimated tox (in percentage) stacked below each other
                  dat <- data.frame(dose=
                                    rep(x@doseGrid, 4L),
                                    group=
                                    rep(1:4, each=length(x@doseGrid)),
                                    linetype=
                                    factor(rep(linetype[c(1, 2, 3, 3)],
                                               each=length(x@doseGrid)),
                                           levels=linetype),
                                    lines=
                                    unlist(x@meanFit) * 100)

                  ## linetypes for the plot
                  lt <- c("True toxicity"=1,
                          "Average estimated toxicity"=1,
                          "95% interval for estimated toxicity"=2)

                  ## colour for the plot
                  col <- c("True toxicity"=1,
                          "Average estimated toxicity"=2,
                          "95% interval for estimated toxicity"=2)

                  ## now create and save the plot
                  thisPlot <- ggplot() +
                      geom_line(aes(x=dose,
                                    y=lines,
                                    group=group,
                                    linetype=linetype,
                                    col=linetype),
                                data=dat)

                  thisPlot <- thisPlot +
                       scale_linetype_manual(values=lt) +
                           scale_colour_manual(values=col) +
                               xlab("Dose level") +
                                   ylab("Probability of DLT [%]")

                  ## add this plot to the bottom
                  ret <-
                      if(moreFromGeneral)
                          gridExtra::arrangeGrob(ret, thisPlot)
                      else
                          thisPlot
              }

              ## then finally plot everything
              ret
          })


##' Plot summaries of the dual-endpoint design simulations
##'
##' This plot method can be applied to \code{\linkS4class{DualSimulationsSummary}}
##' objects in order to summarize them graphically. Possible \code{type} of
##' plots at the moment are those listed in
##' \code{\link{plot,SimulationsSummary,missing-method}} plus:
##' \describe{
##' \item{meanBiomarkerFit}{Plot showing the average fitted dose-biomarker curve across
##' the trials, together with 95\% credible intervals, and comparison with the
##' assumed truth (as specified by the \code{trueBiomarker} argument to
##' \code{\link{summary,DualSimulations-method}})}
##' }
##' You can specify any subset of these in the \code{type} argument.
##'
##' @param x the \code{\linkS4class{DualSimulationsSummary}} object we want
##' to plot from
##' @param y missing
##' @param type the types of plots you want to obtain.
##' @param \dots not used
##' @return A single \code{\link[ggplot2]{ggplot}} object if a single plot is
##' asked for, otherwise a \code{\link{gridExtra}{gTree}} object.
##'
##' @importFrom ggplot2 geom_histogram ggplot aes xlab ylab geom_line
##' scale_linetype_manual scale_colour_manual
##' @importFrom gridExtra arrangeGrob
##' 
##' @example examples/Simulations-method-plot-DualSimulationsSummary.R
##' @export
##' @keywords methods
setMethod("plot",
          signature=
          signature(x="DualSimulationsSummary",
                    y="missing"),
          def=
          function(x,
                   y,
                   type=
                   c("nObs",
                     "doseSelected",
                     "propDLTs",
                     "nAboveTarget",
                     "meanFit",
                     "meanBiomarkerFit"),
                   ...){

              ## which plots should be produced?
              type <- match.arg(type,
                                several.ok=TRUE)
              stopifnot(length(type) > 0L)

              ## substract the specific plot types for dual-endpoint
              ## designs
              typeReduced <- setdiff(type,
                                     "meanBiomarkerFit")

              ## are there more plots from general?
              moreFromGeneral <- (length(typeReduced) > 0)

              ## if so, then produce these plots
              if(moreFromGeneral)
              {
                  ret <- callNextMethod(x=x, y=y, type=typeReduced)
              }

              ## is the meanBiomarkerFit plot requested?
              if("meanBiomarkerFit" %in% type)
              {
                  ## which types of lines do we have?
                  linetype <- c("True biomarker",
                                "Average estimated biomarker",
                                "95% interval for estimated biomarker")

                  ## create the data frame, with
                  ## true biomarker, average estimated biomarker, and 95% (lower, upper)
                  ## estimated biomarker stacked below each other
                  dat <- data.frame(dose=
                                    rep(x@doseGrid, 4L),
                                    group=
                                    rep(1:4, each=length(x@doseGrid)),
                                    linetype=
                                    factor(rep(linetype[c(1, 2, 3, 3)],
                                               each=length(x@doseGrid)),
                                           levels=linetype),
                                    lines=
                                    unlist(x@meanBiomarkerFit))

                  ## linetypes for the plot
                  lt <- c("True biomarker"=1,
                          "Average estimated biomarker"=1,
                          "95% interval for estimated biomarker"=2)

                  ## colour for the plot
                  col <- c("True biomarker"=1,
                          "Average estimated biomarker"=2,
                          "95% interval for estimated biomarker"=2)

                  ## now create and save the plot
                  thisPlot <- ggplot() +
                      geom_line(aes(x=dose,
                                    y=lines,
                                    group=group,
                                    linetype=linetype,
                                    col=linetype),
                                data=dat)

                  thisPlot <- thisPlot +
                       scale_linetype_manual(values=lt) +
                           scale_colour_manual(values=col) +
                               xlab("Dose level") +
                                   ylab("Biomarker level")

                  ## add this plot to the bottom
                  ret <-
                      if(moreFromGeneral)
                          gridExtra::arrangeGrob(ret, thisPlot, heights=c(2/3, 1/3))
                      else
                          thisPlot
              }

              ## then finally plot everything
              ret
          })


## --------------------------------------------------------------------------------------------------------
##' Summarize the simulations, relative to a given truth
##'
##' @param object the \code{\linkS4class{PseudoSimulations}} object we want to
##' summarize
##' @param truth a function which takes as input a dose (vector) and returns the
##' true probability (vector) for toxicity
##' @param targetEndOfTrial the target probability of DLE wanted to achieve at the end of a trial
##' @param targetDuringTrial the target probability of DLE wanted to achieve during a trial
##' 
##' @param \dots Additional arguments can be supplied here for \code{truth}
##' @return an object of class \code{\linkS4class{PseudoSimulationsSummary}}
##' 
##' @example examples/Simulations-method-summarySIMsingle.R
##' @export
##' @keywords methods
setMethod("summary",
          signature=
            signature(object="PseudoSimulations"),
          def=
            function(object,
                     truth,
                     targetEndOfTrial=0.3,
                     targetDuringTrial=0.35,
                     ...){
              ##extract dose grid
              doseGrid <- object@data[[1]]@doseGrid
              
              ##evaluate true DLE at doseGrid
              trueDLE <- truth(doseGrid)
              
              ##Inverse function of the truth function
              inverse = function (f, lower = -100, upper = 100) {
                function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
              }
              
              ##Function to obtain corresponsing dose level given target prob
              TD <- inverse(truth, 0, max(doseGrid))
              
              ##Find the dose corresponding to the target dose during trial
              targetDoseEndOfTrial <- as.numeric(TD(targetEndOfTrial))
              
              
              ##Find the dose corresponding to the target does end of trial
              targetDoseDuringTrial <- as.numeric(TD(targetDuringTrial))
              
              ##what are the levels above target End of Trial?
              xAboveTargetEndOfTrial <- which(trueDLE > targetEndOfTrial)
              
              ##what are the levels above target During Trial?
              xAboveTargetDuringTrial<- which(trueDLE > targetDuringTrial)
              
              
              ##proportion of DLEs in this trial
              propDLE<- sapply(object@data,
                               function(d) {
                                 mean(d@y)
                               })
              ### mean toxicity risk
              meanToxRisk <- sapply(object@data,
                                    function(d){
                                      mean(trueDLE[d@xLevel])
                                    })
              
              ## doses selected for MTD
              doseSelected <- object@doses
              
              ## replace NA by 0
              doseSelected[is.na(doseSelected)] <- 0
              
              ## dose most often selected as MTD
              doseMostSelected <-
                as.numeric(names(which.max(table(doseSelected))))
              xMostSelected <-
                match(doseMostSelected,
                      table=doseGrid)
              
              ## observed toxicity rate at dose most often selected
              ## Note: this does not seem very useful!
              ## Reason: In case of a fine grid, few patients if any
              ## will have been treated at this dose.
              tmp <-
                sapply(object@data,
                       function(d){
                         whichAtThisDose <- which(d@x == doseMostSelected)
                         nAtThisDose <- length(whichAtThisDose)
                         nDLTatThisDose <- sum(d@y[whichAtThisDose])
                         return(c(nAtThisDose=nAtThisDose,
                                  nDLTatThisDose=nDLTatThisDose))
                       })
              
              obsToxRateAtDoseMostSelected <-
                mean(tmp["nDLTatThisDose",]) / mean(tmp["nAtThisDose",])
              
              ## number of patients overall
              nObs <- sapply(object@data,
                             slot,
                             "nObs")
              
              ## number of patients treated above target End of trial
              nAboveTargetEndOfTrial <- sapply(object@data,
                                               function(d){
                                                 sum(d@xLevel %in% xAboveTargetEndOfTrial)
                                               })
              
              ## number of patients treated above target During trial
              nAboveTargetDuringTrial <- sapply(object@data,
                                                function(d){
                                                  sum(d@xLevel %in% xAboveTargetDuringTrial)
                                                })
              
              toxAtDoses <- truth(doseSelected)
              
              
              ## Proportion of trials selecting target TDEndOfTrial and TDDuringTrial
              TDtargetEndOfTrialAtDoseGrid<- doseGrid[max(which(targetDoseEndOfTrial-doseGrid >=0))]
              
              TDtargetDuringTrialAtDoseGrid<- doseGrid[max(which(targetDoseDuringTrial-doseGrid >=0))]
              nsim <- length(object@data)
              
              propAtTargetEndOfTrial <- (length(which(object@doses==TDtargetEndOfTrialAtDoseGrid)))/nsim
              propAtTargetDuringTrial <- (length(which(object@doses==TDtargetDuringTrialAtDoseGrid)))/nsim
              
              ##fitted probDLE at dose most often selected
              ##find names in the fit list (check it is with or without samples)
              FitNames<- sapply(object@fit,names)
              
              
              if ("probDLE" %in% FitNames){
                
                fitAtDoseMostSelected <- sapply(object@fit,
                                                function(f){
                                                  f$probDLE[xMostSelected]
                                                })
                meanFitMatrix <- sapply(object@fit,
                                        "[[",
                                        "probDLE")
                
                meanFit <- list(truth=
                                  truth(doseGrid),
                                average=rowMeans(meanFitMatrix))
              } else {
                
                ## fitted toxicity rate at dose most often selected
                fitAtDoseMostSelected <-
                  sapply(object@fit,
                         function(f){
                           f$middle[xMostSelected]
                         })
                
                ## mean fitted toxicity (average, lower and upper quantiles)
                ## at each dose level
                ## (this is required for plotting)
                meanFitMatrix <- sapply(object@fit,
                                        "[[",
                                        "middle")
                meanFit <- list(truth=
                                  truth(doseGrid),
                                average=
                                  rowMeans(meanFitMatrix),
                                lower=
                                  apply(meanFitMatrix,
                                        1L,
                                        quantile,
                                        0.025),
                                upper=
                                  apply(meanFitMatrix,
                                        1L,
                                        quantile,
                                        0.975))
                
              }
              
              ## give back an object of class GeneralSimulationsSummary,
              ## for which we then define a print / plot method
              ret <- .PseudoSimulationsSummary(
                targetEndOfTrial=targetEndOfTrial,
                targetDoseEndOfTrial=targetDoseEndOfTrial,
                targetDuringTrial=targetDuringTrial,
                targetDoseDuringTrial=targetDoseDuringTrial,
                nsim=length(object@data),
                propDLE=propDLE,
                meanToxRisk=meanToxRisk,
                doseSelected=doseSelected,
                doseMostSelected=doseMostSelected,
                obsToxRateAtDoseMostSelected=obsToxRateAtDoseMostSelected,
                nObs=nObs,
                nAboveTargetEndOfTrial=nAboveTargetEndOfTrial,
                nAboveTargetDuringTrial=nAboveTargetDuringTrial,
                toxAtDosesSelected=toxAtDoses,
                propAtTargetEndOfTrial=propAtTargetEndOfTrial,
                propAtTargetDuringTrial=propAtTargetDuringTrial,
                doseGrid=doseGrid,
                fitAtDoseMostSelected=fitAtDoseMostSelected,
                meanFit=meanFit)
              
              return(ret)
            })
## ========================================================================================================
##' Show the summary of the simulations
##'
##' @param object the \code{\linkS4class{PseudoSimulationsSummary}} object we want
##' to print
##' @return invisibly returns a data frame of the results with one row and
##' appropriate column names
##'
##' @example examples/Simulations-method-showSIMsingle.R
##' @export
##' @keywords methods

setMethod("show",
          signature=
            signature(object="PseudoSimulationsSummary"),
          def=
            function(object){
              
              r <- Report$new(object=object,
                              df=
                                as.data.frame(matrix(nrow=1,
                                                     ncol=0)),
                              dfNames=character())
              cat("Summary of",
                  r$dfSave(object@nsim, "nsim"),
                  "simulations\n\n")
              
              cat("Target prob of DLE End of trial was",
                  r$dfSave(object@targetEndOfTrial * 100,
                           "targetEndOfTrial"),"%\n")
              
              cat("Target dose End of Trial was",
                  r$dfSave(object@targetDoseEndOfTrial,
                           "targetDoseEndOfTrial"),"\n")
              
              cat("Target prob of DLE during trial was",
                  r$dfSave(object@targetDuringTrial * 100,
                           "targetDuringTrial"),"%\n")
              
              cat("Target dose during Trial was",
                  r$dfSave(object@targetDoseDuringTrial,
                           "targetDoseDuringTrial"),"\n")
              
              r$report("nObs",
                       "Number of patients overall",
                       percent=FALSE)
              r$report("nAboveTargetEndOfTrial",
                       "Number of patients treated above target End of Trial",
                       percent=FALSE)
              
              r$report("nAboveTargetDuringTrial",
                       "Number of patients treated above target during trial",
                       percent=FALSE)
              
              r$report("propDLE",
                       "Proportions of DLE in the trials")
              r$report("meanToxRisk",
                       "Mean toxicity risks for the patients")
              r$report("doseSelected",
                       "Doses selected as MTD (TD End of Trial)",
                       percent=FALSE, digits=1)
              r$report("toxAtDosesSelected",
                       "True toxicity at doses selected")
              
              cat("Proportion of trials selecting target End of Trial:",
                  r$dfSave(object@propAtTargetEndOfTrial * 100,
                           "percentAtTarget"),
                  "%\n")
              
              
              cat("Proportion of trials selecting target During Trial:",
                  r$dfSave(object@propAtTargetDuringTrial * 100,
                           "percentAtTarget"),
                  "%\n")
              
              cat("Dose most often selected as MTD (TDEndOfTrial):",
                  r$dfSave(object@doseMostSelected,
                           "doseMostSelected"),
                  "\n")
              cat("Observed toxicity rate at dose most often selected:",
                  r$dfSave(round(object@obsToxRateAtDoseMostSelected * 100),
                           "obsToxRateAtDoseMostSelected"),
                  "%\n")
              r$report("fitAtDoseMostSelected",
                       "Fitted probabilities of DLE at dose most often selected")
              
              
              ## finally assign names to the df
              ## and return it invisibly
              names(r$df) <- r$dfNames
              invisible(r$df)
            })
## -------------------------------------------------------------------------------------------
##' Plot summaries of the pseudo simulations
##'
##' Graphical display of the simulation summary
##'
##' This plot method can be applied to \code{\linkS4class{PseudoSimulationsSummary}}
##' objects in order to summarize them graphically. This can be used when only DLE responses are involved
##' in the simulations. This also applied to results with or without samples generated during the simulations
##'
##' @param x the \code{\linkS4class{PseudoSimulationsSummary}} object we want
##' to plot from
##' @param y missing
##' @param type the types of plots you want to obtain.
##' @param \dots not used
##' @return A single \code{\link[ggplot2]{ggplot}} object if a single plot is
##' asked for, otherwise a \code{\link{gridExtra}{gTree}} object.
##'
##' @importFrom ggplot2 geom_histogram ggplot aes xlab ylab geom_line
##' scale_linetype_manual scale_colour_manual
##' @importFrom gridExtra arrangeGrob
##' 
##' @example examples/Simulations-method-plotSUMsingle.R
##' @export
##' @keywords methods
##' 

setMethod("plot",
          signature=
            signature(x="PseudoSimulationsSummary",
                      y="missing"),
          def=
            function(x,
                     y,
                     type=
                       c("nObs",
                         "doseSelected",
                         "propDLE",
                         "nAboveTargetEndOfTrial",
                         "meanFit"),
                     ...){
              
              ## convenience function to make histograms
              myHist <- function(x, description)
              {
                dat <- data.frame(x=x)
                ggplot() +
                  geom_histogram(aes(x=x, y=100*..density..),
                                 data=dat, binwidth=1, origin=-0.5) +
                  xlab(description)+
                  ylab("Percent")
              }
              
              
              ## which plots should be produced?
              type <- match.arg(type,
                                several.ok=TRUE)
              stopifnot(length(type) > 0L)
              
              ## start the plot list
              plotList <- list()
              plotIndex <- 0L
              
              ## distribution of overall sample size
              if("nObs" %in% type)
              {
                plotList[[plotIndex <- plotIndex + 1L]] <-
                  myHist(x=x@nObs,
                         description="Number of patients in total")
              }
              
              ## distribution of final MTD estimate
              if("doseSelected" %in% type)
              {
                plotList[[plotIndex <- plotIndex + 1L]] <-
                  myHist(x=x@doseSelected,
                         description="MTD estimate")
              }
              
              ## distribution of proportion of DLTs
              if("propDLE" %in% type)
              {
                plotList[[plotIndex <- plotIndex + 1L]] <-
                  myHist(x=x@propDLE * 100,
                         description="Proportion of DLE [%]")
              }
              
              ## distribution of number of patients treated at too much tox
              if("nAboveTargetEndOfTrial" %in% type)
              {
                plotList[[plotIndex <- plotIndex + 1L]] <-
                  myHist(x=x@nAboveTargetEndOfTrial,
                         description="Number of patients above target")
              }
              
              
              ## first combine these small plots
              if(length(plotList))
              {
                ret <-
                  ## if there is only one plot
                  if(identical(length(plotList),
                               1L))
                  {
                    ## just use that
                    plotList[[1L]]
                  } else {
                    ## multiple plots in this case
                    do.call(gridExtra::arrangeGrob,
                            plotList)
                  }
              }
              
              ##the meanFit plot
              
              if ("meanFit" %in% type)
              {## Find if DLE samples are generated in the simulations
                ## by checking if there the lower limits of the 95% Credibility
                ## interval are calculated
                if (!is.null(x@meanFit$lower)) {
                  
                  ## which types of lines do we have?
                  linetype <- c("True toxicity",
                                "Average estimated toxicity",
                                "95% interval for estimated toxicity")
                  ## create the data frame, with
                  ## true tox, average estimated tox, and 95% (lower, upper)
                  ## estimated tox (in percentage) stacked below each other
                  dat <- data.frame(dose=
                                      rep(x@doseGrid, 4L),
                                    group=
                                      rep(1:4, each=length(x@doseGrid)),
                                    linetype=
                                      factor(rep(linetype[c(1, 2, 3, 3)],
                                                 each=length(x@doseGrid)),
                                             levels=linetype),
                                    lines=
                                      unlist(x@meanFit) * 100)
                  
                  ## linetypes for the plot
                  lt <- c("True toxicity"=1,
                          "Average estimated toxicity"=1,
                          "95% interval for estimated toxicity"=2)
                  
                  ## colour for the plot
                  col <- c("True toxicity"=1,
                           "Average estimated toxicity"=2,
                           "95% interval for estimated toxicity"=2)
                  
                  ## now create and save the plot
                  thisPlot <- ggplot() +
                    geom_line(aes(x=dose,
                                  y=lines,
                                  group=group,
                                  linetype=linetype,
                                  col=linetype),
                              data=dat)
                  
                  thisPlot <- thisPlot +
                    scale_linetype_manual(values=lt) +
                    scale_colour_manual(values=col) +
                    xlab("Dose level") +
                    ylab("Probability of DLE [%]")} else {
                      ## which types of lines do we have?
                      linetype <- c("True toxicity",
                                    "Average estimated toxicity")
                      
                      ## create the data frame, with
                      ## true tox, average estimated tox
                      ## estimated tox (in percentage) stacked below each other
                      dat <- data.frame(dose=
                                          rep(x@doseGrid, 2L),
                                        group=
                                          rep(1:2, each=length(x@doseGrid)),
                                        linetype=
                                          factor(rep(linetype[c(1, 2)],
                                                     each=length(x@doseGrid)),
                                                 levels=linetype),
                                        lines=
                                          unlist(x@meanFit) * 100)
                      
                      ## linetypes for the plot
                      lt <- c("True toxicity"=1,
                              "Average estimated toxicity"=1)
                      
                      ## colour for the plot
                      col <- c("True toxicity"=1,
                               "Average estimated toxicity"=2)
                      
                      ## now create and save the plot
                      thisPlot <- ggplot() +
                        geom_line(aes(x=dose,
                                      y=lines,
                                      group=group,
                                      linetype=linetype,
                                      col=linetype),
                                  data=dat)
                      
                      thisPlot <- thisPlot +
                        scale_linetype_manual(values=lt) +
                        scale_colour_manual(values=col) +
                        xlab("Dose level") +
                        ylab("Probability of DLE [%]")}
              }
              
           
              ## then add this plot at the bottom
              ret <- gridExtra::arrangeGrob(ret,thisPlot)
              ret
              
            })
 ## --------------------------------------------------------------------------------------
##' Plot simulations
##'
##' Summarize the simulations with plots
##'
##' This plot method can be applied to \code{\linkS4class{PseudoDualSimulations}}
##' objects in order to summarize them graphically. Possible \code{type}s of
##' plots at the moment are: \describe{ \item{trajectory}{Summary of the
##' trajectory of the simulated trials} \item{dosesTried}{Average proportions of
##' the doses tested in patients} \item{sigma2}{The variance of the efficacy responses}} 
##' You can specify one or both of these in the
##' \code{type} argument.
##'
##' @param x the \code{\linkS4class{PseudoDualSimulations}} object we want
##' to plot from
##' @param y missing
##' @param type the type of plots you want to obtain.
##' @param \dots not used
##' @return A single \code{\link[ggplot2]{ggplot}} object if a single plot is
##' asked for, otherwise a \code{\link{gridExtra}{gTree}} object.
##'
##' @importFrom ggplot2 ggplot geom_step geom_bar aes xlab ylab
##' scale_linetype_manual
##' @importFrom gridExtra arrangeGrob
##' 
##' @example examples/Simulations-method-plotSIMDual.R
##' @export
##' @keywords methods
setMethod("plot",
          signature= 
            signature(x="PseudoDualSimulations",
                      y="missing"),
          def=
            function(x,
                     y,
                     type=
                       c("trajectory",
                         "dosesTried",
                         "sigma2"),
                     ...){
              ## start the plot list
              plotList <- list()
              plotIndex <- 0L
              
              ## which plots should be produced?
              type <- match.arg(type,
                                several.ok=TRUE)
              stopifnot(length(type) > 0L)
              
              ## substract the specific plot types for
              ## dual-endpoint simulation results
              typeReduced <- setdiff(type,
                                     "sigma2")
              
              ## are there more plots from general?
              moreFromGeneral <- (length(typeReduced) > 0)
              
              ## if so, then produce these plots
              if(moreFromGeneral)
              {
                genPlot <- callNextMethod(x=x, y=y, type=typeReduced)
              }
              
              ## now to the specific dual-endpoint plots:
              
              ## Efficacy variance estimates boxplot
              if("sigma2" %in% type)
              {
                ## save the plot
                plotList[[plotIndex <- plotIndex + 1L]] <-
                  qplot(factor(0), y=y, data=data.frame(y=x@sigma2est), geom="boxplot",
                        xlab="", ylab="Efficacy variance estimates") +
                  coord_flip() + scale_x_discrete(breaks=NULL)
              }
              
              
              ## then finally plot everything
              if(identical(length(plotList),
                           0L))
              {
                return(genPlot)
              } else if(identical(length(plotList),
                                  1L))
              {
                ret <- plotList[[1L]]
              } else {
                ret <- do.call(gridExtra::arrangeGrob,
                               plotList)
              }
              
              if(moreFromGeneral)
              {
                ret <- gridExtra::arrangeGrob(genPlot, ret,heights=c(2/3, 1/3))
              }
              
              return(ret)
            })
## ---------------------------------------------------------------------------------
##'
##' This plot method can be applied to \code{\linkS4class{PseudoDualFlexiSimulations}}
##' objects in order to summarize them graphically. Possible \code{type}s of
##' plots at the moment are: \describe{ \item{trajectory}{Summary of the
##' trajectory of the simulated trials} \item{dosesTried}{Average proportions of
##' the doses tested in patients} \item{sigma2}{The variance of the efficacy responses} 
##' \item{sigma2betaW}{The variance of the random walk model}} 
##' You can specify one or both of these in the
##' \code{type} argument.
##'
##' @param x the \code{\linkS4class{PseudoDualFlexiSimulations}} object we want
##' to plot from
##' @param y missing
##' @param type the type of plots you want to obtain.
##' @param \dots not used
##' @return A single \code{\link[ggplot2]{ggplot}} object if a single plot is
##' asked for, otherwise a \code{\link{gridExtra}{gTree}} object.
##'
##' @importFrom ggplot2 ggplot geom_step geom_bar aes xlab ylab
##' scale_linetype_manual
##' @importFrom gridExtra arrangeGrob
##' 
##' @example examples/Simulations-method-plotSIMDualFlexi.R
##' @export
##' @keywords methods
setMethod("plot",
          signature= 
            signature(x="PseudoDualFlexiSimulations",
                      y="missing"),
          def=
            function(x,
                     y,
                     type=
                       c("trajectory",
                         "dosesTried",
                         "sigma2",
                         "sigma2betaW"),
                     ...){
              ## start the plot list
              plotList <- list()
              plotIndex <- 0L
              
              ## which plots should be produced?
              type <- match.arg(type,
                                several.ok=TRUE)
              stopifnot(length(type) > 0L)
              
              ## substract the specific plot types for
              ## dual-endpoint simulation results
              typeReduced <- setdiff(type,"sigma2betaW")
              
              ## are there more plots from general?
              moreFromGeneral <- (length(typeReduced) > 0)
              
              ## if so, then produce these plots
              if(moreFromGeneral)
              {
                genPlot <- callNextMethod(x=x, y=y, type=typeReduced)
              }
              
              ## now to the specific dual-endpoint plots:
              ## random walk model variance estimates boxplot
              
              if("sigma2betaW" %in% type)
              {
                ## save the plot
                plotList[[plotIndex <- plotIndex + 1L]] <-
                  qplot(factor(0), y=y, data=data.frame(y=x@sigma2betaWest), geom="boxplot",
                        xlab="", ylab="Random walk model variance estimates") +
                  coord_flip() + scale_x_discrete(breaks=NULL)
              }
              
              ## then finally plot everything
              if(identical(length(plotList),
                           0L))
              {
                return(genPlot)
              } else if(identical(length(plotList),
                                  1L))
              {
                ret <- plotList[[1L]]
              } else {
                ret <- do.call(gridExtra::arrangeGrob,
                               plotList)
              }
              
              if(moreFromGeneral)
              {
                ret <- gridExtra::arrangeGrob(genPlot, ret,heights=c(2/3, 1/3))
              }
              
              return(ret)
            })

## -----------------------------------------------------------------------------------------
##' Summary for Pseudo Dual responses simulations, relative to a given pseudo DLE and efficacy model 
##' (except the EffFlexi class model)
##' 
##' @param object the \code{\linkS4class{PseudoDualSimulations}} object we want to summarize
##' @param trueDLE a function which takes as input a dose (vector) and returns the true probability (vector)
##' of DLE
##' @param trueEff a function which takes as input a dose (vector) and returns the mean efficacy value(s) (vector).
##' @param targetEndOfTrial the target probability of DLE that are used at the end of a trial. Default at 0.3.
##' @param targetDuringTrial the target probability of DLE that are used during the trial. Default at 0.35.
##' @param \dots Additional arguments can be supplied here for \code{trueDLE} and \code{trueEff}
##' @return an object of class \code{\linkS4class{PseudoDualSimulationsSummary}}
##' 
##' @example examples/Simulations-method-summarySIMDual.R
##' @export
##' @keywords methods
setMethod("summary",
          signature=
            signature(object="PseudoDualSimulations"),
          def=
            function(object,
                     trueDLE,
                     trueEff,
                     targetEndOfTrial=0.3,
                     targetDuringTrial=0.35,
                     ...){
              ##call the parent method
              start <- callNextMethod(object=object,
                                      truth=trueDLE,
                                      targetEndOfTrial=targetEndOfTrial,
                                      targetDuringTrial=targetDuringTrial,
                                      ...)
              doseGrid <- object@data[[1]]@doseGrid
              
              ## ## dose level most often selected as MTD (TDtargetEnd of Trial)
              xMostSelected <-
                match(start@doseMostSelected,
                      table=doseGrid)
              
              ##check if true Eff is a function
              ## check if special case applies
              isTrueEffFx <- is.function(trueEff)
              
              
              if (isTrueEffFx) {
                negtrueGainfn <- function(dose)
                {return(-(trueEff(dose))/(1+(trueDLE(dose)/(1-trueDLE(dose)))))}
                Gstar <- optim(exp(1),negtrueGainfn,method="BFGS")$par
                maxGainValue <- -(optim(exp(1),negtrueGainfn,method="BFGS")$value)
                GstarAtDoseGrid <- doseGrid[max(which(Gstar-doseGrid >= 0 ))]
              } else {
              trueGain <- (trueEff)/(1+(trueDLE(doseGrid)/(1-trueDLE(doseGrid))))
              maxGainValue<-max(trueGain)
              Gstar <- doseGrid[which.max(trueGain)]
              GstarAtDoseGrid <- Gstar
              }
              
              
              
              ##find names in the fit efficacy list (check it is with or without samples)
              FitNames<- sapply(object@fitEff,names)
              if ("ExpEff" %in% FitNames){
                ## fitted efficacy level at dose most often selected
                EffFitAtDoseMostSelected <- sapply(object@fitEff,
                                                   function(f){
                                                     f$ExpEff[xMostSelected]
                                                   })
                meanEffFitMatrix <- sapply(object@fitEff,
                                           "[[",
                                           "ExpEff")
                
                meanEffFit <- list(truth=
                                     trueEff(doseGrid),
                                   average=rowMeans(meanEffFitMatrix))
                
              } else {## fitted efficacy level at dose most often selected
                EffFitAtDoseMostSelected <-
                  sapply(object@fitEff,
                         function(f){
                           f$middle[xMostSelected]
                         })
                
                ## mean fitted  curve (average, lower and upper quantiles)
                ## at each dose level
                ## (this is required for plotting)
                meanEffFitMatrix <- sapply(object@fitEff,
                                           "[[",
                                           "middle")
                
                ## check if special case applies
                
                if (isTrueEffFx) {TRUTHeff<- trueEff(doseGrid)} else {TRUTHeff <- trueEff}
                
                meanEffFit <- list(truth=
                                     TRUTHeff,
                                   average=
                                     rowMeans(meanEffFitMatrix),
                                   lower=
                                     apply(meanEffFitMatrix,
                                           1L,
                                           quantile,
                                           0.025),
                                   upper=
                                     apply(meanEffFitMatrix,
                                           1L,
                                           quantile,
                                           0.975))}
              
              ## give back an object of class PseudoDualSimulationsSummary,
              ## for which we then define a print / plot method
              ret <- .PseudoDualSimulationsSummary(
                start,
                targetGstar=Gstar,
                targetGstarAtDoseGrid=GstarAtDoseGrid,
                EffFitAtDoseMostSelected=EffFitAtDoseMostSelected,
                meanEffFit=meanEffFit)
              
              return(ret)
            })
## --------------------------------------------------------------------------------------------------
##' Summary for Pseudo Dual responses simulations given a pseudo DLE model and the Flexible efficacy model.
##' 
##' @param object the \code{\linkS4class{PseudoDualFlexiSimulations}} object we want to summarize
##' @param trueDLE a function which takes as input a dose (vector) and returns the true probability of DLE (vector)
##' @param trueEff a vector which takes as input the true mean efficacy values at all dose levels (in order)
##' @param targetEndOfTrial the target probability of DLE that are used at the end of a trial. Default at 0.3.
##' @param targetDuringTrial the target probability of DLE that are used during the trial. Default at 0.35.
##' @param \dots Additional arguments can be supplied here for \code{trueDLE} and \code{trueEff}
##' @return an object of class \code{\linkS4class{PseudoDualSimulationsSummary}}
##' 
##' @example examples/Simulations-method-summarySIMDualFlexi.R
##' @export
##' @keywords methods
setMethod("summary",
          signature=
            signature(object="PseudoDualFlexiSimulations"),
          def=
            function(object,
                     trueDLE,
                     trueEff,
                     targetEndOfTrial=0.3,
                     targetDuringTrial=0.35,
                     ...){
              ##call the parent method
              start <- callNextMethod(object=object,
                                      trueDLE=trueDLE,
                                      trueEff=trueEff,
                                      targetEndOfTrial=targetEndOfTrial,
                                      targetDuringTrial=targetDuringTrial,
                                      ...)
              
              
              ## give back an object of class PseudoDualSimulationsSummary,
              ## for which we then define a print / plot method
              ret <- .PseudoDualSimulationsSummary(start)
              
              return(ret)
            })

## ----------------------------------------------------------------------------------------
##' Show the summary of Pseudo Dual simulations summary
##' 
##' @param object the \code{\linkS4class{PseudoDualSimulationsSummary}} object we want to print
##' @return invisibly returns a data frame of the results with one row and appropriate column names
##' 
##' 
##' @example examples/Simulations-method-showSIMDual.R
##' @export
##' @keywords methods
setMethod("show",
signature=
  signature(object="PseudoDualSimulationsSummary"),
def=
  function(object){
    
    ##call the parent method
    df <- callNextMethod(object)
    dfNames <- names(df)
    
    ##start report object
    r <- Report$new(object=object,
                    df=df,
                    dfNames=dfNames)
    
    ##add three reporting lines
    cat("Target Gstar, the dose which gives the maximum gain value was",
        r$dfSave(object@targetGstar,
                     "targetGstar"),"\n")
    cat("Target Gstar at dose Grid was",
        r$dfSave(object@targetGstarAtDoseGrid,
                 "targetGstarAtDoseGrid"),"\n")
            
    r$report("EffFitAtDoseMostSelected",
             "Fitted expected efficacy level at dose most often selected",
              percent=FALSE,
              digits=1)
    ## and return the updated information
    names(r$df) <- r$dfNames
    invisible(r$df)
  })

## --------------------------------------------------------------------------------------------------
##' Plot the summary of Pseudo Dual Simulations summary
##' 
##' This plot method can be applied to \code{\linkS4class{PseudoDualSimulationsSummary}} objects in order
##' to summarize them graphically. Possible \code{type} of plots at the moment are those listed in
##' \code{\link{plot,PseudoSimulationsSummary,missing-method}} plus: 
##' \describe{\item{meanEffFit}{Plot showing the fitted dose-efficacy curve. If no samples are involved, only the
##' average fitted dose-efficacy curve across the trials will be ploted. If samples (DLE and efficacy) are involved, 
##' the average fitted dose-efficacy curve across the trials, together with the 95\% credibility interval; and comparison 
##' with the assumed truth (as specified by the \code{trueEff} argument to 
##' \code{\link{summary,PseudoDualSimulations-method}})}}
##' You can specify any subset of these in the \code{type} argument. 
##' 
##' @param x the \code{\linkS4class{PseudoDualSimulationsSummary}} object we want to plot from
##' @param y missing
##' @param type the types of plots you want to obtain.
##' @param \dots not used
##' @return A single \code{\link[ggplot2]{ggplot}} object if a single plot is
##' asked for, otherwise a \code{\link{gridExtra}{gTree}} object.
##'
##' @importFrom ggplot2 geom_histogram ggplot aes xlab ylab geom_line
##' scale_linetype_manual scale_colour_manual
##' @importFrom gridExtra arrangeGrob
##'
##' @example examples/Simulations-method-plotSUMDual.R
##' @export
##' @keywords methods
setMethod("plot",
          signature=
            signature(x="PseudoDualSimulationsSummary",
                      y="missing"),
          def=
            function(x,
                     y,
                     type=
                       c("nObs",
                         "doseSelected",
                         "propDLE",
                         "nAboveTargetEndOfTrial",
                         "meanFit",
                         "meanEffFit"),
                     ...){
              
              ## which plots should be produced?
              type <- match.arg(type,
                                several.ok=TRUE)
              stopifnot(length(type) > 0L)
              
              ## substract the specific plot types for dual-endpoint
              ## designs
              typeReduced <- setdiff(type,
                                     "meanEffFit")
              
              ## are there more plots from general?
              moreFromGeneral <- (length(typeReduced) > 0)
              
              ## if so, then produce these plots
              if(moreFromGeneral)
              {
                ret <- callNextMethod(x=x, y=y, type=typeReduced)
              }
              
              ## is the meanBiomarkerFit plot requested?
              if("meanEffFit" %in% type)
              { ## Find if Effsamples are generated in the simulations
                ## by checking if there the lower limits of the 95% Credibility
                ## interval are calculated
                if (!is.null(x@meanEffFit$lower)) {
                  ## which types of lines do we have?
                  linetype <- c("True Expected Efficacy",
                                "Average estimated expected efficacy",
                                "95% interval for estimated expected efficacy")
                  
                  ## create the data frame, with
                  ## true biomarker, average estimated expected efficacy, and 95% (lower, upper)
                  ## estimated biomarker stacked below each other
                  dat <- data.frame(dose=
                                      rep(x@doseGrid, 4L),
                                    group=
                                      rep(1:4, each=length(x@doseGrid)),
                                    linetype=
                                      factor(rep(linetype[c(1, 2, 3, 3)],
                                                 each=length(x@doseGrid)),
                                             levels=linetype),
                                    lines=
                                      unlist(x@meanEffFit))
                  
                  ## linetypes for the plot
                  lt <- c("True Expected Efficacy"=1,
                          "Average estimated expected efficacy"=1,
                          "95% interval for estimated expected efficacy"=2)
                  
                  ## colour for the plot
                  col <- c("True Expected Efficacy"=1,
                           "Average estimated expected efficacy"=4,
                           "95% interval for estimated expected efficacy"=4)
                  
                  ## now create and save the plot
                  thisPlot <- ggplot() +
                    geom_line(aes(x=dose,
                                  y=lines,
                                  group=group,
                                  linetype=linetype,
                                  col=linetype),
                              data=dat)
                  
                  thisPlot <- thisPlot +
                    scale_linetype_manual(values=lt) +
                    scale_colour_manual(values=col) +
                    xlab("Dose level") +
                    ylab("Expected Efficacy level")} else {
                      linetype <- c("True Expected Efficacy",
                                    "Average estimated expected efficacy")
                      
                      ## create the data frame, with
                      ## true biomarker, average estimated expected efficacy
                      dat <- data.frame(dose=
                                          rep(x@doseGrid, 2L),
                                        group=
                                          rep(1:2, each=length(x@doseGrid)),
                                        linetype=
                                          factor(rep(linetype[c(1, 2)],
                                                     each=length(x@doseGrid)),
                                                 levels=linetype),
                                        lines=
                                          unlist(x@meanEffFit))
                      
                      ## linetypes for the plot
                      lt <- c("True Expected Efficacy"=1,
                              "Average estimated expected efficacy"=1)
                      
                      ## colour for the plot
                      col <- c("True Expected Efficacy"=1,
                               "Average estimated expected efficacy"=4)
                      
                      ## now create and save the plot
                      thisPlot <- ggplot() +
                        geom_line(aes(x=dose,
                                      y=lines,
                                      group=group,
                                      linetype=linetype,
                                      col=linetype),
                                  data=dat)
                      
                      thisPlot <- thisPlot +
                        scale_linetype_manual(values=lt) +
                        scale_colour_manual(values=col) +
                        xlab("Dose level") +
                        ylab("Expected Efficacy level")
                    }
                
                ## add this plot to the bottom
                ret <-
                  if(moreFromGeneral)
                    gridExtra::arrangeGrob(ret, thisPlot, heights=c(2/3, 1/3))
                else
                  thisPlot
              }
              
              ## then finally plot everything
              ret
            })

 ## ------------------------------------------------------------------------


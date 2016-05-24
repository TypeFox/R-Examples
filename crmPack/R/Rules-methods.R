#####################################################################################
## Author: Daniel Sabanes Bove [sabanesd *a*t* roche *.* com],
##         Wai Yin Yeung [ w *.* yeung1 *a*t* lancaster *.* ac *.* uk]
## Project: Object-oriented implementation of CRM designs
##
## Time-stamp: <[Rules-methods.R] by DSB Die 09/06/2015 21:29>
##
## Description:
## Encapsulate the rule functions in formal methods.
##
## History:
## 07/02/2014   file creation
## 10/07/2015   Adding more classes
###################################################################################

##' @include Model-methods.R
##' @include Samples-class.R
##' @include Rules-class.R
{}

## ============================================================

## --------------------------------------------------
## Find out what is the next best dose
## --------------------------------------------------


##' Find the next best dose
##'
##' Compute the recommended next best dose.
##'
##' This function outputs the next best dose recommendation based on the
##' corresponding rule \code{nextBest}, the posterior \code{samples} from the
##' \code{model} and the underlying \code{data}.
##'
##' @param nextBest The rule, an object of class \code{\linkS4class{NextBest}}
##' @param doselimit The maximum allowed next dose. If this is an empty (length
##' 0) vector, then no dose limit will be applied in the course of dose
##' recommendation calculation, and a corresponding warning is given.
##' @param samples the \code{\linkS4class{Samples}} object
##' @param model The model input, an object of class \code{\linkS4class{Model}}
##' @param data The data input, an object of class \code{\linkS4class{Data}}
##' @param \dots possible additional arguments without method dispatch
##' @return a list with the next best dose (element \code{value})
##' on the grid defined in \code{data}, and a plot depicting this recommendation
##' (element \code{plot})
##'
##' @export
##' @keywords methods
setGeneric("nextBest",
           def=
               function(nextBest, doselimit, samples, model, data, ...){

                   ## there should be no default method,
                   ## therefore just forward to next method!
                   standardGeneric("nextBest")
               },
           valueClass="list")

## --------------------------------------------------
## The MTD method
## --------------------------------------------------

##' @describeIn nextBest Find the next best dose based on the MTD rule
##'
##' @example examples/Rules-method-NextBestMTD.R
##' @importFrom ggplot2 ggplot geom_density xlab ylab xlim aes geom_vline
##' geom_text
setMethod("nextBest",
          signature=
          signature(nextBest="NextBestMTD",
                    doselimit="numeric",
                    samples="Samples",
                    model="Model",
                    data="Data"),
          def=
              function(nextBest, doselimit, samples, model, data, ...){

                  if(identical(length(doselimit), 0L))
                  {
                      warning("doselimit is empty, therefore no dose limit will be applied")
                  }

              ## First, generate the MTD samples.
              mtdSamples <- dose(prob=nextBest@target,
                                 model,
                                 samples)

              ## then derive the next best dose
              mtdEstimate <- nextBest@derive(mtdSamples=mtdSamples)

              ## be sure which doses are ok with respect to maximum
              ## possible dose - if one was specified
              dosesOK <-
                  if(length(doselimit))
                      which(data@doseGrid <= doselimit)
                  else
                      seq_along(data@doseGrid)

              ## but now, round to the next possible grid point
              index <- which.min(abs(data@doseGrid[dosesOK] - mtdEstimate))
              ret <- data@doseGrid[dosesOK][index]

              ## produce plot
              plot1 <- ggplot() +
                      geom_density(data=
                                   data.frame(x=mtdSamples),
                                   aes(x=x),
                                   fill = "grey50", colour = "grey50") +
                          xlab("MTD") + ylab("Posterior density") +
                              xlim(range(data@doseGrid))

              plot1 <- plot1 +
                  geom_vline(xintercept=mtdEstimate, colour="black", lwd=1.1) +
                      geom_text(data=
                                data.frame(x=mtdEstimate),
                                aes(x, 0,
                                    label = "Est", hjust=+1, vjust = +1),
                                colour="black")

              if(length(doselimit))
              {
                  plot1 <- plot1 +
                      geom_vline(xintercept=doselimit, colour="red", lwd=1.1) +
                          geom_text(data=
                                        data.frame(x=doselimit),
                                    aes(x, 0,
                                        label = "Max", hjust = +1, vjust = +1),
                                    colour="red")
              }

              plot1 <- plot1 +
                  geom_vline(xintercept=ret, colour="blue", lwd=1.1) +
                      geom_text(data=
                                data.frame(x=ret),
                                aes(x, 0,
                                    label = "Next", hjust = 0, vjust = +1),
                                colour="blue")

              ## return next best dose and plot
              return(list(value=ret,
                          plot=plot1))
          })

## --------------------------------------------------
## The NCRM method
## --------------------------------------------------

##' @describeIn nextBest Find the next best dose based on the NCRM method
##'
##' @example examples/Rules-method-NextBestNCRM.R
##' @importFrom ggplot2 ggplot geom_bar xlab ylab ylim aes geom_vline
##' geom_hline geom_point
##' @importFrom gridExtra arrangeGrob
setMethod("nextBest",
          signature=
          signature(nextBest="NextBestNCRM",
                    doselimit="numeric",
                    samples="Samples",
                    model="Model",
                    data="Data"),
          def=
              function(nextBest, doselimit, samples, model, data, ...){

                  if(identical(length(doselimit), 0L))
                  {
                      warning("doselimit is empty, therefore no dose limit will be applied")
                  }

              ## first we have to get samples from the dose-tox
              ## curve at the dose grid points.
              probSamples <- matrix(nrow=sampleSize(samples@options),
                                    ncol=data@nGrid)

              ## evaluate the probs, for all samples.
              for(i in seq_len(data@nGrid))
              {
                  ## Now we want to evaluate for the
                  ## following dose:
                  probSamples[, i] <- prob(dose=data@doseGrid[i],
                                           model,
                                           samples)
              }

              ## Now compute probabilities to be in target and
              ## overdose tox interval
              probTarget <-
                  colMeans((probSamples >= nextBest@target[1]) &
                           (probSamples <= nextBest@target[2]))

              probOverdose <-
                  colMeans((probSamples > nextBest@overdose[1]) &
                           (probSamples <= nextBest@overdose[2]))

              ## which doses are eligible after accounting
              ## for maximum possible dose and  discarding overdoses?
              dosesBelowLimit <-
                  if(length(doselimit))
                      (data@doseGrid <= doselimit)
                  else
                      rep(TRUE, length(data@doseGrid))

              dosesOK <- which(dosesBelowLimit &
                               (probOverdose < nextBest@maxOverdoseProb))

              ## check if there are doses that are OK
              if(length(dosesOK))
              {
                  ## what is the recommended dose level?

                  ## if maximum target probability is higher than some numerical
                  ## threshold, then take that level, otherwise stick to the
                  ## maximum level that is OK:
                  doseLevel <-
                      if(max(probTarget[dosesOK]) > 0.05)
                      {
                          which.max(probTarget[dosesOK])
                      } else {
                          which.max(data@doseGrid[dosesOK])
                      }

                  ret <- data@doseGrid[dosesOK][doseLevel]
              } else {
                  ## if none of the doses is OK:
                  doseLevel <- NA
                  ret <- NA
              }

              ## produce plot

              ## first for the target probability
              plot1 <- ggplot() +
                  geom_bar(data=
                           data.frame(x=data@doseGrid,
                                      y=probTarget * 100),
                           aes(x=x, y=y),
                           stat="identity",
                           position="identity",
                           width=1,
                           colour="darkgreen",
                           fill="darkgreen") +
                               xlab("Dose") +
                                   ylab(paste("Target probability [%]")) +
                                       ylim(c(0, 100))

              if(length(doselimit))
              {
                  plot1 <- plot1 +
                      geom_vline(xintercept=doselimit,
                                 lwd=1.1,
                                 lty=2,
                                 colour="black")
              }

              if(length(dosesOK))
              {
                  plot1 <- plot1 +
                      geom_vline(xintercept=data@doseGrid[max(dosesOK)],
                                 lwd=1.1,
                                 lty=2,
                                 colour="red")

                  plot1 <- plot1 +
                      geom_point(data=
                                 data.frame(x=ret,
                                            y=probTarget[dosesOK][doseLevel] *
                                            100 + 0.03),
                                 aes(x=x, y=y),
                                 size=3,
                                 pch=25,
                                 col="red",
                                 bg="red")
              }

              ## second for the overdosing probability
              plot2 <- ggplot() +
                  geom_bar(data=
                           data.frame(x=data@doseGrid,
                                      y=probOverdose * 100),
                           aes(x=x, y=y),
                           stat="identity",
                           position="identity",
                           width=1,
                           colour="red",
                           fill="red") +
                               xlab("Dose") +
                                   ylab("Overdose probability [%]") +
                                       ylim(c(0, 100))

              plot2 <- plot2 +
                  geom_hline(yintercept=nextBest@maxOverdoseProb * 100,
                             lwd=1.1,
                             lty=2,
                             colour="black")

              ## now plot them below each other
              plotJoint <- gridExtra::arrangeGrob(plot1, plot2, nrow=2)

              ## return value and plot
              return(list(value=ret,
                          plot=plotJoint))
          })


##' @describeIn nextBest Find the next best dose based on the NCRM method when
##' two parts trial is used - todo: need an example here for DataParts
##' @example examples/Rules-method-NextBestNCRM-DataParts.R
setMethod("nextBest",
          signature=
          signature(nextBest="NextBestNCRM",
                    doselimit="numeric",
                    samples="Samples",
                    model="Model",
                    data="DataParts"),
          def=
              function(nextBest, doselimit, samples, model, data, ...){

              ## exception when we are in part I or about to start part II!
              if(all(data@part == 1L))
              {
                  ## here we will always propose the highest possible dose
                  ## (assuming that the dose limit came from reasonable
                  ## increments rule, i.e. inrementsRelativeParts)
                  if(identical(length(doselimit), 0L))
                  {
                      stop("doselimit needs to be given for Part I")
                  }

                  return(list(value=doselimit,
                              plot=NULL))
              } else {
                  ## otherwise we will just do the standard thing
                  callNextMethod(nextBest, doselimit, samples, model, data, ...)
              }
          })


## --------------------------------------------------
## The 3+3 method
## --------------------------------------------------

##' @describeIn nextBest Find the next best dose based on the 3+3 method
##' @example examples/Rules-method-NextBestThreePlusThree.R
setMethod("nextBest",
          signature=
          signature(nextBest="NextBestThreePlusThree",
                    doselimit="missing",
                    samples="missing",
                    model="missing",
                    data="Data"),
          def=
          function(nextBest, doselimit, samples, model, data, ...){

              ## split the DLTs into the dose level groups
              dltSplit <- split(data@y,
                                factor(data@x,
                                       levels=data@doseGrid))

              ## number of patients and number of DLTs per dose level group
              nPatients <- sapply(dltSplit, length)
              nDLTs <- sapply(dltSplit, sum)

              ## what was the last dose level tested?
              lastLevel <- tail(data@xLevel, 1)

              ## if there are less than 1/3 DLTs at that level
              if(nDLTs[lastLevel]/nPatients[lastLevel] < 1/3)
              {
                  ## we could escalate, unless this is the highest
                  ## level or the higher level was tried already
                  ## (in which case it was not safe)
                  if((lastLevel == length(data@doseGrid)) ||
                     (nPatients[lastLevel+1] > 0))
                  {
                      nextLevel <- lastLevel
                  } else {
                      nextLevel <- lastLevel + 1
                  }
              } else if(nDLTs[lastLevel]/nPatients[lastLevel] > 1/3) {
                  ## rate here is too high, therefore deescalate
                  nextLevel <- lastLevel - 1
              } else {
                  ## otherwise: rate is 1/3 == 2/6,
                  ## then it depends on the number of patients:
                  ## if more than 3, then deescalate, otherwise stay.
                  nextLevel <-
                      if(nPatients[lastLevel] > 3)
                          lastLevel - 1
                      else
                          lastLevel
              }

              ## do we stop here? only if we have no MTD
              ## or the next level has been tried enough (more than
              ## three patients already)
              stopHere <-
                  if(nextLevel == 0)
                  {
                      TRUE
                  } else {
                      nPatients[nextLevel] > 3
                  }

              ## return value and plot
              return(list(value=
                          if(nextLevel == 0) NA else data@doseGrid[nextLevel],
                          stopHere=stopHere))
          })



## --------------------------------------------------
## The method for the dual endpoint model
## --------------------------------------------------

##' @describeIn nextBest Find the next best dose based on the dual endpoint
##' model
##' @example examples/Rules-method-NextBestDualEndpoint.R
##' @importFrom ggplot2 ggplot geom_bar xlab ylab ylim aes geom_vline
##' geom_hline geom_point
##' @importFrom gridExtra arrangeGrob
setMethod("nextBest",
          signature=
          signature(nextBest="NextBestDualEndpoint",
                    doselimit="numeric",
                    samples="Samples",
                    model="DualEndpoint",
                    data="Data"),
          def=
          function(nextBest, doselimit, samples, model, data, ...){

              if(identical(length(doselimit), 0L))
                  {
                      warning("doselimit is empty, therefore no dose limit will be applied")
                  }

              ## get the biomarker level samples
              ## at the dose grid points.
              biomLevelSamples <- matrix(nrow=sampleSize(samples@options),
                                         ncol=data@nGrid)

              ## evaluate the biomLevels, for all samples.
              for(i in seq_len(data@nGrid))
              {
                  ## Now we want to evaluate for the
                  ## following dose:
                  biomLevelSamples[, i] <- biomLevel(dose=data@doseGrid[i],
                                                     xLevel=i,
                                                     model,
                                                     samples)
              }
              ## biomLevelSamples <- samples@data$betaW


              ## now get samples from the dose-tox
              ## curve at the dose grid points.
              probSamples <- matrix(nrow=sampleSize(samples@options),
                                    ncol=data@nGrid)

              ## evaluate the probs, for all samples.
              for(i in seq_len(data@nGrid))
              {
                  ## Now we want to evaluate for the
                  ## following dose:
                  probSamples[, i] <- prob(dose=data@doseGrid[i],
                                           model,
                                           samples)
              }

              # If there is an 'Emax' parameter, target biomarker level will
              # be relative to 'Emax', otherwise will be relative to the
              # maximum biomarker level achieved in the given dose range.
              if("Emax" %in% names(samples@data)){

                  ## For each sample, look which dose is maximizing the
                  ## simultaneous probability to be in the target biomarker
                  ## range and below overdose toxicity
                  probTarget <- numeric(ncol(biomLevelSamples))
                  probTarget <- sapply(seq(1,ncol(biomLevelSamples)),
                                       function(x){
                                           sum(biomLevelSamples[, x] >= nextBest@target[1]*samples@data$Emax &
                                               biomLevelSamples[, x] <= nextBest@target[2]*samples@data$Emax &
                                               probSamples[, x] <= nextBest@overdose[1]) / nrow(biomLevelSamples)
                                       })
              }else{

                  ## For each sample, look which was the minimum dose giving
                  ## relative target level
                  targetIndex <- apply(biomLevelSamples, 1L,
                                       function(x){
                                           rnx <- range(x)
                                           min(which((x >= nextBest@target[1] * diff(rnx) + rnx[1]) &
                                                     (x <= nextBest@target[2] * diff(rnx) + rnx[1] + 1e-10))
                                              )
                                       })

                  probTarget <- numeric(ncol(biomLevelSamples))
                  tab <- table(targetIndex)
                  probTarget[as.numeric(names(tab))] <- tab
                  probTarget <- probTarget / nrow(biomLevelSamples)
              }

              ## Now compute probabilities to be in
              ## overdose tox interval
              probOverdose <-
                  colMeans((probSamples > nextBest@overdose[1]) &
                           (probSamples <= nextBest@overdose[2]))

              ## which doses are eligible after accounting
              ## for maximum possible dose and discarding overdoses?
              dosesBelowLimit <-
                  if(length(doselimit))
                      (data@doseGrid <= doselimit)
                  else
                      rep(TRUE, length(data@doseGrid))

              dosesOK <- which(dosesBelowLimit &
                               (probOverdose < nextBest@maxOverdoseProb))

              ## check if there are doses that are OK
              if(length(dosesOK))
              {
                  ## what is the recommended dose level?

                  ## if maximum target probability is higher than some numerical
                  ## threshold, then take that level, otherwise stick to the
                  ## maximum level that is OK:
                  doseLevel <-
                      if(max(probTarget[dosesOK]) > 0.05)
                      {
                          which.max(probTarget[dosesOK])
                      } else {
                          which.max(data@doseGrid[dosesOK])
                      }

                  ret <- data@doseGrid[dosesOK][doseLevel]
              } else {
                  ## if none of the doses is OK:
                  doseLevel <- NA
                  ret <- NA
              }

              ## produce plot

              ## first for the target probability
              plot1 <- ggplot() +
                  geom_bar(data=
                           data.frame(x=data@doseGrid,
                                      y=probTarget * 100),
                           aes(x=x, y=y),
                           stat="identity",
                           position="identity",
                           width=1,
                           colour="darkgreen",
                           fill="darkgreen") +
                               xlab("Dose") +
                                   ylab(paste("Target probability [%]")) +
                                       ylim(c(0, 100))

              if(length(doselimit))
              {
                  plot1 <- plot1 +
                      geom_vline(xintercept=doselimit,
                                 lwd=1.1,
                                 lty=2,
                                 colour="black")
              }

              if(length(dosesOK))
              {
                  plot1 <- plot1 +
                      geom_vline(xintercept=data@doseGrid[max(dosesOK)],
                                 lwd=1.1,
                                 lty=2,
                                 colour="red")

                  plot1 <- plot1 +
                      geom_point(data=
                                 data.frame(x=ret,
                                            y=probTarget[dosesOK][doseLevel] *
                                            100 + 0.03),
                                 aes(x=x, y=y),
                                 size=3,
                                 pch=25,
                                 col="red",
                                 bg="red")
              }

              ## second for the overdosing probability
              plot2 <- ggplot() +
                  geom_bar(data=
                           data.frame(x=data@doseGrid,
                                      y=probOverdose * 100),
                           aes(x=x, y=y),
                           stat="identity",
                           position="identity",
                           width=1,
                           colour="red",
                           fill="red") +
                               xlab("Dose") +
                                   ylab("Overdose probability [%]") +
                                       ylim(c(0, 100))

              plot2 <- plot2 +
                  geom_hline(yintercept=nextBest@maxOverdoseProb * 100,
                             lwd=1.1,
                             lty=2,
                             colour="black")

              ## now plot them below each other
              plotJoint <- gridExtra::arrangeGrob(plot1, plot2, nrow=2)

              ## return value and plot
              return(list(value=ret,
                          plot=plotJoint))
          })



## ============================================================


## --------------------------------------------------
## Determine the maximum possible next dose
## --------------------------------------------------

##' Determine the maximum possible next dose
##'
##' Determine the upper limit of the next dose based on the increments rule.
##'
##' This function outputs the maximum possible next dose, based on the
##' corresponding rule \code{increments} and the \code{data}.
##'
##' @param increments The rule, an object of class
##' \code{\linkS4class{Increments}}
##' @param data The data input, an object of class \code{\linkS4class{Data}}
##' @param \dots further arguments
##' @return the maximum possible next dose
##'
##' @export
##' @keywords methods
setGeneric("maxDose",
           def=
           function(increments, data, ...){
               ## there should be no default method,
               ## therefore just forward to next method!
               standardGeneric("maxDose")
           },
           valueClass="numeric")


## --------------------------------------------------
## The maximum allowable relative increments in intervals method
## --------------------------------------------------

##' @describeIn maxDose Determine the maximum possible next dose based on
##' relative increments
##' 
##' @example examples/Rules-method-maxDose-IncrementsRelative.R
setMethod("maxDose",
          signature=
          signature(increments="IncrementsRelative",
                    data="Data"),
          def=
          function(increments, data, ...){
              ## determine what was the last dose
              lastDose <- tail(data@x, 1)

              ## determine in which interval this dose was
              lastInterval <-
                  findInterval(x=lastDose,
                               vec=increments@intervals)

              ## so the maximum next dose is
              ret <-
                  (1 + increments@increments[lastInterval]) *
                      lastDose

              return(ret)
          })


## --------------------------------------------------
## The maximum allowable relative increments, with special rules for
## part 1 and beginning of part 2, method method
## --------------------------------------------------

##' @describeIn maxDose Determine the maximum possible next dose based on
##' relative increments and part 1 and 2
##' @example examples/Rules-method-maxDose-IncrementsRelativeParts.R
setMethod("maxDose",
          signature=
          signature(increments="IncrementsRelativeParts",
                    data="DataParts"),
          def=
          function(increments, data, ...){

              ## determine if there are already cohorts
              ## belonging to part 2:
              alreadyInPart2 <- any(data@part == 2L)

              ## if so, we just call the next higher method
              if(alreadyInPart2)
              {
                  callNextMethod(increments, data, ...)
              } else {
                  ## otherwise we have special rules.

                  ## what dose level (index) has the highest dose
                  ## so far?
                  lastDoseLevel <- match(max(data@x),
                                         data@part1Ladder)

                  ## determine the next maximum dose
                  ret <-
                      if(data@nextPart == 1L)
                      {
                          ## here the next cohort will still be in part 1.
                          ## Therefore we just make one step on the part 1 ladder:
                          data@part1Ladder[lastDoseLevel + 1L]
                      } else {
                          ## the next cohort will start part 2.

                          ## if there was a DLT so far:
                          if(any(data@y == 1L))
                          {
                              data@part1Ladder[lastDoseLevel + increments@dltStart]
                          } else {
                              ## otherwise
                              if(increments@cleanStart > 0)
                              {
                                  ## if we want to start part 2 higher than
                                  ## the last part 1 dose, use usual increments
                                  callNextMethod(increments, data, ...)
                              } else {
                                  ## otherwise
                                  data@part1Ladder[lastDoseLevel + increments@cleanStart]
                              }
                          }
                      }

                  return(ret)
              }
          })


## --------------------------------------------------
## The maximum allowable relative increments in terms of DLTs
## --------------------------------------------------

##' @describeIn maxDose Determine the maximum possible next dose based on
##' relative increments determined by DLTs so far
##' 
##' @example examples/Rules-method-maxDose-IncrementsRelativeDLT.R
setMethod("maxDose",
          signature=
          signature(increments="IncrementsRelativeDLT",
                    data="Data"),
          def=
          function(increments, data, ...){
              ## determine what was the last dose
              lastDose <- tail(data@x, 1)

              ## determine how many DLTs have occurred so far
              dltHappened <- sum(data@y)

              ## determine in which interval this is
              interval <-
                  findInterval(x=dltHappened,
                               vec=increments@DLTintervals)

              ## so the maximum next dose is
              ret <-
                  (1 + increments@increments[interval]) *
                      lastDose

              return(ret)
          })


## ============================================================

## --------------------------------------------------
## "AND" combination of stopping rules
## --------------------------------------------------

##' The method combining two atomic stopping rules
##'
##' @param e1 First \code{\linkS4class{Stopping}} object
##' @param e2 Second \code{\linkS4class{Stopping}} object
##' @return The \code{\linkS4class{StoppingAll}} object
##'
##' @example examples/Rules-method-and-stopping-stopping.R
##' @keywords methods
setMethod("&",
          signature(e1="Stopping",
                    e2="Stopping"),
          def=
          function(e1, e2){
              StoppingAll(list(e1, e2))
          })

##' The method combining a stopping list and an atomic
##'
##' @param e1 \code{\linkS4class{StoppingAll}} object
##' @param e2 \code{\linkS4class{Stopping}} object
##' @return The modified \code{\linkS4class{StoppingAll}} object
##'
##' @example examples/Rules-method-and-stoppingAll-stopping.R
##' @keywords methods
setMethod("&",
          signature(e1="StoppingAll",
                    e2="Stopping"),
          def=
          function(e1, e2){
              e1@stopList <- c(e1@stopList,
                               e2)
              return(e1)
          })

##' The method combining an atomic and a stopping list
##'
##' @param e1 \code{\linkS4class{Stopping}} object
##' @param e2 \code{\linkS4class{StoppingAll}} object
##' @return The modified \code{\linkS4class{StoppingAll}} object
##'
##' @example examples/Rules-method-and-stopping-stoppingAll.R
##' @keywords methods
setMethod("&",
          signature(e1="Stopping",
                    e2="StoppingAll"),
          def=
          function(e1, e2){
              e2@stopList <- c(e1,
                               e2@stopList)
              return(e2)
          })

## --------------------------------------------------
## "OR" combination of stopping rules
## --------------------------------------------------

##' The method combining two atomic stopping rules
##'
##' @param e1 First \code{\linkS4class{Stopping}} object
##' @param e2 Second \code{\linkS4class{Stopping}} object
##' @return The \code{\linkS4class{StoppingAny}} object
##'
##' @aliases |,Stopping,Stopping-method
##' @name or-Stopping-Stopping
##' @example examples/Rules-method-or-stopping-stopping.R
##' @keywords methods
setMethod("|",
          signature(e1="Stopping",
                    e2="Stopping"),
          def=
          function(e1, e2){
              StoppingAny(list(e1, e2))
          })

##' The method combining a stopping list and an atomic
##'
##' @param e1 \code{\linkS4class{StoppingAny}} object
##' @param e2 \code{\linkS4class{Stopping}} object
##' @return The modified \code{\linkS4class{StoppingAny}} object
##'
##' @aliases |,StoppingAny,Stopping-method
##' @name or-Stopping-StoppingAny
##' @example examples/Rules-method-or-stoppingAny-stopping.R
##' @keywords methods
setMethod("|",
          signature(e1="StoppingAny",
                    e2="Stopping"),
          def=
          function(e1, e2){
              e1@stopList <- c(e1@stopList,
                               e2)
              return(e1)
          })

##' The method combining an atomic and a stopping list
##'
##' @param e1 \code{\linkS4class{Stopping}} object
##' @param e2 \code{\linkS4class{StoppingAny}} object
##' @return The modified \code{\linkS4class{StoppingAny}} object
##'
##' @aliases |,Stopping,StoppingAny-method
##' @name or-StoppingAny-Stopping
##' @example examples/Rules-method-or-stopping-stoppingAny.R
##' @keywords methods
setMethod("|",
          signature(e1="Stopping",
                    e2="StoppingAny"),
          def=
          function(e1, e2){
              e2@stopList <- c(e1,
                               e2@stopList)
              return(e2)
          })



## --------------------------------------------------
## Stop the trial?
## --------------------------------------------------

##' Stop the trial?
##'
##' This function returns whether to stop the trial.
##'
##' @param stopping The rule, an object of class
##' \code{\linkS4class{Stopping}}
##' @param dose the recommended next best dose
##' @param samples the \code{\linkS4class{Samples}} object
##' @param model The model input, an object of class \code{\linkS4class{Model}}
##' @param data The data input, an object of class \code{\linkS4class{Data}}
##' @param \dots additional arguments
##'
##' @return logical value: \code{TRUE} if the trial can be stopped, \code{FALSE}
##' otherwise. It should have an attribute \code{message} which gives the reason
##' for the decision.
##'
##' @export
##' @example examples/Rules-method-CombiningStoppingRulesAndOr.R
##' @keywords methods
setGeneric("stopTrial",
           def=
           function(stopping, dose, samples, model, data, ...){
               ## if the recommended next dose is NA,
               ## stop in any case.
               if(is.na(dose))
               {
                   return(structure(TRUE,
                                    message="Recommended next best dose is NA"))
               }

               ## there should be no default method,
               ## therefore just forward to next method!
               standardGeneric("stopTrial")
           },
           valueClass="logical")


## --------------------------------------------------
## Stopping based on multiple stopping rules
## --------------------------------------------------

##' @describeIn stopTrial Stop based on multiple stopping rules
##' @example examples/Rules-method-stopTrial-StoppingList.R
setMethod("stopTrial",
          signature=
          signature(stopping="StoppingList",
                    dose="ANY",
                    samples="ANY",
                    model="ANY",
                    data="ANY"),
          def=
          function(stopping, dose, samples, model, data, ...){
              ## evaluate the individual stopping rules
              ## in the list
              individualResults <-
                if(missing(samples))
                {
                  lapply(stopping@stopList,
                         stopTrial,
                         dose=dose,
                         model=model,
                         data=data,
                         ...)
                } else {
                  lapply(stopping@stopList,
                         stopTrial,
                         dose=dose,
                         samples=samples,
                         model=model,
                         data=data,
                         ...)
                }

              ## summarize to obtain overall result
              overallResult <- stopping@summary(as.logical(individualResults))

              ## retrieve individual text messages,
              ## but let them in the list structure
              overallText <- lapply(individualResults, attr, "message")

              return(structure(overallResult,
                               message=overallText))
          })

## --------------------------------------------------
## Stopping based on fulfillment of all multiple stopping rules
## --------------------------------------------------

##' @describeIn stopTrial Stop based on fulfillment of all multiple stopping
##' rules
##' 
##' @example examples/Rules-method-stopTrial-StoppingAll.R
setMethod("stopTrial",
          signature=
          signature(stopping="StoppingAll",
                    dose="ANY",
                    samples="ANY",
                    model="ANY",
                    data="ANY"),
          def=
          function(stopping, dose, samples, model, data, ...){
              ## evaluate the individual stopping rules
              ## in the list
              individualResults <-
                if(missing(samples))
                {
                  lapply(stopping@stopList,
                         stopTrial,
                         dose=dose,
                         model=model,
                         data=data,
                         ...)
                } else {
                  lapply(stopping@stopList,
                         stopTrial,
                         dose=dose,
                         samples=samples,
                         model=model,
                         data=data,
                         ...)
                }

              ## summarize to obtain overall result
              overallResult <- all(as.logical(individualResults))

              ## retrieve individual text messages,
              ## but let them in the list structure
              overallText <- lapply(individualResults, attr, "message")

              return(structure(overallResult,
                               message=overallText))
          })


## --------------------------------------------------
## Stopping based on fulfillment of any stopping rule
## --------------------------------------------------

##' @describeIn stopTrial Stop based on fulfillment of any stopping rule
##' 
##' @example examples/Rules-method-stopTrial-StoppingAny.R
setMethod("stopTrial",
          signature=
          signature(stopping="StoppingAny",
                    dose="ANY",
                    samples="ANY",
                    model="ANY",
                    data="ANY"),
          def=
          function(stopping, dose, samples, model, data, ...){
              ## evaluate the individual stopping rules
              ## in the list
              individualResults <-
                if(missing(samples))
                {
                  lapply(stopping@stopList,
                         stopTrial,
                         dose=dose,
                         model=model,
                         data=data,
                         ...)
                } else {
                  lapply(stopping@stopList,
                         stopTrial,
                         dose=dose,
                         samples=samples,
                         model=model,
                         data=data,
                         ...)
                }

              ## summarize to obtain overall result
              overallResult <- any(as.logical(individualResults))

              ## retrieve individual text messages,
              ## but let them in the list structure
              overallText <- lapply(individualResults, attr, "message")

              return(structure(overallResult,
                               message=overallText))
          })




## --------------------------------------------------
## Stopping based on number of cohorts near to next best dose
## --------------------------------------------------

##' @describeIn stopTrial Stop based on number of cohorts near to next best dose
##' 
##' @example examples/Rules-method-stopTrial-StoppingCohortsNearDose.R
setMethod("stopTrial",
          signature=
          signature(stopping="StoppingCohortsNearDose",
                    dose="numeric",
                    samples="ANY",
                    model="ANY",
                    data="Data"),
          def=
          function(stopping, dose, samples, model, data, ...){
              ## determine the range where the cohorts must lie in
              lower <- (100 - stopping@percentage) / 100 * dose
              upper <- (100 + stopping@percentage) / 100 * dose

              ## which patients lie there?
              indexPatients <- which((data@x >= lower) & (data@x <= upper))

              ## how many cohorts?
              nCohorts <- length(unique(data@cohort[indexPatients]))

              ## so can we stop?
              doStop <- nCohorts >= stopping@nCohorts

              ## generate message
              text <- paste(nCohorts,
                            " cohorts lie within ",
                            stopping@percentage,
                            "% of the next best dose ",
                            dose,
                            ". This ",
                            ifelse(doStop, "reached", "is below"),
                            " the required ",
                            stopping@nCohorts,
                            " cohorts",
                            sep="")

              ## return both
              return(structure(doStop,
                               message=text))
          })


## -------------------------------------------------------------
## Stopping based on number of patients near to next best dose
## -------------------------------------------------------------

##' @describeIn stopTrial Stop based on number of patients near to next best
##' dose
##' 
##' @example examples/Rules-method-stopTrial-StoppingPatientsNearDose.R
setMethod("stopTrial",
          signature=
          signature(stopping="StoppingPatientsNearDose",
                    dose="numeric",
                    samples="ANY",
                    model="ANY",
                    data="Data"),
          def=
          function(stopping, dose, samples, model, data, ...){
              ## determine the range where the cohorts must lie in
              lower <- (100 - stopping@percentage) / 100 * dose
              upper <- (100 + stopping@percentage) / 100 * dose

              ## how many patients lie there?
              nPatients <- sum((data@x >= lower) & (data@x <= upper))

              ## so can we stop?
              doStop <- nPatients >= stopping@nPatients

              ## generate message
              text <- paste(nPatients,
                            " patients lie within ",
                            stopping@percentage,
                            "% of the next best dose ",
                            dose,
                            ". This ",
                            ifelse(doStop, "reached", "is below"),
                            " the required ",
                            stopping@nPatients,
                            " patients",
                            sep="")

              ## return both
              return(structure(doStop,
                               message=text))
          })

## --------------------------------------------------
## Stopping based on minimum number of cohorts
## --------------------------------------------------

##' @describeIn stopTrial Stop based on minimum number of cohorts
##' 
##' @example examples/Rules-method-stopTrial-StoppingMinCohorts.R
setMethod("stopTrial",
          signature=
          signature(stopping="StoppingMinCohorts",
                    dose="ANY",
                    samples="ANY",
                    model="ANY",
                    data="Data"),
          def=
          function(stopping, dose, samples, model, data, ...){
              ## determine number of cohorts
              nCohorts <- length(unique(data@cohort))

              ## so can we stop?
              doStop <- nCohorts >= stopping@nCohorts

              ## generate message
              text <-
                  paste("Number of cohorts is",
                        nCohorts,
                        "and thus",
                        ifelse(doStop, "reached", "below"),
                        "the prespecified minimum number",
                        stopping@nCohorts)

              ## return both
              return(structure(doStop,
                               message=text))
          })

## --------------------------------------------------
## Stopping based on minimum number of patients
## --------------------------------------------------

##' @describeIn stopTrial Stop based on minimum number of patients
##' 
##' @example examples/Rules-method-stopTrial-StoppingMinPatients.R
setMethod("stopTrial",
          signature=
          signature(stopping="StoppingMinPatients",
                    dose="ANY",
                    samples="ANY",
                    model="ANY",
                    data="Data"),
          def=
          function(stopping, dose, samples, model, data, ...){
              ## so can we stop?
              doStop <- data@nObs >= stopping@nPatients

              ## generate message
              text <-
                  paste("Number of patients is",
                        data@nObs,
                        "and thus",
                        ifelse(doStop, "reached", "below"),
                        "the prespecified minimum number",
                        stopping@nPatients)

              ## return both
              return(structure(doStop,
                               message=text))
          })


## --------------------------------------------------
## Stopping based on probability of target tox interval
## --------------------------------------------------

##' @describeIn stopTrial Stop based on probability of target tox interval
##' 
##' @example examples/Rules-method-stopTrial-StoppingTargetProb.R
setMethod("stopTrial",
          signature=
          signature(stopping="StoppingTargetProb",
                    dose="numeric",
                    samples="Samples",
                    model="Model",
                    data="ANY"),
          def=
          function(stopping, dose, samples, model, data, ...){
              ## first we have to get samples from the dose-tox
              ## curve at the dose.
              probSamples <- prob(dose=dose,
                                  model,
                                  samples)

              ## Now compute probability to be in target interval
              probTarget <-
                  mean((probSamples >= stopping@target[1]) &
                       (probSamples <= stopping@target[2]))

              ## so can we stop?
              doStop <- probTarget >= stopping@prob

              ## generate message
              text <-
                  paste("Probability for target toxicity is",
                        round(probTarget * 100),
                        "% for dose",
                        dose,
                        "and thus",
                        ifelse(doStop, "above", "below"),
                        "the required",
                        round(stopping@prob * 100),
                        "%")

              ## return both
              return(structure(doStop,
                               message=text))
          })


## --------------------------------------------------
## Stopping based on MTD distribution
## --------------------------------------------------

##' @describeIn stopTrial Stop based on MTD distribution
##' 
##' @example examples/Rules-method-stopTrial-StoppingMTDdistribution.R
setMethod("stopTrial",
          signature=
          signature(stopping="StoppingMTDdistribution",
                    dose="numeric",
                    samples="Samples",
                    model="Model",
                    data="ANY"),
          def=
          function(stopping, dose, samples, model, data, ...){
              ## First, generate the MTD samples.

              ## add prior data and samples to the
              ## function environment so that they
              ## can be used.
              mtdSamples <- dose(prob=stopping@target,
                                 model,
                                 samples)

              ## what is the absolute threshold?
              absThresh <- stopping@thresh * dose

              ## what is the probability to be above this dose?
              prob <- mean(mtdSamples > absThresh)

              ## so can we stop?
              doStop <- prob >= stopping@prob

              ## generate message
              text <-
                  paste("Probability of MTD above",
                        round(stopping@thresh * 100),
                        "% of current dose",
                        dose,
                        "is",
                        round(prob * 100),
                        "% and thus",
                        ifelse(doStop, "above", "below"),
                        "the required",
                        round(stopping@prob * 100),
                        "%")

              ## return both
              return(structure(doStop,
                               message=text))
          })


## --------------------------------------------------
## Stopping based on probability of targeting biomarker
## --------------------------------------------------

##' @describeIn stopTrial Stop based on probability of targeting biomarker
##' 
##' @example examples/Rules-method-stopTrial-StoppingTargetBiomarker.R
setMethod("stopTrial",
          signature=
          signature(stopping="StoppingTargetBiomarker",
                    dose="numeric",
                    samples="Samples",
                    model="DualEndpoint",
                    data="ANY"),
          def=
          function(stopping, dose, samples, model, data, ...){
              ## compute the target biomarker prob at this dose

              ## get the biomarker level samples
              ## at the dose grid points.
              biomLevelSamples <- matrix(nrow=sampleSize(samples@options),
                                         ncol=data@nGrid)

              ## evaluate the biomLevels, for all samples.
              for(i in seq_len(data@nGrid))
              {
                  ## Now we want to evaluate for the
                  ## following dose:
                  biomLevelSamples[, i] <- biomLevel(dose=data@doseGrid[i],
                                                     xLevel=i,
                                                     model,
                                                     samples)
              }

              ## If there is an 'Emax' parameter, target biomarker level will
              ## be relative to 'Emax', otherwise will be relative to the
              ## maximum biomarker level achieved in the given dose range.
              if("Emax" %in% names(samples@data)){

                  ## For each sample, look which dose is maximizing the
                  ## simultaneous probability to be in the target biomarker
                  ## range and below overdose toxicity
                  probTarget <- numeric(ncol(biomLevelSamples))
                  probTarget <- sapply(seq(1,ncol(biomLevelSamples)),
                                       function(x){
                                           sum(biomLevelSamples[, x] >= stopping@target[1]*samples@data$Emax &
                                               biomLevelSamples[, x] <= stopping@target[2]*samples@data$Emax &
                                               probSamples[, x] <= stopping@overdose[1]) / nrow(biomLevelSamples)
                                       })
              }else{

                  ## For each sample, look which was the minimum dose giving
                  ## relative target level
                  targetIndex <- apply(biomLevelSamples, 1L,
                                       function(x){
                                           rnx <- range(x)
                                           min(which((x >= stopping@target[1] * diff(rnx) + rnx[1]) &
                                                     (x <= stopping@target[2] * diff(rnx) + rnx[1] + 1e-10))
                                              )
                                       })

                  probTarget <- numeric(ncol(biomLevelSamples))
                  tab <- table(targetIndex)
                  probTarget[as.numeric(names(tab))] <- tab
                  probTarget <- probTarget / nrow(biomLevelSamples)
              }

              ## so for this dose we have:
              probTarget <- probTarget[which(data@doseGrid == dose)]

              ## so can we stop?
              doStop <- probTarget >= stopping@prob

              ## generate message
              text <-
                  paste("Probability for target biomarker is",
                        round(probTarget * 100),
                        "% for dose",
                        dose,
                        "and thus",
                        ifelse(doStop, "above", "below"),
                        "the required",
                        round(stopping@prob * 100),
                        "%")

              ## return both
              return(structure(doStop,
                               message=text))
          })

## --------------------------------------------------
## Stopping when the highest dose is reached
## --------------------------------------------------

##' @describeIn stopTrial Stop when the highest dose is reached
##' 
##' @example examples/Rules-method-stopTrial-StoppingHighestDose.R
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingHighestDose",
                      dose="numeric",
                      samples="ANY",
                      model="ANY",
                      data="Data"),
          def=
            function(stopping, dose, samples, model, data, ...){
              isHighestDose <- (dose == data@doseGrid[data@nGrid])
              return(structure(isHighestDose,
                               message=
                                 paste("Next best dose is", dose, "and thus",
                                       ifelse(isHighestDose, "the",
                                              "not the"),
                                       "highest dose")))
            })

## ============================================================

## --------------------------------------------------
## "MAX" combination of cohort size rules
## --------------------------------------------------

##' "MAX" combination of cohort size rules
##'
##' This function combines cohort size rules by taking
##' the maximum of all sizes.
##'
##' @param \dots Objects of class \code{\linkS4class{CohortSize}}
##' @return the combination as an object of class
##' \code{\linkS4class{CohortSizeMax}}
##'
##' @seealso \code{\link{minSize}}
##' @export
##' @keywords methods
setGeneric("maxSize",
           def=
           function(...){
               ## there should be no default method,
               ## therefore just forward to next method!
               standardGeneric("maxSize")
           },
           valueClass="CohortSizeMax")

##' @describeIn maxSize The method combining cohort size rules by taking maximum
##' @example examples/Rules-method-maxSize.R
setMethod("maxSize",
          "CohortSize",
          def=
          function(...){
              CohortSizeMax(list(...))
          })

## --------------------------------------------------
## "MIN" combination of cohort size rules
## --------------------------------------------------

##' "MIN" combination of cohort size rules
##'
##' This function combines cohort size rules by taking
##' the minimum of all sizes.
##'
##' @param \dots Objects of class \code{\linkS4class{CohortSize}}
##' @return the combination as an object of class
##' \code{\linkS4class{CohortSizeMin}}
##'
##' @seealso \code{\link{maxSize}}
##' @export
##' @keywords methods
setGeneric("minSize",
           def=
           function(...){
               ## there should be no default method,
               ## therefore just forward to next method!
               standardGeneric("minSize")
           },
           valueClass="CohortSizeMin")

##' @describeIn minSize The method combining cohort size rules by taking minimum
##' @example examples/Rules-method-minSize.R
setMethod("minSize",
          "CohortSize",
          def=
          function(...){
              CohortSizeMin(list(...))
          })


## --------------------------------------------------
## Determine the size of the next cohort
## --------------------------------------------------

##' Determine the size of the next cohort
##'
##' This function determines the size of the next cohort.
##'
##' @param cohortSize The rule, an object of class
##' \code{\linkS4class{CohortSize}}
##' @param dose the next dose
##' @param data The data input, an object of class \code{\linkS4class{Data}}
##' @param \dots additional arguments
##'
##' @return the size as integer value
##'
##' @export
##' @keywords methods
setGeneric("size",
           def=
           function(cohortSize, dose, data, ...){
               ## if the recommended next dose is NA,
               ## don't check and return 0
               if(is.na(dose))
               {
                   return(0L)
               }

               ## there should be no default method,
               ## therefore just forward to next method!
               standardGeneric("size")
           },
           valueClass="integer")

## --------------------------------------------------
## The dose range method
## --------------------------------------------------

##' @describeIn size Determine the cohort size based on the range into which the
##' next dose falls into
##' 
##' @example examples/Rules-method-size-CohortSizeRange.R
setMethod("size",
          signature=
          signature(cohortSize="CohortSizeRange",
                    dose="ANY",
                    data="Data"),
          def=
          function(cohortSize, dose, data, ...){

              ## determine in which interval the next dose is
              interval <-
                  findInterval(x=dose,
                               vec=cohortSize@intervals)

              ## so the cohort size is
              ret <- cohortSize@cohortSize[interval]

              return(ret)
          })

## --------------------------------------------------
## The DLT range method
## --------------------------------------------------

##' @describeIn size Determine the cohort size based on the number of DLTs so
##' far
##' 
##' @example examples/Rules-method-size-CohortSizeDLT.R
setMethod("size",
          signature=
          signature(cohortSize="CohortSizeDLT",
                    dose="ANY",
                    data="Data"),
          def=
          function(cohortSize, dose, data, ...){

              ## determine how many DLTs have occurred so far
              dltHappened <- sum(data@y)

              ## determine in which interval this is
              interval <-
                  findInterval(x=dltHappened,
                               vec=cohortSize@DLTintervals)

              ## so the cohort size is
              ret <- cohortSize@cohortSize[interval]

              return(ret)
          })

## --------------------------------------------------
## Size based on maximum of multiple cohort size rules
## --------------------------------------------------

##' @describeIn size Size based on maximum of multiple cohort size rules
##' @example examples/Rules-method-size-CohortSizeMax.R
setMethod("size",
          signature=
          signature(cohortSize="CohortSizeMax",
                    dose="ANY",
                    data="Data"),
          def=
          function(cohortSize, dose, data, ...){
              ## evaluate the individual cohort size rules
              ## in the list
              individualResults <-
                  sapply(cohortSize@cohortSizeList,
                         size,
                         dose=dose,
                         data=data,
                         ...)

              ## summarize to obtain overall result
              overallResult <- max(individualResults)

              return(overallResult)
          })

## --------------------------------------------------
## Size based on minimum of multiple cohort size rules
## --------------------------------------------------

##' @describeIn size Size based on minimum of multiple cohort size rules
##' @example examples/Rules-method-size-CohortSizeMin.R
setMethod("size",
          signature=
          signature(cohortSize="CohortSizeMin",
                    dose="ANY",
                    data="Data"),
          def=
          function(cohortSize, dose, data, ...){
              ## evaluate the individual cohort size rules
              ## in the list
              individualResults <-
                  sapply(cohortSize@cohortSizeList,
                         size,
                         dose=dose,
                         data=data,
                         ...)

              ## summarize to obtain overall result
              overallResult <- min(individualResults)

              return(overallResult)
          })

## --------------------------------------------------
## Constant cohort size
## --------------------------------------------------

##' @describeIn size Constant cohort size
##' @example examples/Rules-method-size-CohortSizeConst.R
setMethod("size",
          signature=
          signature(cohortSize="CohortSizeConst",
                    dose="ANY",
                    data="Data"),
          def=
          function(cohortSize, dose, data, ...){
              return(cohortSize@size)
          })


## --------------------------------------------------
## Cohort size based on the parts
## --------------------------------------------------

##' @describeIn size Cohort size based on the parts
##' @example examples/Rules-method-size-CohortSizeParts.R
setMethod("size",
          signature=
          signature(cohortSize="CohortSizeParts",
                    dose="ANY",
                    data="DataParts"),
          def=
          function(cohortSize, dose, data, ...){
              return(cohortSize@sizes[data@nextPart])
          })


## ================================================================
## The nextBest method based only on DLE responses with samples
## ================================================================
##' @describeIn nextBest Find the next best dose based on the 'NextBestTDsamples'
##' class rule. This a method based only on the DLE responses and for 
##' \code{\linkS4class{LogisticIndepBeta}} model class object involving DLE samples
##' 
##' @importFrom ggplot2 ggplot geom_density xlab ylab xlim aes geom_vline
##' geom_text
##' 
##' @example examples/Rules-method-nextbest_TDsamples.R
##' 
##' @export
##' @keywords methods
setMethod("nextBest",
          signature=
            signature(nextBest="NextBestTDsamples",
                      doselimit="numeric",
                      samples="Samples",
                      model="LogisticIndepBeta",
                      data="Data"),
          def=
            function(nextBest, doselimit, samples, model, data, ...){
              ## First, generate the TDtargetDuringTrial (TDtarget During a Trial) and 
              ## TDtargetEndOfTrial (TDtarget at the EndOfTrial) samples.
              
              TDtargetDuringTrialSamples <- dose(prob=nextBest@targetDuringTrial,
                                                 model,
                                                 samples)
              
              TDtargetEndOfTrialSamples <- dose(prob=nextBest@targetEndOfTrial,
                                                model,
                                                samples)
              
              ## then derive the prior/posterior mean of the above two samples
              
              TDtargetDuringTrialEstimate <- nextBest@derive(TDsamples=TDtargetDuringTrialSamples)
              TDtargetEndOfTrialEstimate <- nextBest@derive(TDsamples=TDtargetEndOfTrialSamples)
              
              ## be sure which doses are ok with respect to maximum
              ## possible dose
              dosesOK <- which(data@doseGrid <= doselimit)
              
              ##Find the index of next dose in the doseGrid
              ##next dose is the dose level closest below the TDtargetDuringTrialEstimate
              index <- suppressWarnings(max(which((signif(TDtargetDuringTrialEstimate,digits=4) - data@doseGrid[dosesOK]) >= 0)))
              ret <- data@doseGrid[dosesOK][index]
              
              
              ##Find the dose level (in doseGrid) closest below the TDtargetEndOfTrialEstimate
              index1 <- suppressWarnings(max(which((signif(TDtargetEndOfTrialEstimate,digits=4) - data@doseGrid[dosesOK]) >= 0)))
              ret1 <- data@doseGrid[dosesOK][index1]
              
              
              ## produce plot
              plot1 <- ggplot() +
                geom_density(data=
                               data.frame(x=TDtargetDuringTrialSamples),
                             aes(x=x),
                             fill = "grey50", colour = "grey50") 
              
              
              plot1 <- plot1 +
                geom_density(data=
                               data.frame(x=TDtargetEndOfTrialSamples),
                             aes(x=x),
                             fill = "grey50", colour = "violet") +
                xlab("TD") + ylab("Posterior density") +
                xlim(range(data@doseGrid))
              
              
              plot1 <- plot1+
                geom_vline(xintercept=TDtargetDuringTrialEstimate, colour="orange", lwd=1.1) +
                annotate("text",label=paste(paste("TD",nextBest@targetDuringTrial*100),"Estimate"),
                         x=TDtargetDuringTrialEstimate,y=0,hjust=-0.1, vjust = -20,size=5,colour="orange")
              
              
              plot1 <- plot1+
                geom_vline(xintercept=TDtargetEndOfTrialEstimate, colour="violet", lwd=1.1) +
                annotate("text",label=paste(paste("TD",nextBest@targetEndOfTrial*100),"Estimate"),
                         x=TDtargetEndOfTrialEstimate,y=0,hjust=-0.1, vjust = -25,size=5,colour="violet")
              
              
              
              if (doselimit > max(data@doseGrid)){maxdoselimit<-max(data@doseGrid)} else {maxdoselimit <-doselimit}
              
              plot1 <- plot1 +
                geom_vline(xintercept=maxdoselimit, colour="red", lwd=1.1) +
                geom_text(data=
                            data.frame(x=maxdoselimit),
                          aes(x, 0,
                              label = "Max", hjust = +1, vjust = -35),
                          colour="red")
              
              plot1 <- plot1 +
                geom_vline(xintercept=ret, colour="blue", lwd=1.1) +
                geom_text(data=
                            data.frame(x=ret),
                          aes(x, 0,
                              label = "Next", hjust = 0.1, vjust = -30),
                          colour="blue")
              
              ## return next best dose and plot
              return(list(nextdose=ret,
                          targetDuringTrial=nextBest@targetDuringTrial,
                          TDtargetDuringTrialEstimate=TDtargetDuringTrialEstimate,
                          targetEndOfTrial=nextBest@targetEndOfTrial,
                          TDtargetEndOfTrialEstimate=TDtargetEndOfTrialEstimate,
                          
                          TDtargetEndOfTrialAtDoseGrid=ret1,
                          plot=plot1))
            })
## -------------------------------------------------------------------------------
## The nextBest method based only on DLE responses without samples
## -----------------------------------------------------------------------------
##' @describeIn nextBest Find the next best dose based on the 'NextBestTD'
##' class rule. This a method based only on the DLE responses and for 
##' \code{\linkS4class{LogisticIndepBeta}} model class object without DLE samples
##' 
##' @param SIM internal command to notify if this method is used within simulations. Default as FALSE 
##' @importFrom ggplot2 ggplot  xlab ylab xlim aes geom_vline
##' geom_text
##' 
##' @example examples/Rules-method-nextbest_TD.R
##' 
##' @export
##' @keywords methods
setMethod("nextBest",
          signature=
            signature(nextBest="NextBestTD",
                      doselimit="numeric",
                      samples="missing",
                      model="LogisticIndepBeta",
                      data="Data"),
          def=
            function(nextBest, doselimit, model, data, SIM=FALSE,...){
              ##Find the target prob During Trial
              targetDuringTrial<- nextBest@targetDuringTrial
              ##Find the target prob End of Trial
              targetEndOfTrial<- nextBest@targetEndOfTrial
              
              mylabel<-targetDuringTrial*100
              label2 <-targetEndOfTrial*100
              
              ## Find the TD30 Estimate and TD(target) Estimate
              
              TDtargetEndOfTrialEstimate <- dose(prob=targetEndOfTrial,
                                                 model)
              
              TDEfourdg<-signif(TDtargetEndOfTrialEstimate,digits=4)
              
              
              TDtargetDuringTrialEstimate<-dose(prob=targetDuringTrial,
                                                model)
              
              TDDfourdg<-signif(TDtargetDuringTrialEstimate,digits=4)
              
              probDLE=prob(dose=data@doseGrid,
                           model=model)
              
              ## be sure which doses are ok with respect to maximum
              ## possible dose
              dosesOK <- which(data@doseGrid <= doselimit)
              
              ##Find the index of next dose in the doseGrid
              ##next dose is the dose level closest below the TDtargetEstimate
              
              index <- suppressWarnings(max(which((TDDfourdg- data@doseGrid[dosesOK]) >= 0)))
              ret <- data@doseGrid[dosesOK][index]
              
              
              ##Find the dose level (in doseGrid) closest below the TD30Estimate
              index <- suppressWarnings(max(which((TDEfourdg - data@doseGrid[dosesOK]) >= 0)))
              retTDE <- data@doseGrid[dosesOK][index]
              
              plotData <- data.frame(dose=data@doseGrid,
                                     probDLE=prob(dose=data@doseGrid,
                                                  model=model))
              
              ##make the plot
              gdata <- with(plotData,
                            data.frame(x=dose,
                                       y=probDLE,
                                       group=rep("Estimated DLE",each=nrow(plotData)),
                                       Type=factor(rep("Estimated DLE",nrow(plotData)),levels="Estimated DLE")))
              
              plot1 <- ggplot(data=gdata, aes(x=x,y=y), group=group) +
                xlab("Dose Levels")+
                ylab(paste("Probability of DLE")) + ylim(c(0,1)) + xlim(c(0,max(data@doseGrid))) +
                geom_line(colour=I("red"), size=1.5)
              
              
              if ((TDDfourdg < min(data@doseGrid))|(TDDfourdg > max(data@doseGrid))) {
                if (SIM==FALSE){
                plot1<-plot1
                print(paste(paste("TD",targetDuringTrial*100),paste("=",paste(TDtargetDuringTrialEstimate," not within dose Grid"))))
                } else {plot1<- plot1}
               } else {plot1 <- plot1+
                  geom_point(data=data.frame(x=TDtargetDuringTrialEstimate,y=targetDuringTrial),aes(x=x,y=y),colour="orange", shape=15, size=8) +
                  annotate("text",label=paste(paste("TD",mylabel),"Estimate"),x=TDtargetDuringTrialEstimate+1,y=targetDuringTrial-0.2,size=5,colour="orange")}
              
              
              if ((TDEfourdg < min(data@doseGrid))|(TDEfourdg > max(data@doseGrid))) {
                if (SIM==FALSE){
                plot1<-plot1
                print(paste(paste("TD",targetEndOfTrial*100),paste("=",paste(TDtargetEndOfTrialEstimate," not within dose Grid"))))
                } else {plot1 <- plot1}
              } else {plot1 <- plot1+
                  geom_point(data=data.frame(x=TDtargetEndOfTrialEstimate,y=0.3),aes(x=x,y=y),colour="violet", shape=16, size=8) +
                  annotate("text",label=paste(paste("TD",label2),"Estimate"),x=TDtargetEndOfTrialEstimate+1,y=targetEndOfTrial-0.1,size=5,colour="violet")}
              
              if (doselimit > max(data@doseGrid)) {maxdoselimit <- max(data@doseGrid)} else {maxdoselimit<-doselimit}
              
              plot1 <- plot1 +
                geom_vline(xintercept=maxdoselimit, colour="brown", lwd=1.1) +
                geom_text(data=
                            data.frame(x=maxdoselimit),
                          aes(x, 0,
                              label = "Max", hjust = +1, vjust = -30),
                          colour="brown")
              
              plot1 <-plot1 +
                geom_vline(xintercept=ret, colour="purple", lwd=1.1) +
                geom_text(data=
                            data.frame(x=ret),
                          aes(x, 0,
                              label = "Next", hjust = 0, vjust = -30),
                          colour="purple")
              
              
              ## return next best dose and plot
              return(list(nextdose=ret,
                          targetDuringTrial=targetDuringTrial,
                          TDtargetDuringTrialEstimate=TDtargetDuringTrialEstimate,
                          targetEndOfTrial=targetEndOfTrial,
                          TDtargetEndOfTrialEstimate=TDtargetEndOfTrialEstimate,
                          TDtargetEndOfTrialatdoseGrid=retTDE,
                          plot=plot1))
            })


## ------------------------------------------------------------------------------------
## the nextBest method based on DLE and efficacy responses without DLE and efficacy samples
## -------------------------------------------------------------------------- ----------
##' @describeIn nextBest for slots \code{nextBest},\code{doselimit}, \code{data} and \code{SIM}. This is 
##' a function to find the next best dose based on the 'NextBestMaxGain'
##' class rule. This a method based on the DLE responses and efficacy responses without DLE and 
##' efficacy samples. 
##' 
##' @param Effmodel the efficacy model of \code{\linkS4class{ModelEff}} class object
##'
##' @importFrom ggplot2 ggplot xlab ylab xlim aes geom_vline
##' geom_text
##' 
##' @example examples/Rules-method-nextbest_MaxGain.R
##' 
##' @export
##' @keywords methods
setMethod("nextBest",
          signature=
            signature(nextBest="NextBestMaxGain",
                      doselimit="numeric",
                      samples="missing",
                      model="ModelTox",
                      data="DataDual"),
          def=
            function(nextBest,doselimit,model,data,Effmodel,SIM=FALSE,...){
              
              stopifnot(is(Effmodel, "ModelEff"))
              
              DuringTrialtargetprob <- nextBest@DLEDuringTrialtarget
              EndOfTrialtargetprob <- nextBest@DLEEndOfTrialtarget
              
              ## Find the TDtarget Estimate for During Trial and End of trial
              
              
              TDtargetEndOfTrialEstimate <- dose(prob=EndOfTrialtargetprob,model=model)
              
              
              TDtargetDuringTrialEstimate<-dose(prob=DuringTrialtargetprob,model=model)
              
              ##Get all prob of DLE at all dose levels
              probDLE=prob(dose=data@doseGrid,
                           model=model)
              
              ##Define gain function
              Gainfun<-function(DOSE){
                -gain(DOSE,DLEmodel=model,Effmodel=Effmodel)
              }
              ##Find the dose which gives the maximum gain
              Gstar<-(optim(min(data@doseGrid),Gainfun, method = "L-BFGS-B", lower=min(data@doseGrid),upper=max(data@doseGrid))$par)
              ##Find the maximum gain value
              
              MaxGain<--(optim(min(data@doseGrid),Gainfun,method = "L-BFGS-B", lower=min(data@doseGrid),upper=max(data@doseGrid))$value)
              ## be sure which doses are ok with respect to maximum
              ## possible dose
              
              dosesOK <- which(data@doseGrid <= doselimit)
              
              ##FIND the next dose which is the minimum between TDtargetDuringTrial and Gstar
              nextdose<-min(TDtargetDuringTrialEstimate,Gstar)
              
              ##Find the dose level in doseGrid closest below nextdose
              
              index <- suppressWarnings(max(which((signif(nextdose,digits=4) - data@doseGrid[dosesOK]) >= 0)))
              
              
              ret <- data@doseGrid[dosesOK][index]
              
              ##Find the dose level in doseGrid closest below TDtargetEndOfTrial
              
              indexE <- suppressWarnings(max(which((signif(TDtargetEndOfTrialEstimate,digits=4) - data@doseGrid[dosesOK]) >= 0)))
              
              
              retE <- data@doseGrid[indexE]
              
              ##Find the dose level in doseGrid closest below TDtargetDuringTrial
              
              indexD <- suppressWarnings(max(which((signif(TDtargetDuringTrialEstimate,digits=4) - data@doseGrid[dosesOK]) >= 0)))
              
              
              retD <- data@doseGrid[indexD]
              
              ##Find the dose level in doseGrid closest below Gstar
              
              Gstarindex <- suppressWarnings(max(which((signif(Gstar,digits=4) - data@doseGrid[dosesOK]) >= 0)))
              
              
              Gstarret <- data@doseGrid[Gstarindex]
              
              
              plotData<-data.frame(dose=rep(data@doseGrid,3),
                                   values=c(prob(dose=data@doseGrid,
                                                 model=model),
                                            ExpEff(dose=data@doseGrid,
                                                   model=Effmodel),
                                            gain(dose=data@doseGrid,
                                                 DLEmodel=model,
                                                 Effmodel=Effmodel)))
              gdata<-with(plotData,
                          data.frame(x=dose,
                                     y=values,
                                     group=c(rep("p(DLE)",length(data@doseGrid)),
                                             rep("Expected Efficacy",length(data@doseGrid)),
                                             rep("Gain",length(data@doseGrid))),
                                     Type=factor("Estimate",levels="Estimate")
                                     
                          ))
              
              plot1 <- ggplot(data=gdata, aes(x=x,y=y))+geom_line(aes(group=group,color=group),size=1.5)+
                ggplot2:::scale_colour_manual(name="curves",values=c("blue","green3","red"))+
                xlab("Dose Level")+ xlim(c(0,max(data@doseGrid)))+
                ylab(paste("Values")) + ylim(c(min(gdata$y),max(gdata$y)))
              
              
              
              if ((signif(TDtargetEndOfTrialEstimate,4) < min(data@doseGrid))|(signif(TDtargetEndOfTrialEstimate,4) > max(data@doseGrid))) {
                if (SIM==FALSE){
                plot1<-plot1
                print(paste(paste("Estimated TD",EndOfTrialtargetprob*100),paste("=",paste(TDtargetEndOfTrialEstimate," not within dose Grid"))))
                } else {plot1 <- plot1} 
              } else {
                  plot1 <-plot1 + geom_point(data=data.frame(x=TDtargetEndOfTrialEstimate,y=EndOfTrialtargetprob),aes(x=x,y=y),colour="violet", shape=16, size=8) +
                    annotate("text",label=paste(paste("TD",EndOfTrialtargetprob*100),"Estimate"),x=TDtargetEndOfTrialEstimate-3,y=0.2,size=5,colour="violet")}
              
              
              
              if ((signif(Gstar,4) < min(data@doseGrid))|(signif(Gstar,4) > max(data@doseGrid))) {
                if (SIM==FALSE){
                plot1<-plot1
                print(paste("Estimated Gstar=",paste(Gstar," not within dose Grid")))} else {plot1 <- plot1} 
              } else {plot1 <- plot1 + 
                  geom_point(data=data.frame(x=Gstar,y=MaxGain),aes(x=x,y=y),colour="green3", shape=17, size=8) +
                  annotate("text",label="Max Gain Estimate",x=Gstar,y=MaxGain-0.1,size=5,colour="green3")}
              
              
              mylabel=format(DuringTrialtargetprob,digits=2)
              
              
              
              if ((signif(TDtargetDuringTrialEstimate,4) < min(data@doseGrid))|(signif(TDtargetDuringTrialEstimate,4) > max(data@doseGrid))) {
                if (SIM==FALSE){
                plot1<-plot1
                print(paste(paste("Estimated TD",DuringTrialtargetprob*100),paste("=",paste(TDtargetDuringTrialEstimate," not within dose Grid"))))
              } else {plot1 <- plot1}
            }  else {
                  plot1 <- plot1+
                    geom_point(data=data.frame(x=signif(TDtargetDuringTrialEstimate,4),y=DuringTrialtargetprob),aes(x=x,y=y),colour="orange", shape=15, size=8) +
                    annotate("text",label=paste(paste("TD",DuringTrialtargetprob*100),"Estimate"),x=TDtargetDuringTrialEstimate+25,
                             y=DuringTrialtargetprob+0.01,size=5,colour="orange")
                }
              
              
              if (doselimit > max(data@doseGrid)) {maxdoselimit <- max(data@doseGrid)} else {maxdoselimit<-doselimit}
              
              plot1 <- plot1 +
                geom_vline(xintercept=maxdoselimit, colour="brown", lwd=1.1) +
                annotate("text",label="Max",x=maxdoselimit-2,y=max(gdata$y),size=5,colour="brown")
              
              
              plot1 <-plot1 +
                geom_vline(xintercept=ret, colour="purple", lwd=1.1) +
                annotate("text",label="Next", x=ret+1, y=max(gdata$y)-0.05,size=5,color="purple")
              
              
              ## return next best dose and plot
              return(list(nextdose=ret,
                          DLEDuringTrialtarget=DuringTrialtargetprob,
                          TDtargetDuringTrialEstimate=TDtargetDuringTrialEstimate,
                          TDtargetDuringTrialAtDoseGrid=retD,
                          DLEEndOfTrialtarget=EndOfTrialtargetprob,
                          TDtargetEndOfTrialEstimate=TDtargetEndOfTrialEstimate,
                          TDtargetEndOfTrialAtDoseGrid=retE,
                          GstarEstimate=Gstar,
                          GstarAtDoseGrid=Gstarret,
                          plot=plot1))
            })
## =====================================================================================
## the nextBest method based on DLE and efficacy responses with DLE and efficacy samples
## -------------------------------------------------------------------------- ----------
##' @describeIn nextBest for slots \code{nextBest},\code{doselimit}, \code{data} and \code{SIM}. This is 
##' a function to find the next best dose based on the 'NextBestMaxGainSamples'
##' class rule. This a method based on the DLE responses and efficacy responses with DLE and 
##' efficacy samples. Effmodel must be of class \code{\linkS4class{Effloglog}} or  
##' \code{\linkS4class{EffFlexi}}.
##' @param Effsamples the efficacy samples of \code{\linkS4class{Samples}} class object
##' 
##' @importFrom ggplot2 ggplot geom_histogram xlab ylab xlim aes geom_vline
##' geom_text
##' @example examples/Rules-method-nextbest_MaxGainSamples.R
setMethod("nextBest",
          signature=
            signature(nextBest="NextBestMaxGainSamples",
                      doselimit="numeric",
                      samples="Samples",
                      model="ModelTox",
                      data="DataDual"),
          
          def=
            function(nextBest,doselimit,samples,model,data,Effmodel,Effsamples,SIM=FALSE,...){
              
              if(is(Effmodel, "Effloglog")) 
                {
              
              ##first get the probDLE samples
              points <- data@doseGrid
              probDLESamples <- matrix(nrow=sampleSize(samples@options),
                                       ncol=length(points))
              
              ## evaluate the probs, for all gain samples.
              for(i in seq_along(points))
              {
                ## Now we want to evaluate for the
                ## following dose:
                probDLESamples[, i] <- prob(dose=points[i],
                                            model=model,
                                            samples=samples)
              }
              probDLE <- apply(probDLESamples,2,FUN=nextBest@TDderive)
              
              DuringTrialtargetprob <- nextBest@DLEDuringTrialtarget
              EndOfTrialtargetprob <- nextBest@DLEEndOfTrialtarget
              
              ## Find the TDtarget samples for During Trial and End of trial
              
              
              TDtargetEndOfTrialSamples <- dose(prob=EndOfTrialtargetprob,
                                                model=model,
                                                samples=samples)
              
              
              TDtargetDuringTrialSamples<-dose(prob=DuringTrialtargetprob,
                                               model=model,
                                               samples=samples)
              
              
              
              ## Find the TDtarget Estimate for During Trial and End of trial
              
              
              TDtargetEndOfTrialEstimate <- nextBest@TDderive(TDtargetEndOfTrialSamples)
              ## Ensure the estimate is within dose range
              #TDtargetEndOfTrialEstimate <- min(TDtargetEndOfTrialEstimate,max(data@doseGrid))
              
              
              TDtargetDuringTrialEstimate<-nextBest@TDderive(TDtargetDuringTrialSamples)
              
              ## Ensure the estimate is within dose range
              #TDtargetDuringTrialEstimate <- min(TDtargetDuringTrialEstimate,max(data@doseGrid))
              
              
              ##we have to get samples from the gain values at all dose levels
              
              ExpEffSamples  <- matrix(nrow=sampleSize(Effsamples@options),
                                       ncol=length(points))
              
              ## evaluate the probs, for all gain samples.
              for(i in seq_along(points))
              {
                ## Now we want to evaluate for the
                ## following dose:
                ExpEffSamples[, i] <- ExpEff(dose=points[i],
                                             model=Effmodel,
                                             samples=Effsamples)
              }
              
              
              ExpEff <- apply(ExpEffSamples,2,FUN=nextBest@Gstarderive)
              
              
              GainSamples <- matrix(nrow=sampleSize(samples@options),
                                    ncol=length(points))
              
              ## evaluate the probs, for all gain samples.
              for(i in seq_along(points))
              {
                ## Now we want to evaluate for the
                ## following dose:
                GainSamples[, i] <- gain(dose=points[i],
                                         DLEmodel=model,
                                         DLEsamples=samples, 
                                         Effmodel=Effmodel,
                                         Effsamples=Effsamples)
              }
              
              ##Find the maximum gain value samples
              MaxGainSamples <- apply(GainSamples,1,max)
              
              ##Obtain Gstar samples, samples for the dose level which gives the maximum gain value
              IndexG <- apply(GainSamples,1,which.max)
              GstarSamples <- data@doseGrid[IndexG]
              
              ##Obtain the Gstar estimate which is the 50th percentile of the Gstar samples
              Gstar <- nextBest@Gstarderive(GstarSamples)
              ##Ensure the estimate is within dose range
              
              #Gstar <- min(Gstar,max(data@doseGrid))
              
              gainvalues <- apply(GainSamples,2,FUN=nextBest@Gstarderive)
              
              dosesOK <- which(data@doseGrid <= doselimit)
              
              ##FIND the next dose which is the minimum between TDtargetDuringTrial and Gstar
              nextdose<-min(TDtargetDuringTrialEstimate,Gstar)
              
              ##Find the dose level in doseGrid closest below nextdose
              
              index <- suppressWarnings(max(which((signif(nextdose,digits=4) - data@doseGrid[dosesOK]) >= 0)))
              
              
              ret <- data@doseGrid[dosesOK][index]
              
              ##Find the dose level in doseGrid closest below TDtargetEndOfTrial
              
              indexE <- suppressWarnings(max(which((signif(TDtargetEndOfTrialEstimate,digits=4) - data@doseGrid[dosesOK]) >= 0)))
              
              
              retE <- data@doseGrid[indexE]
              
              ##Find the dose level in doseGrid closest below TDtargetDuringTrial
              
              indexD <- suppressWarnings(max(which((signif(TDtargetDuringTrialEstimate,digits=4) - data@doseGrid[dosesOK]) >= 0)))
              
              
              retD <- data@doseGrid[indexD]
              
              ##Find the dose level in doseGrid closest below Gstar
              
              Gstarindex <- suppressWarnings(max(which((signif(Gstar,digits=4) - data@doseGrid[dosesOK]) >= 0)))
              
              
              Gstarret <- data@doseGrid[Gstarindex]
              
              
              
              plotData<-data.frame(dose=rep(data@doseGrid,3),
                                   values=c(probDLE,
                                            ExpEff,
                                            gainvalues))
              
              gdata<-with(plotData,
                          data.frame(x=dose,
                                     y=values,
                                     group=c(rep("p(DLE)",length(data@doseGrid)),
                                             rep("Expected Efficacy",length(data@doseGrid)),
                                             rep("Gain",length(data@doseGrid))),
                                     Type=factor("Estimate",levels="Estimate")
                                     
                          ))
              
              
              ## produce plot
              plot1 <- ggplot() +
                geom_histogram(data=
                                 data.frame(x=GstarSamples),
                               aes(x=x),
                               fill = "darkgreen", colour = "green3", binwidth = 25) +
                xlab("Gstar")+ xlim(c(0,max(data@doseGrid)))+
                ylab("Posterior density")
              
              if (signif(TDtargetDuringTrialEstimate,4) < min(data@doseGrid)|signif(TDtargetDuringTrialEstimate,4) > max(data@doseGrid)) {
                if (SIM==FALSE){
                plot1<-plot1 
                print(paste(paste("Estimated TD",DuringTrialtargetprob*100),paste("=",paste(TDtargetDuringTrialEstimate," not within dose Grid"))))
                } else {plot1 <- plot1}
              } else {
                  plot1 <- plot1+
                    geom_vline(xintercept=TDtargetDuringTrialEstimate, colour="orange", lwd=1.1) +
                    annotate("text",label=paste(paste("TD",DuringTrialtargetprob*100),"Estimate"),
                             x=TDtargetDuringTrialEstimate,y=0,hjust=-0.1, vjust = -20,size=5,colour="orange")}
              
              if (signif(TDtargetEndOfTrialEstimate,4) < min(data@doseGrid)|signif(TDtargetEndOfTrialEstimate,4) > max(data@doseGrid)) {
                if (SIM==FALSE){
                plot1<-plot1 
                print(paste(paste("Estimated TD",EndOfTrialtargetprob*100),paste("=",paste(TDtargetEndOfTrialEstimate," not within dose Grid"))))
                } else {plot1 <- plot1}
              } else {
                  plot1 <- plot1+
                    geom_vline(xintercept=TDtargetEndOfTrialEstimate, colour="violet", lwd=1.1) +
                    annotate("text",label=paste(paste("TD",EndOfTrialtargetprob*100),"Estimate"),
                             x=TDtargetEndOfTrialEstimate,y=0,hjust=-0.1, vjust = -25,size=5,colour="violet")}
              
              if (signif(Gstar,4) < min(data@doseGrid)|signif(Gstar,4) > max(data@doseGrid)) {
                if (SIM==FALSE){
                plot1<-plot1
                print(paste("Estimated Gstar=",paste(Gstar," not within dose Grid")))
                } else {plot1 <- plot1}
              } else {           
                  plot1 <- plot1+
                    geom_vline(xintercept=Gstar, colour="green", lwd=1.1) +
                    annotate("text",label=" Gstar Estimate",
                             x=Gstar,y=0,hjust=-0.1, vjust = -25,size=5,colour="green")}
              
              
              if (doselimit > max(data@doseGrid)){maxdoselimit<-max(data@doseGrid)} else {maxdoselimit <-doselimit}
              
              plot1 <- plot1 +
                geom_vline(xintercept=maxdoselimit, colour="red", lwd=1.1) +
                geom_text(data=
                            data.frame(x=maxdoselimit),
                          aes(x, 0,
                              label = "Max", hjust = +1, vjust = -35),
                          colour="red")
              
              plot1 <- plot1 +
                geom_vline(xintercept=ret, colour="blue", lwd=1.1) +
                geom_text(data=
                            data.frame(x=ret),
                          aes(x, 0,
                              label = "Next", hjust = 0.1, vjust = -30),
                          colour="blue")
              
              
              
              ## return next best dose and plot
              return(list(nextdose=ret,
                          DLEDuringTrialtarget=DuringTrialtargetprob,
                          TDtargetDuringTrialEstimate=TDtargetDuringTrialEstimate,
                          TDtargetDuringTrialAtDoseGrid=retD,
                          DLEEndOfTrialtarget=EndOfTrialtargetprob,
                          TDtargetEndOfTrialEstimate=TDtargetEndOfTrialEstimate,
                          TDtargetEndOfTrialAtDoseGrid=retE,
                          GstarEstimate=Gstar,
                          GstarAtDoseGrid=Gstarret,
                          plot=plot1))
              } else if(is(Effmodel, "EffFlexi")) {
              
                ##first get the probDLE samples
                points <- data@doseGrid
                probDLESamples <- matrix(nrow=sampleSize(samples@options),
                                         ncol=length(points))
                
                ## evaluate the probs, for all gain samples.
                for(i in seq_along(points))
                {
                  ## Now we want to evaluate for the
                  ## following dose:
                  probDLESamples[, i] <- prob(dose=points[i],
                                              model=model,
                                              samples=samples)
                }
                probDLE <- apply(probDLESamples,2,FUN=nextBest@TDderive)
                
                DuringTrialtargetprob <- nextBest@DLEDuringTrialtarget
                EndOfTrialtargetprob <- nextBest@DLEEndOfTrialtarget
                
                ## Find the TDtarget samples for During Trial and End of trial
                
                
                TDtargetEndOfTrialSamples <- dose(prob=EndOfTrialtargetprob,
                                                  model=model,
                                                  samples=samples)
                
                
                TDtargetDuringTrialSamples<-dose(prob=DuringTrialtargetprob,
                                                 model=model,
                                                 samples=samples)
                
                
                
                ## Find the TDtarget Estimate for During Trial and End of trial
                
                
                TDtargetEndOfTrialEstimate <- nextBest@TDderive(TDtargetEndOfTrialSamples)
                ## Ensure the estimate is within dose range
                #TDtargetEndOfTrialEstimate <- min(TDtargetEndOfTrialEstimate,max(data@doseGrid))
                
                
                TDtargetDuringTrialEstimate<-nextBest@TDderive(TDtargetDuringTrialSamples)
                
                ## Ensure the estimate is within dose range
                #TDtargetDuringTrialEstimate <- min(TDtargetDuringTrialEstimate,max(data@doseGrid))
                
                
                ##we have to get samples from the gain values at all dose levels
                
                ExpEffsamples <- Effsamples@data$ExpEff
                
                ExpEff <- apply(ExpEffsamples,2,FUN=nextBest@Gstarderive)
                
                
                GainSamples <- matrix(nrow=sampleSize(samples@options),
                                      ncol=length(points))
                
                ## evaluate the probs, for all gain samples.
                for(i in seq_along(points))
                {
                  ## Now we want to evaluate for the
                  ## following dose:
                  GainSamples[, i] <- gain(dose=points[i],
                                           DLEmodel=model,
                                           DLEsamples=samples, 
                                           Effmodel=Effmodel,
                                           Effsamples=Effsamples)
                }
                
                ##Find the maximum gain value samples
                MaxGainSamples <- apply(GainSamples,1,max)
                
                ##Obtain Gstar samples, samples for the dose level which gives the maximum gain value
                IndexG <- apply(GainSamples,1,which.max)
                GstarSamples <- data@doseGrid[IndexG]
                
                ##Obtain the Gstar estimate which is the 50th percentile of the Gstar samples
                Gstar <- nextBest@Gstarderive(GstarSamples)
                ##Ensure the estimate is within dose range
                
                #Gstar <- min(Gstar,max(data@doseGrid))
                
                gainvalues <- apply(GainSamples,2,FUN=nextBest@Gstarderive)
                
                dosesOK <- which(data@doseGrid <= doselimit)
                
                ##FIND the next dose which is the minimum between TDtargetDuringTrial and Gstar
                nextdose<-min(TDtargetDuringTrialEstimate,Gstar)
                
                ##Find the dose level in doseGrid closest below nextdose
                
                index <- suppressWarnings(max(which((signif(nextdose,digits=4) - data@doseGrid[dosesOK]) >= 0)))
                
                
                ret <- data@doseGrid[dosesOK][index]
                
                ##Find the dose level in doseGrid closest below TDtargetEndOfTrial
                
                indexE <- suppressWarnings(max(which((signif(TDtargetEndOfTrialEstimate,digits=4) - data@doseGrid[dosesOK]) >= 0)))
                
                
                retE <- data@doseGrid[indexE]
                
                ##Find the dose level in doseGrid closest below TDtargetDuringTrial
                
                indexD <- suppressWarnings(max(which((signif(TDtargetDuringTrialEstimate,digits=4) - data@doseGrid[dosesOK]) >= 0)))
                
                
                retD <- data@doseGrid[indexD]
                
                ##Find the dose level in doseGrid closest below Gstar
                
                Gstarindex <- suppressWarnings(max(which((signif(Gstar,digits=4) - data@doseGrid[dosesOK]) >= 0)))
                
                
                Gstarret <- data@doseGrid[Gstarindex]
                
                
                
                plotData<-data.frame(dose=rep(data@doseGrid,3),
                                     values=c(probDLE,
                                              ExpEff,
                                              gainvalues))
                
                gdata<-with(plotData,
                            data.frame(x=dose,
                                       y=values,
                                       group=c(rep("p(DLE)",length(data@doseGrid)),
                                               rep("Expected Efficacy",length(data@doseGrid)),
                                               rep("Gain",length(data@doseGrid))),
                                       Type=factor("Estimate",levels="Estimate")
                                       
                            ))
                
                
                ## produce plot
                plot1 <- ggplot() +
                  geom_histogram(data=
                                   data.frame(x=GstarSamples),
                                 aes(x=x),
                                 fill = "darkgreen", colour = "green3", binwidth = 25) +
                  xlab("Gstar")+ xlim(c(0,max(data@doseGrid)))+
                  ylab("Posterior density")
                
                if (signif(TDtargetDuringTrialEstimate,4) < min(data@doseGrid)|signif(TDtargetDuringTrialEstimate,4) > max(data@doseGrid)) {
                  if (SIM==FALSE){
                  plot1<-plot1
                  print(paste(paste("Estimated TD",DuringTrialtargetprob*100),paste("=",paste(TDtargetDuringTrialEstimate," not within dose Grid"))))
                  } else {plo1 <- plot1}
                } else {
                    plot1 <- plot1+
                      geom_vline(xintercept=TDtargetDuringTrialEstimate, colour="orange", lwd=1.1) +
                      annotate("text",label=paste(paste("TD",DuringTrialtargetprob*100),"Estimate"),
                               x=TDtargetDuringTrialEstimate,y=0,hjust=-0.1, vjust = -20,size=5,colour="orange")}
                
                if (signif(TDtargetEndOfTrialEstimate,4) < min(data@doseGrid)|signif(TDtargetEndOfTrialEstimate,4) > max(data@doseGrid)) {
                  if (SIM==FALSE){
                  plot1<-plot1 
                  print(paste(paste("Estimated TD",EndOfTrialtargetprob*100),paste("=",paste(TDtargetEndOfTrialEstimate," not within dose Grid"))))
                  } else {plot1 <- plot1}
                } else {
                    plot1 <- plot1+
                      geom_vline(xintercept=TDtargetEndOfTrialEstimate, colour="violet", lwd=1.1) +
                      annotate("text",label=paste(paste("TD",EndOfTrialtargetprob*100),"Estimate"),
                               x=TDtargetEndOfTrialEstimate,y=0,hjust=-0.1, vjust = -25,size=5,colour="violet")}
                
                if (signif(Gstar,4) < min(data@doseGrid)|signif(Gstar,4) > max(data@doseGrid)) {
                  if (SIM==FALSE){
                  plot1<-plot1
                  print(paste("Estimated Gstar=",paste(Gstar," not within dose Grid")))
                  } else {plot1 <- plot1}
                } else {           
                    plot1 <- plot1+
                      geom_vline(xintercept=Gstar, colour="green", lwd=1.1) +
                      annotate("text",label=" Gstar Estimate",
                               x=Gstar,y=0,hjust=+0.6, vjust = -25,size=5,colour="green")}
                
                
                if (doselimit > max(data@doseGrid)){maxdoselimit<-max(data@doseGrid)} else {maxdoselimit <-doselimit}
                
                plot1 <- plot1 +
                  geom_vline(xintercept=maxdoselimit, colour="red", lwd=1.1) +
                  geom_text(data=
                              data.frame(x=maxdoselimit),
                            aes(x, 0,
                                label = "Max", hjust = +1, vjust = -35),
                            colour="red")
                
                plot1 <- plot1 +
                  geom_vline(xintercept=ret, colour="blue", lwd=1.1) +
                  geom_text(data=
                              data.frame(x=ret),
                            aes(x, 0,
                                label = "Next", hjust = 0.1, vjust = -30),
                            colour="blue")
                
                
                
                ## return next best dose and plot
                return(list(nextdose=ret,
                            DLEDuringTrialtarget=DuringTrialtargetprob,
                            TDtargetDuringTrialEstimate=TDtargetDuringTrialEstimate,
                            TDtargetDuringTrialAtDoseGrid=retD,
                            DLEEndOfTrialtarget=EndOfTrialtargetprob,
                            TDtargetEndEstimate=TDtargetEndOfTrialEstimate,
                            TDtargetEndOfTrialAtDoseGrid=retE,
                            GstarEstimate=Gstar,
                            GstarAtDoseGrid=Gstarret,
                            plot=plot1))
                
            } else stop("Effmodel needs to be of class Effloglog or EffFlexi")
              
              })
 
## ------------------------------------------------------------------------------------------------
## Stopping based on a target ratio of the upper to the lower 95% credibility interval
## ------------------------------------------------------------------------------------------------
##' @describeIn stopTrial Stop based on 'StoppingTDCIRatio' class when 
##' reaching the target ratio of the upper to the lower 95% credibility 
##' interval of the estimate (TDtargetEndOfTrial). This is a stopping rule which incorporate only 
##' DLE responses and DLE samples are given
##' 
##' @example examples/Rules-method-stopTrialCITDsamples.R
##' 
##' @export
##' @keywords methods
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingTDCIRatio",
                      dose="ANY",
                      samples="Samples",
                      model="ModelTox",
                      data="ANY"),
          def=
            function(stopping,dose,samples,model,data,...){
              
              targetEndOfTrial <- stopping@targetEndOfTrial
              ##check id targetEndOfTrial is a probability
              stopifnot(is.probability(targetEndOfTrial))
              
              ## find the TDtarget End of Trial samples
              TDtargetEndOfTrialSamples <- dose(prob=targetEndOfTrial,
                                                model=model,
                                                samples=samples)
              
              ##Find the upper and lower limit of the 95% credibility interval
              CI <- quantile(TDtargetEndOfTrialSamples, prob=c(0.025,0.975))
              
              ##The ratio of the upper to the lower 95% credibility interval
              ratio <- as.numeric(CI[2]/CI[1])
              
              
              ##so can we stop? 
              doStop <- ratio <= stopping@targetRatio
              ##generate messgae
              text <- paste("95% CI is (",CI[1], "," , CI[2], "), Ratio =",ratio, "is " , ifelse(doStop,"is less than or equal to","greater than"),
                            "targetRatio =", stopping@targetRatio)
              ##return both
              return(structure(doStop,
                               messgae=text))
            })

## ----------------------------------------------------------------------------------------------
## Stopping based on a target ratio of the upper to the lower 95% credibility interval
## ------------------------------------------------------------------------------------------------
##' @describeIn stopTrial Stop based on 'StoppingTDCIRatio' class
##' when reaching the target ratio of the upper to the lower 95% credibility 
##' interval of the estimate (TDtargetEndOfTrial). This is a stopping rule which incorporate only 
##' DLE responses and no DLE samples are involved
##' @example examples/Rules-method-stopTrialCITD.R
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingTDCIRatio",
                      dose="ANY",
                      samples="missing",
                      model="ModelTox",
                      data="ANY"),
          def=
            function(stopping,dose,model,data,...){
              targetEndOfTrial <- stopping@targetEndOfTrial
              
              ##check if targetEndOfTrial is a probability
              stopifnot(is.probability(targetEndOfTrial))
              
              ## find the TDtarget End of Trial
              TDtargetEndOfTrial <- dose(prob=targetEndOfTrial,
                                         model=model)
              
              ##Find the variance of the log of the TDtargetEndOfTrial(eta)
              M1 <- matrix(c(-1/(model@phi2), - (log(targetEndOfTrial/(1-targetEndOfTrial))-model@phi1)/(model@phi2)^2),1,2)
              M2 <- model@Pcov
              
              varEta <- M1%*%M2%*%t(M1)
              
              ##Find the upper and lower limit of the 95% credibility interval
              CI <- c()
              CI[2] <- exp(log(TDtargetEndOfTrial) + 1.96* sqrt(varEta))
              CI[1] <- exp(log(TDtargetEndOfTrial) - 1.96* sqrt(varEta))
              
              ##The ratio of the upper to the lower 95% credibility interval
              ratio <- as.numeric(CI[2]/CI[1])
              
              
              ##so can we stop? 
              doStop <- ratio <= stopping@targetRatio
              ##generate messgae
              text <- paste("95% CI is (",CI[1], "," , CI[2], "), Ratio =", ratio, "is " , ifelse(doStop,"is less than or equal to","greater than"),
                            "targetRatio =", stopping@targetRatio)
              ##return both
              return(structure(doStop,
                               messgae=text))
            })

## --------------------------------------------------------------------------------------------------
## Stopping based on a target ratio of the upper to the lower 95% credibility interval
## ------------------------------------------------------------------------------------------------
##' @describeIn stopTrial Stop based on reaching the target ratio of the upper to the lower 95% credibility 
##' interval of the estimate (the minimum of Gstar and TDtargetEndOfTrial). This is a stopping rule which 
##' incorporate DLE and efficacy responses and DLE and efficacy samples are also used.
##' 
##' @param TDderive the function which derives from the input, a vector of the posterior samples called 
##' \code{TDsamples} of the dose
##' which has the probability of the occurrence of DLE equals to either the targetDuringTrial or
##' targetEndOfTrial, the final next best TDtargetDuringTrial (the dose with probability of the 
##' occurrence of DLE equals to the targetDuringTrial)and TDtargetEndOfTrial estimate.
##' @param Effmodel the efficacy model of \code{\linkS4class{ModelEff}} class object
##' @param Effsamples the efficacy samples of \code{\linkS4class{Samples}} class object
##' @param Gstarderive the function which derives from the input, a vector of the posterior Gstar (the dose
##' which gives the maximum gain value) samples 
##' called \code{Gstarsamples}, the final next best Gstar estimate.
##' 
##' @example examples/Rules-method-stopTrialCIMaxGainSamples.R
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingGstarCIRatio",
                      dose="ANY",
                      samples="Samples",
                      model="ModelTox",
                      data="DataDual"),
          def=
            function(stopping,dose,samples,model,data,TDderive, Effmodel,Effsamples,Gstarderive,...){
              targetEndOfTrial <- stopping@targetEndOfTrial
              
              ##checks
              stopifnot(is.probability(targetEndOfTrial))
              stopifnot(is(Effmodel,"ModelEff"))
              stopifnot(is(Effsamples,"Samples"))
              stopifnot(is.function(TDderive))
              stopifnot(is.function(Gstarderive))
              
              ## find the TDtarget End of Trial samples
              TDtargetEndOfTrialSamples <- dose(prob=targetEndOfTrial,
                                                model=model,
                                                samples=samples)
              ##Find the TDtarget End of trial estimate
              TDtargetEndOfTrialEstimate <- TDderive(TDtargetEndOfTrialSamples)
              
              ##Find the gain value samples then the GstarSamples
              points <- data@doseGrid
              
              GainSamples <- matrix(nrow=sampleSize(samples@options),
                                    ncol=length(points))
              
              ## evaluate the probs, for all gain samples.
              for(i in seq_along(points))
              {
                ## Now we want to evaluate for the
                ## following dose:
                GainSamples[, i] <- gain(dose=points[i],
                                         DLEmodel=model,
                                         DLEsamples=samples, 
                                         Effmodel=Effmodel,
                                         Effsamples=Effsamples)
              }
              
              ##Find the maximum gain value samples
              MaxGainSamples <- apply(GainSamples,1,max)
              
              ##Obtain Gstar samples, samples for the dose level which gives the maximum gain value
              IndexG <- apply(GainSamples,1,which.max)
              GstarSamples <- data@doseGrid[IndexG]
              
              ##Find the Gstar estimate
              
              Gstar <- Gstarderive(GstarSamples)
              
              ## Find which is smaller (TDtargetEndOfTrialEstimate or Gstar)
              
              if (TDtargetEndOfTrialEstimate <= Gstar){
                
                ##Find the upper and lower limit of the 95% credibility interval
                CI <- quantile(TDtargetEndOfTrialSamples, prob=c(0.025,0.975)) 
                chooseTD <- TRUE} else {
                  
                  CI <- quantile(GstarSamples, prob=c(0.025,0.975))
                  chooseTD <- FALSE}
              
              ##The ratio of the upper to the lower 95% credibility interval
              ratio <- as.numeric(CI[2]/CI[1]) 
              
              
              
              ##so can we stop? 
              doStop <- ratio <= stopping@targetRatio
              ##generate messgae
              text <- paste(ifelse(chooseTD,"TD30 estimate","Gstar estimate"), "is smaller with 95% CI (", CI[1], ",", CI[2], 
                            ") and its ratio =", 
                            ratio, "is " , ifelse(doStop,"is less than or equal to","greater than"),
                            "targetRatio =", stopping@targetRatio)
              ##return both
              return(structure(doStop,
                               messgae=text))
            })

## -----------------------------------------------------------------------------------------------
## Stopping based on a target ratio of the upper to the lower 95% credibility interval
## ------------------------------------------------------------------------------------------------
##' @describeIn stopTrial Stop based on reaching the target ratio of the upper to the lower 95% credibility 
##' interval of the estimate (the minimum of Gstar and TDtargetEndOfTrial). This is a stopping rule which 
##' incorporate DLE and efficacy responses without DLE and efficacy samples involved.
##' @example examples/Rules-method-stopTrialCIMaxGain.R
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingGstarCIRatio",
                      dose="ANY",
                      samples="missing",
                      model="ModelTox",
                      data="DataDual"),
          def=
            function(stopping,dose,model,data,Effmodel,...){
              
              targetEndOfTrial <- stopping@targetEndOfTrial
              
              ##checks
              stopifnot(is.probability(targetEndOfTrial))
              stopifnot(is(Effmodel,"ModelEff"))
              
              
              ## find the TDtarget End of Trial
              TDtargetEndOfTrial <- dose(prob=targetEndOfTrial,
                                         model=model)
              
              ##Find the dose with maximum gain value
              Gainfun<-function(DOSE){
                -gain(DOSE,DLEmodel=model,Effmodel=Effmodel)
              }
              Gstar<-(optim(min(data@doseGrid),Gainfun,method = "L-BFGS-B",lower=min(data@doseGrid),upper=max(data@doseGrid))$par)
              MaxGain<--(optim(min(data@doseGrid),Gainfun,method = "L-BFGS-B",lower=min(data@doseGrid),upper=max(data@doseGrid))$value)
              logGstar <- log(Gstar)
              
              
              if (Gstar <= TDtargetEndOfTrial){
                chooseTD <- FALSE
                ##From paper
                
                meanEffGstar <- Effmodel@theta1+Effmodel@theta2*log(logGstar)
                
                denom <- (model@phi2)*(meanEffGstar)*(1+logGstar*model@phi2)
                
                dgphi1 <- -(meanEffGstar*logGstar*model@phi2-Effmodel@theta2)/denom
                
                dgphi2 <- -((meanEffGstar)*logGstar+meanEffGstar*(logGstar)^2*model@phi2-Effmodel@theta2*logGstar)/denom
                
                dgtheta1 <- -(logGstar*model@phi2)/denom
                
                dgtheta2 <- -(logGstar*exp(model@phi1+model@phi2*logGstar)*model@phi2*log(logGstar)-1-exp(model@phi1+model@phi2*logGstar))/denom
                
                #DLEPRO <- exp(model@phi1+model@phi2*logGstar)
                
                #dgphi1 <- Effmodel@theta2*DLEPRO - logGstar*model@phi2*meanEffGstar*DLEPRO
                
                #dgphi2 <- logGstar*DLEPRO *(Effmodel@theta2-(meanEffGstar)+model@phi2)
                
                #dgtheta1 <- -logGstar*DLEPRO*model@phi2
                
                #dgtheta2 <- 1+DLEPRO-logGstar*DLEPRO*model@phi2*log(logGstar)
                
                deltaG <- matrix(c(dgphi1,dgphi2,dgtheta1,dgtheta2),4,1)
                
                
                ##Find the variance of the log Gstar
                ##First find the covariance matrix of all the parameters, phi1, phi2, theta1 and theta2
                ## such that phi1 and phi2 and independent of theta1 and theta2
                emptyMatrix <- matrix(0,2,2)
                covBETA <-  cbind(rbind(model@Pcov,emptyMatrix),rbind(emptyMatrix,Effmodel@Pcov))
                varlogGstar <- t(deltaG)%*%covBETA%*%deltaG
                
                
                
                ##Find the upper and lower limit of the 95% credibility interval
                CI <- c()
                CI[2] <- exp(logGstar + 1.96* sqrt(varlogGstar))
                CI[1] <- exp(logGstar - 1.96* sqrt(varlogGstar))
                
                ##The ratio of the upper to the lower 95% credibility interval
                ratio <- as.numeric(CI[2]/CI[1]) } else {
                  chooseTD <- TRUE 
                  
                  ##Find the variance of the log of the TDtargetEndOfTrial(eta)
                  M1 <- matrix(c(-1/(model@phi2), - (log(targetEndOfTrial/(1-targetEndOfTrial))-model@phi1)/(model@phi2)^2),1,2)
                  M2 <- model@Pcov
                  
                  varEta <- M1%*%M2%*%t(M1)
                  
                  ##Find the upper and lower limit of the 95% credibility interval
                  CI <- c()
                  CI[2] <- exp(log(TDtargetEndOfTrial) + 1.96* sqrt(varEta))
                  CI[1] <- exp(log(TDtargetEndOfTrial) - 1.96* sqrt(varEta))
                  
                  ##The ratio of the upper to the lower 95% credibility interval
                  ratio <- as.numeric(CI[2]/CI[1])
                } 
              
              
              ##so can we stop? 
              doStop <- ratio <= stopping@targetRatio
              ##generate messgae
              text <- paste(ifelse(chooseTD,"TD30 estimate","Gstar estimate"), "is smaller with 95% CI (", CI[1], ",", CI[2], 
                            ") and its ratio =", ratio, "is " , ifelse(doStop,"is less than or equal to","greater than"),
                            "targetRatio =", stopping@targetRatio)
              ##return both
              return(structure(doStop,
                               message=text))
            })


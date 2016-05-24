################################################################################
################################################################################
## author Till Junge <till.junge@altermail.ch>                                ##
##                                                                            ##
## Copyright (c) UNIL (Universite de Lausanne)                                ##
## NCCR - LIVES (National Centre of Competence in Research "LIVES -           ##
## Overcoming vulnerability: life course perspectives",                       ##
## <http://www.lives-nccr.ch/>)                                               ##
##                                                                            ##
## spacom is free software: you can redistribute it and/or modify it under    ##
## the terms of the GNU General Public License as published by the Free       ##
## Software Foundation, either version 2 of the License or any later version. ##
##                                                                            ##
## spacom is distributed in the hope that it will be useful, but WITHOUT ANY  ##
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  ##
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      ##
## details, see <http://www.gnu.org/licenses/>.                               ##
################################################################################
################################################################################


## provides a ResampleMLSpawExactObject needed for the
## analysis and performs all consistency checks
MakeResampleMLSpawExactObject <-function(individual.level.data,
                                         context.id,
                                         formula,
                                         precise.data,
                                         confidence.intervals,
                                         nb.resamples,
                                         individual.sample.seed,
                                         verbose,
                                         obj=NULL) {

  ## create an empty ResampleMLSpawExactObject to fill, using the
  ## parent class MLSpawExactObject
  if (is.null(obj)){
    obj <- new("ResampleMLSpawExactObject")
  }
  obj <-
    MakeMLSpawExactObject(individual.level.data=individual.level.data,
                          context.id=context.id,
                          formula=formula,
                          precise.data=precise.data,
                          verbose=verbose,
                          obj)

  ## makes sure the number of resamples is a positive integer
  obj@nb.resamples <- checkNbResamples(nb.resamples)

  ## make sure the sample seed is ok
  obj@individual.sample.seed <- checkSeed(individual.sample.seed)

  ## make sure the confidence intervals are correctly defined and contain the
  ## median
  obj@percentiles <-
    checkConfidenceIntervals(confidence.intervals, obj@nb.resamples)


  return(obj)
}

## ## Definition of an output object for describeResampledIndividualContext
printSraweObject <- function(x, is.print=TRUE) {
  cat("New model result \n")
  if (!("contextual.sample.seed" %in% slotNames(x))) {
    cat("Stratified resampling for precise contextual values\n")
  } else {
    cat("Stratified resampling for contextual variables aggregated with error")
  }
  cat(paste("\nNumber of resamples: ", x@nb.resamples, "\n"))
  cat("\nLinear mixed model fit by REML:\n")
  if (is.print) {
    print(x@model.fit)
  } else {
    show(x@model.fit)
  }
  cat("\nFixed effects:\n")
  if (is.print) {
    print(x@fixed)
  } else {
    show(x@fixed)
  }
  cat("\nRandom effects:\n")
  if (is.print) {
    print(x@random.var)
  } else {
    show(x@random.var)
  }
  cat("\nStandardised fixed effects:\n")
  if (TRUE==all.equal(dim(x@betas), c(0,0))) {
    cat("None\n")
  } else {
    if (is.print) {
      print(x@betas)
    } else {
      show(x@betas)
    }
  }
  cat("\nslots for direct acces to data:\n  usage: <object.name>@<slot.name>\n")
  for (name in slotNames(x)) {
    if (is.print) {
      print(name)
    } else {
      show(name)
    }
  }
}

setMethod("print", signature="ResampleMLSpawOutput",
          definition=printSraweObject)
setMethod("show", signature="ResampleMLSpawOutput",
          definition=function(object) { printSraweObject(object, is.print=FALSE)})



PerformResampleMLSpawExact <- function(obj, ...){
  ## creation of an output data to bundle all the return data
  output.obj <- new("ResampleMLSpawOutput")
  output.obj@nb.resamples <- obj@nb.resamples

  helper.functions <-
    randomEffectsClosure(nrow(obj@precise.data),
                         obj@nb.resamples)
  getRanefs <- helper.functions$getRanefs
  fillMatrices <- helper.functions$fillMatrices
  getTmpList <- helper.functions$getTmpList
  getBetaResults <- helper.functions$getBetaResults

  ## custom stateful iterator object which yields a different resample everytime
  ## it is used in the loop
  individual.iterator <-
    sampleIterator(obj@individual.sample.seed,
                   obj@individual.level.data[[obj@context.id]],
                   obj@nb.resamples)

  output.obj@individual.sample.seed <- obj@individual.sample.seed

  ## loops through all resamples of the context data and aggregates them at the
  ## upper level. The results are stored in a list of matrices, each containing
  ## the resampled data for one variable
  start.time <- proc.time()[3]

  ## These two variables are only declared so that they won't appear to be global
  ## variables without visible binding during R CMD check
  individual.slice <- NULL
  column <- NULL
  analyses <-
    foreach(individual.slice=individual.iterator, column=1:obj@nb.resamples,
            .inorder=TRUE) %do% {
              sml.obj <-
                new("MLSpawExactObject",
                    individual.level.data =
                      obj@individual.level.data[individual.slice,],
                    precise.data = obj@precise.data,
                    context.id = obj@context.id,
                    formula = obj@formula,
                    nb.analyses = obj@nb.analyses)

              lme <- PerformMLSpawExact(sml.obj, ...)
              fillMatrices(lme@lme, lme@beta, column)
              elapsed.time <- proc.time()[3] - start.time
              if (obj@verbose){
                cat("\rcomputed step ", column, " of ", obj@nb.resamples,
                    ". ETA = ", as.integer((obj@nb.resamples/column-1)*elapsed.time))
              }
              NULL
            }


  cat("\rspacom done                                                        \n")
  tmp.list <- getTmpList()

  output.obj@ranefs <- getRanefs()
  for (name in names(tmp.list)) {
    ## compute mean (mandatory)
    mu <- rowMeans(tmp.list[[name]])
    frame <- data.frame(mean=mu)
    ## compute standard deviation (mandatory)
    frame$sd <- apply(tmp.list[[name]], 1, sd)
    ## compute percentiles
    perc.frame <- t(apply(tmp.list[[name]], 1, quantile, probs=obj@percentiles))
    slot(output.obj, name) <- cbind(frame, perc.frame)
  }

  return(output.obj)
}


ResampleMLSpawExact <-function(individual.level.data,
                               context.id,
                               formula,
                               precise.data,
                               confidence.intervals=c(.95),
                               nb.resamples=1000,
                               individual.sample.seed=NULL,
                               verbose=TRUE,
                               ...) {
  obj <-
    MakeResampleMLSpawExactObject(individual.level.data = individual.level.data,
                                  context.id = context.id,
                                  formula = formula,
                                  precise.data = precise.data,
                                  confidence.intervals = confidence.intervals,
                                  nb.resamples = nb.resamples,
                                  individual.sample.seed=individual.sample.seed,
                                  verbose=verbose)
  model <- PerformResampleMLSpawExact(obj,...)
}







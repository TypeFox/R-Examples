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


## provides a ResampleMLSpawAggregateObject needed for the analysis and
## performs all consistency checks
MakeResampleMLSpawAggregateObject <-function(individual.level.data,
                                             context.id,
                                             formula,
                                             aggregates,
                                             precise.data=NULL,
                                             confidence.intervals,
                                             individual.sample.seed=NULL,
                                             verbose) {
  ## check the aggregates input object
  if (!is(aggregates, "SpawAggregateOutput")) {
    stop("The input parameter 'aggregates' has to be of type ",
         "'SpawAggregateOutput'. You can create such an object using the ",
         "function 'SpawAggregate'")
  }

  ## get the number of samples from the SpawAggregate output
  nb.resamples <- length(aggregates)

  ## create an empty ResampleMLSpawAggregateObject to fill, using the
  ## parent class MLSpawExactObject

  obj <-
    MakeResampleMLSpawExactObject(individual.level.data=individual.level.data,
                                  context.id=context.id,
                                  formula=formula,
                                  precise.data=precise.data,
                                  confidence.intervals,
                                  nb.resamples=nb.resamples,
                                  individual.sample.seed=individual.sample.seed,
                                  verbose,
                                  obj=new("ResampleMLSpawAggregateObject"))

  obj@aggregates <- aggregates
  return(obj)
}

ResampleMLSpawAggregate <-function(individual.level.data,
                                   context.id,
                                   formula,
                                   aggregates,
                                   precise.data=NULL,
                                   confidence.intervals=c(.95),
                                   individual.sample.seed=NULL,
                                   verbose=TRUE,
                                   ...) {
  obj <-
    MakeResampleMLSpawAggregateObject(individual.level.data =
                                      individual.level.data,
                                      context.id=context.id,
                                      formula=formula,
                                      aggregates=aggregates,
                                      precise.data=precise.data,
                                      confidence.intervals=confidence.intervals,
                                      individual.sample.seed=
                                      individual.sample.seed,
                                      verbose=verbose)
  model <- PerformResampleMLSpawAggregate(obj,...)
}


PerformResampleMLSpawAggregate <- function(obj, ...){
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
    foreach(
      individual.slice=individual.iterator,
      column=1:obj@nb.resamples,
      .inorder=TRUE) %do% {
        ## need to merge precise data with the aggregates bootstrap sample
        precise.data <- merge(obj@precise.data,
                              obj@aggregates[column],
                              by=obj@context.id)
        sml.obj <-
          new("MLSpawExactObject",
              individual.level.data =
                obj@individual.level.data[individual.slice,],
              precise.data = precise.data,
              context.id = obj@context.id,
              formula = obj@formula,
              nb.analyses = obj@nb.analyses)

        lme <- PerformMLSpawExact(sml.obj, ...)
        fillMatrices(lme@lme, lme@beta, column)
        elapsed.time <- proc.time()[3]-start.time
        if (obj@verbose) {
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

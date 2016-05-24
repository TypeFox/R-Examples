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


makeMoranObject <-
  function(ml.spaw.obj,
           distance.matrix,
           bandwidths,
           kernel,
           confidence.intervals,
           verbose) {
  obj <- new("MoranObject")
  if (is(ml.spaw.obj, "ResampleMLSpawOutput")) {
    obj@ranefs <- ml.spaw.obj@ranefs
  } else if (is(ml.spaw.obj, "MLSpawExactOutput")) {
    obj@ranefs <- as.matrix(ranef(ml.spaw.obj)[[1]])
  } else if (is(ml.spaw.obj, "matrix")) {
    obj@ranefs <- ml.spaw.obj
  } else {
    stop("The parameter 'ml.spaw.obj' has to be either an output object from ",
         "any of spacom's MLSpaw methods * or a matrix. You specified ",
         "an object of class '", class(ml.spaw.obj), "'.\n",
         "  1) MLSpawExact\n  2) ResampleMLSpawExact\n  ",
         "3) ResampleMLSpawAggregate")
  }

  obj@nb.resamples <- ncol(obj@ranefs)

  obj@weights.object <- makeWeightsObject(distance.matrix, kernel, moran=TRUE)

  obj@bandwidths <- checkBandwidths(bandwidths)
  obj@nb.moron <- length(bandwidths)

  ## make sure the confidence intervals are correctly defined and contain the
  ## median
  if(obj@nb.resamples > 1) {
    obj@percentiles <-
      checkConfidenceIntervals(confidence.intervals, obj@nb.resamples)
  }

  ## check verbose flag
  obj@verbose <- check.flag(verbose, "verbose")

  return(obj)
}


MLSpawResidMoran <-
  function(ml.spaw.obj,
           distance.matrix,
           bandwidths,
           kernel=NULL,
           confidence.intervals=c(.95),
           verbose=TRUE) {
    obj <- makeMoranObject(
      ml.spaw.obj=ml.spaw.obj,
      distance.matrix=distance.matrix,
      bandwidths=bandwidths,
      kernel=kernel,
      confidence.intervals=confidence.intervals,
      verbose=verbose)

  tmp.matrix <- matrix(nrow=obj@nb.moron, ncol=obj@nb.resamples)
  row.names(tmp.matrix) <- lapply(bandwidths,
                                  FUN=function(x){paste('h.',x,sep='')})
  for (row in 1:obj@nb.moron) {
    message("Computing I for bandwidth ", obj@bandwidths[row])
    weight <- performWeights(obj@weights.object, obj@bandwidths[row])
    weight <- apply(weight, 2, "/", rowSums(weight))
    ordered.ids <- order(colnames(weight))
    reverse.order <- order(ordered.ids)
    weight <- mat2listw(weight)
    start.time <- proc.time()[3]
    reordered.ranefs <- as.matrix(obj@ranefs[reverse.order,])
    for (column in 1:obj@nb.resamples) {
      a <- reordered.ranefs[ , column]
      tmp.matrix[row, column] <-
        as.numeric(lm.morantest(lm(a~1, data=as.data.frame(reordered.ranefs)),
                                listw=weight)$estimate[1])
      elapsed.time <- proc.time()[3] - start.time
      if (obj@verbose) {
        cat("\rcomputed step ", column, " of ", obj@nb.resamples,
            ". ETA = ", as.integer((obj@nb.resamples/column-1)*elapsed.time))
      }
    }
    cat(
      "\rmoran done                                                         \n")
  }

  if (ncol(tmp.matrix) > 1) {
  ## compute mean (mandatory)
    mu <- rowMeans(tmp.matrix)
    frame <- data.frame(mean=mu)
     ## compute standard deviation (mandatory)
    frame$sd <- apply(tmp.matrix, 1, sd)

    ## compute percentiles
    perc.frame <- t(apply(tmp.matrix, 1, quantile, probs=obj@percentiles))
    return (cbind(frame, perc.frame))
  } else {
    frame <- data.frame(tmp.matrix)
    names(frame) <- "I"
    return(frame)
  }
}

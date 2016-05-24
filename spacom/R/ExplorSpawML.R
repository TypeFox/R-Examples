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

MakeExploreMLSpawObject <-
  function(individual.level.data,
           contextual.name,
           contextual.data,
           precise.data,
           context.id,
           formula,
           distance.matrix,
           multilevel.bandwidths,
           design.weight.names=NULL,
           aggregation.function="mean",
           kernel=NULL,
           additional.args=NULL,
           verbose=NULL)
{
  obj <- new("ExploreSpawML")

  ## make sure individual.level.data is a data.frame
  obj@individual.level.data <- checkType(individual.level.data, c("data.frame"),
                                         "individual.level.data")

  ## make sure contextual.data is either NULL or a data.frame
  obj@contextual.data <- checkType(obj=contextual.data,
                                   classes=c("NULL", "data.frame"),
                                   name="contextual.data")

  obj@precise.data <- checkType(obj=precise.data,
                                   classes=c("NULL", "data.frame"),
                                   name="precise.data")
  contextual.variables <- names(GetContext(obj))
  if (!is.null(precise.data)) {
    contextual.variables <- c(contextual.variables,
                              names(precise.data))
  }
  obj <- checkContextualNames(obj, contextual.name, contextual.variables)

  obj@design.weight.name <-
    checkExistence(design.weight.names,
                   contextual.variables)


  ## make sure context.id is a name in contextual.data
  obj@context.id <-
    checkContextId(
      context.id=context.id,
      individual.level.data.names=names(obj@individual.level.data),
      precise.data.names=contextual.variables,
      skip.individual.check=FALSE)

  ## extract number of upper level units
  obj@nb.area <- length(levels(as.factor(GetContext(obj)[[context.id]])))

  ## make sure the formula isn't a string anymore
  obj@formula <- checkFormula(formula)

  ## check the distance matrix for consistency
  if (is(distance.matrix, "data.frame")) {
    distance.matrix <- as.matrix(distance.matrix)
  }
  if (!is(distance.matrix, "matrix") && !is(distance.matrix, "Matrix")) {
    stop("The distance matrix has to be of class 'matrix' or 'Matrix'. You ",
         "specified an object of class '", class(distance.matrix), "'.")
  }

  if (!nrow(distance.matrix) == ncol(distance.matrix)) {
    stop("The distance matrix has to be square of size ", obj@nb.area,
         " (number of areas, you specified a matrix of size ",
         nrow(distance.matrix),"x", ncol(distance.matrix))
  }
  if (!nrow(distance.matrix)==obj@nb.area) {
    stop("The distance matrix has to be square of size ", obj@nb.area,
         " (number of areas), you specified a square matrix of size ",
         nrow(distance.matrix))
  }
  obj@distance.matrix <- Matrix(distance.matrix)

  ## check multilevel.bandwidths for consistency
  obj@multilevel.bandwidths <- checkBandwidths(multilevel.bandwidths)


  ## check aggregation.function for consistency
  aggregation.function <- checkType(aggregation.function, "character",
                                    "aggregation.function")

  obj@aggregation.function <- checkAggregationFunctions(aggregation.function,
                                                        1,
                                                        obj@contextual.names
                                                        )

  ## deal with kernel function
  obj@kernel <- checkKernel(kernel)

  ## make sure the additional.args are ok
  obj@ additional.args <- additional.args

  ## deal with verbose flag
  obj@verbose <- check.flag(verbose, "verbose")

  return(obj)
}

performExploreMLSpaw <- function(obj){
  ## prepare a weights object to compute weight.matrices on the fly
  weight.object <- new("weightsObject",
                       distance.matrix=obj@distance.matrix,
                       kernel=obj@kernel,
                       moran=FALSE)
  coded.name <- obj@contextual.names[[1]]
  name <- substr(coded.name, 1, nchar(coded.name)-5)

  output.list <- list()
  for (bandwidth in obj@multilevel.bandwidths) {
    ## compute the weight matrix
    weight.matrix <- performWeights(weight.object, bandwidth)

    ## if precise.data is not NULL, we're in exact mode
    if (!is.null(obj@precise.data)) {
      context <- SpawExact(precise.data=obj@precise.data,
                           context.id=obj@context.id,
                           contextual.names=name,
                           contextual.weight.matrices=weight.matrix)
    } else {
      context <- SpawAggregate(
        contextual.data=obj@contextual.data,
        context.id=obj@context.id,
        contextual.names=name,
        contextual.weight.matrices=weight.matrix,
        nb.resamples=0,
        aggregation.functions=obj@aggregation.function,
        design.weight.names=obj@design.weight.name,
        additional.args=obj@additional.args,
        verbose=obj@verbose)
    }

    ## screw the pooch
    formula <-
      as.formula(paste(as.character(obj@formula[2]), '~',
                       as.character(obj@formula[3]), '+',
                       name, ".1", sep=""))
    model <- MLSpawExact(individual.level.data=obj@individual.level.data,
                         context.id=obj@context.id,
                         formula=formula,
                         precise.data=context,
                         verbose=obj@verbose)
    output.list[[paste("bandwidth = ", bandwidth)]] <- model
  }

  return(output.list)
}

ExploreMLSpawExact <-
  function(individual.level.data,
           contextual.name,
           context.id,
           formula,
           distance.matrix,
           multilevel.bandwidths,
           precise.data,
           kernel=NULL,
           verbose=TRUE)
{
obj <-
  MakeExploreMLSpawObject(
    individual.level.data=individual.level.data,
    contextual.name=contextual.name,
    contextual.data=NULL,
    precise.data=precise.data,
    context.id=context.id,
    formula=formula,
    distance.matrix=distance.matrix,
    multilevel.bandwidths=multilevel.bandwidths,
    design.weight.names=NULL,
    aggregation.function="weighted.mean",
    kernel=kernel,
    additional.args=NULL,
    verbose=verbose)
  output.obj <- performExploreMLSpaw(obj)
}

ExploreMLSpawAggregate <-
  function(individual.level.data,
           contextual.name,
           contextual.data,
           context.id,
           formula,
           distance.matrix,
           multilevel.bandwidths,
           design.weight.names=NULL,
           aggregation.function="weighted.mean",
           kernel=NULL,
           additional.args=NULL,
           verbose=TRUE)
{
obj <-
  MakeExploreMLSpawObject(
    individual.level.data=individual.level.data,
    contextual.name=contextual.name,
    contextual.data=contextual.data,
    precise.data=NULL,
    context.id=context.id,
    formula=formula,
    distance.matrix=distance.matrix,
    multilevel.bandwidths=multilevel.bandwidths,
    design.weight.names=design.weight.names,
    aggregation.function=aggregation.function,
    kernel=kernel,
    additional.args=additional.args,
    verbose=verbose)

  output.obj <- performExploreMLSpaw(obj)
}

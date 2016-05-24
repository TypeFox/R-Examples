################################################################################
################################################################################
## author Till Junge <till.junge@altermail.ch>                                ##
## author Mathieu Cossutta <mcossutta@gmail.com>                              ##
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


## provides a SpawExactObject with basic consistency checks
MakeSpawExactObject <- function(contextual.data,
                                context.id,
                                contextual.names,
                                contextual.weight.matrices,
                                population.weight.names,
                                verbose,
                                obj=NULL){
  ## create an empty SpawExactObject
  pedantic.check <- FALSE
  if (is.null(obj)){
    obj <- new("SpawExactObject")
    pedantic.check <- TRUE
  }

  ## make sure contextual.data is a data.frame
  obj@contextual.data <- checkType(contextual.data, c("data.frame"),
                                   "contextual.data")

  ## extract number of upper level units
  context.ids <- unique(obj@contextual.data[[context.id]])
  obj@nb.area <- length(context.ids)
  if (pedantic.check && nrow(obj@contextual.data) != obj@nb.area) {
    stop("The contextual data has to have as many rows as there are areas, ",
         obj@nb.area,". ",
         "You specified a frame with ", nrow(obj@contextual.data), " rows.")}

  ## make sure context.id is a name in contextual.data
  obj@context.id <-
    checkContextId(context.id=context.id,
                   individual.level.data.names=NULL,
                   precise.data.names=names(obj@contextual.data),
                   skip.individual.check=TRUE)

  ## make sure the contextual names are in a list
  obj <- checkContextualNames(obj, contextual.names, names(contextual.data))

  ## make sure the weight matrix list is suitable
  weight.list <-
    checkAllWeightMatrices(contextual.weight.matrices, obj@contextual.names,
                           obj@nb.area, obj@nb.analyses)
  obj@contextual.weight.matrices <- weight.list[[1]]
  matrix.size <- weight.list[[2]]

  ## if there are contextual.weight.matrices, we need to make sure that
  ## contex.id is either numbered from 1:nb.area, or that it corresponds to
  ## the colnames of the matrix
  if (is.null(contextual.weight.matrices)){
    column.names <- context.ids
  } else {
    column.names <- checkSameColNames(obj@contextual.weight.matrices)
  }
  if (is.null(column.names)) {
    stop("your weight matrix needs column names ")
    column.names <- as.character(1:matrix.size)
  }
  if (!all(context.ids %in% column.names)) {
    stop("The column names in the weight matrices do not correspond to the ",
         "context identifiers. Your context ids must all appear in the column ",
         "names of the weight matrix\n")
  }
  if (!all(column.names %in% context.ids)) {
    slice <- column.names %in% context.ids
    warning(
      "The following context id(s) appear in the weight matrix but not in the ",
      "data: '", paste(column.names[!slice], collapse="', '"), "'. Please ",
      "make sure this is what you want")
  }
  num.id <- 1:length(column.names)
  names(num.id) <- column.names
  obj@numeric.context.id <- 'NUMERIC_AREA_ID'
  obj@contextual.data[[obj@numeric.context.id]] <-
    num.id[as.character(obj@contextual.data[[obj@context.id]])]
  class(obj@unique.context.ids) <- class(column.names[1])
  counter = 1
  for (name in column.names) {
    if (name %in% context.ids) {
      obj@unique.context.ids[counter] <- name
      counter <- counter+1
    }
  }

  ## make sure the population.weight.names are in
  obj@population.weight.names <-
    checkAllDesignWeights(population.weight.names,
                          contextual.names=obj@contextual.names,
                          context.variables=names(obj@contextual.data),
                          nb.analyses=obj@nb.analyses,
                          level="population")
  
  ## check verbose flag
  obj@verbose <- check.flag(verbose, "verbose")
  return(obj)
}

performSpawExact <- function(obj){
  ## prepare a zero-filled named list to keep count
  count.list <- list()
  decoded.names <- lapply(obj@contextual.names, function(x)
                          {return(substr(x, 1, nchar(x)-5))})
  for(name in decoded.names) {
    count.list[[name]] <- 0}



  ## the order of the columns in the matrix do not necessarily correspond to
  ## order of the areas in the data. Hence we need to look out for it
  context.order <- order(obj@contextual.data[[obj@context.id]])

  ## prepare the merge dataframe into which the new weighted contextual data
  ## are loaded
  merge.data <- data.frame(obj@contextual.data[[obj@context.id]][context.order])
  names(merge.data) <- obj@context.id
  ## loop through the contextual variables to be analysed

  for (i in 1:obj@nb.analyses){
    coded.name = obj@contextual.names[[i]]
    len <- nchar(coded.name)
    name <- substr(coded.name, 1, len-5)

    context <- as.matrix(obj@contextual.data[[name]][context.order], ncol=1)
    population.weight.name <- obj@population.weight.names[[coded.name]]
    if (!is.null(population.weight.name)) {
      weight <- as.matrix(obj@contextual.data[[population.weight.name]][context.order],
                          ncol)
      context <- context * weight/mean(weight)
    }
    spatial.weight = obj@contextual.weight.matrices[[coded.name]]
    if (!is.null(spatial.weight)) {
      spatial.order <- order(colnames(spatial.weight))
      spatial.weight <- spatial.weight[spatial.order, spatial.order]
      tryCatch( {
          context = as.matrix(spatial.weight %*% context/rowSums(spatial.weight))
      }, error=function(e) {
          print("weight class = ")
          print(class(spatial.weight))
          print("rowSums(weight):");
          print(rowSums(spatial.weight));
          stop("fucked up when trying to multiply the weight matrix of dim  ", dim(spatial.weight), "with the context ", context)})
    }
    ## store the aggregated context in merge.data for later merge with
    ## the rest, compute the appropriate renaming of the contextual variables
    count.list[[name]] <- count.list[[name]] + 1
    new.name <- paste(name, ".", count.list[[name]], sep='')

    merge.data[[new.name]] <- context

  }
  return(merge.data)
}

SpawExact <- function(precise.data,
                      context.id,
                      contextual.names,
                      contextual.weight.matrices,
                      population.weight.names=NULL){
  obj <-MakeSpawExactObject(contextual.data=precise.data,
                            context.id,
                            contextual.names,
                            contextual.weight.matrices,
                            population.weight.names,
                            verbose=FALSE)
  output <- performSpawExact(obj)
  return(output)
}



## provides a SpawExactObject with basic consistency checks
MakeSpawAggregateObject <- function(contextual.data,
                                    context.id,
                                    contextual.names,
                                    contextual.weight.matrices,
                                    nb.resamples,
                                    aggregation.functions,
                                    confidence.intervals,
                                    design.weight.names,
                                    sample.seed,
                                    additional.args,
                                    verbose){
  ## use the parent class' maker function
  obj <-
  MakeSpawExactObject(contextual.data = contextual.data,
                      context.id = context.id,
                      contextual.names = contextual.names,
                      contextual.weight.matrices = contextual.weight.matrices,
                      population.weight.names=NULL,
                      verbose,
                      obj=new("SpawAggregateObject"))
  ##
  obj@nb.resamples <- as.integer(nb.resamples)
  if (obj@nb.resamples == 1) {
    stop("One resample makes no sense. If you just want to aggregate the ",
         "data, that's nb.resamples = 0.\n")
  }

  ## make sure the contextual names are in a list
  obj <- checkContextualNames(obj, contextual.names, names(contextual.data))


  ## make sure the aggregation function list is suitable
  obj@aggregation.functions <- checkAggregationFunctions(aggregation.functions,
                                                         obj@nb.analyses,
                                                         obj@aggregation.names)

  ## make sure the confidence intervals are correctly defined and contain the
  ## median
  obj@percentiles <-
    checkConfidenceIntervals(confidence.intervals, obj@nb.resamples)

  ## make sure the design.weight.names are suitable
  obj@design.weight.names <-
    checkAllDesignWeights(design.weight.names,
                          obj@contextual.names,
                          names(obj@contextual.data),
                          obj@nb.analyses)

  ## make sure the sample seed is ok
  if (nb.resamples != 0) {
    obj@sample.seed <- checkSeed(sample.seed)
  } else {
    obj@sample.seed <- as.integer(0)
  }


  ## make sure the additional.args are ok
  obj@ additional.args <-
    getHash(obj@contextual.names,
            checkNestedArguments(args=additional.args,
                                 context.variables=names(obj@contextual.data),
                                 description="additional argument",
                                 max.level=2))
  return(obj)
}

print.descr <- function(object) {
  cat("aggregated data with descriptives:\n")
  print(object@frames)
  cat("\nslots for direct access to data:\n  usage: <object.name>@<slot.name>\n")
  is.print <- TRUE
  for (name in slotNames(object)) {
    if (is.print) {
      print(name)
    } else {
      show(name)
    }
  }
}
setMethod("show", signature="SpawAggregateOutput",definition=print.descr)

mergeMethod <- function(x, y, ...) {
  nb.area.x <- dim(x@aggregated.samples[[1]])[1]
  nb.area.y <- dim(y@aggregated.samples[[1]])[1]
  nb.samples.x <- dim(x@aggregated.samples[[1]])[2]
  nb.samples.y <- dim(y@aggregated.samples[[1]])[2]

  if (nb.area.x != nb.area.y) {
    stop(
      "The two objects you're trying to merge do not have the same number of ",
      "contextual units")
  }
  if (nb.samples.x != nb.samples.y) {
    stop(
      "The two objects you're trying to merge do not have the same number of ",
      "bootstrap samples")
  }
  obj <- new("SpawAggregateOutput",
             seed=NULL,
             aggregated.samples=c(x@aggregated.samples, y@aggregated.samples),
             frames=c(x@frames, y@frames))
}
setGeneric('merge')
setMethod("merge", signature("SpawAggregateOutput", "SpawAggregateOutput"), definition=mergeMethod)



performSpawAggregate <- function(obj, exploratory=FALSE) {
  ## creation of an output data to bundle all the return data
  output.obj <- new("SpawAggregateOutput")

  ## temporary array object to gather the results
  aggregation.results <-
    array(dim=c(obj@nb.area, max(obj@nb.resamples, 1), obj@nb.analyses))

  ## inserts the analyses at the right position in the results
  fillAggregationResults <- function(frame, index){
    for (column in 2:ncol(frame)) {
      aggregation.results[,index, column-1] <<- frame[[column]]
    }
  }

  ## custom stateful iterator object which yields a different resample every
  ## time it is used in the loop
  ## if obj@nb.resamples == 0, then the data is unchanged

  iterator <- sampleIterator(obj@sample.seed,
                             obj@contextual.data[[obj@numeric.context.id]],
                             obj@nb.resamples)
  output.obj@seed <- obj@sample.seed

  ## loops through all resamples of the context data and aggregates them at the
  ## upper level. The results are stored in a list of matrices, each containing
  ## the resampled data for one variable
  gc(FALSE)
  start.time <- proc.time()[3]

  ## These two variables are only declared so that they won't appear to be global
  ## variables without visible binding during R CMD check
  slice <- NULL
  column <- NULL
  
  analyses <-
    foreach(slice=iterator,
            column=1:max(obj@nb.resamples,1),
            .inorder=TRUE) %do% {
    frame <-
      performAggregation(contextual.names = obj@contextual.names,
                         context.id = obj@numeric.context.id,
                         nb.area = obj@nb.area,
                         nb.analyses=obj@nb.analyses,
                         design.weight.names = obj@design.weight.names,
                         contextual.data = obj@contextual.data[slice,],
                         aggregation.functions = obj@aggregation.functions,
                         contextual.weight.matrices =
                           obj@contextual.weight.matrices,
                         additional.args = obj@additional.args)
    fillAggregationResults(frame, column)
    if (obj@verbose) {
        elapsed.time <- proc.time()[3] - start.time
        cat("\rcomputed step ", column, " of ", obj@nb.resamples,
            ". ETA = ", as.integer((obj@nb.resamples/column-1)*elapsed.time))
    }
    NULL
  }
  gc(FALSE)
  cat("\rdescription done                                                   \n")
  if (exploratory) {
    return(matrix(aggregation.results, nrow=obj@nb.area, ncol=obj@nb.resamples))
  }
  variable.names <- names(frame)[1:obj@nb.analyses+1]

  for (i in 1:length(variable.names)) {
    name <- variable.names[i]
    ## add the computed matrices to the output object
    output.obj@aggregated.samples[[name]] <-
      matrix(aggregation.results[,,i],
             nrow=obj@nb.area)
    ## here, loop through the upper level aggregation functions (mean, sd,
    ## median, percentiles)
    frame = data.frame(obj@unique.context.ids)
    names(frame) <- obj@context.id

    if (obj@nb.resamples > 0) {
      ## add mean(mandatory)
      frame[['mean']] <-
        rowMeans(aggregation.results[,,i])
      ## add standard deviation(mandatory)
      frame[['sd']] <-
        apply(aggregation.results[,,i],1,sd)
      ## add percentiles
      perc.frame <- t(apply(aggregation.results[,,i], 1, quantile, probs=obj@percentiles))
      output.obj@frames[[name]] <- cbind(frame, perc.frame)
    } else {
      frame[[name]] <- aggregation.results[, , i]
      output.obj@frames[[name]] <- frame
    }
  }
  if (obj@nb.resamples > 0) {
    return(output.obj)
  } else {
    return.frame <- data.frame(obj@unique.context.ids)
    names(return.frame) <- obj@context.id
    for (frame in output.obj@frames) {
      for (name in names(frame)[-1]) {
        return.frame[[name]] <- frame[[name]]
      }
    }
    return(return.frame)
  }
}

SpawAggregate <- function(contextual.data,
                          context.id,
                          contextual.names,
                          contextual.weight.matrices=NULL,
                          nb.resamples=1000,
                          aggregation.functions='weighted.mean',
                          confidence.intervals=.95,
                          design.weight.names=NULL,
                          sample.seed=NULL,
                          additional.args=NULL,
                          verbose=TRUE){
  input.obj <- MakeSpawAggregateObject(contextual.data,
                                       context.id,
                                       contextual.names,
                                       contextual.weight.matrices,
                                       nb.resamples,
                                       aggregation.functions,
                                       confidence.intervals,
                                       design.weight.names,
                                       sample.seed,
                                       additional.args,
                                       verbose)
  return(performSpawAggregate(input.obj))
}

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

lengthMethod <- function(x) {
  return(dim(x@aggregated.samples[[1]])[2])
}
setMethod("length", signature=c("SpawAggregateOutput"),
          definition=lengthMethod)

getSampleMethod <- function(object, sample.num) {
  if (sample.num < 1 || sample.num > length(object)) {
    stop("Out of bounds. You tried to access the ", sample.num, "th sample of ",
         "a SpawAggregateOutput which contains ", length(object), " samples.")
  }

  frame <- data.frame(object@frames[[1]][[1]])
  names(frame) <- names(object@frames[[1]])[1]
  for (frame.id in 1:length(object@frames)) {
    name <- names(object@frames)[frame.id]
    frame[[name]] <- object@aggregated.samples[[frame.id]][,sample.num]
  }
  return(frame)
}

setMethod("names",
          signature=c("SpawAggregateOutput"),
          definition=function(x){return (names(x@frames))})

names.assigment.SpawAggregateOutput <- function(x, value){
  names(x@frames) <- value
  names(x@aggregated.samples) <- value
  return(x)
}
setMethod("names<-",
          signature=c("SpawAggregateOutput"),
          definition=names.assigment.SpawAggregateOutput
        )

getSample <- function(object, sample.num){}
setGeneric("getSample")
setMethod("getSample", signature=c("SpawAggregateOutput"),
          definition=getSampleMethod)

setMethod("[",
          signature=c("SpawAggregateOutput"),
          definition=function(x, i, j, ..., drop) { return(getSample(x,i))})

## define a method which yields the contextual data, if it exists, or the
## individual level data in the opposite case
GetContextMethod <- function(object) {
  if (is.null(object@ contextual.data)) {
    object@ individual.level.data
  } else {
    object@ contextual.data
  }
}

GetContext <- function(object){}
setGeneric("GetContext")
setMethod("GetContext", signature=c("SpawExactObject"),
          definition=function(object) {return(object@contextual.data)})
setMethod("GetContext", signature=c("SpawAggregateObject"),
          definition=function(object) {return(object@contextual.data)})
setMethod("GetContext", signature=c("ExploreSpawML"),
          definition=GetContextMethod)

getHash <- function(keys, data) {
  indices <- list()

  if (!is(data, "list")){
    for (key in keys) {
      indices[[key]] <- 1
    }
    ret.data <- list()
    if (is.null(data)){
      data <- "RETARDEDNULL"
    }
    ret.data[[1]] <- data
    if (is.null(data)) {
      ret.data[[2]] <- 1
    }
    return(new("MultiKeyHash", indices=indices, data=ret.data))
  }
  for (i in seq(length(keys))) {
    indices[[keys[[i]]]] <- i
    if(is.null(data[[i]])) {
      data[[i]] <- "RETARDEDNULL"
    }
  }
  return (new("MultiKeyHash", indices=indices, data=data))
}

slc <- function(x, i, j, ..., drop) {
  item <- x@data[[x@indices[[i]]]]
  if(is(item, "character") && item == "RETARDEDNULL") {
    return (NULL)
  }
  return(item)
}
setMethod("[",
          signature=c("MultiKeyHash"),
          definition=slc)
setMethod("[[",
          signature=c("MultiKeyHash"),
          definition=slc)

len <- function(x) {
  return(length(x@indices))
}

setMethod("length",
          signature=c("MultiKeyHash"),
          definition=len)

assgnmt <- function(x, i, j, value) {
  if (is.null(value)) {
    value <- "RETARDEDNULL"
  }
  index <- length(x@data)+1
  x@data[[index]] <- value
  x@indices[[i]] <- index
  return(x)
}
setMethod(`[<-`,
          signature=c("MultiKeyHash"),
          definition=assgnmt)

setMethod(`[[<-`,
          signature=c("MultiKeyHash"),
          definition=assgnmt)



## replace the replicate function but for a list output
makeList <- function(nb.elem, object){
  ret.list <- list()
  for (i in seq (nb.elem)){
    ret.list[[i]] <- object
  }
  return(ret.list)
}

## check that a flag is TRUE or FALSE
check.flag <- function(flag, name = "unnamed") {
  if (!is(flag, "logical")){
    stop ("The flag '", name, "' needs to be of class 'logical', not '",
          class(flag), "'")
  }
  return(flag)
}

## checks whether a variable name is NULL or denotes a column of obj
checkExistence <- function(name,
                           context.variables,
                           description="design.weight") {
  if (is.null(name)) {
    return(name)
  }
  else if (is(name,"character")) {
    if (name %in% context.variables) {
      return(name)
    } else {
      stop("The ", paste(description), " name '", name,
           "' is not in the column names {", paste(context.variables, collapse=", "), "}\n")
    }
  } else {
    stop("The ", paste(description), " names have to be of type 'character' or 'NULL'",
         " or a list thereof. You gave an object of type ", class(name),
         "\n")
  }
}


## Checks whether an object is a nested list of existing variable names
## according to checkExistence
checkNestedArguments <- function(args,
                                 context.variables,
                                 description,
                                 max.level) {

  ## recursive check of nested list above the root level
  checkNestedList <- function(arglist, level) {
    if (level == 0) {
      stop ("reached bottom level of nested list!\n")
    }
    for (arg in arglist) {
      if (is(arg, "list")) {
        checkNestedList(arg, level-1)
      } else {
        checkExistence(arg, context.variables, description)
      }
    }
  }
  if (!is(args, "list")) {
    args <- checkExistence(args,
                           context.variables,
                           description)
  } else {
    checkNestedList(args, max.level)
  }
  return(args)
}

## checks whether all variable names chosen for individual weights are suitable
checkAllDesignWeights <- function(design.weight.names, contextual.names,
                                  context.variables,
                                  nb.analyses,
                                  level="individual"){
  if (!is(design.weight.names, "list")) {
    design.weight.names <- checkExistence(design.weight.names,
                                              context.variables,
                                              level)
    return(getHash(contextual.names, design.weight.names))
  } else {
    len.mat <- length(design.weight.names)
    if (len.mat != nb.analyses) {
      stop("You specified ", len.mat, level, " weight variables and ",
           nb.analyses, " contextual names. ",
           "You have to specify either 1 ", level, " weight name (which may be ",
           "NULL) or as many as you specified contextual names.")}
    tmp.list <- list()
    for (i in seq(len.mat)) {
      tmp.list[[i]] <- checkExistence(design.weight.names[[i]],
                                      context.variables,
                                      level)}
    return(getHash(contextual.names, tmp.list))
  }
}






## checks dimensions of square matrix
checkWeightMatrix <- function(matrix, nb.area, name){
  nb.rows <- nrow(matrix)
  if (is.null(matrix)){}
  else if (is(matrix, "matrix") || is(matrix, "Matrix")){
    if (nb.rows!=ncol(matrix) || nb.rows<nb.area){
      stop(cat("The weight matrix given for", name, "is", nb.rows, "x",
               ncol(matrix), "but should be square and at least", nb.area,
               "x", nb.area, '\n'))}
  } else {
    stop(cat("The weight matrices have to be of type 'matrix', 'Matrix' or ",
             "NULL. You gave an object of type", class(matrix), ".\n"))
  }
  return(nb.rows)
}

## checks whether all all weight matrices are suitable
checkAllWeightMatrices <- function(contextual.weight.matrices,
                                   contextual.names,
                                   nb.area, nb.analyses){
  nb.colnames <- 0
  if (!is(contextual.weight.matrices,"list")) {
    nb.colnames <-
      checkWeightMatrix(contextual.weight.matrices,
                        nb.area,
                        contextual.names[[1]])
    return(list(getHash(contextual.names, contextual.weight.matrices),
                nb.colnames))
  } else {
    len.mat <- length(contextual.weight.matrices)
    tmp.list <- list()
    if (len.mat != nb.analyses){
      stop("You specified ", len.mat, " weight matrice(s) and ",
           nb.analyses, " contextual name(s). ",
           "You have to specify either 1 weight (which may be NULL) or as ",
           "many as you specified contextual names")}
    for (i in seq(len.mat)){
      nb.colnames <- max(checkWeightMatrix(contextual.weight.matrices[[i]],
                                           nb.area,
                                           contextual.names[[i]]),
                         nb.colnames)
      tmp.list[[i]] <- contextual.weight.matrices[[i]]
    }
    length(tmp.list) <- len.mat
    return(list(getHash(contextual.names, tmp.list), nb.colnames))
  }
}

## Make sure all weight matrix have the same colnames (if any)
checkSameColNames <- function(matrices) {
  names <- NULL
  update.names <- function(matrix) {
    if (is(matrix, "character") && matrix == "RETARDEDNULL") {
    } else {
      if (is.null(names)) {
        names <<- colnames(matrix)
      } else {
        if (!all.equal(names, colnames(matrix))) {
          stop("All weight matrices need to have the same columns names")
        }
      }
    }
  }
  for (matrix in matrices@data) {
    update.names(matrix)
  }
  return(names)
}

## makes sure that the contextual names are in a list
checkContextualNames <- function(obj, contextual.names, contextual.variables){

  if (!is(contextual.names, "list")) {
    contextual.names <- as.list(contextual.names)}
  for (name in as.list(contextual.names)) {
    if (!name %in% contextual.variables) {
      stop("The variable '", name, "' is not in the context data")}}

  obj@nb.aggregations <- as.integer(0)
  obj@nb.precise.weightings <- as.integer(0)

  ## make an array of individualised contextual names to avoid name clashes
  counter <- 0
  coder <- function(x){
    counter <<- counter + 1
    return(paste(x,sprintf("__%03d",counter), sep=""))
  }
  obj@contextual.names <- lapply(contextual.names, coder)
  obj@nb.analyses <- length(obj@contextual.names)

  ## split the contextual names into names for aggregation and names for
  ## precise contexts
  if (length(obj@contextual.names)>0) {
    for (i in 1:length(obj@contextual.names)) {
      name <- contextual.names[[i]]
      if (name %in% names(GetContext(obj))) {
        obj@nb.aggregations <- obj@nb.aggregations + as.integer(1)
        obj@aggregation.names[[obj@nb.aggregations]] <- obj@contextual.names[[i]]
      } else if (!is.null(obj@precise.data)) {
        if (name %in% names(obj@precise.data)) {
          obj@nb.precise.weightings <- obj@nb.precise.weightings + as.integer(1)
          obj@precise.names[[obj@nb.precise.weightings]] <- obj@contextual.names[[i]]
        }
      }
    }
  }
  return(obj)
}

## Simple checkNames

checkNames <- function(contextual.names,names)
{

  for (name in contextual.names) {
    if (!name %in% names) {
      stop("The variable '", name, "' is not in the context data")}}
  return(contextual.names)
}


## makes sure context ids exist in the right dataframes
checkContextId <- function(context.id,
                           individual.level.data.names,
                           precise.data.names=NULL,
                           skip.individual.check=FALSE){
  if (!skip.individual.check && !context.id %in% individual.level.data.names){
    stop("The value given for context.id '", context.id, "' is not in the ",
         "names of individual.level.data")
  }
  if (!is.null(precise.data.names)){
    if(!context.id %in% precise.data.names){
      stop("The value given for context.id '", context.id, "' is not in the ",
           "names of precise.data")
    }
  }
  return(context.id)
}

checkContexts <- function(contextual.contexts,
                          individual.contexts) {
  a.not.in.b <- function(a, b) {
    slice <- !(a %in% b)
    return(a[slice])
  }
  dupl.contextual <- duplicated(contextual.contexts)
  if (any(dupl.contextual)) {
    stop("Context ids may appear only once in aggregated data.\n",
         "The following id(s) appeared more than once: '",
         paste(contextual.contexts[dupl.contextual], collapse="', '"), "'.\n")
  }
  if (is.null(individual.contexts)) {
    return()
  }
  unique.individual <- unique(individual.contexts)
  ind.not.in.cont <- a.not.in.b(unique.individual, contextual.contexts)
  cont.not.in.ind <- a.not.in.b(contextual.contexts, unique.individual)
  if (length(ind.not.in.cont)!=0) {
    warning(
      "The following context id(s) appear in the individual level data, but ",
      "not in the aggregated data: '",
      paste(ind.not.in.cont, collapse="', '"), "'.\n",
      "This is almost certainly not what you wanted.\n")
  }
  if (length(cont.not.in.ind)!=0) {
    warning(
      "The following context id(s) appear in the aggregated data, but ",
      "not in the individual level data: '",
      paste(ind.not.in.cont, collapse="', '"), "'.\n",
      "Please make sure this is what you wanted.\n")
  }
}

## makes sure a object is either one of defined types





checkType <- function(obj, classes, name){
  cls <- class(obj)
  if (!cls %in% classes){
    stop(name, " has to be of either of the types {",
         paste(classes, collapse=", "), "}. You gave an object of class ", cls)
  }
  return(obj)
}

## makes sure the aggregation funciton list is suitable
checkAggregationFunctions <- function(aggregation.functions,
                                      nb.analyses,
                                      aggregation.names){
  aggregation.functions <- checkType(aggregation.functions,
                                     c("character", "list"),
                                     "aggregation.functions")
  if (!class(aggregation.functions) == "list") {
    aggregation.functions <- as.list(replicate(nb.analyses, aggregation.functions))
  } else {
    len.mat <- length(aggregation.functions)
    if (len.mat != 1 && len.mat != nb.analyses) {
      stop("You specified ", len.mat, " aggregation functions and ",
           nb.analyses, " contextual names. ",
           "You have to specify either 1 function (which by default is 'mean') or as many ",
           "as you specified contextual names")}
    for (i in seq(len.mat)){
      checkType(aggregation.functions[[i]], c("character"),
                paste("aggregation.functions[[", i, "]]", sep=""))
    }

    if (len.mat == 1){
      aggregation.functions <- replicate(nb.analyses, aggregation.functions)
    }
  }
  names(aggregation.functions) <- aggregation.names
  return(aggregation.functions)
}

## makes sure that the object given as formula can be coerced into a formula
checkFormula <- function(formula){
  tryCatch(formula <- as.formula(formula),
           error=function(er) {
             stop(cat('what you gave as formula "', formula,
                      '" could not be coerced into a formula. Please give ',
                      'a character string or a formula\n'))})
  return(formula)
}

## makes sure that the number of resamples is positive and integer
checkNbResamples <- function(nb.resamples) {
  if (!is(nb.resamples, "numeric")) {
    stop("nb.resamples has to be a numeric (integer)  value. You specified ",
         "an object of class ", class(nb.resamples))
  }
  nb.resamples <- as.integer(nb.resamples)
  if (nb.resamples < 1) {
    stop("nb.resamples has to be larger than 1. You specified ",
         nb.resamples)
  }
  return(nb.resamples)
}

## makes sure that the confidence intervals are correctly defined and contain
## the median
checkConfidenceIntervals <- function(confidence.intervals, nb.resamples){
  tryCatch(percentiles <- as.numeric(confidence.intervals),
           error=function(er){
             stop("what you gave as percentiles '", confidence.intervals,
                  "' could not ",
                  "be coerced into numeric. Please give a numeric vector of ",
                  "percentiles")})
  len.confidence.intervals <- length(confidence.intervals)
  len.percentiles = 2*len.confidence.intervals + 1
  percentiles <- rep(0, len.percentiles)
  names <- rep('', len.percentiles)
  ## insure that the median is part of the percentiles
  for (i in seq(len.confidence.intervals)) {
    if (confidence.intervals[i] >= 1 || confidence.intervals[i] <= 0) {
      stop("Confidence interval ", i, " = ", confidence.intervals[i], " ",
           "is outside the interval I = {x|0<x<1}. Make sure all confidence ",
           "intervals are within I")}
    exclude <- 1 - confidence.intervals[i]
    if ((exclude*nb.resamples < 50) && (nb.resamples > 0)){
      warning("Warning, the confidence interval ", i, " = ",
              confidence.intervals[i], " excludes only ",
              ceiling(exclude*nb.resamples), " individuals. It may be ",
              "unreliable.")
    }
    percentiles[2*i-1] <- .5*exclude
    names[2*i-1] <- paste("CI_", confidence.intervals[i], "_lower", sep='')
    percentiles[2*i] <- 1-percentiles[2*i-1]
    names[2*i] <- paste("CI_", confidence.intervals[i], "_upper", sep='')
  }
  percentiles[len.percentiles] <- .5
  names[len.percentiles] <- "median"
  names(percentiles) <- names
  return(percentiles)
}

## makes sure the kernel function is a function (yeah, I heard it too)
checkKernel <- function(kernel) {
  if (is.null(kernel)) {
    kernel <- function(distance.matrix, h) {
      return(Matrix(.5^((distance.matrix/h)^2), sparse=TRUE))
    }
  } else if (!is(kernel, "function")) {
    stop("The kernel function has to be of class 'function'. You specified a ",
         "object of class '", class(kernel), "'")
  }

  return(kernel)
}


##makes sure that specified bandwidths are suitable
checkBandwidths <- function(bandwidths){
  if (!is(bandwidths, "numeric")) {
    stop("The bandwidths have to be a vector of numeric. You specified an ",
         "object of class '", class(bandwidths), "'.")
  }
  if (any(bandwidths<0)) {
    stop("only positive bandwidths are allowed. You specified ",
         paste(bandwidths, collapse=", "))
  }
  return(bandwidths)
}


## makes sure a sample seed is a valid random seed or a list of samples
checkSeed <- function(sample.seed, nb.resamples){

  runif(1) # single call to runif to assure the existence of .Random.seed
  current.seed <- .Random.seed
  if (is(sample.seed, "numeric")){
    if (!length(sample.seed) ==1){
      tryCatch(.Random.seed <- sample.seed, runif(1),
               error=function(er) {
                 stop("You specified a numeric of length ", length(sample.seed), " for ",
                      "the sample seed :", sample.seed, "\n Please specify one of the ",
                      "following instead:\na) a single numeric\nb) a .Random.seed",
                      "\nc) NULL\n",
                      "If you do not understand this error message, you probably want",
                      "to specify either NULL (if you do not to reproduce the same",
                      "samples) or a single digit number (for reproducibility)")})
      return(sample.seed)
    } else {
      set.seed(sample.seed)
      return(.Random.seed)
    }
  } else if (is.null(sample.seed)) {
    return(.Random.seed)
  } else {
    stop("You specified a ", class(sample.seed), " of length ",
         length(sample.seed), " for ",
         "the sample seed :", sample.seed, "\n Please specify one of the ",
         "following instead:\na) a single numeric\nb) a .Random.seed",
         "\nc) NULL\n",
         "If you do not understand this error message, you probably want",
         "to specify either NULL (if you do not to reproduce the same",
         "samples) or a single digit number (for reproducibility)")
  }
}

computeWeights <- function(design.weight.name,
                           design.weights,
                           area.ids,
                           weight.matrix,
                           area) {
  if (!is.null(weight.matrix)) {
    weight = weight.matrix[area, area.ids]
  } else {
    weight <- as.numeric(area.ids==area)
  }
  if (!is.null(design.weight.name)) {
    weight = weight * design.weights
  }
  return(weight)
}

performAggregation <- function(contextual.names,
                               context.id,
                               nb.area,
                               nb.analyses,
                               design.weight.names,
                               contextual.data,
                               aggregation.functions,
                               contextual.weight.matrices,
                               additional.args,
                               formula.str=NULL){
  ## prepare a zero-filled named list to keep count
  count.list <- list()
  decoded.names <- lapply(contextual.names, function(x)
                          {return(substr(x, 1, nchar(x)-5))})
  for(name in decoded.names) {
    count.list[[name]] <- 0}

  ## prepare the merge dataframe into which the new weighted contextual data
  ## are loaded
  merge.data <- data.frame(1:nb.area)
  names(merge.data) <- context.id
  ## loop through the contextual variables to be analysed

  if (nb.analyses > 0) {
    for (i in 1:nb.analyses) {
      coded.name <- contextual.names[[i]]
      len <- nchar(coded.name)
      name <- substr(coded.name, 1, len-5)

      aggregated.context <- matrix(nrow=nb.area, ncol=1)
      for (area in 1:nb.area) {
        design.weight.name <- design.weight.names[[coded.name]]
        additional.arg.name <- additional.args[[coded.name]]

        arguments=list()
        arguments[[name]] <- contextual.data[[name]]
        ## compute combined spatial and design weights
        weights <- computeWeights(design.weight.name,
                                  contextual.data[[design.weight.name]],
                                  contextual.data[[context.id]],
                                  contextual.weight.matrices[[coded.name]],
                                  area)
        ## if either spatial or design weights are specified, we need to
        ## add them to the arguments
        if (!is.null(design.weight.name) ||
            !is.null(contextual.weight.matrices[[coded.name]])) {
          arguments[["weights"]] <- weights
          ## additional arguments are acceptable only if weights are used
          for (argname in additional.arg.name) {
            arguments[[argname]] <- contextual.data[[argname]]
          }
        } else {
          arguments[[name]] <- arguments[[name]][weights==1]
          ## now we need to be careful with the weights, which may legitimately
          ## be NULL but still be followed by additional arguments
          if (length(additional.arg.name) != 0) {
            arguments[["weights"]] <- weights[weights==1]
            for (argname in additional.arg.name) {
              arguments[[argname]] <- contextual.data[[argname]][weights==1]
            }
          }
        }
        ## store names in case we need them for an error message
        arg.names <- names(arguments)
        ## and then drop them so the function call happens by position
        names(arguments) <- NULL
        error.handler <- function(er) {
          stop(
            "Spacom encountered the following error when aggregating:\n'",
            er, "'.\nThis is probably not a spacom error, but a problem with ",
            "the aggregation function and the arguments you provided.\n",
            "You've specified that the function '",
            aggregation.functions[[coded.name]],
            "' be called using the column(s) '",
            paste(arg.names, collapse="', '"),
            "' as arguments. Make sure this what you intended.")
        }
        tryCatch(aggregated.context[area,1] <-
                 do.call(aggregation.functions[[coded.name]], arguments),
                 error=error.handler)
      }

      ## store the aggregated context in merge.data for later merge with
      ## the rest, compute the appropriate renaming of the contextual variables
      count.list[[name]] <- count.list[[name]] + 1
      new.name <- paste(name, ".", count.list[[name]], sep='')

      merge.data[[new.name]] <- aggregated.context
      ## add the contextual variable to the formula
      if (!is.null(formula.str)) {
        formula.str <- paste(formula.str, " + ", new.name, sep=)
      }
    }
  } else {
    merge.data <- NULL
  }

  if (is.null(formula.str)) {
    return(merge.data)
  } else {
    return(list(merge.data, formula.str))
  }
}

## resamples contextual data at the individual level as used for stratified
## resampling used in descrw and in srawe
resample <- function(area.list){
  sampler <- function(group){
    size <- length(group)
    row.slice <- sample(1:size, size=size, replace=TRUE)
    return(group[row.slice])
  }
  split.contexts = split(1:length(area.list), area.list)
  resampled.groups <- lapply(split.contexts, sampler)
  ## cont.dat$lino=1:nrow(cont.dat)
  ## return(do.call(rbind,
  ##                lapply(split(cont.dat, area.list),
  ##                       function(x) x[sample(1:nrow(x), size=nrow(x), replace=TRUE),])))
  return(do.call(c, resampled.groups))
}



## a stateful iterator which yields either presampled samples or creates them
## on the fly
## if the number of resamples is 0, then the orignal data order is preserved
sampleIterator <- function(sample.seed, area.list = NULL,
                           nb.resamples = NULL)
{
  if (is.null(area.list)) {
    return(iter(sample.seed))
  } else if (nb.resamples == 0) {
    index <- 0
    nextEl <- function() {
      if (index > 0) {
        stop("StopIteration")
      }
      index <<- index + 1
      return(1:length(area.list))
    }
  } else {
    .Random.seed <<- sample.seed
    ran.seed <- sample.seed
    index <- 1
    nextEl <- function() {
      if (index <= nb.resamples) {
        index <<- index + 1
      } else {
        stop("StopIteration")
      }
      current.seed <- .Random.seed
      .Random.seed <<- ran.seed
      sample <- resample(area.list)
      ran.seed <<- .Random.seed
      .Random.seed <<- current.seed
      return(sample)
    }
  }
  obj <- list(nextElem=nextEl)
  class(obj) <- c("random.sample.iterator", "abstractiter", "iter")
  return(obj)
}


wt.gini <- function(data, weights=rep(1, length(data))){
  s = which(is.finite(data * weights))
  weights = weights[s]
  data = data[s]

  odata <- order(data)
  data <- data[odata]
  weights <- weights[odata]/sum(weights)
  p <- cumsum(weights)
  nu <- cumsum(weights*data)
  n <- length(nu)
  nu <- nu / nu[n]
  sum(nu[-1]*p[-n]) - sum(nu[-n]*p[-1])
}

wt.gini.group <- function(data, weights=rep(1, length(data)), groups) {
  ## page 104 in Horizontalinequalities and conflict
  x.bar <- weighted.mean(data, weights)
  unique.groups <- unique(groups)
  nb.groups <- length(unique.groups)

  x.bar.group <- numeric(nb.groups)
  population.shares <- numeric(nb.groups)
  ##compute weighted group means
  for (i in 1:nb.groups) {
    group <- unique.groups[i]
    group.data <- data[groups==group]
    group.weights <- weights[groups==group]
    x.bar.group[i] <- weighted.mean(group.data, group.weights)
    population.shares[i] <- sum(group.weights)
  }

  ## gini
  gini <- wt.gini(x.bar.group, population.shares)
  return(gini)
}


wt.GRI <- function(data, weights= rep(1, length(data)), groups) {
    unique.groups <- unique(groups)
    nb.groups <- length(unique.groups)
    nb.binoms <- nb.groups*(nb.groups-1)/2

    r <- numeric(nb.groups)
    for (i in 1:nb.groups) {
        group <- unique.groups[i]
        group.data <- data[groups==group]
        group.weights <- weights[groups==group]
        r[i] <- weighted.mean(group.data, group.weights)
    }
    one.minus.GRI <- 0
    for (i in 1:(nb.groups-1)) {
        for (j in (i+1):nb.groups) {
            one.minus.GRI <- one.minus.GRI + abs(r[i]-r[j])
        }
    }
    max.diff = 0
    for (i in 1:nb.groups) {
        max.diff = max.diff + floor(i/2)
    }
    m <- max.diff/nb.binoms
    GRI <- 1-one.minus.GRI/(nb.binoms*m)
    return(GRI)
}
## Gini when x takes a finite number of value

wt.gini.categ <- function(data, weights=rep(1, length(data)))
  {
    weights<- weights/sum(weights)
    Newframe <- aggregate(list(poids=weights),by=list(value=data),FUN=sum)
    1-sum(Newframe$poids^2)
  }

wt.var <- function (data, weights=rep(1, length(data)))
{
  s = which(is.finite(data + weights))
  weights = weights[s]
  data = data[s]
  databar = weighted.mean(data, weights)
  return(sum(weights * (data - databar)^2) *
         (sum(weights)/(sum(weights)^2 - sum(weights^2))))
}



wt.sd <- function (data, weights=rep(1, length(data)))
{
  return(sqrt(wt.var(data, weights)))
}


wt.RS <- function(data, weights=rep(1, length(data)))
{
  s = which(is.finite(data))
  weights = weights[s]
  data = data[s]
  weightsotal <- sum(weights)
  databar <- sum(weights * data)/weightsotal
  n <- sum( abs((data - databar)/databar)*weights)/(2*weightsotal)
  return(n)
}


wt.Atkinson <- function(data, weights=rep(1, length(data)))
{

  s = which(is.finite(data * weights))
  weights = weights[s]
  data = data[s]
  databar <- weighted.mean(data,weights)
  datalogbar <- weighted.mean(log(data),weights)


  parameter <- 0.5
  if(parameter==1)
    A <- 1 - (exp(datalogbar)/databar)
  else
    {
      data <- sum(weights*(data/databar)^(1-parameter))/sum(weights)
      A <- 1 - data^(1/(1-parameter))
    }
  return(A)
}


wt.Theil <- function(data, weights=rep(1, length(data)))
{
  s = which(is.finite(data * weights))
  weights = weights[s]
  data = data[s]

  t = which(!(data==0))
  weights = weights[t]
  data = data[t]

  databar <- weighted.mean(data,weights)
  d <- data/databar*log(data/databar)
  Th <-weighted.mean(d,weights)

  return(Th)
}


##This closure returns the functions fillMatrices, getRanefs, getTmpList
randomEffectsClosure <- function(nb.area, nb.resamples) {
  ## temporary data structures to gather the results
  aic.results <- matrix(nrow=5, ncol=nb.resamples)
  row.names(aic.results) <- c('AIC', 'BIC', 'logLik', 'deviance', 'REMLdev')
  fixed.effect.results <- NULL # will have to be filled at first occurence
  random.var.results <- NULL # will have to be filled at first occurence
  ranefs <- matrix(nrow=nb.area, ncol=nb.resamples)
  beta.results <- NULL # will have to be filled at first occurence

  n <- 0
  nn <- 0

  prepareRandomEffects <- function(lme.obj, column) {
    vc <- lme4::VarCorr(lme.obj)
    if (is.null(random.var.results)) {
      n <<- nrow(vc[[1]])
      nn <<- n*(n-1)/2

      colnames <- character((1+n)*n+1) # to store 2 sd, (n*n-n)/2+n cov, (n*n-n)/2 corr and sc

      colnames[1] <- 'sd.e'
      rand.names <-row.names(vc[[1]])
      colnames[2:(n+1)] <-
        lapply(rand.names,
               FUN=function(x){return(paste("VAR.", x, sep=""))})
      colnames[(n+2):(2*n+1)] <-
        lapply(rand.names,
               FUN=function(x){return(paste("SD.", x, sep=""))})

      if (n>1) {
        ind <- 2*n+1
        for (i in 2:n) {
          for (j in 1:(i-1)) {
            ind <- ind + 1
            colnames[ind] <-
              paste( "COV.", rand.names[i], ".", rand.names[j], sep="")
            colnames[ind+nn] <-
              paste("CORR.", rand.names[i], ".", rand.names[j], sep="")
          }
        }
      }

      random.var.results <<- matrix(nrow=(1+n)*n+1, ncol=nb.resamples)
      row.names(random.var.results) <<- colnames
    }

    random.var.results[1, column] <<- attr(vc, "sc")
    mat <- vc[[1]]
    random.var.results[2:(n+1), column] <<- diag(mat)
    random.var.results[(n+2):(2*n+1), column] <<- attr(mat,"stddev")

    if (n>1) {
      ind <- 2*n+1
      for (i in 2:n) {
        for (j in 1:(i-1)) {
          ind <- ind + 1
          random.var.results[ind, column] <<- mat[i,j]
          random.var.results[ind+nn, column] <<- attr(mat, "correlation")[i,j]
        }
      }
    }
  }

  ## inserts the analyses at the right positions in the result matrices
  fillMatrices <- function(lme.obj, beta, index) {
    if (is.null(fixed.effect.results)) {
      ##create the fixed.effect.results matrix of the right dims and name rows
      fix <- lme4::fixef(lme.obj)
      fixed.effect.results <<- matrix(nrow=length(fix),
                                      ncol=nb.resamples)
      row.names(fixed.effect.results) <<- names(fix)

      ## TODO check whether we can make beta conditional when used with empty model
      if (length(beta)>0) {
        beta.results <<- matrix(nrow=length(beta), ncol=nb.resamples)
      }
    }
    fixed.effect.results[, index] <<- lme4::fixef(lme.obj)
    prepareRandomEffects(lme.obj, index)

    aic.results[1, index] <<- AIC(lme.obj)
    aic.results[2, index] <<- BIC(lme.obj)
    aic.results[3, index] <<- as.numeric(logLik(lme.obj))
    aic.results[4, index] <<- deviance(lme.obj, REML=FALSE)
    aic.results[5, index] <<- deviance(lme.obj, REML=TRUE)
    ## TODO check whether we can make beta conditional when used with empty model
    if (!is.null(beta.results)){
      beta.results[, index] <<- beta
      rownames(beta.results) <<- names(beta)
    }


    ranefs[ ,index] <<- lme4::ranef(lme.obj)[[1]][,1]
  }

  getRanefs <- function(){
    return(ranefs)
  }

  getTmpList <- function(){
    tmp.list <- list(fixed=fixed.effect.results,
                model.fit=aic.results,
                random.var=random.var.results)
    if (!is.null(beta.results)) {
      tmp.list[["betas"]] <- beta.results
    }

    return(tmp.list)
  }

  getBetaResults <- function(){
    return(beta.results)
  }

  return(list(fillMatrices=fillMatrices,
              getRanefs=getRanefs,
              getTmpList=getTmpList,
              getBetaResults=getBetaResults))
}

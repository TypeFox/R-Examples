###############################################################################
#    rsubgroup package R classes
# 
#    This file is part of the rsubgroup package.
#    Copyright (C) 2011-2014 by Martin Atzmueller
#    
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
    
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Contact: Martin Atzmueller (martin@atzmueller.net)
###############################################################################

library(foreign)
library(rJava)

setGeneric(".CreateARFFProvider",
    function(source, name, ...) {
      standardGeneric(".CreateARFFProvider")
    }
)

setMethod(".CreateARFFProvider", signature(source = "data.frame", name = "character"),
    function(source, name, ...) {
      # Creates a dataset provider (converting the dataframe)
      arff.vector <- NULL # This prevents the R CHECK NOTE:
      # 'No visible bindingfor global variable Note in R CMD check'
      
      con <- textConnection("arff.vector", open = "w", local = TRUE)
      write.arff(source, con)
      flush(con)
      close(con)
      rm(con)
      arff <- paste(arff.vector, "", collapse="\n")
      provider <- .jnew("org/vikamine/kernel/xpdl/ARFFAsStringDatasetProvider", arff, name)
      rm(arff.vector)
      return(provider)
    }
)

setMethod(".CreateARFFProvider", signature(source = "character", name = "character"),
    function(source, name, ...) {
      # Creates a dataset provider given a file name
      provider <- .jnew("org/vikamine/kernel/xpdl/FileDatasetProvider", source)
      return(provider)
    }
)

.CreateOntologyForData <- function(provider, dataset) {
  # Creates the ontology object for the respective dataset
  ontology <- .jcall(provider, "Lorg/vikamine/kernel/data/Ontology;","getDataset", dataset)
  return(ontology)
}

.CreateSimpleSDTask <- function(ontology, target) {
  # Creates a simple subgroup discovery task
  .FreeMemory()
  simpleTask <- new(J("org/vikamine/kernel/subgroup/search/SDSimpleTask"), ontology)
  if (!is.null(target$value)) {
    selector <- new(J("org/vikamine/kernel/subgroup/selectors/DefaultSGSelector"), ontology, target$attribute, target$value)
    target <- new(J("org/vikamine/kernel/subgroup/target/SelectorTarget"), selector)
  } else {
    attribute <- J(ontology, "getAttribute", target$attribute)
    target <- new(J("org/vikamine/kernel/subgroup/target/NumericTarget"), attribute)
  }
  J(simpleTask, "setTarget", target)
  return(simpleTask)
}

CreateSDTask <- function(source, target, config = new("SDTaskConfig")) {
  # Creates a subgroup discovery task
  #
  # Args:
  #   source: A data source, i.e., dataframe or file (name)
  #   target: The target variable
  #   config: A SDTaskConfig
  #
  # Returns:
  #   A subgroup discovery task
  .FreeMemory()
  provider <- .CreateARFFProvider(source, "data")
  ontology <- .CreateOntologyForData(provider, "data")
  task <- .CreateSimpleSDTask(ontology, target)
  J(task, "setQualityFunction", config@qf)
  J(task, "setSDMethod", config@method)
  J(task, "setMaxSGCount", as.integer(config@k))
  J(task, "setMinQualityLimit", as.double(config@minqual))
  J(task, "setMinSubgroupSize", as.double(config@minsize))
  J(task, "setMaxSGDSize", as.integer(config@maxlen))
  J(task, "setSuppressStrictlyIrrelevantSubgroups", config@relfilter)
  J(task, "setIgnoreDefaultValues", config@nodefaults)
  if (config@postfilter != "") {
    J(task, "setPostFilter", config@postfilter)
  }
  if (is.null(config@attributes)) {
    attributesArrayObject <- .GetAllAttributesAsJArray(ontology = ontology)
    J(task, "setAttributes", attributesArrayObject)
  } else if ((!is.null(config@attributes)) && (length(config@attributes) > 0)) {
    J(task, "setAttributes", .jarray(config@attributes))  
  } else {
    J(task, "setAttributes", .jarray(character(0))) 
  }
  return(task)
}

as.target <- function(attribute=NULL, value=NULL) {
  # Creates a target variable object given attribute and value (for nominals)
  #
  # Args:
  #   attribute: The respective attribute
  #   value: For nominals, the respective value; for numeric NULL
  #
  # Returns:
  #   A target object representation
  if (!is.null(attribute) && !is.null(value))
    return(list(attribute=attribute, value=value))
  else if (!is.null(attribute))
    return(list(attribute=attribute))
  else
    return(NULL)
}

.GetParameters <- function(task, sg) {
  target <- J(task, "getTarget")
  if (J(target, "isBoolean")) {
    size <- J(J(sg, "getStatistics"), "getSubgroupSize")
    p <- J(J(sg, "getStatistics"), "getP")
    p0 <- J(J(sg, "getStatistics"), "getP0")
    return(list(p = p, p0 = p0, size = size))
  } else if (J(target, "isNumeric")) {
    size <- J(J(sg, "getStatistics"), "getSubgroupSize")
    mean <- J(J(sg, "getStatistics"), "getSGMean")
    popMean <- J(J(sg, "getStatistics"), "getPopulationMean")
    return(list(mean = mean, populationMean = popMean, size = size))
  } else {
    stop("Unknown target")
  }
}

.ConvertDescription <- function(sgDescription) {
  # Internal function for converting a (Java) SGDescription consisting
  # of a set of selection expressions into a character vector of strings
  # representing these
  # Args:
  #   sgDescription: A (Java) SGDescription object
  #
  # Returns:
  #   A character vector
  sgSelectorArray <-
      J("org.vikamine.kernel.subgroup.search.SDSimpleTask")$getSimpleDescription(sgDescription)
  return(as.character(sgSelectorArray))
}

DiscoverSubgroupsByTask <- function(task, as.df = FALSE) {
  # Internal function for setting up and performing subgroup discovery
  # Args:
  #   task: A subgroup discovery task
  #
  # Returns:
  #   A list of subgroup patterns
  sgSet <- J(task, "performSubgroupDiscovery")
  sgList <- J(sgSet, "toSortedList", FALSE)
  sgArray <- .jevalArray(J(sgList, "toArray"))
  
  patterns = list()
  for (sg in sgArray) {
    #description <- as.character(J(J(sg, "getSGDescription"), "getDescription"))
    sgDescription <- J(sg, "getSGDescription")
    description <- .ConvertDescription(sgDescription)
    quality <- J(sg, "getQuality")
    size <- J(J(sg, "getStatistics"), "getSubgroupSize")
    parameters = .GetParameters(task, sg)    
    pattern <- new("Pattern", description=description, quality=quality, size=size, parameters=parameters)
    patterns = append(patterns, pattern)
  }
  
  if (as.df) {
    dataFrameRules <- ToDataFrame(patterns)
    return(dataFrameRules)
  } else {
    return(patterns)
  }
}

DiscoverSubgroups <- function(source, target, config=new("SDTaskConfig"), as.df=FALSE) {
  # Performs subgroup discovery according to target and config on data
  #
  # Args:
  #   data: A dataframe
  #   target: A target variable (constructed by as.target)
  #   config: a SDTaskConfig configuration for the algorithm
  #
  # Returns:
  #   A list of subgroup patterns
  task <- CreateSDTask(source, target, config)
  result <- DiscoverSubgroupsByTask(task, as.df)
  return(result)
}


.FormatDoubleSignificantDigits <- function(double, ndigits=2) {
  # Internal function: Prints double according to
  if (is.numeric(ndigits)) {
    sprintf(paste("%.", ndigits, "f", sep=""), double)
  } else {
    double
  }
}



ToDataFrame <- function(patterns, ndigits=2) {
  # Transforms a list/vector of patterns into a dataframe
  #
  # Args:
  #   patterns: List of patterns
  #   ndigits: Number of significant digits for floats
  #
  # Returns:
  #   The dataframe containing the pattern information
  isNumeric = FALSE  
  descriptions <- list()
  length(descriptions) = length(patterns)
  qualities <-list()
  length(qualities) = length(patterns)
  sizes <- list()
  length(sizes) <- length(patterns)
  ps <- list()
  
  i = 1
  for (pattern in patterns) {
    descriptions[i] = paste(pattern@description, collapse=", ")
    qualities[i] = .FormatDoubleSignificantDigits(pattern@quality, ndigits)
    sizes[i] = pattern@size
    if (!is.null(pattern@parameters$mean)) {
      ps[i] = .FormatDoubleSignificantDigits(pattern@parameters$mean, ndigits)
      isNumeric = TRUE
    } else {
      ps[i] = .FormatDoubleSignificantDigits(pattern@parameters$p, ndigits)
      isNumeric = FALSE
    }
    i = i + 1
  }
  if (isNumeric) {
    dataframe <- data.frame(
        quality=as.vector(qualities, "numeric"),
        mean=as.vector(ps, "numeric"), 
        size=as.vector(sizes, "numeric"),
        description=as.vector(descriptions, "character"))
  } else {
    dataframe <- data.frame(
        quality=as.vector(qualities, "numeric"),
        p=as.vector(ps, "numeric"), 
        size=as.vector(sizes, "numeric"),
        description=as.vector(descriptions, "character"))
  }
  return(dataframe)
}

.FreeMemory <- function(...) {
  # Call the R garbage collection
  # Then call Java garbage collection
  gc(...)
  .jcall("java/lang/System", method = "gc")
  invisible()
}
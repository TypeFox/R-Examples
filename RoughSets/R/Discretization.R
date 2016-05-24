#############################################################################
#
#  This file is a part of the R package "RoughSets".
#
#  Author: Lala Septem Riza and Andrzej Janusz
#  Supervisors: Chris Cornelis, Francisco Herrera, Dominik Slezak and Jose Manuel Benitez
#  Copyright (c):
#       DiCITS Lab, Sci2s group, DECSAI, University of Granada and
#       Institute of Mathematics, University of Warsaw
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' It is a wrapper function for all discretization methods based on RST.
#' It provides an interface that allows users to use the discretization methods easily.
#'
#' The discretization is used to convert numeric attributes into nominal ones
#' in an information system. It is usually a preliminary step for the most of methods based on the rough set theory, which
#' need nominal attributes, for exemple, to compute the indiscernibility relation.
#'
#' Output of this function is an object of a class \code{Discretization} which
#' contains cut values. The function \code{\link{SF.applyDecTable}} can be used
#' to generate a new (discretized) decision table from the computed cuts. Type of
#' all attributes in the resulting table will be changed into nominal (i.e. ordered factors).
#'
#' All implemented supervised discretization methods need a nominal decision attribute.
#' Furthermore, especially for the method type \code{"global.discernibility"}, all conditional attributes
#' must be numeric. A different method needs to be chosen in a case when
#' a data set contains attributes of mixed types (numeric and nominal).
#'
#' @title The wrapper function for discretization methods
#' @author Andrzej Janusz
#'
#' @param decision.table an object inheriting from the \code{"DecisionTable"} class, which represents a decision system.
#'        See \code{\link{SF.asDecisionTable}}.
#' @param type.method a character representing a discretization method to be used in the computations.
#'        Currently it can be one of the following methods:
#'         \itemize{
#'           \item \code{"global.discernibility"}: See \code{\link{D.global.discernibility.heuristic.RST}}.
#'           \item \code{"local.discernibility"}: See \code{\link{D.local.discernibility.heuristic.RST}}.
#'           \item \code{"unsupervised.intervals"}: See \code{\link{D.discretize.equal.intervals.RST}}.
#'           \item \code{"unsupervised.quantiles"}: See \code{\link{D.discretize.quantiles.RST}}.
#'        }
#' @param  ... parameters that are passed to the discretization methods. See the manual of particular functions.
#'
#' @return An object of a class \code{"Discretization"} which stores cuts for each conditional attribute. It contains the following components:
#' \itemize{
#'   \item \code{cut.values}: a list representing cut values for each of numeric attributes. NULL value means that
#'         no cut was selected for a given attribute.
#'   \item \code{type.method}: the type of method which is used to define cut values.
#'   \item \code{type.task}: the type of task which is \code{"discretization"}.
#'   \item \code{model}: the type of model which is \code{"RST"}.
#' }
#'
#' @seealso \code{\link{BC.LU.approximation.RST}}, \code{\link{FS.reduct.computation}}, \code{\link{SF.applyDecTable}}.
#'
#' @examples
#' #################################################################
#' ## Example: Determine cut values and generate new decision table
#' #################################################################
#' data(RoughSetData)
#' wine.data <- RoughSetData$wine.dt
#' cut.values1 <- D.discretization.RST(wine.data,
#'                                    type.method = "unsupervised.quantiles",
#'                                    nOfIntervals = 3)
#'
#' ## generate a new decision table
#' wine.discretized1 <- SF.applyDecTable(wine.data, cut.values1)
#' dim(wine.discretized1)
#' lapply(wine.discretized1, unique)
#'
#' cut.values2 <- D.discretization.RST(wine.data,
#'                                     type.method = "global.discernibility")
#'
#' wine.discretized2 <- SF.applyDecTable(wine.data, cut.values2)
#' dim(wine.discretized2)
#' lapply(wine.discretized2, unique)
#'
#' @export
D.discretization.RST <- function(decision.table, type.method = "unsupervised.quantiles", ...){
	if (!(type.method %in% c("global.discernibility", "local.discernibility",
                           "unsupervised.intervals", "unsupervised.quantiles"))) {
    stop("Unrecognized discretization type.")
	}

	if (!inherits(decision.table, "DecisionTable")) {
		stop("Provided data should inherit from the \'DecisionTable\' class.")
	}

	if (is.null(attr(decision.table, "decision.attr"))) {
		decisionIdx = ncol(decision.table) + 1
	} else decisionIdx = attr(decision.table, "decision.attr")

	if (all(attr(decision.table, "nominal.attrs")[-decisionIdx])) {
		stop("All the conditional attributes are already nominal.")
	} else {
		if(any(attr(decision.table, "nominal.attrs")[-decisionIdx]) && type.method == "global.discernibility") {
			stop("This discretization method is not implemented for decision tables with mixed attribute types.")
		}
	}

	cut.values = switch(type.method,
    	                unsupervised.quantiles = D.discretize.quantiles.RST(decision.table, ...),
    	                unsupervised.intervals = D.discretize.equal.intervals.RST(decision.table, ...),
                      global.discernibility = D.global.discernibility.heuristic.RST(decision.table, ...),
	                    local.discernibility = D.local.discernibility.heuristic.RST(decision.table, ...) )

	return(cut.values)
}

#' This function implements unsupervised discretization into intervals containing similar number of instances ("quantile-based").
#'
#' This approach belongs to a class of unsupervised discretization methods
#' since it does not consider the class labels. Each numeric attribute is divided in \code{k} intervals which contain approximately
#' the same number of data instances (objects).
#' Detailed information regarding this method can be found in (Dougherty et al, 1995).
#'
#' It should be noted that the output of this function is an object of a class \code{"Discretization"}
#' which contains the cut values.
#' The function \code{\link{SF.applyDecTable}} has to be used in order to generate the new (discretized) decision table.
#'
#' @title The quantile-based discretization
#' @author Andrzej Janusz
#'
#' @param decision.table an object inheriting from the \code{"DecisionTable"} class, which represents a decision system.
#'        See \code{\link{SF.asDecisionTable}}.
#' @param nOfIntervals a positive integer giving the number of intervals.
#'
#' @return An object of a class \code{"Discretization"} which stores cuts for each conditional attribute.
#'         See \code{\link{D.discretization.RST}}.
#'
#' @seealso \code{\link{D.discretize.equal.intervals.RST}}, \code{\link{D.global.discernibility.heuristic.RST}},
#'          \code{\link{D.local.discernibility.heuristic.RST}}, \code{\link{SF.applyDecTable}}.
#'          A wrapper function for all available discretization methods: \code{\link{D.discretization.RST}}
#'
#' @references
#' J. Dougherty, R. Kohavi, and M. Sahami, "Supervised and Unsupervised Discretization of Continuous Features",
#' In A. Prieditis & S. J. Russell, eds. Work. Morgan Kaufmann, p. 194-202 (1995).
#'
#' @examples
#' #################################################################
#' ## Example: Determine cut values and generate new decision table
#' #################################################################
#' data(RoughSetData)
#' wine.data <- RoughSetData$wine.dt
#' cut.values <- D.discretize.quantiles.RST(wine.data, nOfIntervals = 5)
#'
#' ## generate a new decision table
#' wine.discretized <- SF.applyDecTable(wine.data, cut.values)
#' dim(wine.discretized)
#' lapply(wine.discretized, unique)
#'
#' @export
D.discretize.quantiles.RST <- function(decision.table, nOfIntervals = 4) {

  nominalAttrs = attr(decision.table, "nominal.attrs")
  if (!is.null(attr(decision.table, "decision.attr"))) {
    nominalAttrs = nominalAttrs[-attr(decision.table, "decision.attr")]
    decision.table = decision.table[-attr(decision.table, "decision.attr")]
  }

  if (nOfIntervals < 1 || nOfIntervals > nrow(decision.table)) {
    stop("Wrong value of the parameter nOfIntervals.")
  }

  if(sum(nominalAttrs) > 0) {
    cutsList = list()
    cutsList[1:ncol(decision.table)] = list(numeric())
    cutsList[!nominalAttrs] = lapply(decision.table[!nominalAttrs], discretize.quantiles, n = nOfIntervals)
  } else cutsList = lapply(decision.table, discretize.quantiles, n = nOfIntervals)

	cutsList = list(cut.values = cutsList, type.method = "unsupervised.quantiles",
                type.task = "discretization", model = "RST")
	cutsList = ObjectFactory(cutsList, classname = "Discretization")
	return(cutsList)
}

#' This function implements unsupervised discretization into intervals of equal size.
#'
#' This approach belongs to a class of unsupervised discretization methods
#' since it does not consider the class labels. Each numeric attribute is divided in \code{k} intervals of equal length.
#' Detailed information regarding this method can be found in (Dougherty et al, 1995).
#'
#' It should be noted that the output of this function is an object of a class \code{"Discretization"}
#' which contains the cut values.
#' The function \code{\link{SF.applyDecTable}} has to be used in order to generate the new (discretized) decision table.
#'
#' @title Unsupervised discretization into intervals of equal length.
#' @author Andrzej Janusz
#'
#' @param decision.table an object inheriting from the \code{"DecisionTable"} class, which represents a decision system.
#'        See \code{\link{SF.asDecisionTable}}.
#' @param nOfIntervals a positive integer giving the number of intervals.
#'
#' @return An object of a class \code{"Discretization"} which stores cuts for each conditional attribute.
#'         See \code{\link{D.discretization.RST}}.
#'
#' @seealso \code{\link{D.discretize.quantiles.RST}}, \code{\link{D.global.discernibility.heuristic.RST}},
#'          \code{\link{D.local.discernibility.heuristic.RST}}, \code{\link{SF.applyDecTable}}.
#'          A wrapper function for all available discretization methods: \code{\link{D.discretization.RST}}
#'
#' @references
#' J. Dougherty, R. Kohavi, and M. Sahami, "Supervised and Unsupervised Discretization of Continuous Features",
#' In A. Prieditis & S. J. Russell, eds. Work. Morgan Kaufmann, p. 194-202 (1995).
#'
#' @examples
#' #################################################################
#' ## Example: Determine cut values and generate new decision table
#' #################################################################
#' data(RoughSetData)
#' wine.data <- RoughSetData$wine.dt
#' cut.values <- D.discretize.equal.intervals.RST(wine.data, nOfIntervals = 3)
#'
#' ## generate a new decision table
#' wine.discretized <- SF.applyDecTable(wine.data, cut.values)
#' dim(wine.discretized)
#' lapply(wine.discretized, unique)
#'
#' @export
D.discretize.equal.intervals.RST <- function(decision.table, nOfIntervals = 4) {

  nominalAttrs = attr(decision.table, "nominal.attrs")
  if (!is.null(attr(decision.table, "decision.attr"))) {
    nominalAttrs = nominalAttrs[-attr(decision.table, "decision.attr")]
		decision.table = decision.table[-attr(decision.table, "decision.attr")]
	}

  if (nOfIntervals < 1 || nOfIntervals > nrow(decision.table)) {
    stop("Wrong value of the parameter nOfIntervals.")
  }

  if(sum(nominalAttrs) > 0) {
	  cutsList = list()
	  cutsList[1:ncol(decision.table)] = list(numeric())
    cutsList[!nominalAttrs] = lapply(decision.table[!nominalAttrs], discretize.equal.intervals, n = nOfIntervals)
  } else cutsList = lapply(decision.table, discretize.equal.intervals, n = nOfIntervals)

	cutsList = list(cut.values = cutsList, type.method = "unsupervised.intervals",
                type.task = "discretization", model = "RST")
	cutsList = ObjectFactory(cutsList, classname = "Discretization")
	return(cutsList)
}

#' It is a function used for computing globally semi-optimal cuts using the maximum discernibility heuristic.
#'
#' A complete description of the implemented algorithm can be found in (Nguyen, 2001).
#'
#' It should be noted that the output of this function is an object of a class \code{"Discretization"}
#' which contains the cut values.
#' The function \code{\link{SF.applyDecTable}} has to be used in order to generate the new (discretized) decision table.
#'
#' @title Supervised discretization based on the maximum discernibility heuristic
#' @author Andrzej Janusz
#'
#' @param decision.table an object inheriting from the \code{"DecisionTable"} class, which represents a decision system.
#'        See \code{\link{SF.asDecisionTable}}.
#'        It should be noted that for this particular method all conditional attributes
#'        must be numeric.
#' @param maxNOfCuts a positive integer indicating the maximum number of allowed cuts.
#' @param attrSampleSize an integer between 1 and the number of conditional attributes (the default). It indicates
#'        the attribute sample size for the Monte Carlo selection of candidating cuts.
#' @param cutCandidatesList an optional list containing candidates for optimal cut values.
#'        By default the candidating cuts are determined automatically.
#' @param discFunction a function used for computation of cuts. Currently only one implementation of maximu discernibility heuristic
#'        is available (the default). However, this parameter can be used to integrate custom implementations of
#'        discretization functions with the \code{RoughSets} package.
#' @param ... additional parameters to the \code{discFunction} (currently unsupported).
#'
#' @seealso \code{\link{D.discretize.quantiles.RST}}, \code{\link{D.discretize.equal.intervals.RST}},
#'          \code{\link{D.local.discernibility.heuristic.RST}} and \code{\link{SF.applyDecTable}}.
#'          A wrapper function for all available discretization methods: \code{\link{D.discretization.RST}}
#'
#' @return An object of a class \code{"Discretization"} which stores cuts for each conditional attribute.
#'         See \code{\link{D.discretization.RST}}.
#'
#' @references
#' S. H. Nguyen, "On Efficient Handling of Continuous Attributes in Large Data Bases",
#' Fundamenta Informaticae, vol. 48, p. 61 - 81 (2001).
#'
#' Jan G. Bazan, Hung Son Nguyen, Sinh Hoa Nguyen, Piotr Synak, and Jakub Wroblewski,
#' "Rough Set Algorithms in Classification Problem", Chapter 2
#' In: L. Polkowski, S. Tsumoto and T.Y. Lin (eds.): Rough Set Methods and Applications
#' Physica-Verlag, Heidelberg, New York, p. 49 - 88 ( 2000).
#'
#' @examples
#' #################################################################
#' ## Example: Determine cut values and generate new decision table
#' #################################################################
#' data(RoughSetData)
#' wine.data <- RoughSetData$wine.dt
#' cut.values <- D.global.discernibility.heuristic.RST(wine.data)
#'
#' ## generate a new decision table:
#' wine.discretized <- SF.applyDecTable(wine.data, cut.values)
#' dim(wine.discretized)
#' lapply(wine.discretized, unique)
#'
#' ## remove attributes with only one possible value:
#' to.rm.idx <- which(sapply(lapply(wine.discretized, unique), function(x) length(x) == 1))
#' to.rm.idx
#' wine.discretized.reduced <- wine.discretized[-to.rm.idx]
#' dim(wine.discretized.reduced)
#'
#' ## check whether the attributes in the reduced data are a super-reduct of the original data:
#' colnames(wine.discretized.reduced)
#' class.idx <- which(colnames(wine.discretized.reduced) == "class")
#' sum(duplicated(wine.discretized.reduced)) == sum(duplicated(wine.discretized.reduced[-class.idx]))
#' ## yes it is
#'
#' @export
D.global.discernibility.heuristic.RST <- function(decision.table, maxNOfCuts = 2*ncol(decision.table),
                                                 attrSampleSize = ncol(decision.table)-1,
                                                 cutCandidatesList = NULL,
                                                 discFunction = global.discernibility, ...)  {

	if (!is.null(attr(decision.table, "decision.attr"))) {
		infoSystem = decision.table[-attr(decision.table, "decision.attr")]
		decisionAttr = factor(decision.table[[attr(decision.table, "decision.attr")]])
	} else {
    stop("A decision attribute is not indicated.")
	}

  if (maxNOfCuts < 1) {
    stop("Wrong value of maxNOfCuts. It should be a positive integer.")
  }

	if (attrSampleSize < 1 || attrSampleSize > ncol(decision.table) - 1) {
	  stop("Wrong value of attrSampleSize. It should be a positive integer, not larger than the number of conditional attributes.")
	}

	if (is.null(cutCandidatesList)) {
		cutCandidatesList = lapply(infoSystem, chooseCutCandidates, decisionAttr)
	} else {
    if (length(cutCandidatesList) != ncol(decision.table) - 1) {
      stop("Wrong length of a list containing candidate cuts. Its length should equal the number of conditional attributes.")
    }
	}

	candidatesCounterVec = sapply(cutCandidatesList, length)
	nonNullIdx = which(sapply(cutCandidatesList, function(x) return(length(x) > 0)))

	cutsList = list()
	cutsList[1:ncol(infoSystem)] = list(numeric())

	tmpCutsList = discFunction(as.list(infoSystem)[nonNullIdx], cutCandidatesList[nonNullIdx],
                             decVec = decisionAttr, nOfCuts = maxNOfCuts,
                             nAttrs = attrSampleSize, ...)
	cutsList[nonNullIdx] = tmpCutsList


	names(cutsList) = colnames(infoSystem)
	cutsList = list(cut.values = cutsList, type.method = "global.discernibility",
                  type.task = "discretization", model = "RST")
	cutsList = ObjectFactory(cutsList, classname = "Discretization")

	return(cutsList)
}


#' It is a function used for computing locally semi-optimal cuts using the local discernibility heuristic.
#'
#' A local (univariate) version of the algorithm described in (Nguyen, 2001) and (Bazan et al., 2000).
#'
#' The output of this function is an object of a class \code{"Discretization"}
#' which contains cut values.
#' The function \code{\link{SF.applyDecTable}} has to be used in order to generate the new (discretized) decision table.
#'
#' @title Supervised discretization based on the local discernibility heuristic
#' @author Andrzej Janusz
#'
#' @param decision.table an object inheriting from the \code{"DecisionTable"} class, which represents a decision system.
#'        See \code{\link{SF.asDecisionTable}}.
#'        It should be noted that for this particular method all conditional attributes
#'        must be numeric.
#' @param maxNOfCuts a positive integer indicating the maximum number of allowed cuts on a single attribute.
#' @param cutCandidatesList an optional list containing candidates for optimal cut values.
#'        By default the candidating cuts are determined automatically.
#' @param discFunction a function used for computation of cuts. Currently only one implementation of the local discernibility heuristic
#'        is available (the default). However, this parameter can be used to integrate custom implementations of
#'        discretization functions with the \code{RoughSets} package.
#'
#' @seealso \code{\link{D.discretize.quantiles.RST}}, \code{\link{D.discretize.equal.intervals.RST}},
#'          \code{\link{D.global.discernibility.heuristic.RST}} and \code{\link{SF.applyDecTable}}.
#'          A wrapper function for all available discretization methods: \code{\link{D.discretization.RST}}
#'
#' @return An object of a class \code{"Discretization"} which stores cuts for each conditional attribute.
#'         See \code{\link{D.discretization.RST}}.
#'
#' @references
#' S. H. Nguyen, "On Efficient Handling of Continuous Attributes in Large Data Bases",
#' Fundamenta Informaticae, vol. 48, p. 61 - 81 (2001).
#'
#' Jan G. Bazan, Hung Son Nguyen, Sinh Hoa Nguyen, Piotr Synak, and Jakub Wroblewski,
#' "Rough Set Algorithms in Classification Problem", Chapter 2
#' In: L. Polkowski, S. Tsumoto and T.Y. Lin (eds.): Rough Set Methods and Applications
#' Physica-Verlag, Heidelberg, New York, p. 49 - 88 ( 2000).
#'
#' @examples
#' #################################################################
#' ## Example: Determine cut values and generate new decision table
#' #################################################################
#' data(RoughSetData)
#' wine.data <- RoughSetData$wine.dt
#' cut.values <- D.local.discernibility.heuristic.RST(wine.data)
#'
#' ## generate a new decision table:
#' wine.discretized <- SF.applyDecTable(wine.data, cut.values)
#' dim(wine.discretized)
#' lapply(wine.discretized, unique)
#'
#' @export
D.local.discernibility.heuristic.RST <- function(decision.table, maxNOfCuts = 2,
                                                 cutCandidatesList = NULL,
                                                 discFunction = local.discernibility)  {

  if (!is.null(attr(decision.table, "decision.attr"))) {
    decisionIdx = attr(decision.table, "decision.attr")
    infoSystem = decision.table[-decisionIdx]
    decisionAttr = factor(decision.table[[decisionIdx]])
  } else {
    stop("A decision attribute is not indicated.")
  }

  if (maxNOfCuts < 1) {
    stop("Wrong value of maxNOfCuts. It should be a positive integer.")
  }

  if (is.null(cutCandidatesList)) {
    if(any(attr(decision.table, "nominal.attrs")[-decisionIdx])) {
      cutCandidatesList = list()
      cutCandidatesList[1:ncol(infoSystem)] = list(numeric())
      cutCandidatesList[!(attr(decision.table, "nominal.attrs")[-decisionIdx])] =
        lapply(infoSystem[!(attr(decision.table, "nominal.attrs")[-decisionIdx])],
               chooseCutCandidates,
               decisionAttr)
    } else {
      cutCandidatesList = lapply(infoSystem, chooseCutCandidates, decisionAttr)
    }
  } else {
    if (length(cutCandidatesList) != ncol(decision.table) - 1) {
      stop("Wrong length of a list containing candidate cuts. Its length should equal the number of conditional attributes.")
    }
  }

  candidatesCounterVec = sapply(cutCandidatesList, length)
  nonNullIdx = which(sapply(cutCandidatesList, function(x) return(length(x) > 0)))

  cutsList = list()
  cutsList[1:ncol(infoSystem)] = list(numeric())

  tmpCutsList = mapply(discFunction,
                       as.list(infoSystem)[nonNullIdx], cutCandidatesList[nonNullIdx],
                       MoreArgs = list(decVec = decisionAttr, nOfCuts = maxNOfCuts,
                                       nDecisions = length(levels(decisionAttr))),
                       SIMPLIFY = FALSE)
  cutsList[nonNullIdx] = tmpCutsList


  names(cutsList) = colnames(infoSystem)
  cutsList = list(cut.values = cutsList, type.method = "local.discernibility",
                  type.task = "discretization", model = "RST")
  cutsList = ObjectFactory(cutsList, classname = "Discretization")

  return(cutsList)
}


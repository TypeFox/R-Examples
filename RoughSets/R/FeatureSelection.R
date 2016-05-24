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
#' This function is a wrapper for computing different types of decision reducts
#' and approximate decision reducts.
#'
#' The implemented methods include the following approaches:
#' \itemize{
#' \item \code{"greedy.heuristic"}: a greedy heuristic method for computation of decision reducts (or approximate decision reducts) based on RST.
#'                                  See \code{\link{FS.greedy.heuristic.reduct.RST}}.
#'
#' \item \code{"DAAR.heuristic"}: Dynamically Adapted Approximate Reduct heuristic, which is a modification of the greedy heuristic with a random probe test to avoid inclusion of irrelevant attributes to the reduct.
#'                               See \code{\link{FS.DAAR.heuristic.RST}}.
#'
#' \item \code{"nearOpt.fvprs"}: the near-optimal reduction algorithm based on FRST.
#'                               See \code{\link{FS.nearOpt.fvprs.FRST}}.
#'
#' \item \code{"permutation.heuristic"}: a permutation-based elimination heuristic for computation of decision reducts based on RST.
#'                                       See \code{\link{FS.permutation.heuristic.reduct.RST}}.
#' }
#' Those methods can be selected by setting the parameter \code{method}.
#' Additionally, \code{\link{SF.applyDecTable}} has been provided to generate a new decision table.
#'
#' @title The reduct computation methods based on RST and FRST
#' @author Andrzej Janusz
#'
#' @param decision.table  an object of a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}.
#' @param method  a character representing the type of computation method to use. See in Section \code{Details}.
#' @param ...  other parameters. See the parameters of \code{\link{FS.greedy.heuristic.reduct.RST}}, \code{\link{FS.DAAR.heuristic.RST}},
#'        \code{\link{FS.nearOpt.fvprs.FRST}} and \code{\link{FS.permutation.heuristic.reduct.RST}}.
#'
#' @return An object of a class \code{"FeatureSubset"}. See \code{\link{FS.greedy.heuristic.reduct.RST}},
#' \code{\link{FS.DAAR.heuristic.RST}}, \code{\link{FS.permutation.heuristic.reduct.RST}} or
#' \code{\link{FS.nearOpt.fvprs.FRST}} for more details.
#'
#' @seealso \code{\link{D.discretization.RST}}, \code{\link{BC.LU.approximation.RST}}
#'
#' @examples
#' ##############################################################
#' ## Example 1: generate reduct and new decision table
#' ## using RST and FRST
#' ##############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## generate a single reduct using RST
#' reduct.1 <- FS.reduct.computation(decision.table, method = "greedy.heuristic")
#'
#' ## generate a single reduct using FRST
#' reduct.2 <- FS.reduct.computation(decision.table, method = "nearOpt.fvprs")
#'
#' ## generate a new decision table using reduct.1
#' new.decTable.1 <- SF.applyDecTable(decision.table, reduct.1)
#'
#' ## generate new decision table using reduct.2
#' new.decTable.2 <- SF.applyDecTable(decision.table, reduct.2)
#'
#' @export
FS.reduct.computation <- function(decision.table, method = "greedy.heuristic", ...){

	if (!(method %in% c("greedy.heuristic", "DAAR.heuristic", "permutation.heuristic", "nearOpt.fvprs"))) {
		stop("Unrecognized attribute reduction method.")
	}

	if (!inherits(decision.table, "DecisionTable")) {
		stop("Provided data should inherit from the \'DecisionTable\' class.")
	}

	nominal.att <- attr(decision.table, "nominal.attrs")
	if (!all(nominal.att) && !(method %in% "nearOpt.fvprs")) stop("Discretize attribures before computing RST reducts.")

	if (is.null(attr(decision.table, "decision.attr"))) stop("A decision attribute is not indicated.")
	else decIdx = attr(decision.table, "decision.attr")

	if (is.null(attr(decision.table, "desc.attrs"))) stop("No description of attribute values.")
	else desc.attrs = attr(decision.table, "desc.attrs")

	## call the chosen method
	reduct = switch(method,
                  greedy.heuristic = FS.greedy.heuristic.reduct.RST(decision.table,
                                                                    attrDescriptions = desc.attrs,
                                                                    decisionIdx = decIdx, ...),
                  DAAR.heuristic = FS.DAAR.heuristic.RST(decision.table,
                                                         attrDescriptions = desc.attrs,
                                                         decisionIdx = decIdx, ...),
                  permutation.heuristic = FS.permutation.heuristic.reduct.RST(decision.table,
                                                                              decisionIdx = decIdx, ...),
                  nearOpt.fvprs = FS.nearOpt.fvprs.FRST(decision.table, ...) )

	return(reduct)
}


#' It is a function implementing the permutation heuristic approach based on RST.
#'
#' Basically there are two steps in this algorithm which are
#' \itemize{
#' \item generating feature subset as a superreduct: In this step, we choose a subset of attributes that
#'       discern all object from different decision classes. It is done by adding consecutive attributes
#'       in an order defined by a permutation of attribute indices. The permutation can be random
#'       or it can be explicitly given (by the parameter \code{permutation}).
#' \item iterative elimination of attributes from the set obtained in the previous step.
#'       It is done in the reverse order to that, defined by the permutation.
#' }
#' More details regarding this algorithm can be found in (Janusz and Slezak, 2012).
#'
#' Additionally, \code{\link{SF.applyDecTable}} has been provided to generate new decision table.
#'
#' @title The permutation heuristic algorithm for computation of a decision reduct
#' @author Andrzej Janusz
#'
#' @param decision.table  an object of a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}.
#' @param permutation  a logical value, an integer vector or\code{NULL} (the default). If an integer vector with a length
#'        equal the cardinality of the conditional attribute set of the decision table is given (it must contain a permutation of
#'        integers from 1:(ncol(decision.table) - 1) ), then it will define the elimination order. Otherwise, if \code{permutation}
#'        is \code{NULL} or \code{TRUE} a random permutation will be generated. In the case when \code{permutation} is FALSE, the
#'        elimination will be performed in the order of attributes in the decision system.
#' @param decisionIdx  an index of the decision attribute. The default value is the last column of a decision table.
#'
#' @return A class \code{"FeatureSubset"} that contains the following components:
#' \itemize{
#' \item \code{reduct}: a list representing a single reduct. In this case, it could be a superreduct or just a subset of features.
#' \item \code{type.method}: a string representing the type of method which is \code{"permutation.heuristic"}.
#' \item \code{type.task}: a string showing the type of task which is \code{"feature selection"}.
#' \item \code{model}: a string representing the type of model. In this case, it is \code{"RST"} which means rough set theory.
#' \item \code{epsilon}: the approximation threshold.
#' }
#'
#' @seealso \code{\link{FS.quickreduct.RST}} and \code{\link{FS.reduct.computation}}.
#'
#' @references
#' A. Janusz and D. Ślęzak, "Utilization of Attribute Clustering Methods for Scalable Computation of Reducts from High-Dimensional Data"
#'										Proceedings of Federated Conference on Computer Science and Information Systems - FedCSIS, p. 295 - 302 (2012).
#' @examples
#' ###################################################
#' ## Example 1: Generate reduct and new decision table
#' ###################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## generate single reduct
#' res.1 <- FS.permutation.heuristic.reduct.RST(decision.table,
#'                                              permutation = NULL,
#'                                              decisionIdx = 5)
#' print(res.1)
#'
#' res.2 <- FS.permutation.heuristic.reduct.RST(decision.table,
#'                                              permutation = 4:1,
#'                                              decisionIdx = 5)
#' print(res.2)
#'
#' ## generate new decision table according to the reduct
#' new.decTable <- SF.applyDecTable(decision.table, res.1)
#' @export
FS.permutation.heuristic.reduct.RST <- function(decision.table,
                                                permutation = NULL,
                                                decisionIdx = ncol(decision.table)){

	## get the data
	names.attrs = colnames(decision.table)

	if (is.null(permutation) || (length(permutation) == 1 && permutation == TRUE)) {
		## shuffle the attributes
		permutation = sample((1:ncol(decision.table))[-decisionIdx])
	}	else {
    if (length(permutation) == 1 && permutation == FALSE) {
		## without shuffling
		  permutation = (1:ncol(decision.table))[-decisionIdx]
    } else {
      if (length(permutation) != ncol(decision.table) - 1 ||
            !all((1:ncol(decision.table))[-decisionIdx] %in% permutation)) {
        stop("The permutation does not match the data.")
      }
    }
	}

	## Initialization
	endFlag = FALSE
	iteration = 1

	## check consistency/certainty which is objects included in lower approximations
	## iteration refers to the number of selected variables to be superreducts
	while (!(sum(duplicated(decision.table[permutation[1:iteration]])) == sum(duplicated(decision.table[c(permutation[1:iteration],decisionIdx)])))) {
		iteration = iteration + 1
		if (iteration > length(permutation))
			stop("Inconsistent decision table - this method is currently implemented only for consistent data.\n\t\tTry a different reduct computation algorithm.")
	}

	## elimating process to get reducts
	redIdxs = permutation[1:iteration]
	if (iteration == 1) endFlag = TRUE
	while (!endFlag) {
		if (sum(duplicated(decision.table[redIdxs[-iteration]])) == sum(duplicated(decision.table[c(redIdxs[-iteration],decisionIdx)]))) {
			redIdxs = redIdxs[-iteration]
		}
		iteration = iteration - 1
		if (iteration == 0 | length(redIdxs) == 1) {
		  endFlag = TRUE
    }
	}

	## get reduct
	reduct <- redIdxs[order(redIdxs)]
	names(reduct) <- names.attrs[reduct]

	## construct class
	mod <- list(reduct = reduct, type.method = "permutation.heuristic",
	            epsilon = 0,
	            type.task = "feature selection", model = "RST")

	class.mod <- ObjectFactory(mod, classname = "FeatureSubset")
	return(class.mod)
}

#' This function implements a greedy heuristic algorithm for computing decision reducts
#' (or approximate decision reducts) based on RST.
#'
#' In this implementation, we provided some attribute subset quality measures which can be
#' passed to the algorithm by the parameter \code{qualityF}. Those measures
#' guide the computations in the search for a decision/approximated reduct. They are used to
#' assess amount of information gained after addition of an attribute. For example,
#' \code{X.entropy} corresponds to the information gain measure.
#'
#' Additionally, this function can use the value of \code{epsilon} parameter in order to compute
#' \eqn{\epsilon}-approximate reducts. The \eqn{\epsilon}-approximate can be defined as an
#' irreducable subset of attributes \code{B}, such that:
#'
#' \eqn{Quality_{\mathcal{A}}(B) \ge (1 - \epsilon)Quality_{\mathcal{A}}(A)},
#'
#' where \eqn{Quality_{\mathcal{A}}(B)} is the value of a quality measure (see possible values
#' of the parameter \code{qualityF}) for an attribute subset \eqn{B} in decision table \eqn{\mathcal{A}}
#' and \eqn{\epsilon} is a numeric value between 0 and 1 expressing the approximation threshold.
#' A lot of monographs provide comprehensive explanations about this topics, for example
#' (Janusz and Stawicki, 2011; Slezak, 2002; Wroblewski, 2001) which are used as the references of this function.
#'
#' Finally, this implementation allows to restrain the computational complexity of greedy
#' searching for decision reducts by setting the value of the parameter \code{nAttrs}. If this
#' parameter is set to a positive integer, the Monte Carlo method of selecting candidating
#' attributes will be used in each iteration of the algorithm.
#'
#' @title The greedy heuristic algorithm for computing decision reducts and approximate decision reducts
#' @author Andrzej Janusz
#
#' @param decision.table an object of a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}.
#' @param attrDescriptions a list containing possible values of attributes (columns) in \code{decision.table}. It usually corresponds to \code{attr(decision.table, "desc.attrs")}.
#' @param decisionIdx an integer value representing an index of the decision attribute.
#' @param qualityF a function used for computation of the quality of attribute subsets.
#'        Currently, the following functions are included:
#'        \itemize{
#'        \item \code{X.entropy}: See \code{\link{X.entropy}}.
#'        \item \code{X.gini}: See \code{\link{X.gini}}.
#'        \item \code{X.nOfConflicts}: See \code{\link{X.nOfConflicts}}.
#'        }
#' @param nAttrs an integer between 1 and the number of conditional attributes. It indicates
#'        the attribute sample size for the Monte Carlo selection of candidating attributes.
#'        If set to \code{NULL} (default) all attributes are used and the algorithm changes
#'        to a standard greedy method for computation of decision reducts.
#' @param epsilon a numeric value between [0, 1) representing an approximate threshold. It
#'        indicates whether to compute approximate reducts or not. If it equals 0 (the default)
#'        a standard decision reduct is computed.
#' @param inconsistentDecisionTable logical indicating whether the decision table is suspected
#'        to be inconsistent or \code{NULL} (the default) which indicated that a test should
#'        be made to determine the data consistency.
#'
#' @return A class \code{"FeatureSubset"} that contains the following components:
#' \itemize{
#' \item \code{reduct}: a list representing a single reduct. In this case, it could be a superreduct or just a subset of features.
#' \item \code{type.method}: a string representing the type of method which is \code{"greedy.heuristic"}.
#' \item \code{type.task}: a string showing the type of task which is \code{"feature selection"}.
#' \item \code{model}: a string representing the type of model. In this case, it is \code{"RST"} which means rough set theory.
#' \item \code{epsilon}: the approximation threshold.
#' }
#'
#' @seealso \code{\link{FS.DAAR.heuristic.RST}} and \code{\link{FS.reduct.computation}}.
#'
#' @references
#' A. Janusz and S. Stawicki, "Applications of Approximate Reducts to the Feature Selection Problem",
#' Proceedings of International Conference on Rough Sets and Knowledge Technology ({RSKT}), vol. 6954, p. 45 - 50 (2011).
#'
#' D. Ślęzak, "Approximate Entropy Reducts", Fundamenta Informaticae, vol. 53, no. 3 - 4, p. 365 - 390 (2002).
#'
#' J. Wróblewski, "Ensembles of Classifiers Based on Approximate Reducts", Fundamenta Informaticae, vol. 47, no. 3 - 4, p. 351 - 360 (2001).
#'
#' @examples
#' ###################################################
#' ## Example 1: Evaluate reduct and generate
#' ##            new decision table
#' ###################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## evaluate a single reduct
#' res.1 <- FS.greedy.heuristic.reduct.RST(decision.table, qualityF = X.entropy,
#'                                         epsilon = 0.0)
#'
#' ## generate a new decision table corresponding to the reduct
#' new.decTable <- SF.applyDecTable(decision.table, res.1)
#' @export
FS.greedy.heuristic.reduct.RST <- function(decision.table,
                                           attrDescriptions = attr(decision.table, "desc.attrs"),
                                           decisionIdx = attr(decision.table, "decision.attr"),
                                           qualityF = X.gini, nAttrs = NULL,
                                           epsilon = 0.0, inconsistentDecisionTable = NULL)  {
  toRmVec = decisionIdx
  attrIdxVec = (1:ncol(decision.table))[-toRmVec]

  if (!is.null(nAttrs)) {
    if (nAttrs == 0 || nAttrs > ncol(decision.table) - 1) {
      stop("There is something wrong with data (too little attributes?) or the parameter nAttrs has a wrong value.")
    }
    tmpAttrSub = sample(attrIdxVec, min(nAttrs, length(attrIdxVec)))
  }  else tmpAttrSub = attrIdxVec

  if (epsilon >= 1 || epsilon < 0) {
    stop("Wrong value of the parameter epsilon. It must be within [0,1) interval.")
  }

  INDrelation = list(1:nrow(decision.table))
  INDsizes = nrow(decision.table)
  decisionChaos = compute_chaos(INDrelation, decision.table[[decisionIdx]],
                                attrDescriptions[[decisionIdx]])
  decisionChaos = qualityF(decisionChaos[[1]])
  attrScoresVec = mapply(qualityGain, decision.table[tmpAttrSub], attrDescriptions[tmpAttrSub],
                         MoreArgs = list(decisionVec = decision.table[[decisionIdx]],
                                         uniqueDecisions = attrDescriptions[[decisionIdx]],
                                         INDclassesList = INDrelation,
                                         INDclassesSizes = INDsizes,
                                         baseChaos = decisionChaos,
                                         chaosFunction = qualityF),
                         SIMPLIFY = TRUE, USE.NAMES = FALSE)
  tmpBestIdx = which.max(attrScoresVec)

  selectedAttrIdxVec  = tmpAttrSub[tmpBestIdx]
  INDrelation = compute_indiscernibility(INDrelation,
                                         as.character(decision.table[[selectedAttrIdxVec]]),
                                         attrDescriptions[[selectedAttrIdxVec]])
  attrIdxVec = (1:ncol(decision.table))[-c(selectedAttrIdxVec, toRmVec)]

  if(is.null(inconsistentDecisionTable)) {
    if(sum(duplicated(decision.table)) == sum(duplicated(decision.table[-decisionIdx]))) {
      inconsistentDecisionTable = FALSE
    } else {
      inconsistentDecisionTable = TRUE
    }
  }

  endFlag = FALSE
  iteration = 1
  if(inconsistentDecisionTable) {
    tmpIND = split(1:nrow(decision.table),
                   do.call(paste, decision.table[-decisionIdx]), drop = TRUE)
    totalChaos = compute_chaos(tmpIND,
                               as.character(decision.table[[decisionIdx]]),
                               attrDescriptions[[decisionIdx]])
    totalChaos = sum(sapply(totalChaos, qualityF) * sapply(tmpIND, length) / nrow(decision.table))
    rm(tmpIND)
  } else totalChaos = 0
  totalDependencyInData = decisionChaos - totalChaos - 10^(-16)
  approxThereshold = (1 - epsilon)*totalDependencyInData

  while (!endFlag) {
    contingencyTabs = lapply(INDrelation,
                             function(x,y) table(y[x]),
                             decision.table[[decisionIdx]])
    chaosVec = sapply(contingencyTabs, qualityF)
    if(any(chaosVec == 0)) {
      tmpIdx = which(chaosVec == 0)
      INDrelation = INDrelation[-tmpIdx]
      contingencyTabs = contingencyTabs[-tmpIdx]
      chaosVec = chaosVec[-tmpIdx]
      rm(tmpIdx)
    }
    sumsVec = sapply(contingencyTabs, sum)
    if(length(INDrelation) > 0) tmpChaos = sum(chaosVec*(sumsVec/nrow(decision.table)))
    else tmpChaos = 0
    tmpDependencyInData = decisionChaos - tmpChaos
    rm(contingencyTabs, sumsVec, chaosVec)

    if (approxThereshold <= tmpDependencyInData) {
      endFlag = TRUE
    }	else {
      if (!is.null(nAttrs)) {
        tmpAttrSub = sample(attrIdxVec, min(nAttrs, length(attrIdxVec)))
      }	else tmpAttrSub = attrIdxVec
      INDsizes = sapply(INDrelation, length)
      attrScoresVec = mapply(qualityGain, decision.table[tmpAttrSub], attrDescriptions[tmpAttrSub],
                             MoreArgs = list(decisionVec = decision.table[[decisionIdx]],
                                             uniqueDecisions = attrDescriptions[[decisionIdx]],
                                             INDclasses = INDrelation,
                                             INDclassesSizes = INDsizes,
                                             baseChaos = tmpChaos,
                                             chaosFunction = qualityF),
                             SIMPLIFY = TRUE, USE.NAMES = FALSE)
      tmpBestIdx = which.max(attrScoresVec)
      selectedAttrIdxVec[iteration + 1] = tmpAttrSub[tmpBestIdx]
      INDrelation = compute_indiscernibility(INDrelation,
                                             as.character(decision.table[[tmpAttrSub[tmpBestIdx]]]),
                                             attrDescriptions[[tmpAttrSub[tmpBestIdx]]])
      if(length(INDrelation) == 0) endFlag = TRUE

      attrIdxVec = (1:ncol(decision.table))[-c(selectedAttrIdxVec, toRmVec)]
      iteration = iteration + 1
    }
  }

  if (iteration > 1) {
    endFlag = FALSE
    iteration = iteration - 1
    while (!endFlag) {
      clsContingencyTab = as.matrix(table(do.call(paste, decision.table[selectedAttrIdxVec[-iteration]]), decision.table[[decisionIdx]]))
      tmpChaos = sum(apply(clsContingencyTab, 1, qualityF)*(rowSums(clsContingencyTab)/nrow(decision.table)))
      tmpDependencyInData = decisionChaos - tmpChaos
      if (approxThereshold <= tmpDependencyInData) {
        selectedAttrIdxVec = selectedAttrIdxVec[-iteration]
      }
      iteration = iteration - 1
      if (iteration == 0) endFlag = TRUE
    }
  }

  reduct <- selectedAttrIdxVec[order(selectedAttrIdxVec)]
  names(reduct) <- colnames(decision.table)[reduct]
  mod <- list(reduct = reduct, type.method = "greedy.heuristic",
              epsilon = epsilon,
              type.task = "feature selection", model = "RST")

  class.mod <- ObjectFactory(mod, classname = "FeatureSubset")
  return(class.mod)
}

#' This function implements the Dynamically Adjusted Approximate Reducts heuristic (DAAR)
#' for feature selection based on RST. The algorithm modifies the greedy approach to selecting
#' attributes by introducing an additional stop condition. The algorithm stops when a random
#' probe (permutation) test fails to reject a hypothesis that the selected attribute introduces
#' illusionary dependency in data (in a context of previously selected attributes).
#'
#' As in the case of \code{\link{FS.greedy.heuristic.reduct.RST}} the implementation can use
#' different attribute subset quality functions (parameter \code{qualityF}) and Monte Carlo
#' generation of candidating attributes (parameter \code{nAttrs}).
#'
#' @title The DAAR heuristic for computation of decision reducts
#' @author Andrzej Janusz
#
#' @param decision.table an object of a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}.
#' @param attrDescriptions a list containing possible values of attributes (columns) in \
#'        code{decision.table}. It usually corresponds to \code{attr(decision.table, "desc.attrs")}.
#' @param decisionIdx an integer value representing an index of the decision attribute.
#' @param qualityF a function used for computation of the quality of attribute subsets.
#'        Currently, the following functions are included:
#'        \itemize{
#'        \item \code{X.entropy}: See \code{\link{X.entropy}}.
#'        \item \code{X.gini}: See \code{\link{X.gini}}.
#'        \item \code{X.nOfConflicts}: See \code{\link{X.nOfConflicts}}.
#'        }
#' @param nAttrs an integer between 1 and the number of conditional attributes. It indicates
#'        the attribute sample size for the Monte Carlo selection of candidating attributes.
#'        If set to \code{NULL} (default) all attributes are used and the algorithm changes
#'        to a standard greedy method for computation of decision reducts.
#' @param allowedRandomness a threshold for attribute relevance. Computations will be terminated
#'        when the relevance of a selected attribute fall below this threshold.
#' @param nOfProbes a number of random probes used for estimating the attribute relevance
#'        (see the references).
#' @param permsWithinINDclasses a logical value indicating whether the permutation test
#'        should be conducted within indescernibility classes.
#' @param inconsistentDecisionTable logical indicating whether the decision table is suspected
#'        to be inconsistent or \code{NULL} (the default) which indicated that a test should
#'        be made to determine the data consistency.
#'
#' @return A class \code{"FeatureSubset"} that contains the following components:
#' \itemize{
#' \item \code{reduct}: a list representing a single reduct. In this case, it could be a superreduct or just a subset of features.
#' \item \code{type.method}: a string representing the type of method which is \code{"greedy.heuristic"}.
#' \item \code{type.task}: a string showing the type of task which is \code{"feature selection"}.
#' \item \code{model}: a string representing the type of model. In this case, it is \code{"RST"} which means rough set theory.
#' \item \code{relevanceProbabilities}: an intiger vector with estimated relevances of selected attributes.
#' \item \code{epsilon}: a value between 0 and 1 representing the estimated approximation threshold.
#' }
#'
#' @seealso \code{\link{FS.greedy.heuristic.reduct.RST}} and \code{\link{FS.reduct.computation}}.
#'
#' @references
#' A. Janusz and S. Stawicki, "Applications of Approximate Reducts to the Feature Selection Problem",
#' Proceedings of International Conference on Rough Sets and Knowledge Technology ({RSKT}), vol. 6954, p. 45 - 50 (2011).
#'
#' A. Janusz and D. Ślęzak, "Random Probes in Computation and Assessment of Approximate Reducts",
#' Proceedings of {RSEISP} 2014, Springer, LNCS vol. 8537: p. 53 - 64 (2014).
#'
#' @examples
#' ###################################################
#' ## Example 1: Evaluate reduct and generate
#' ##            new decision table
#' ###################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## evaluate a single reduct
#' res.1 <- FS.DAAR.heuristic.RST(decision.table)
#'
#' ## generate a new decision table corresponding to the reduct
#' new.decTable <- SF.applyDecTable(decision.table, res.1)
#' @export
FS.DAAR.heuristic.RST = function(decision.table,
                                 attrDescriptions = attr(decision.table, "desc.attrs"),
                                 decisionIdx = attr(decision.table, "decision.attr"),
                                 qualityF = X.gini, nAttrs = NULL,
                                 allowedRandomness = 1/ncol(decision.table),
                                 nOfProbes = ncol(decision.table),
                                 permsWithinINDclasses = FALSE,
                                 inconsistentDecisionTable = NULL)
{
  toRmVec = decisionIdx
  attrIdxVec = (1:ncol(decision.table))[-toRmVec]
  relevanceProbVec = numeric()
  if (!is.null(nAttrs)) {
    if (nAttrs == 0 || nAttrs > ncol(decision.table) - 1) {
      stop("There is something wrong with data (too little attributes?) or the parameter nAttrs has a wrong value.")
    }
    tmpAttrSub = sample(attrIdxVec, min(nAttrs, length(attrIdxVec)))
  }
  else tmpAttrSub = attrIdxVec
  if (allowedRandomness >= 1 || allowedRandomness < 0) {
    stop("Wrong value of the allowedRandomness parameter. It must be within [0,1) interval.")
  }

  INDrelation = list(1:nrow(decision.table))
  INDsizes = nrow(decision.table)
  decisionChaos = compute_chaos(INDrelation, decision.table[[decisionIdx]],
                                attrDescriptions[[decisionIdx]])
  decisionChaos = qualityF(decisionChaos[[1]])
  attrScoresVec = mapply(qualityGain, decision.table[tmpAttrSub], attrDescriptions[tmpAttrSub],
                         MoreArgs = list(decisionVec = decision.table[[decisionIdx]],
                                         uniqueDecisions = attrDescriptions[[decisionIdx]],
                                         INDclasses = INDrelation,
                                         INDclassesSizes = INDsizes,
                                         baseChaos = decisionChaos,
                                         chaosFunction = qualityF),
                         SIMPLIFY = TRUE, USE.NAMES = FALSE)
  tmpBestIdx = which.max(attrScoresVec)
  relevanceProbVec = computeRelevanceProb(INDrelation, INDsizes, decision.table[[tmpAttrSub[tmpBestIdx]]],
                                          uniqueValues = attrDescriptions[[tmpAttrSub[tmpBestIdx]]],
                                          attrScore = attrScoresVec[tmpBestIdx],
                                          decisionVec = decision.table[[decisionIdx]],
                                          uniqueDecisions = attrDescriptions[[decisionIdx]],
                                          baseChaos = decisionChaos,
                                          qualityF = qualityF,
                                          nOfProbes = nOfProbes,
                                          withinINDclasses = permsWithinINDclasses)

  selectedAttrIdxVec = tmpAttrSub[tmpBestIdx]
  INDrelation = compute_indiscernibility(INDrelation,
                                         as.character(decision.table[[selectedAttrIdxVec]]),
                                         attrDescriptions[[selectedAttrIdxVec]])
  attrIdxVec = (1:ncol(decision.table))[-c(selectedAttrIdxVec, toRmVec)]

  if(is.null(inconsistentDecisionTable)) {
    if(sum(duplicated(decision.table)) == sum(duplicated(decision.table[-decisionIdx]))) {
      inconsistentDecisionTable = FALSE
    } else {
      inconsistentDecisionTable = TRUE
    }
  }

  endFlag = FALSE
  iteration = 1
  if(inconsistentDecisionTable) {
    tmpIND = split(1:nrow(decision.table),
                   do.call(paste, decision.table[-decisionIdx]), drop = TRUE)
    totalChaos = compute_chaos(tmpIND,
                               as.character(decision.table[[decisionIdx]]),
                               attrDescriptions[[decisionIdx]])
    totalChaos = sum(sapply(totalChaos, qualityF) * sapply(tmpIND, length) / nrow(decision.table))
    rm(tmpIND)
  } else totalChaos = 0
  totalDependencyInData = decisionChaos - totalChaos - 10^(-16)

  while (!endFlag) {
    contingencyTabs = lapply(INDrelation,
                             function(x,y) table(y[x]),
                             decision.table[[decisionIdx]])
    chaosVec = sapply(contingencyTabs, qualityF)
    if(any(chaosVec == 0)) {
      tmpIdx = which(chaosVec == 0)
      INDrelation = INDrelation[-tmpIdx]
      contingencyTabs = contingencyTabs[-tmpIdx]
      chaosVec = chaosVec[-tmpIdx]
      rm(tmpIdx)
    }
    if(length(INDrelation) > 0) {
      sumsVec = sapply(contingencyTabs, sum)
      tmpChaos = sum(chaosVec*(sumsVec/nrow(decision.table)))
      tmpDependencyInData = decisionChaos - tmpChaos
      rm(contingencyTabs, sumsVec, chaosVec)

      if (!is.null(nAttrs)) {
        tmpAttrSub = sample(attrIdxVec, min(nAttrs, length(attrIdxVec)))
      }
      else tmpAttrSub = attrIdxVec

      INDsizes = sapply(INDrelation, length)
      attrScoresVec = mapply(qualityGain, decision.table[tmpAttrSub], attrDescriptions[tmpAttrSub],
                             MoreArgs = list(decisionVec = decision.table[[decisionIdx]],
                                             uniqueDecisions = attrDescriptions[[decisionIdx]],
                                             INDclasses = INDrelation,
                                             INDclassesSizes = INDsizes,
                                             baseChaos = tmpChaos,
                                             chaosFunction = qualityF),
                             SIMPLIFY = TRUE, USE.NAMES = FALSE)
      tmpBestIdx = which.max(attrScoresVec)
      tmpProbeP = computeRelevanceProb(INDrelation, INDsizes, decision.table[[tmpAttrSub[tmpBestIdx]]],
                                       uniqueValues = attrDescriptions[[tmpAttrSub[tmpBestIdx]]],
                                       attrScore = attrScoresVec[tmpBestIdx],
                                       decisionVec = decision.table[[decisionIdx]],
                                       uniqueDecisions = attrDescriptions[[decisionIdx]],
                                       baseChaos = tmpChaos,
                                       qualityF = qualityF,
                                       nOfProbes = nOfProbes,
                                       withinINDclasses = permsWithinINDclasses)
      if(tmpProbeP > (1 - allowedRandomness)) {
        selectedAttrIdxVec[iteration + 1] = tmpAttrSub[tmpBestIdx]
        INDrelation = compute_indiscernibility(INDrelation,
                                               as.character(decision.table[[tmpAttrSub[tmpBestIdx]]]),
                                               attrDescriptions[[tmpAttrSub[tmpBestIdx]]])
        if(length(INDrelation) == 0) {
          endFlag = TRUE
          approxThereshold = tmpDependencyInData
        }

        attrIdxVec = (1:ncol(decision.table))[-c(selectedAttrIdxVec, toRmVec)]
        iteration = iteration + 1
        relevanceProbVec = c(relevanceProbVec, tmpProbeP)
      } else {
        endFlag = TRUE
        approxThereshold = tmpDependencyInData
      }
      rm(tmpAttrSub, attrScoresVec, tmpBestIdx)
    } else {
      endFlag = TRUE
      approxThereshold = decisionChaos
    }
  }

  if (iteration > 1) {
    endFlag = FALSE
    iteration = iteration - 1
    while (!endFlag) {
      clsContingencyTab = as.matrix(table(do.call(paste,
                                                  decision.table[selectedAttrIdxVec[-iteration]]),
                                          decision.table[[decisionIdx]]))
      tmpChaos = sum(apply(clsContingencyTab, 1, qualityF) *
                       (rowSums(clsContingencyTab)/nrow(decision.table)))
      tmpDependencyInData = decisionChaos - tmpChaos
      if (approxThereshold <= tmpDependencyInData) {
        selectedAttrIdxVec = selectedAttrIdxVec[-iteration]
      }
      iteration = iteration - 1
      if (iteration == 0)
        endFlag = TRUE
    }
  }

  attrsOrder <- order(selectedAttrIdxVec)
  reduct <- selectedAttrIdxVec[attrsOrder]
  relevanceProbVec <- relevanceProbVec[attrsOrder]
  names(reduct) <- colnames(decision.table)[reduct]
  mod <- list(reduct = reduct, type.method = "DAAR.heuristic",
              relevanceProbabilities = relevanceProbVec,
              epsilon = 1 -approxThereshold/totalDependencyInData,
              type.task = "feature selection", model = "RST")
  class.mod <- ObjectFactory(mod, classname = "FeatureSubset")
  return(class.mod)
}

#' This function is a wrapper for computing different types of decision superreducts
#' (i.e. attribute subsets which do not lose any information regarding the decisions
#' but are not require to be irreducable).
#'
#' Currently, there are implemented three methods that can be used with this function:
#' \itemize{
#' \item \code{"greedy.heuristic.superreduct"}: it is a greedy heuristic method which employs several quality measures from RST.
#'            See \code{\link{FS.greedy.heuristic.superreduct.RST}}.
#' \item \code{"quickreduct.frst"}: it is a feature selection function based on the fuzzy QuickReduct algorithm on FRST.
#'            See \code{\link{FS.quickreduct.FRST}}.
#' \item \code{"quickreduct.rst"}: it is a feature selection function based on the RST QuickReduct algorithm.
#'            See \code{\link{FS.quickreduct.RST}}.
#' }
#' These methods can be selected by assigning an appropriate value of the parameter \code{method}.
#' Additionally, \code{\link{SF.applyDecTable}} is provided to generate the new decision table.
#'
#' @title The superreduct computation based on RST and FRST
#' @author Andrzej Janusz
#'
#' @param decision.table an object of a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}.
#' @param method a character representing the type of a method to use for computations. See in Section \code{Details}.
#' @param ... other parameters corresponding to the chosen \code{method}.
#'
#' @return A class \code{"FeatureSubset"}.
#'
#' @seealso \code{\link{FS.greedy.heuristic.superreduct.RST}}, \code{\link{FS.quickreduct.RST}}, \code{\link{FS.quickreduct.FRST}}.
#'
#' @examples
#' ###############################################################
#' ## Example 1: generate reduct and new decision table using RST
#' ###############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## generate single superreduct
#' res.1 <- FS.feature.subset.computation(decision.table,
#'                                        method = "quickreduct.rst")
#'
#' ## generate new decision table according to the reduct
#' new.decTable <- SF.applyDecTable(decision.table, res.1)
#'
#' ###############################################################
#' ## Example 2: generate reduct and new decision table using FRST
#' ###############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$housing7.dt
#'
#' ## generate single superreduct
#' res.2 <- FS.feature.subset.computation(decision.table,
#'                                        method = "quickreduct.frst")
#'
#' ## generate new decision table according to the reduct
#' new.decTable <- SF.applyDecTable(decision.table, res.2)
#' @export
FS.feature.subset.computation <- function(decision.table, method = "greedy.heuristic.superreduct", ...) {

	if (!(method %in% c("greedy.heuristic.superreduct", "quickreduct.rst",
                      "quickreduct.frst"))) {
		stop("Unrecognized attribute reduction method.")
	}

	if (!inherits(decision.table, "DecisionTable")) {
		stop("Provided data should inherit from the \'DecisionTable\' class.")
	}

	nominal.att <- attr(decision.table, "nominal.attrs")
	if (!all(nominal.att) && !(method %in% "quickreduct.frst")) stop("Discretize attribures before computing RST reducts.")

	if (is.null(attr(decision.table, "decision.attr"))) stop("A decision attribute is not indicated.")
	else decIdx = attr(decision.table, "decision.attr")

	if (is.null(attr(decision.table, "desc.attrs"))) stop("No description of attribute values.")
	else desc.attrs = attr(decision.table, "desc.attrs")

	superreduct = switch(method,
                       greedy.heuristic.superreduct = FS.greedy.heuristic.superreduct.RST(decision.table,
                                                                                          attrDescriptions = desc.attrs,
                                                                                          decisionIdx = decIdx, ...),
                       quickreduct.rst = FS.quickreduct.RST(decision.table),
                       quickreduct.frst = FS.quickreduct.FRST(decision.table) )

	return(superreduct)
}

#' This is a function for implementing the QuickReduct algorithm for feature selection based
#' on RST proposed by (Shen and Chouchoulas, 2000). The algorithm produces only one feature subset that could be a superreduct.
#'
#' This algorithm considers the dependency degree (see \code{\link{A.Introduction-RoughSets}})
#' of the addition of each attribute to the current reduct candidate. Then the best candidate will be chosen.
#' This process continues until the dependency of the subset equals to the dependency of the full dataset.
#'
#' Additionally, in \code{control} parameter, we provide one component which is
#' \code{randomize}. It has a boolean value: \code{TRUE} or \code{FALSE} that means we want to perform
#' quickreduct by evaluating attributes randomly or all attributes in decision table.
#'
#' It should be noted that this function does not give the new decision table directly.
#' The other additional function called \code{\link{SF.applyDecTable}} is used to produce new decision table based on
#' information about the reduct from this function.
#'
#' @title QuickReduct algorithm based on RST
#' @author Lala Septem Riza
#'
#' @param decision.table an object of a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}.
#' @param control other parameters. It contains the following component:
#'        \itemize{
#'        \item \code{randomize}: it has a boolean value. For the detailed description, see in Section \code{Details}.
#'              The default value is \code{FALSE}.
#'        }
#'
#' @return A class \code{"FeatureSubset"} that contains the following components:
#' \itemize{
#' \item \code{reduct}: a list representing single reduct. In this case, it could be super reduct or just subset of feature.
#' \item \code{type.method}: a string representing a type of method which is \code{"quickreduct"}.
#' \item \code{type.task}: a string showing type of task which is \code{"feature selection"}.
#' \item \code{model}: a string representing a type of model. In this case, it is \code{"RST"} which means rough set theory.
#' }
#'
#' @seealso \code{\link{FS.quickreduct.FRST}}
#'
#' @references
#' Q. Shen and A. Chouchoulas, "A Modular Approach to Generating Fuzzy Rules with Reduced Attributes for the Monitoring of Complex Systems",
#' Engineering Applications of Artificial Intelligence, vol. 13, p. 263 - 278 (2000).
#'
#' @examples
#' ###################################################
#' ## Example 1: Evaluate reduct and generate
#' ##            new decision table
#' ###################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## evaluate single reduct
#' res.1 <- FS.quickreduct.RST(decision.table)
#'
#' ## generate new decision table according to the reduct
#' new.decTable <- SF.applyDecTable(decision.table, res.1)
#'
#' @export
FS.quickreduct.RST <- function(decision.table, control = list()){

	names.attrs <-  names(attr(decision.table, "desc.attrs"))
	super.reduct <- quickreduct.alg(decision.table, type.method = "quickreduct.RST", control = control)
    names(super.reduct) = names.attrs[super.reduct]
	mod <- list(reduct = super.reduct, type.method = "quickreduct",
	            type.task = "feature selection", model = "RST")

	class.mod <- ObjectFactory(mod, classname = "FeatureSubset")
	return(class.mod)
}

#' It is used to get a feature subset (superreduct) based on the greedy heuristic algorithm
#' employing some quality measurements. Regarding the quality measurements, the detailed description can be seen in \code{\link{FS.greedy.heuristic.reduct.RST}}.
#'
#' @title The greedy heuristic method for determining superreduct based on RST
#' @author Andrzej Janusz
#'
#' @param decision.table an object of a \code{"DecisionTable"} class representing a decision table.
#'        See \code{\link{SF.asDecisionTable}}.
#' @param attrDescriptions a list containing possible values of attributes (columns)
#'        in \code{decision.table}. It usually corresponds to \code{attr(decision.table, "desc.attrs")}.
#' @param decisionIdx a integer value representing an index of decision attribute.
#' @param qualityF a function for calculating a quality of an attribute subset.
#'        See \code{\link{FS.greedy.heuristic.reduct.RST}}.
#' @param nAttrs an integer between 1 and the number of conditional attributes. It indicates
#'        the attribute sample size for the Monte Carlo selection of candidating attributes.
#'        If set to \code{NULL} (default) all attributes are used and the algorithm changes
#'        to a standard greedy method for computation of decision reducts.
#' @param inconsistentDecisionTable logical indicating whether the decision table is suspected
#'        to be inconsistent or \code{NULL} (the default) which indicated that a test should
#'        be made to determine the data consistency.
#'
#' @return A class \code{"FeatureSubset"} that contains the following components:
#' \itemize{
#' \item \code{reduct}: a list representing a single reduct. In this case, it could be a superreduct or just a subset of features.
#' \item \code{type.method}: a string representing the type of method which is \code{"greedy.heuristic.superreduct"}.
#' \item \code{type.task}: a string showing the type of task which is \code{"feature selection"}.
#' \item \code{model}: a string representing the type of model. In this case, it is \code{"RST"} which means rough set theory.
#' }
#'
#' @seealso \code{\link{FS.quickreduct.RST}} and \code{\link{FS.feature.subset.computation}}.
#'
#' @references
#' A. Janusz and S. Stawicki, "Applications of Approximate Reducts to the Feature Selection Problem",
#' Proceedings of International Conference on Rough Sets and Knowledge Technology ({RSKT}), vol. 6954, p. 45 - 50 (2011).
#'
#' D. Ślęzak, "Approximate Entropy Reducts", Fundamenta Informaticae, vol. 53, no. 3 - 4, p. 365 - 390 (2002).
#'
#' J. Wroblewski, "Ensembles of Classifiers Based on Approximate Reducts", Fundamenta Informaticae, vol. 47, no. 3 - 4, p. 351 - 360 (2001).
#'
#' @examples
#' ###################################################
#' ## Example 1: Evaluate reduct and generate
#' ##            new decision table
#' ###################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## evaluate single reduct
#' res.1 <- FS.greedy.heuristic.superreduct.RST(decision.table, qualityF = X.nOfConflicts)
#' print(res.1)
#'
#' ## generate new decision table according to the reduct
#' new.decTable <- SF.applyDecTable(decision.table, res.1)
#' @export
FS.greedy.heuristic.superreduct.RST <- function(decision.table,
                                                attrDescriptions = attr(decision.table, "desc.attrs"),
                                                decisionIdx = attr(decision.table, "decision.attr"),
                                                qualityF = X.gini, nAttrs = NULL,
                                                inconsistentDecisionTable = NULL)  {
  toRmVec = decisionIdx
  attrIdxVec = (1:ncol(decision.table))[-toRmVec]

  if (!is.null(nAttrs)) {
    if(nAttrs == 0 || nAttrs > ncol(decision.table) - 1) {
      stop("There is something wrong with data (too little attributes?) or the parameter nAttrs has a wrong value.")
    }
    tmpAttrSub = sample(attrIdxVec, min(nAttrs, length(attrIdxVec)))
  }  else tmpAttrSub = attrIdxVec

  INDrelation = list(1:nrow(decision.table))
  INDsizes = nrow(decision.table)
  decisionChaos = compute_chaos(INDrelation, decision.table[[decisionIdx]],
                                attrDescriptions[[decisionIdx]])
  decisionChaos = qualityF(decisionChaos[[1]])
  attrScoresVec = mapply(qualityGain, decision.table[tmpAttrSub], attrDescriptions[tmpAttrSub],
                         MoreArgs = list(decisionVec = decision.table[[decisionIdx]],
                                         uniqueDecisions = attrDescriptions[[decisionIdx]],
                                         INDclasses = INDrelation,
                                         INDclassesSizes = INDsizes,
                                         baseChaos = decisionChaos,
                                         chaosFunction = qualityF),
                         SIMPLIFY = TRUE, USE.NAMES = FALSE)
  tmpBestIdx = which.max(attrScoresVec)

  selectedAttrIdxVec  = tmpAttrSub[tmpBestIdx]
  INDrelation = compute_indiscernibility(INDrelation,
                                         as.character(decision.table[[selectedAttrIdxVec]]),
                                         attrDescriptions[[selectedAttrIdxVec]])
  attrIdxVec = (1:ncol(decision.table))[-c(selectedAttrIdxVec, toRmVec)]

  if(is.null(inconsistentDecisionTable)) {
    if(sum(duplicated(decision.table)) == sum(duplicated(decision.table[-decisionIdx]))) {
      inconsistentDecisionTable = FALSE
    } else {
      inconsistentDecisionTable = TRUE
    }
  }

  endFlag = F
  iteration = 1
  if(inconsistentDecisionTable) {
    tmpIND = split(1:nrow(decision.table),
                   do.call(paste, decision.table[-decisionIdx]), drop = TRUE)
    totalChaos = compute_chaos(tmpIND,
                               as.character(decision.table[[decisionIdx]]),
                               attrDescriptions[[decisionIdx]])
    totalChaos = sum(sapply(totalChaos, qualityF) * sapply(tmpIND, length) / nrow(decision.table))
    rm(tmpIND)
  } else totalChaos = 0
  totalDependencyInData = decisionChaos - totalChaos - 10^(-16)

  while (!endFlag) {
    contingencyTabs = lapply(INDrelation,
                             function(x,y) table(y[x]),
                             decision.table[[decisionIdx]])
    chaosVec = sapply(contingencyTabs, qualityF)
    if(any(chaosVec == 0)) {
      tmpIdx = which(chaosVec == 0)
      INDrelation = INDrelation[-tmpIdx]
      contingencyTabs = contingencyTabs[-tmpIdx]
      chaosVec = chaosVec[-tmpIdx]
      rm(tmpIdx)
    }
    sumsVec = sapply(contingencyTabs, sum)
    if(length(INDrelation) > 0) tmpChaos = sum(chaosVec*(sumsVec/nrow(decision.table)))
    else tmpChaos = 0
    tmpDependencyInData = decisionChaos - tmpChaos
    rm(contingencyTabs, sumsVec, chaosVec)
    if (totalDependencyInData <= tmpDependencyInData) endFlag = TRUE
    else {
      if (!is.null(nAttrs)) {
        tmpAttrSub = sample(attrIdxVec, min(nAttrs, length(attrIdxVec)))
      }  else tmpAttrSub = attrIdxVec
      INDsizes = sapply(INDrelation, length)
      attrScoresVec = mapply(qualityGain, decision.table[tmpAttrSub], attrDescriptions[tmpAttrSub],
                             MoreArgs = list(decisionVec = decision.table[[decisionIdx]],
                                             uniqueDecisions = attrDescriptions[[decisionIdx]],
                                             INDclasses = INDrelation,
                                             INDclassesSizes = INDsizes,
                                             baseChaos = tmpChaos,
                                             chaosFunction = qualityF),
                             SIMPLIFY = TRUE, USE.NAMES = FALSE)
      tmpBestIdx = which.max(attrScoresVec)
      selectedAttrIdxVec[iteration + 1] = tmpAttrSub[tmpBestIdx]
      INDrelation = compute_indiscernibility(INDrelation,
                                             as.character(decision.table[[tmpAttrSub[tmpBestIdx]]]),
                                             attrDescriptions[[tmpAttrSub[tmpBestIdx]]])
      if(length(INDrelation) == 0) endFlag = TRUE

      attrIdxVec = (1:ncol(decision.table))[-c(selectedAttrIdxVec, toRmVec)]
      iteration = iteration + 1
    }
  }

  super.reduct <- selectedAttrIdxVec[order(selectedAttrIdxVec)]
  names(super.reduct) = colnames(decision.table)[super.reduct]
  mod <- list(reduct = super.reduct, type.method = "greedy.heuristic.superreduct",
              type.task = "feature selection", model = "RST")

  class.mod <- ObjectFactory(mod, classname = "FeatureSubset")
  return(class.mod)
}


#' It is a function implementing the fuzzy QuickReduct algorithm
#' for feature selection based on FRST.
#' The fuzzy QuickReduct is a modification of QuickReduct based on RST (see \code{\link{FS.quickreduct.RST}}).
#'
#' In this function, we provide an algorithm proposed by
#'(Jensen and Shen, 2002) which is fuzzy QuickReduct. Then, the algorithm has been modified by (Bhatt and Gopal, 2005) to improve stopping criteria.
#' This function is aimed to implement both algorithms. These algorithms can be executed by assigning the parameter \code{type.QR}
#' with \code{"fuzzy.QR"} and \code{"modified.QR"} for fuzzy quickreduct and modified fuzzy quickreduct
#' algorithms, respectively. Additionally, in the \code{control} parameter, we provide one component which is
#' \code{randomize} having boolean values: \code{TRUE} or \code{FALSE}. \code{randomize = TRUE} means that
#' we evaluate some (or not all) attributes randomly along iteration. It will be useful if we have a large number of attributes
#' in a decision table.
#'
#' In this function, we have considered many approaches of the lower and upper approximations.
#' The following list shows considered methods and their descriptions. Additionally, those approaches can be executed by
#' assigning the following value to the parameter \code{type.method}.
#' \itemize{
#'    \item \code{"fuzzy.dependency"}: It is based on the degree of dependency using the implication/t-norm model approximation (Jensen and Shen, 2009).
#'                The detailed concepts about this approximation have been explained in \code{\link{B.Introduction-FuzzyRoughSets}}
#'                and
#'
#'                \code{\link{BC.LU.approximation.FRST}}.
#'    \item \code{"fuzzy.boundary.reg"}: It is based on the fuzzy boundary region proposed by (Jensen and Shen, 2009).
#'                This algorithm introduced the usage of the total uncertainty degree \eqn{\lambda_B(Q)}
#'                for all concepts of feature subset \eqn{B} and decision attribute \eqn{Q}.
#'                The total uncertainty degree is used as a parameter to select appropriate features.
#'   \item \code{"vqrs"}: It is based on vaquely quantified rough set (VQRS)
#'                proposed by (Cornelis and Jensen, 2008). See also \code{\link{BC.LU.approximation.FRST}}.
#'   \item \code{"owa"}: Based on ordered weighted average (OWA) based fuzzy rough set, (Cornelis et al, 2010) proposed
#'                the degree of dependency as a parameter employed in the algorithm to select appropriate features. The explanation
#'                about lower and upper approximations based on OWA can be found in \code{\link{BC.LU.approximation.FRST}}.
#'   \item \code{"rfrs"}: It is based on degree of dependency that is obtained by performing
#'                the robust fuzzy rough sets proposed by (Hu et al, 2012).
#'                The detailed concepts about this approximation have been explained in \code{\link{BC.LU.approximation.FRST}}.
#   \item \code{"sfrs"}: It is based on degree of dependency that is obtained by performing
#                soft fuzzy rough sets (SFRS) proposed by (Hu et al, 2010).
#                The detailed concepts about this approximation have been explained in \code{\link{BC.LU.approximation.FRST}}.
#                Additionally, it should be noted that this method is good in selecting features on
#                dataset containing continuous conditional features only.
#'   \item \code{"min.positive.reg"}: Based on measure introduced in (Cornelis et al, 2010) which considers the most problematic element in
#'              the positive region, defined using the implicator/t-norm model.
#'   \item \code{"fvprs"}: It is based on degree of dependency proposed by (Zhao et al, 2009).
#'                The degree is obtained by using fuzzy lower approximation based on
#'                fuzzy variable precision rough set model.
#'   \item \code{"fuzzy.discernibility"}: This approach attempts to combine the the decision-relative discernibility matrix
#'               and the fuzzy QuickReduct algorithm. (Jensen and Shen, 2009) introduced a measurement which is the degree of satisfaction to select the attributes.
#'   \item \code{"beta.pfrs"}: Based on \eqn{\beta}-precision fuzzy rough sets (\eqn{\beta}-PFRS) proposed by (Salido and Murakami, 2003),
#'                the degree of dependency as a parameter employed in the algorithm to select appropriate features. The explanation
#'                about lower and upper approximations based on \eqn{\beta}-PFRS can be found in \code{\link{BC.LU.approximation.FRST}}.
#' }
#'
#' It should be noted that the parameter \code{type.method} is related to parameter \code{control}.
#' In other words, we only set the components in the \code{control} parameter that related to the chosen type of method.
#' The following is a list showing the components of \code{control} needed by each type of methods.
#' \itemize{
#' \item \code{type.method = "fuzzy.dependency"}:
#'
#' \code{control <- list(t.implicator, type.relation, type.aggregation)}
#'
#' \item \code{type.method = "fuzzy.boundary.reg"}:
#'
#' \code{control <- list(t.implicator, type.relation, type.aggregation)}
#'
#' \item \code{type.method = "vqrs"}:
#'
#' \code{control <- list(alpha, q.some, q.most, type.aggregation)}
#'
#' \item \code{type.method = "owa"}:
#'
#' \code{control <- list(t.implicator, type.relation, m.owa, type.aggregation)}
#'
# \item \code{type.method = "sfrs"}:
#
# \code{control <- list(penalty.fact, type.aggregation)}
#
#' \item \code{type.method = "rfrs"}:
#'
#' \code{control <- list(t.implicator, type.relation, type.rfrs,}
#'
#'                 \code{k.rfrs, type.aggregation)}
#'
#' \item \code{type.method = "min.positive.reg"}:
#'
#' \code{control <- list(alpha, t.implicator, type.relation, type.aggregation)}
#'
#' \item \code{type.method = "fuzzy.discernibility"}:
#'
#' \code{control <- list(alpha, t.implicator, type.relation, type.aggregation)}
#'
#' \item \code{type.method = "fvprs"}:
#'
#' \code{control <- list(alpha.precision, t.implicator, type.relation, type.aggregation)}
#'
#' \item \code{type.method = "beta.pfrs"}:
#'
#' \code{control <- list(t.implicator, type.relation, beta.quasi, type.aggregation)}
#' }
#' The descriptions of each component can be seen in the documentation of the \code{control} parameter.
#'
#' It should be noted that this function does not give the new decision table directly.
#' An additional function called \code{\link{SF.applyDecTable}} is used to produce new decision table based on
#' information about the reduct from this function. See Section \code{Examples}.
#'
#' @title The fuzzy QuickReduct algorithm based on FRST
#' @author Lala Septem Riza
#'
#' @param decision.table an object of a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}.
#' @param type.method a string representing the type of methods.
#'         The complete description can be found in Section \code{Details}.
#' @param type.QR a string expressing the type of QuickReduct algorithm which is one of the two following algorithms:
#'         \itemize{
#'         		\item \code{"fuzzy.QR"}: it is the original fuzzy rough QuickReduct algorithm based on (Jensen and Shen, 2002).
#'         		\item \code{"modified.QR"}: it is the modified QuickReduct algorithm based on (Bhatt and Gopal, 2005).
#'         }
#' @param control a list of other parameters as follows.
#'         \itemize{
#'          \item \code{type.aggregation}: a type of aggregation operator. See \code{\link{BC.IND.relation.FRST}}.
#'          \item \code{t.implicator}: a type of implicator function. See \code{\link{BC.LU.approximation.FRST}}.
#'                The default value is \code{"lukasiewicz"}.
#'          \item \code{type.relation}: a type of indiscernibility relation. See \code{\link{BC.IND.relation.FRST}}.
#'                 The default value is \code{type.relation = c("tolerance", "eq.3")}.
#'          \item \code{alpha}: a real number between 0 and 1 expressing a threshold value or stopping criterion.
#'                 The following methods use the parameter: \code{"vqrs"},
#'
#'                 \code{"min.positive.reg"}, and \code{"fuzzy.discernibility"}.
#'                 The default value is 0.95.
#'          \item \code{alpha.precision}: a real number between 0 and 1 expressing variable precision (\eqn{\alpha}) for \code{"fvprs"}.
#'                 See \code{\link{BC.LU.approximation.FRST}}. The default value is 0.05.
#'          \item \code{q.some}: a pair of numeric values for the alpha and beta parameter of VQRS for the quantifier \code{some}.
#'                 The default value is \code{q.some = c(0.1, 0.6)}.
#'
#'                 See \code{\link{BC.LU.approximation.FRST}}.
#'          \item \code{q.most}: a pair of numeric values for the alpha and beta parameter of VQRS for the quantifier \code{most}.
#'                 The default value is \code{q.most = c(0.2, 1)}.
#'
#'                 See \code{\link{BC.LU.approximation.FRST}}.
#'          \item \code{m.owa}: a numeric value to define the parameter in OWA. The default value is the mean number of objects.
#          \item \code{penalty.fact}: a real number representing penalty factor on soft fuzzy rough sets.
#                 The default value is 0.8. See \code{\link{BC.LU.approximation.FRST}}.
#'          \item \code{type.rfrs}: a type of robust fuzzy rough sets.
#'
#'                 The default is \code{type.rfrs = "k.trimmed.min")}.
#'
#'                 See \code{\link{BC.LU.approximation.FRST}}.
#'          \item \code{k.rfrs}: a value between 0 and length of data representing index of considered data.
#'                 The default is \code{k.rfrs = round(0.5*nrow(decision.table))}.
#'                See \code{\link{BC.LU.approximation.FRST}}.
#'          \item \code{beta.quasi}: a number between 0 and 1 representing \eqn{\beta}-precision t-norms and t-conorms.
#'                    The default value is 0.05.
#'          \item \code{randomize}: a boolean value to define whether selecting attributes randomly or not. For more detail,
#'                see in Section \code{Details}. The default value is \code{FALSE}.
#'          }
#'         It should be noted that instead of supplying all the above parameters, we only set
#'         those parameters needed by the considered method. See in Section \code{Details}.
#'         Also, we provide some examples to illustrate how the parameters are used.
#'
#' @seealso \code{\link{FS.quickreduct.RST}} and \code{\link{FS.feature.subset.computation}}.
#' @return A class \code{"FeatureSubset"} that contains the following components:
#' \itemize{
#' \item \code{reduct}: a list representing a single reduct. In this case, it could be a superreduct or just a subset of feature.
#' \item \code{type.method}: a string representing the type of method.
#' \item \code{type.task}: a string showing the type of task which is \code{"feature selection"}.
#' \item \code{model}: a string representing the type of model. In this case, it is \code{"FRST"} which means fuzzy rough set theory.
#' }
#' @references
#' C. Cornelis, G. Hurtado Martin, R. Jensen, and D. Slezak,
#' "Feature Selection with Fuzzy Decision Reducts", Information Sciences, vol. 180, no. 2, p. 209 - 224 (2010).
#'
#' C. Cornelis, N. Verbiest, and R. Jensen, "Ordered Weighted Average Based Fuzzy Rough Sets",
#' Proceedings of the 5th International Conference on Rough Sets and Knowledge Technology (RSKT 2010),
#' p. 78 - 85 (2010).
#'
#' C. Cornelis and R. Jensen, "A Noise-tolerant Approach to Fuzzy-rough Feature Selection",
#' Proceedings of the 2008 IEEE International Conference on Fuzzy Systems (FUZZ-IEEE 2008),
#' p. 1598 - 1605 (2008).
#'
#' J. M. F. Salido and S. Murakami, "Rough Set Analysis of a General Type of Fuzzy Data
#' Using Transitive Aggregations of Fuzzy Similarity Relations",
#' Fuzzy Sets Syst., vol. 139, p. 635 - 660 (2003).
#'
# Q. Hu, S. An, and D. Yu, "Soft Fuzzy Rough Sets for Robust Feature Evaluation and Selection",
# Information Sciences, vol. 180, p. 4384 - 4400 (2010).
#
#' Q. Hu, L. Zhang, S. An, D. Zhang, and D. Yu, "On Robust Fuzzy Rough Set Models",
#' IEEE Trans. on Fuzzy Systems, vol. 20, no. 4, p. 636 - 651 (2012).
#'
#' R. B. Bhatt and M. Gopal, "On Fuzzy-rough Sets Approach to Feature Selection",
#' Pattern Recognition Letters, vol. 26, no. 7, p. 965 - 975 (2005).
#'
#' R. Jensen and Q. Shen, "Fuzzy-rough Sets for Descriptive Dimensionality Reduction",
#' In: Proceedings of IEEE International Conference on Fuzzy System, FUZZ-IEEE, p. 29 - 34 (2002).
#'
#' R. Jensen and Q. Shen, "New Approaches to Fuzzy-rough Feature Selection",
#' IEEE Transactions on Fuzzy Systems, vol. 17, no. 4, p. 824 - 838 (2009).
#'
#' S. Y. Zhao, E. C. C. Tsang, and D. G. Chen,
#' "The Model of Fuzzy Variable Precision Rough Sets",
#' IEEE Trans. Fuzzy Systems, vol. 17, no. 2,
#' p. 451 - 467 (2009).
#'
#' @examples
#' ##########################################################
#' ## Example 1: Dataset containing nominal values on all attributes
#' ##########################################################
#'
#' data(RoughSetData)
#' decision.table <- RoughSetData$housing7.dt
#'
#' ########## using fuzzy lower approximation ##############
#' control <- list(t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"),
#'                type.aggregation = c("t.tnorm", "lukasiewicz"))
#' reduct.1 <- FS.quickreduct.FRST(decision.table, type.method = "fuzzy.dependency",
#'                             type.QR = "fuzzy.QR", control = control)
#'
#' ########## using fuzzy boundary region ##############
#' \dontrun{control <- list(t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"),
#'                 type.aggregation = c("t.tnorm", "lukasiewicz"))
#' reduct.2 <- FS.quickreduct.FRST(decision.table, type.method = "fuzzy.boundary.reg",
#'                             type.QR = "fuzzy.QR", control = control)
#'
#' ########## using vaguely quantified rough sets (VQRS) #########
#' control <- list(alpha = 0.9, q.some = c(0.1, 0.6), q.most = c(0.2, 1),
#'                 type.aggregation = c("t.tnorm", "lukasiewicz"))
#' reduct.3 <- FS.quickreduct.FRST(decision.table, type.method = "vqrs",
#'                             type.QR = "fuzzy.QR", control = control)
#'
#' ########## ordered weighted average (OWA) #########
#' control <- list(t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"),
#'                 m.owa = 3, type.aggregation = c("t.tnorm","lukasiewicz"))
#' reduct.4 <- FS.quickreduct.FRST(decision.table, type.method = "owa",
#'                             type.QR = "fuzzy.QR", control = control)
#'
#' ########## robust fuzzy rough sets (RFRS) #########
#' control <- list(t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"),
#'                type.rfrs = "k.trimmed.min", type.aggregation = c("t.tnorm", "lukasiewicz"),
#'                k.rfrs = 0)
#' reduct.5 <- FS.quickreduct.FRST(decision.table, type.method = "rfrs",
#'                             type.QR = "fuzzy.QR", control = control)
#'
#' ########## using min positive region (delta) ###########
#' control <- list(alpha = 1, t.implicator = "lukasiewicz",
#'                 type.relation = c("tolerance", "eq.1"), type.aggregation =
#'                                 c("t.tnorm", "lukasiewicz"))
#' reduct.6 <- FS.quickreduct.FRST(decision.table, type.method = "min.positive.reg",
#'                             type.QR = "fuzzy.QR", control = control)
#'
#' ########## using FVPRS approximation ##############
#' control <- list(alpha.precision = 0.05, t.implicator = "lukasiewicz",
#'                type.aggregation = c("t.tnorm", "lukasiewicz"),
#'                type.relation = c("tolerance", "eq.1"))
#' reduct.7 <- FS.quickreduct.FRST(decision.table, type.method = "fvprs",
#'                             type.QR = "fuzzy.QR", control = control)
#'
# ########## using SFRS approximation ##############
# control <- list(penalty.fact = 1, type.aggregation = c("t.tnorm", "lukasiewicz"))
# reduct.8 <- FS.quickreduct.FRST(decision.table, type.method = "sfrs",
#                             type.QR = "fuzzy.QR", control = control)
#
#' ########## using beta.PFRS approximation ##############
#' control <- list(t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"),
#'                 beta.quasi = 0.05, type.aggregation = c("t.tnorm", "lukasiewicz"))
#' reduct.8 <- FS.quickreduct.FRST(decision.table, type.method = "beta.pfrs",
#'                             type.QR = "fuzzy.QR", control = control)
#'
#' ########## using fuzzy discernibility matrix ##############
#' control <- list(alpha = 1, type.relation = c("tolerance", "eq.1"),
#'                type.aggregation = c("t.tnorm", "lukasiewicz"),
#'                 t.implicator = "lukasiewicz")
#' reduct.9 <- FS.quickreduct.FRST(decision.table, type.method = "fuzzy.discernibility",
#'                             type.QR = "fuzzy.QR", control = control)}
#'
#' ##########################################################
#' ## Example 2: Dataset containing nominal and continuous values
#' ## In this case, we only provide one method but others work in
#' ## the same way.
#' ## In this example, we will show how to get the
#' ## new decision table as well
#' ##########################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ########## using fuzzy lower approximation ##############
#' control <- list(type.aggregation = c("t.tnorm", "lukasiewicz"),
#'                t.implicator = "lukasiewicz", type.relation = c("tolerance", "eq.1"))
#' reduct.1 <- FS.quickreduct.FRST(decision.table, type.method = "fuzzy.dependency",
#'                             type.QR = "fuzzy.QR", control = control)
#'
#' ## get new decision table based on reduct
#' new.decTable <- SF.applyDecTable(decision.table, reduct.1)
#'
#' @export
FS.quickreduct.FRST <- function(decision.table, type.method = "fuzzy.dependency", type.QR = "fuzzy.QR", control = list()) {

	## execute quickreduct algorithm
	super.reduct <- quickreduct.alg(decision.table, type.method, type.QR, control)

	## construct FeatureSubset class
	names.attrs <-  names(attr(decision.table, "desc.attrs"))
    names(super.reduct) = names.attrs[super.reduct]
	mod <- list(reduct = super.reduct, type.method = type.method, type.task = "feature selection", model = "FRST")
	class.mod <- ObjectFactory(mod, classname = "FeatureSubset")

	return(class.mod)
}

#' This is a function implementing the near-optimal reduction algorithm by employing
#' fuzzy variable precision rough sets (FVPRS) for feature selection
#' based on FRST proposed by (Zhao et al, 2009).
#'
#' The near-optimal algorithm is an algorithm to find one reduct only rather than all reducts. It modifies the \eqn{\alpha}-reduction based on
#' discernibility matrix by using a heuristic algorithm. To get basic knowledge about discernibility matrix,
#' users can refers to the \code{"alpha.red"} discernibility type in \code{\link{BC.discernibility.mat.FRST}} .
#'
#' It should be noted that this function does not give the new decision table directly.
#' The other additional function called \code{\link{SF.applyDecTable}} is used to produce the new decision table based on
#' information about the reduct from this function.
#'
#' @title The near-optimal reduction algorithm based on fuzzy rough set theory
#' @author Lala Septem Riza
#'
#' @param decision.table  an object of a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}.
#'                        In this case, the decision attribute must be nominal.
#' @param alpha.precision a numeric value representing variable precision of FVPRS.
#'
#'        See \code{\link{BC.LU.approximation.FRST}}.
#' @seealso \code{\link{BC.discernibility.mat.FRST}}
#' @return A class \code{"FeatureSubset"} that contains the following components:
#' \itemize{
#' \item \code{reduct}: a list representing a single reduct. In this case, it could be a superreduct or just a subset of features.
#' \item \code{type.method}: a string representing the type of method which is \code{"near.optimal.fvprs"}.
#' \item \code{type.task}: a string showing the type of task which is \code{"feature selection"}.
#' \item \code{model}: a string representing the type of model. In this case, it is \code{"FRST"} which means fuzzy rough set theory.
#' }
#' @references
#' S. Zhao, E. C. C. Tsang, and D. Chen, "The Model of Fuzzy Variable Precision Rough Sets",
#' IEEE Trans. on Fuzzy Systems, vol. 17, no. 2, p. 451 - 467 (2009).
#'
#' @examples
#' #########################################################
#' ## Example 1: Hiring dataset containing 8 objects with 5 attributes
#' #########################################################
#' \dontrun{data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## get reduct as FeatureSubset class
#' reduct.1 <- FS.nearOpt.fvprs.FRST(decision.table)
#'
#' ## get new decision table according to the reduct
#' new.decTable <- SF.applyDecTable(decision.table, reduct.1)}
#'
#' #########################################################
#' ## Example 2: Pima dataset containing 7 objects with 9 attributes
#' #########################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$pima7.dt
#'
#' ## get reduct
#' reduct.2 <- FS.nearOpt.fvprs.FRST(decision.table)
#'
#' ## get new decision table according to the reduct
#' new.decTable <- SF.applyDecTable(decision.table, reduct.2)
#' @export
FS.nearOpt.fvprs.FRST <- function(decision.table, alpha.precision = 0.05) {

	if (is.null(attr(decision.table, "decision.attr"))){
		stop("A decision attribute is not indicated.")
	}

	if (attr(decision.table, "nominal.attrs")[attr(decision.table, "decision.attr")] == FALSE){
		stop("The decision attribute must be nominal values")
	}

	## get the data
	objects <- decision.table
	nominal.att <- attr(decision.table, "nominal.attrs")
	desc.attrs <- attr(decision.table, "desc.attrs")
	num.att <- ncol(objects)
	num.object <- nrow(objects)
	names.attr <- t(colnames(objects))
	type.discernibility = "alpha.red"
	num.red = "near.opt"
	type.relation = c("tolerance", "eq.1")
	t.implicator = "lukasiewicz"
	type.LU = "implicator.tnorm"
	type.aggregation <- c("t.tnorm", "lukasiewicz")

	## build decision-relative discernibility matrix
	res.temp <- build.discMatrix.FRST(decision.table, type.discernibility, num.red, alpha.precision, type.relation,
		                               t.implicator, type.LU, type.aggregation = type.aggregation)

	disc.mat = res.temp$disc.mat
	disc.list = res.temp$disc.list

	## perform near-optimal reduction using FVPRS
	## delete redundant of discernibility matrix
	disc.list <- unique(disc.list)

	## delete blank characters
	disc.list.temp <- c()
	j <- 1
	for (i in 1 : length(disc.list)){
		if (!(is.character(disc.list[[i]]) & length(disc.list[[i]]) == 0) && !is.na(disc.list[[i]])){
			disc.list.temp[j] <- list(disc.list[[i]])
			j <- j + 1
		}
	}
	disc.list <- disc.list.temp

	if (length(disc.list) > 1){
		## get core
		core.red <- disc.list[which(lapply(disc.list, length) == 1)]

		if (length(core.red) != 0){
			## delete element == core (for each core)
			disc.list <- disc.list[which(lapply(disc.list, length) > 1)]

			## delete element containing core
			for (i in 1 : length(core.red)){
				disc.list <- disc.list[which(lapply(disc.list, function(x){core.red[i] %in% x}) == FALSE)]
			}
		}		else {
			core.red <- c()
		}
		ii <- 1
		while (ii <= length(disc.list)){

			## get the max of freq.
			new.red <- names(which.max(table(unlist(disc.list))))

			## add new.red into reducts
			core.red <- c(core.red, new.red)

			## delete element containing new.red
			disc.list <- disc.list[which(lapply(disc.list, function(x){new.red %in% x}) == FALSE)]

			ii <- ii + 1
		}
		## delete redundant elements
		reduct <- unlist(unique(core.red))
        reduct <- which(as.character(names.attr) %in% reduct)

		## construct FeatureSubset class
		mod <- list(reduct = reduct, type.method = "near.optimal.fvprs", type.task = "feature selection", model = "FRST")
		class.mod <- ObjectFactory(mod, classname = "FeatureSubset")

		return(class.mod)

	}	else if (length(disc.list) == 1){

		reduct <- unlist(sort(core.red))
		reduct = which(as.character(names.attr) %in% reduct)
        colnames(reduct) = names.attr[reduct]

		## construct FeatureSubset class
		mod <- list(reduct = reduct, type.method = "near.optimal.fvprs", type.task = "feature selection", model = "FRST")
		class.mod <- ObjectFactory(mod, classname = "FeatureSubset")

		return(class.mod)
	}
}

#' A wrapper function used for generating all decision reducts of a decision system. The reducts
#' are obtained from a discernibility matrix which can be computed using methods based on RST
#' and FRST. Therefore, it should be noted that before calling the function, we need to
#' compute a discernibility matrix using \code{\link{BC.discernibility.mat.RST}} or
#' \code{\link{BC.discernibility.mat.FRST}}.
#'
#' @title A function for computing all decision reducts of a decision system
#' @author Andrzej Janusz
#'
#' @param discernibilityMatrix an \code{"DiscernibilityMatrix"} object representing
#' a discernibility matrix of a decision system.
#'
#' @return An object of a class \code{"ReductSet"}.
#'
#' @seealso \code{\link{BC.discernibility.mat.RST}}, \code{\link{BC.discernibility.mat.FRST}}.
#'
#' @examples
#' ########################################################
#' ## Example 1: Generate all reducts and
#' ##            a new decision table using RST
#' ########################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## build the decision-relation discernibility matrix
#' res.2 <- BC.discernibility.mat.RST(decision.table, range.object = NULL)
#'
#' ## generate all reducts
#' reduct <- FS.all.reducts.computation(res.2)
#'
#' ## generate new decision table
#' new.decTable <- SF.applyDecTable(decision.table, reduct, control = list(indx.reduct = 1))
#'
#' ##############################################################
#' ## Example 2: Generate all reducts and
#' ##            a new decision table using FRST
#' ##############################################################
#' \dontrun{data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## build the decision-relation discernibility matrix
#' control.1 <- list(type.relation = c("crisp"),
#'                 type.aggregation = c("crisp"),
#'                 t.implicator = "lukasiewicz", type.LU = "implicator.tnorm")
#' res.1 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "standard.red",
#'                                     control = control.1)
#'
#' ## generate single reduct
#' reduct <- FS.all.reducts.computation(res.1)
#'
#' ## generate new decision table
#' new.decTable <- SF.applyDecTable(decision.table, reduct, control = list(indx.reduct = 1))}
#' @export
FS.all.reducts.computation <- function(discernibilityMatrix) {

  if (!inherits(discernibilityMatrix, "DiscernibilityMatrix")) {
    stop("The argument is not in the class of \'DiscernibilityMatrix\' objects.")
  }

  reductSet <- convertCNFtoDNF(discernibilityMatrix$disc.list)

  core <- computeCore(reductSet)

  reductSet <- list(decision.reduct = reductSet,
                    core = core,
                    discernibility.type = discernibilityMatrix$type.discernibility,
                    type.task = "computation of all reducts",
                    type.model = discernibilityMatrix$type.model)
  reductSet <- ObjectFactory(reductSet, classname = "ReductSet")

  return(reductSet)
}

#' It is a function for computing one reduct from a discernibility matrix - it can use
#' the greedy heuristic or a randomized (Monte Carlo) search.
#'
#' @title Computing one reduct from a discernibility matrix
#' @author Andrzej Janusz
#'
#' @param discernibilityMatrix a \code{"DiscernibilityMatrix"} class representing the discernibility matrix of RST and FRST.
#' @param greedy a boolean value indicating whether the greedy heuristic or a randomized search should be used in computations.
#' @param power a numeric representing a parameter of the randomized search heuristic.
#'
#' @return A class \code{"ReductSet"}.
#'
#' @seealso \code{\link{BC.discernibility.mat.RST}} and \code{\link{BC.discernibility.mat.FRST}}.
#'
#' @references
#' Jan G. Bazan, Hung Son Nguyen, Sinh Hoa Nguyen, Piotr Synak, and Jakub Wroblewski,
#' "Rough Set Algorithms in Classification Problem", Chapter 2
#' In: L. Polkowski, S. Tsumoto and T.Y. Lin (eds.): Rough Set Methods and Applications
#' Physica-Verlag, Heidelberg, New York, p. 49 - 88 ( 2000).
#'
#' @examples
#' ########################################################
#' ## Example 1: Generate one reducts and
#' ##            a new decision table using RST
#' ########################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## build the decision-relation discernibility matrix
#' res.1 <- BC.discernibility.mat.RST(decision.table, range.object = NULL)
#'
#' ## generate all reducts
#' reduct <- FS.one.reduct.computation(res.1)
#'
#' ## generate new decision table
#' new.decTable <- SF.applyDecTable(decision.table, reduct, control = list(indx.reduct = 1))
#'
#' ##############################################################
#' ## Example 2: Generate one reducts and
#' ##            a new decision table using FRST
#' ##############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## build the decision-relation discernibility matrix
#' control <- list(type.relation = c("crisp"),
#'                 type.aggregation = c("crisp"),
#'                 t.implicator = "lukasiewicz", type.LU = "implicator.tnorm")
#' res.2 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "standard.red",
#'                                     control = control)
#'
#' ## generate single reduct
#' reduct <- FS.one.reduct.computation(res.2)
#'
#' ## generate new decision table
#' new.decTable <- SF.applyDecTable(decision.table, reduct, control = list(indx.reduct = 1))
#' @export
FS.one.reduct.computation <- function(discernibilityMatrix, greedy = TRUE, power = 1) {

  if(!inherits(discernibilityMatrix, "DiscernibilityMatrix")) {
    stop("The argument is not in the class of \'DiscernibilityMatrix\' objects.")
  }

  CNFclauses = discernibilityMatrix$disc.list
  clauseLengths = sapply(CNFclauses, length)

  reduct = character()
  tmpIdx = (clauseLengths > 0)
  attrInvertedIdx = list()

  while(any(tmpIdx)) {
    attributeCounts = table(unlist(CNFclauses[tmpIdx], use.names = FALSE))

    if(greedy) bestAttr = which.max(attributeCounts)
    else bestAttr = sample(length(attributeCounts), 1, prob = attributeCounts^power)

    tmpName = names(attributeCounts)[bestAttr]
    reduct = append(reduct, tmpName)
    attrInvertedIdx = append(attrInvertedIdx, list(sapply(CNFclauses, function(clause) tmpName %in% clause)))
    tmpIdx = tmpIdx & !attrInvertedIdx[[length(attrInvertedIdx)]]
  }
  rm(tmpIdx, tmpName, attributeCounts, bestAttr)

  if(length(reduct) > 1) {
    attrInvertedIdx = do.call(cbind, attrInvertedIdx)
    tmpSums = rowSums(attrInvertedIdx)
    for(i in (ncol(attrInvertedIdx) - 1):1) {
      if(all(tmpSums > 1) || !(any(attrInvertedIdx[tmpSums == 1,i]))) {
        reduct = reduct[-i]
        attrInvertedIdx = attrInvertedIdx[ , -i, drop = FALSE]
        tmpSums = rowSums(attrInvertedIdx)
      }
    }
  }

  reduct <- list(decision.reduct = list(reduct),
                 core = NULL,
                 discernibility.type = discernibilityMatrix$type.discernibility,
                 type.task = "computation of one reduct from a discernibility matrix",
                 type.model = discernibilityMatrix$type.model)
  reduct <- ObjectFactory(reduct, classname = "ReductSet")
  return(reduct)
}




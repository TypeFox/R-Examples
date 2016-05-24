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
## it is used to represent rules in order to get the standard format
# @param rules original form of rules
# @param type.method a type of method
# @param decision.table a decision table
# @param t.similarity a type of similarity equation
# @param t.tnorm a type of t-tnorm
std.rules <- function(rules, type.method = "RI.hybridFS.FRST", decision.table, t.similarity = NULL, t.tnorm = NULL){
	## get data
	objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	decision.attr <- attr(decision.table, "decision.attr")
	consequent.attr <- names(desc.attrs)[decision.attr]

	if (any(type.method == c("RI.hybridFS.FRST", "RI.GFRS.FRST"))){
		all.attr = colnames(objects)

		if (nominal.att[length(nominal.att)] == TRUE){
			type.task <- c("classification")
		}
		else {
			type.task <- c("regression")
		}

		## check duplication attributes in rules for collecting their objects
		duplicated.pt <- duplicated(rules$attributes)
		new.obj <- c()
		ii <- 1

		## collect the objects that have the same attributes in the rules
		for (i in 1 : length(duplicated.pt)){
			if (duplicated.pt[i] == FALSE){
				new.obj <- append(new.obj, list(rules$objects[[i]]))
				ii <- ii + 1
			}
			else {
				new.obj[[ii - 1]] <- c(new.obj[[ii - 1]], rules$objects[[i]])
			}
		}

		if (type.method == "RI.hybridFS.FRST"){
			rules <- list(attributes = rules$attributes[which(duplicated.pt == FALSE)], objects = new.obj)
			attrs.rules <- unique(unlist(rules$attributes))
			rules <- ch.rules(rules, type.method = type.method, decision.table)

			## calculate variance and range of data
			attrs.p <- c(which(colnames(objects) %in% attrs.rules))
			res.varRange <- cal.var.range(objects, attrs.p, nominal.att)
			variance.data <- res.varRange$variance.data
			rownames(variance.data) <- NULL
			colnames(variance.data) <- colnames(objects)[attrs.p]
			range.data <- res.varRange$range.data
			colnames(range.data) <- colnames(objects)[attrs.p]

			mod <- list(rules = rules, type.model = "FRST", type.method = type.method, type.task = type.task,
						t.similarity = t.similarity, t.tnorm = t.tnorm, variance.data = variance.data, range.data = range.data,
						antecedent.attr = colnames(objects)[attrs.p], consequent.attr = consequent.attr, nominal.att = nominal.att[attrs.p])
			class.mod <- ObjectFactory(mod, classname = "RuleSetFRST")
			return(class.mod)
		}
		else if (type.method == "RI.GFRS.FRST"){
			rules <- list(attributes = rules$attributes[which(duplicated.pt == FALSE)], objects = new.obj,
			            threshold = rules$threshold[sort(unlist(new.obj))])
			attrs.rules <- unique(unlist(rules$attributes))
			rules <- ch.rules(rules, type.method = type.method, decision.table)

			## calculate variance and range of data
			attrs.p <- c(which(colnames(objects) %in% attrs.rules))
			res.varRange <- cal.var.range(objects, attrs.p, nominal.att)
			variance.data <- res.varRange$variance.data
			rownames(variance.data) <- NULL
			colnames(variance.data) <- colnames(objects)[attrs.p]
			range.data <- res.varRange$range.data
			colnames(range.data) <- colnames(objects)[attrs.p]

			mod <- list(rules = rules, type.model = "FRST", type.method = type.method,
			            type.task = type.task, t.similarity = t.similarity, t.tnorm = t.tnorm, variance.data = variance.data, range.data = range.data,
						antecedent.attr = colnames(objects)[attrs.p], consequent.attr = consequent.attr, nominal.att = nominal.att[attrs.p])
			class.mod <- ObjectFactory(mod, classname = "RuleSetFRST")
			return(class.mod)
		}
	}
}

## it is used to change rules in order to get the standard format
# @param rules original form of rules
# @param type.method a type of method
# @param decision.table
ch.rules <- function(rules, type.method = "RI.hybridFS.FRST", decision.table){
	tra.objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	rule <- c()

	for (i in 1 : length(rules$attributes)){
		exit <- FALSE
		ii <- 1
		while (exit == FALSE) {
			## get value for antecedent and consequent of each rule
			anteVal <- tra.objects[rules$objects[[i]][ii], c(rules$attributes[[i]]), drop = FALSE]
			conqVal <- tra.objects[rules$objects[[i]][ii], ncol(tra.objects), drop = FALSE]
			if (is.null(rule)){
				tmpRule <- cbind(anteVal, conqVal)
				rownames(tmpRule) <- "value"
				rule <- list(tmpRule)
			}
			else {
				tmpRule <- cbind(anteVal, conqVal)
				rownames(tmpRule) <- "value"
				rule <- append(rule, list(tmpRule))
			}

			## checking stopping criteria
			if (ii == length(rules$objects[[i]])){
				exit <- TRUE
			}
			else {
				ii <- ii + 1
			}
		}
	}

	if (type.method == "RI.GFRS.FRST"){
		return(list(rules = rule, threshold = rules$threshold))
	}
	else {
		return(list(rules = rule))
	}
}

#an auxiliary function for computing laplace estimate of rule's confidence
laplaceEstimate <- function(rule, dataS, clsVec, uniqueCls, suppIdx = NULL) {
	if (is.null(suppIdx)) {
		if (length(rule$idx) > 1) {
			dataS <- do.call(paste, dataS[rule$idx])
		}
		else {
			dataS <- as.character(dataS[[rule$idx]])
		}

		tmpValues = paste(rule$values, collapse = " ")
		suppIdx = which(dataS == tmpValues)
	}

	clsFreqs = table(clsVec[suppIdx])
	maxIdx = which.max(clsFreqs)
	nOfCorrect = clsFreqs[maxIdx]

	rule$consequent = names(clsFreqs)[maxIdx]
	rule$support = as.integer(suppIdx)
	rule$laplace = (nOfCorrect + 1)/(length(suppIdx) + length(uniqueCls))
	return(rule)
}

# An auxiliary function that performs the best-first voting strategy.
#
# @title the best-first voting strategy function
# @param object a class "RuleSetRST".
# @param ruleValues values of rules.
# @param consequents values on the consequent part.
# @param majorityCls  a value representing majority class of decision attribute.
# @param ... other parameters.
# @return predicted values.
bestFirst <- function(object, ruleValues, consequents, majorityCls, ...) {
	matchIdx <- which(ruleValues == object)
	if (length(matchIdx) > 0) {
		prediction <- consequents[matchIdx[1]]
	}
	else {
		prediction <- majorityCls
	}

	return(prediction)
}

# An auxiliary function that performs the first-win voting strategy.
#
# @title the first-win voting strategy
# @param a data vector representing the object for which predictions are to be made
# @param object a class "RuleSetRST"
# @return a predicted value
firstWin <- function(object, ruleSet) {
  rulesCount = length(ruleSet)
  i = 1
  endFlag = FALSE
  prediction = NULL

  while(!endFlag) {
    if(all(object[ruleSet[[i]]$idx] == ruleSet[[i]]$values))  {
      endFlag = T
      prediction = ruleSet[[i]]$consequent
    } else  {
      if(i >= rulesCount) endFlag = TRUE
      i = i + 1
    }
  }

  if(is.null(prediction)) prediction = attr(ruleSet, "majorityCls")
  return(prediction)
}

# An auxiliary function that classifies a data object using custom voting strategy.
rulesVoting <- function(object, ruleSet, votingMethod = X.ruleStrength, ...) {

  rulesIdx = which(sapply(ruleSet,
                          function(rule, obj) all(obj[rule$idx] == rule$values),
                          object))

  prediction = NULL
  voteCounter = rep(0, length(attr(ruleSet, "uniqueCls")))
  names(voteCounter) = attr(ruleSet, "uniqueCls")

  if(length(rulesIdx) > 0) {
    for(i in rulesIdx) {
      voteCounter[ruleSet[[i]]$consequent] = voteCounter[ruleSet[[i]]$consequent] + votingMethod(ruleSet[[i]], ...)
    }
    tmpPreds = names(voteCounter)[which(voteCounter == max(voteCounter))]
    if(length(tmpPreds) == 1) {
      prediction = tmpPreds
    } else {
      prediction = names(attr(ruleSet, "clsProbs")[tmpPreds])[which.max(attr(ruleSet, "clsProbs")[tmpPreds])]
    }
  }

  if(is.null(prediction)) prediction = attr(ruleSet, "majorityCls")

  return(prediction)
}

#' A function returning a weight of rule's vote understood as strength of the rule.
#' It is defined as a product of a cardinality of a support of a rule and the length of this rule.
#'
#' @title Rule voting by strength of the rule
#' @author Andrzej Janusz
#'
#' @param rule a decision rule, i.e. element of a "RuleSetRST" object
#'
#' @return a numerical weight of the vote
#'
#' @seealso Other currently available voting methods are: \code{\link{X.laplace}}, \code{\link{X.rulesCounting}}.
#'
#' @export
X.ruleStrength <- function(rule) {
	return(length(rule$support) * length(rule$idx))
}

#' A function returning a weight of rule's vote understood as the Laplace estimate of its confidence.
#'
#' @title Rule voting by the Laplace estimate
#' @author Andrzej Janusz
#'
#' @param rule a decision rule, i.e. element of an "RuleSetRST" object
#'
#' @return a numerical weight of the vote
#'
#' @seealso Other currently available voting methods are: \code{\link{X.ruleStrength}}, \code{\link{X.rulesCounting}}.
#'
#' @export
X.laplace <- function(rule) {
	return(rule$laplace)
}

#' A function returning an equal vote's weight for every rule. It corresponds to voting by counting the
#' matching rules.
#'
#' @title Rule voting by counting matching rules
#' @author Andrzej Janusz
#'
#' @param rule a decision rule, i.e. element of an "RuleSetRST" object
#'
#' @return a numerical weight of the vote
#'
#' @seealso Other currently available voting methods are: \code{\link{X.ruleStrength}}, \code{\link{X.laplace}}.
#'
#' @export
X.rulesCounting = function(rule) {
	return(1)
}

# It is an auxiliary function for greedy covering the data with decision rules using CN2 algorithm.
findBestCN2Rule = function(dataSet, clsVec, uniqueCls, K = 10)  {
  descriptorsList = lapply(dataSet, function(x) sort(unique(x)))

  descriptorCandidates = list()
  for(i in 1:length(descriptorsList)) {
    descriptorCandidates = c(descriptorCandidates,
                             lapply(descriptorsList[[i]],
                                    function(v, x) return(list(idx = x, values = v)), i))
  }

  ruleCandidates = lapply(descriptorCandidates, laplaceEstimate, dataSet, clsVec, uniqueCls)

  endFlag = F
  ruleScores = sapply(ruleCandidates, function(x) return(x$laplace))
  bestIdx = which.max(ruleScores)
  bestRule = ruleCandidates[[bestIdx]]
  ruleCandidates = ruleCandidates[order(ruleScores, decreasing = T)[1:K]]

  while(!endFlag) {

    ruleCandidates = addDescriptorCN2(ruleCandidates, descriptorCandidates)
    ruleCandidates = lapply(ruleCandidates, laplaceEstimate, dataSet, clsVec, uniqueCls)
    ruleScores = sapply(ruleCandidates, function(x) return(x$laplace))

    if(bestRule$laplace < max(ruleScores)) {
      bestIdx = which.max(ruleScores)
      bestRule = ruleCandidates[[bestIdx]]
      ruleCandidates = ruleCandidates[order(ruleScores, decreasing = T)[1:K]]
		}
		else endFlag = T

  }
  return(bestRule)
}

# Adding a descriptor to a rule in CN2 algorithm
addDescriptorCN2 <- function(rulesList, descCandidates) {
  candidates = list()
  for(i in 1:length(rulesList)) {
    candidates = c(candidates,
                   lapply(descCandidates,
                          function(descriptor, rule) return(list(idx = c(rule$idx, descriptor$idx),
                                                                 values = c(as.character(rule$values), as.character(descriptor$values)))),
                          rulesList[[i]]))
  }

  return(candidates)
}

# Computation of a covering of a lower approximation of a concept by decision rules using the LEM2 algorithm.
# It is an auxiliary function.
computeLEM2covering <- function(concept, attributeValuePairs, decisionValues, uniqueCls) {

  if(length(concept) == 0) stop("Empty lower approximation of a decision class.")
  uncoveredConcept = concept
  rules = list()

  while(length(uncoveredConcept) > 0) {
    support = 1:length(decisionValues)
    rules[[length(rules) + 1]] = list(idx = integer(0), values = character(), support = support)
    tmpAttributeValuePairs = attributeValuePairs
    selectedAttributeValuePairs = list()

    while(!all(support %in% concept)) {
      tmpSupp = intersect(support, uncoveredConcept)
      correctLaplaceVec = sapply(tmpAttributeValuePairs,
                                 function(rule) {
                                   nOfCorrect = length(intersect(rule$support, tmpSupp));
                                   if(nOfCorrect > 0) laplace = nOfCorrect + rule$laplace  #(nOfCorrect + 1)/(length(tmpSupp) + length(uniqueCls))   #(length(rule$support) + length(uniqueCls))
                                   else laplace = 0;
                                   return(laplace)
                                 })
      tmpIdx = which.max(correctLaplaceVec)
      tmpRule = tmpAttributeValuePairs[[tmpIdx]]
      selectedAttributeValuePairs[[length(selectedAttributeValuePairs) + 1]] = tmpRule
      tmpAttributeValuePairs = tmpAttributeValuePairs[-tmpIdx]
      support = intersect(support, tmpRule$support)
    }
    rm(tmpSupp, tmpIdx, tmpRule, tmpAttributeValuePairs)

    toRmIdx = integer(0)
    if(length(selectedAttributeValuePairs) > 1) {
      for(i in 1:length(selectedAttributeValuePairs)) {
        suppList = lapply(selectedAttributeValuePairs[-c(toRmIdx, i)], function(x) x$support)
        if(length(suppList) > 0) {
          tmpSupport = Reduce(intersect, suppList)
        } else {
          tmpSupport = -1
        }
        if(all(tmpSupport %in% concept)) toRmIdx = c(toRmIdx, i)
      }
    }

    if(length(toRmIdx) > 0) {
      selectedAttributeValuePairs = selectedAttributeValuePairs[-toRmIdx]
      suppList = lapply(selectedAttributeValuePairs, function(x) x$support)
      support = Reduce(intersect, suppList)
    }

    idxVec = sapply(selectedAttributeValuePairs, function(x) x$idx)
    valuesVec = sapply(selectedAttributeValuePairs, function(x) x$values)
    rules[[length(rules)]]$support = support
    rules[[length(rules)]]$values = valuesVec
    rules[[length(rules)]]$idx = idxVec

    uncoveredConcept = setdiff(uncoveredConcept, support)
  }

  toRmIdx = integer(0)
  for(i in 1:length(rules)) {
    suppList = lapply(rules[-c(toRmIdx, i)], function(x) x$support)
    tmpSupport = Reduce(intersect, suppList)
    if(all(concept %in% tmpSupport)) toRmIdx = c(toRmIdx, i)
  }

  if(length(toRmIdx) > 0) {
    rules = rules[-toRmIdx]
  }

  return(rules)
}

# Computation of a covering of a lower approximation of a concept by decision rules using the AQ algorithm.
# It is an auxiliary function.
computeAQcovering <- function(concept, attributeValuePairs, dataTab, epsilon = 0.05, K = 2) {

  if(length(concept) == 0) stop("Empty lower approximation of a decision class.")
  uncoveredConcept = concept
  coverCounts = rep(0, length(concept))
  names(coverCounts) = as.character(concept)
  rules = list()

  while(length(uncoveredConcept) > 0) {
    seedIdx = sample(uncoveredConcept, 1)
    selectedAttributeValuePairs = mapply(function(avps, v) avps[[as.character(v)]],
                                         attributeValuePairs, dataTab[seedIdx,],
                                         SIMPLIFY = FALSE)

    suppList = lapply(selectedAttributeValuePairs, function(x) x$support)
    support = Reduce(intersect, suppList)
    rules[[length(rules) + 1]] = list(idx = integer(0), values = character(), support = support)

    attrOrdering = sample(1:length(selectedAttributeValuePairs))
    suppList = suppList[attrOrdering]
    selectedAttributeValuePairs = selectedAttributeValuePairs[attrOrdering]

    for(i in length(attrOrdering):1)  {
      tmpSupport = Reduce(intersect, suppList[-i])
      if(length(tmpSupport) > 0 && sum(tmpSupport %in% concept)/length(tmpSupport) >= 1-epsilon) {
        support = tmpSupport
        suppList = suppList[-i]
        selectedAttributeValuePairs = selectedAttributeValuePairs[-i]
      }
    }
    rm(tmpSupport, attrOrdering, suppList)

    idxVec = sapply(selectedAttributeValuePairs, function(x) x$idx)
    valuesVec = sapply(selectedAttributeValuePairs, function(x) x$values)
    rules[[length(rules)]]$support = support
    rules[[length(rules)]]$values = valuesVec
    rules[[length(rules)]]$idx = idxVec

    coveredConcept = intersect(support, concept)
    coverCounts[as.character(coveredConcept)] = coverCounts[as.character(coveredConcept)] + 1
    uncoveredConcept = concept[which(coverCounts < K)]
  }

  for(i in length(rules):1) {
    tmpCoveredConcept = intersect(concept, rules[[i]]$support)
    tmpCoverCounts = coverCounts
    tmpCoverCounts[as.character(tmpCoveredConcept)] = tmpCoverCounts[as.character(tmpCoveredConcept)] - 1
    if(all(tmpCoverCounts >= K)) {
      rules = rules[-i]
      coverCounts = tmpCoverCounts
    }
  }
  rm(tmpCoveredConcept, tmpCoverCounts)

  return(rules)
}


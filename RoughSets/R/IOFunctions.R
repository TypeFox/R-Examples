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
#' This function can be used to import data sets from files and then construct a \code{DecisionTable} object. It uses
#' \code{\link{read.table}} function from \code{base} R.
#'
#' The data should be in a tabular format containing rows and columns, where every row represents
#' an object/instance, while columns represent attributes of the objects.
#'
#' @title Reading tabular data from files.
#' @author Andrzej Janusz
#'
#' @param filename a path with a file name.
#' @param decision.attr an integer indicating an index of the decision attribute. See \code{\link{SF.asDecisionTable}}.
#' @param indx.nominal an integer vector with indices of attributes which should be considered as nominal.
#'        See \code{\link{SF.asDecisionTable}}.
#' @param ... additional parameters which are passed to the \code{read.table} function. See \code{\link{read.table}}.
#'
#' @return An object of the \code{"DecisionTable"} class. See \code{\link{SF.asDecisionTable}}.
#'
#' @examples
#' #############################################################
#' ## Example 1: data set saved in a file
#' #############################################################
#' ## Let us assume we have the following data which has been already saved to the file "tes.dat"
#' data <- data.frame(c(0.12, 0.23, 0.24), c(1,3,2), c(10, 12, 18), c("a", "a", "b"), c(1, 1, 2))
#' \dontrun{write.table(data, file = "tes.dat", row.names = FALSE, col.names = FALSE,
#'                     fileEncoding ="")}
#'
#' ## Then we would generate decision table from tes.dat file.
#' ## in this case, we want to define that second and third attributes are nominal and continuous,
#' ## respectively.
#' \dontrun{decision.table <- SF.read.DecisionTable(filename = "tes.dat", decision.attr = 5,
#'                   indx.nominal = c(2, 5), sep= " ", col.names = c("v1", "v2", "v3", "v4", "o1"))}
#'
#' @export
SF.read.DecisionTable <- function(filename, decision.attr = NULL, indx.nominal = NULL, ...) {

  if(is.null(decision.attr)) {
    warning("A decision attribute is not indicated - the data will be treated as an information system.")
  }

  dataset = utils::read.table(file = filename, ...)
  decision.table = SF.asDecisionTable(dataset, decision.attr = decision.attr, indx.nominal = indx.nominal)
  return(decision.table)
}

#' This function converts \code{data.frames} into \code{DecisionTable} objects. This is a standard data representation
#' in the \code{RoughSets} package.
#'
#' An object of the \code{"DecisionTable"} class adds a few attributes to a standard data.frame:
#' \itemize{
#'   \item \code{desc.attrs}: a list containing the names of attributes and their range/possible symbolic values.
#'                            There are two kinds of representation in this parameters which depend on whether the attributes are
#'                            nominal or numeric, for example:
#'                            \itemize{
#'                              \item nominal attribute: \code{a = c(1,2,3)} means that the attribute \code{a} has values 1, 2, and 3.
#'                              \item numeric attribute: \code{a = c(10, 20)} means that the attribute \code{a} has values between 10 and 20.
#'                            }
#'   \item \code{nominal.attrs}: a logical vector whose length equals the number of columns in the data. In this vector \code{TRUE} values
#'                            indicate that the corresponding attribute is a nominal. For example:
#'                            \code{nominal.attrs = c(FALSE, TRUE, FALSE)} means that the first and third attributes
#'                            are numeric and the second one is nominal.
#'  \item \code{decision.attr}: a numeric value representing the index of the decision attribute. It is necessary to define
#'                            the index of the decision attribute in order to construct a proper decision system. If the value
#'                            of \code{decision.attr} is NULL, the constructed object will correspond to an information system.
#'                            It is strongly recommended to place the decision attribute as the last data column.
#' }
#' \code{"DecisionTable"} objects allow to use all methods of standard data.frame objects.
#' The function \code{\link{SF.read.DecisionTable}} can be used to import data from a file and then construct \code{DecisionTable} object.
#'
#' @title Converting a data.frame into a \code{DecisionTable} object
#' @author Andrzej Janusz
#'
#' @param dataset data.frame that contains objects/instances and attributes/features in its rows and columns, respectively.
#'        See in Section \code{Details}.
#' @param decision.attr an integer value representing the index position of the decision attribute. If this parameter is ignored, then
#'        the function will treat the data as an information system or newdata/test data. In other words,
#'        it is necessary to define the index of the decision attribute in order to construct a decision table (e.g. a training data set).
#' @param indx.nominal a logical vector indicating nominal attributes in the data.
#'         If this parameter is not given, then the function will use a heuristic to guess which of the attributes are nominal.
#'         The following rules will be applied used:
#' 			\itemize{
#' 			\item an attribute contains character values or factors: it will be recognized as a nominal attribute.
#' 			\item an attribute contains integer or numeric values: it will be recognized as a numeric attribute.
#' 			\item indx.nominal: the indicated attributes will be considered as nominal.
#' 			}
#'
#' @return An object of the \code{"DecisionTable"} class.
#'
#' @seealso \code{\link{SF.read.DecisionTable}}, \code{\link{SF.applyDecTable}}.
#'
#' @examples
#' ################################################################
#' ## Example : converting from datasets in data.frame
#' ##            into decision table
#' ################################################################
#' ## Let use iris which is available in R be dataset
#' decision.table <- SF.asDecisionTable(dataset = iris, decision.attr = 5,
#'                   indx.nominal = 5)
#' @export
SF.asDecisionTable <- function(dataset, decision.attr = NULL, indx.nominal = NULL) {

  nominal.attrs = rep(FALSE, ncol(dataset))
  if (!is.null(indx.nominal)) {
		nominal.attrs[indx.nominal] = TRUE
  }

  class.list = lapply(dataset, class)
  nominal.attrs[sapply(class.list, function(x) any(x %in% c("factor", "character")))] = TRUE

  if(any(nominal.attrs)) {
    dataset[nominal.attrs] = lapply(dataset[nominal.attrs], factor)
  }

  ## construct desc.attrs as a description of attributes
  desc.attrs = list()
  desc.attrs[1:ncol(dataset)] = numeric(1)
  indx.nominal = which(nominal.attrs)
  if(length(indx.nominal) > 0) {
		desc.attrs[indx.nominal] = lapply(dataset[indx.nominal], levels)
    if(sum(nominal.attrs) != ncol(dataset)) {
      desc.attrs[which(!nominal.attrs)] = lapply(dataset[which(!nominal.attrs)], range, na.rm = TRUE)
    }
	}
	else {
    desc.attrs = lapply(dataset, range)
  }

  ## construct the class "DecisionTable"
  names(desc.attrs) <- colnames(dataset)
  decision.table = data.frame(dataset)
  attr(decision.table, "nominal.attrs") = nominal.attrs
  attr(decision.table, "desc.attrs") = desc.attrs
  attr(decision.table, "decision.attr") = decision.attr
  decision.table = ObjectFactory(decision.table, "DecisionTable")

  return(decision.table)
}

#' This function enables the output of a summary of the rule induction methods.
#'
#' @title The summary function of rules based on FRST
#' @author Lala Septem Riza
#'
#' @param object a \code{"RuleSetFRST"} object. See \code{\link{RI.hybridFS.FRST}} and \code{\link{RI.GFRS.FRST}}.
#' @param ... the other parameters.
#' @return a description that contains the following information:
#' \itemize{
#' \item The type of the considered model.
#' \item The type of the considered method.
#' \item The type of the considered task.
#' \item The type of similarity.
#' \item The type of triangular norm.
#' \item The names of attributes and their type (whether nominal or not).
#' \item The interval of the data.
#' \item the variance values of the data.
#' \item The rules. Every rule constitutes two parts which are IF and THEN parts.
#'       For example, \code{"IF pres is around 90 and preg is around 8 THEN class is 2"}.
#'       See \code{\link{RI.GFRS.FRST}}.
#' }
#' @examples
#' ###########################################################
#' ## Example 1: Regression problem
#' ###########################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$housing7.dt
#'
#' control <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), type.relation =
#'                 c("tolerance", "eq.3"), t.implicator = "lukasiewicz")
#' res.1 <- RI.hybridFS.FRST(decision.table, control)
#'
#' summary(res.1)
#' ###########################################################
#' ## Example 2: Classification problem
#' ##############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$pima7.dt
#'
#' control <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), type.relation =
#'                 c("tolerance", "eq.3"), t.implicator = "lukasiewicz")
#' res.2 <- RI.hybridFS.FRST(decision.table, control)
#'
#' summary(res.2)
#' @export
#' @method summary RuleSetFRST
summary.RuleSetFRST <- function(object, ...){

 if(!inherits(object, "RuleSetFRST")) stop("not a legitimate object in this package")
 cat("The type of the considered model: ", "\n")
 print(object$type.model)
 cat("The type of the considered method: ", "\n")
 print(object$type.method)
 cat("The type of the considered task: ", "\n")
 print(object$type.task)
 cat("The type of similarity: ", "\n")
 print(object$t.similarity)
 cat("The type of triangular norm: ", "\n")
 print(object$t.tnorm)
 cat("The names of attributes and their type (whether nominal or not): ", "\n")
 temp = matrix(object$nominal.att, nrow = 1)
 colnames(temp) <- matrix(object$antecedent.attr, nrow = 1)
 print(temp)
 cat("The interval of the data: ", "\n")
 temp <- colnames(object$variance.data)
 colnames(object$range.data) <- temp
 print(object$range.data)
 cat("The variance values of the data: ", "\n")
 print(object$variance.data)
 cat("The rules : ", "\n")
 rules <- toStr.rules(object$rules$rules, object$type.task, object$nominal.att)
 print(rules)

 invisible(object)
}

#' This function enables the output of a summary of the rule induction methods.
#'
#' @title The summary function of rules based on RST
#' @author Lala Septem Riza and Andrzej Janusz
#'
#' @param object a \code{"RuleSetRST"} object. See \code{\link{RI.indiscernibilityBasedRules.RST}}.
#' @param ... the other parameters.
#' @return a description that contains the following information:
#' \itemize{
#' \item The type of the considered model.
#' \item The type of the considered method.
#' \item The type of the considered task.
#' \item The rules. Every rule constitutes two parts which are IF and THEN parts.
#'       For example, \code{"IF pres is around 90 and preg is around 8 THEN class is 2; (support=4;laplace=0.67)"}.
#' }
#' @examples
#' ###########################################################
#' ## Example : Classification problem
#' ###########################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## determine feature subset/reduct
#' reduct <- FS.permutation.heuristic.reduct.RST(decision.table,  permutation = NULL)
#'
#' rules <- RI.indiscernibilityBasedRules.RST(decision.table, reduct)
#'
#' summary(rules)
#' @export
#' @method summary RuleSetRST
summary.RuleSetRST <- function(object, ...){

 if(!inherits(object, "RuleSetRST")) stop("not a legitimate object in this package")
 cat("The type of the considered model: ", "\n")
 print("RST")
 cat("The type of the considered method: ", "\n")
 print(attr(object, "method"))
 print(object)

 invisible(object)
}


#' A print method for RuleSetRST objects.
#'
#' @title The print function for RST rule sets
#' @author Andrzej Janusz
#'
#' @param x a \code{"RuleSetRST"} object. See \code{\link{RI.LEM2Rules.RST}}.
#' @param howMany an integer giving the number of rules to be printed.
#'        The default is minimum from 10 and the total number of rules in the set.
#' @param ... the other parameters.
#' @return prints its argument and returns it invisibly
#' @examples
#' ###########################################################
#' ## Example : Printing of a decision rule set problem
#' ###########################################################
#' data(RoughSetData)
#' hiring.data <- RoughSetData$hiring.dt
#'
#' rules <- RI.LEM2Rules.RST(hiring.data)
#'
#' rules             # all rules are printed
#' print(rules, 2)   # only the first two rules are printed
#'
#' # printing a subset of rules
#' rules[2:3]
#' @export
#' @method print RuleSetRST
print.RuleSetRST <- function(x, howMany = min(10, length(x)), ...){

  if(!inherits(x, "RuleSetRST")) stop("not a legitimate object in this package")
  if(howMany < 1 || howMany > length(x)) stop("wrong parameter value")
  if(length(x) > 1) {
    cat("A set consisting of ", length(x), " rules:\n")
  } else {
    if(length(x) == 1) {
      cat("A set consisting of 1 rule:\n")
    } else stop("Empty rule set")
  }
  xTmp <- x[1:min(howMany,length(x))]
  rules <- sapply(xTmp, convertRuleIntoCharacter, attr(x, 'colnames'))
  rules <- mapply(function(x, n) paste(n, ". ", x, sep = ""), rules, 1:length(rules), SIMPLIFY = FALSE)
  lapply(rules, function(x) cat(x, "\n"))
  if(length(x) > howMany) {
    if(length(x) - howMany > 1) cat(paste0('... and ', length(x) - howMany, ' other rules.\n'))
    else cat('... and 1 other rule.\n')
  }

  invisible(x)
}

#' A function for converting a set of rules into their character representation.
#'
#' @title The \code{as.character} method for RST rule sets
#' @author Andrzej Janusz
#'
#' @param x a \code{"RuleSetRST"} object. See \code{\link{RI.LEM2Rules.RST}}.
#' @param ... the other parameters.
#' @return Converts rules from a set into their character representation.
#' @examples
#' ###########################################################
#' ## Example : Converting a set of decision rules
#' ###########################################################
#' data(RoughSetData)
#' hiring.data <- RoughSetData$hiring.dt
#'
#' rules <- RI.LEM2Rules.RST(hiring.data)
#'
#' as.character(rules)
#' @export
#' @method as.character RuleSetRST
as.character.RuleSetRST = function(x, ...) {

  if(!inherits(x, "RuleSetRST")) stop("not a legitimate object in this package")
  colNames <- attr(x, "colnames")
  rules <- sapply(x, convertRuleIntoCharacter, colNames)
  rules <- sub('\n\t\t', ' ', rules)
  rules
}


# auxiliary function used in print and as.character methods of RuleSetRST objects
convertRuleIntoCharacter = function(x, colNames) {
  desc <- paste(colNames[x$idx[1]], x$values[1], sep = " is ")
  if(length(x$values) > 1) {
    for (j in 2 : length(x$values)){
      temp <- paste(colNames[x$idx[j]], x$values[j], sep = " is ")
      desc <- paste(desc, temp, sep = " and ")
    }
  }
  cons <- paste(attr(x, "dec.attr"), paste(x$consequent, ";\n\t\t(supportSize=",
                                           length(x$support), "; ", "laplace=",
                                           x$laplace,")", sep=""), sep = c(" is "))
  rule <- paste("IF", desc, "THEN", cons)
  rule
}

#' Subsetting a set of decision rules.
#'
#' @title The \code{[.} method for \code{"RuleSetRST"} objects
#' @author Andrzej Janusz
#'
#' @param x a \code{"RuleSetRST"} object from which to extract rules(s) or in which to replace rules(s).
#'        See \code{\link{RI.LEM2Rules.RST}}.
#' @param i integer indices specifying elements to be extracted or replaced.
#' @param ... the other parameters.
#' @return A subset of rules.
#' @examples
#' ###########################################################
#' ## Example : Subsetting a set of decision rules
#' ###########################################################
#' data(RoughSetData)
#' hiring.data <- RoughSetData$hiring.dt
#'
#' rules <- RI.LEM2Rules.RST(hiring.data)
#'
#' rules
#'
#' # taking a subset of rules
#' rules[1:3]
#' rules[c(TRUE,FALSE,TRUE,FALSE)]
#'
#' # replacing a subset of rules
#' rules2 <- rules
#' rules2[c(2,4)] <- rules[c(1,3)]
#' rules2
#' @export
#' @aliases Extract.RuleSetRST
#' @method [ RuleSetRST
"[.RuleSetRST" = function (x, i, ...) {
  tmp <- attributes(x)
  x <- as.list(x)[i]
  attributes(x) <- tmp
  x
}

#' A function for converting a set of rules into a list.
#'
#' @title The \code{as.list} method for RST rule sets
#' @author Andrzej Janusz
#'
#' @param x a \code{"RuleSetRST"} object. See \code{\link{RI.LEM2Rules.RST}}.
#' @param ... the other parameters.
#' @return Converts rules from a set into a list.
#' @examples
#' ###########################################################
#' ## Example : Converting a set of decision rules
#' ###########################################################
#' data(RoughSetData)
#' hiring.data <- RoughSetData$hiring.dt
#'
#' rules <- RI.LEM2Rules.RST(hiring.data)
#'
#' as.list(rules)
#' @export
#' @method as.list RuleSetRST
as.list.RuleSetRST = function(x, ...) {
  class(x) = 'list'
  x
}


#' This is a print method for FeatureSubset objects.
#'
#' @title The print method of FeatureSubset objects
#' @author Andrzej Janusz
#'
#' @param x an object inheriting from \code{"FeatureSubset"} class. See \code{\link{FS.reduct.computation}}.
#' @param ...  parameters passes to other functions (currently omitted).
#' @return Prints its argument and returns it invisibly.
#' @examples
#' ###########################################################
#' ## Example : Computation of a decision reduct
#' ###########################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' res.1 <- FS.reduct.computation(decision.table)
#' print(res.1)
#' @export
#' @method print FeatureSubset
print.FeatureSubset <- function(x, ...){

  if(!inherits(x, "FeatureSubset")) stop("not a legitimate FeatureSubset object")
  cat("A feature subset consisting of", length(x$reduct), " attributes:\n")
  cat(paste(names(x$reduct), collapse = ", "), "\n", sep = "")

  invisible(x)
}


#' This function enables the output of a summary of the indiscernibility relation functions.
#'
#' @title The summary function for an indiscernibility relation
#' @author Lala Septem Riza
#'
#' @param object a \code{"IndiscernibilityRelation"} object. See \code{\link{BC.IND.relation.FRST}}
#'
#'        and \code{\link{BC.IND.relation.RST}}.
#' @param ... the other parameters.
#' @return a description that contains the following information. For FRST model:
#' @examples
#' ###########################################################
#' ## Example 1: Dataset containing nominal values for
#' ## all attributes.
#' ###########################################################
#' ## Decision table is represented as data frame
#' dt.ex1 <- data.frame(c(1,0,2,1,1,2,2,0), c(0, 1,0, 1,0,2,1,1),
#'                         c(2,1,0,0,2,0,1,1), c(2,1,1,2,0,1,1,0), c(0,2,1,2,1,1,2,1))
#' colnames(dt.ex1) <- c("aa", "bb", "cc", "dd", "ee")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 5)
#'
#' ## In this case, we only consider the second and third attributes.
#' attributes <- c(2, 3)
#'
#' #### calculate fuzzy indiscernibility relation ####
#' ## in this case, we are using "crisp" as a type of relation and type of aggregation
#' control.ind <- list(type.relation = c("crisp"), type.aggregation = c("crisp"))
#' IND <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.ind)
#'
#' summary(IND)
#' @export
#' @method summary IndiscernibilityRelation
summary.IndiscernibilityRelation <- function(object, ...){

 if(!inherits(object, "IndiscernibilityRelation")) stop("not a legitimate object in this package")
 cat("The name of model: ", object$type.model, "\n")
 if (object$type.model == "FRST"){
	cat("The type of aggregation: ", "\n")
	print(object$type.aggregation)
 }
 cat("The name of relation: ", "\n")
 print(object$type.relation)
 cat("The matrix of indiscernibility relation: ", "\n")
 print(object$IND.relation)

  invisible(object)
}


#' This function enables the output of a summary of the lower and upper approximations.
#'
#' @title The summary function of lower and upper approximations based on RST and FRST
#' @author Lala Septem Riza
#'
#' @param object a \code{"LowerUpperApproximation"} object. See \code{\link{BC.LU.approximation.FRST}} and \code{\link{BC.LU.approximation.RST}}.
#' @param ... the other parameters.
#' @examples
#' #######################################
#' ## Example: Using simple data set
#' #######################################
#' dt.ex1 <- data.frame(c(1,0,2,1,1,2,2,0), c(0, 1,0, 1,0,2,1,1),
#'                         c(2,1,0,0,2,0,1,1), c(2,1,1,2,0,1,1,0), c(0,2,1,2,1,1,2,1))
#' colnames(dt.ex1) <- c("aa", "bb", "cc", "dd", "ee")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 5,
#'                                      indx.nominal = c(1:5))
#'
#' P <- c(2,3)
#'
#' ####### Compute indiscernibility relation #######
#' IND <- BC.IND.relation.RST(decision.table, feature.set = P)
#'
#' ####### Compute lower and upper approximation #####
#' roughset <- BC.LU.approximation.RST(decision.table, IND)
#'
#' summary(roughset)
#' @export
#' @method summary LowerUpperApproximation
summary.LowerUpperApproximation <- function(object, ...){

 if(!inherits(object, "LowerUpperApproximation")) stop("not a legitimate object in this package")
 cat("The name of model: ", object$type.model, "\n")
 cat("The model of lower/upper approximations: ", object$type.LU, "\n")
 cat("The lower approximation: ", "\n")
 if (object$type.model == c("FRST")){
 print(object$fuzzy.lower)
 }
 else {
	for (i in 1:length(object$lower.approximation)){
		names(object$lower.approximation[[i]]) <- NULL
	}
	print(object$lower.approximation)
 }
 cat("The upper approximation: ", "\n")
 if (object$type.model == c("FRST")){
	print(object$fuzzy.upper)
 }
 else {
	for (i in 1:length(object$upper.approximation)){
		names(object$upper.approximation[[i]]) <- NULL
	}
	print(object$upper.approximation)
 }

  invisible(object)
}

#' This function enables the output of a summary of the positive region and degree of dependency.
#'
#' @title The summary function of positive region based on RST and FRST
#' @author Lala Septem Riza
#'
#' @param object a \code{"PositiveRegion"} object. See \code{\link{BC.positive.reg.FRST}} and \code{\link{BC.positive.reg.RST}}.
#' @param ... the other parameters.
#' @examples
#' dt.ex1 <- data.frame(c(1,0,2,1,1,2,2,0), c(0, 1,0, 1,0,2,1,1),
#'                         c(2,1,0,0,2,0,1,1), c(2,1,1,2,0,1,1,0), c(0,2,1,2,1,1,2,1))
#' colnames(dt.ex1) <- c("aa", "bb", "cc", "dd", "ee")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 5,
#'                                     indx.nominal = c(1:5))
#'
#' ## in this case, we consider second and third attributes only
#' P <- c(2,3)
#'
#' ####### Perform indiscernibility relation #######
#' IND <- BC.IND.relation.RST(decision.table, feature.set = P)
#'
#' ####### Perform lower and upper approximations #####
#' roughset <- BC.LU.approximation.RST(decision.table, IND)
#'
#' ####### Determine the positive region ######
#' region <- BC.positive.reg.RST(decision.table, roughset)
#'
#' summary(region)
#' @export
#' @method summary PositiveRegion
summary.PositiveRegion <- function(object, ...){
 if(!inherits(object, "PositiveRegion")) stop("not a legitimate object in this package")
 cat("The name of model: ", object$type.model, "\n")
 cat("The positive region: ", "\n")
 if (object$type.model == "FRST"){
	print(object$positive.freg)
 }
 else print(object$positive.reg)
 cat("The degree of dependency: ", "\n")
 print(object$degree.dependency)

  invisible(object)
}

# It is used to build class of rough set and fuzzy rough set theories. Currently, its implementation is very basic and
# does no argument checking, as it is only used internally.
#
# @title The object factory for RoughSets objects
# @param mod a list containing all the attributes for the object
# @param classname a class name
# @return an object of type \code{RoughSet}
# @aliases RoughSets-object
ObjectFactory <- function(mod, classname){
	class(mod) <- unique(c(classname, class(mod)))
	return(mod)
}

#' It is used to apply a particular object/model for obtaining a new decision table. In other words, in order to use the function,
#' the models, which are objects of missing value completion, feature selection, instance selection, or
#' discretization, have been calculated previously .
#'
#' @title Apply for obtaining a new decision table
#' @author Lala Septem Riza and Andrzej Janusz
#' @param decision.table a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}.
#' @param object a class resulting from feature selection (e.g., \code{\link{FS.reduct.computation}}), discretization (e.g., \code{\link{D.discretization.RST}}),
#'               instance selection functions
#'              (e.g., \code{\link{IS.FRIS.FRST}}), and missing value completion (e.g., \code{\link{MV.missingValueCompletion}}).
#' @param control a list of other parameters which are \code{indx.reduct} representing an index of the chosen decision reduct. It is only considered when
#'               we calculate all reducts using \code{\link{FS.all.reducts.computation}}. The default value is that the first reduct will be chosen.
#' @return A new decision table. Especially for the new decision table resulting from discretization, we
#'         obtain a different representation. Values are expressed in intervals instead of labels. For example,
#'         \eqn{a_1 = [-Inf, 1.35]} refers to the value \eqn{a_1} has a value in that range.
#' @examples
#' #############################################################
#' ## Example 1: The feature selection in RST
#' ## using quickreduct
#' #############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## generate reducts
#' red.1 <- FS.quickreduct.RST(decision.table)
#'
#' new.decTable <- SF.applyDecTable(decision.table, red.1)
#'
#' #############################################################
#' ## Example 2: The feature selection in FRST
#' ## using fuzzy.QR (fuzzy quickreduct)
#' #############################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## fuzzy quickreduct using fuzzy lower approximation
#' control <- list(decision.attr = c(5), t.implicator = "lukasiewicz",
#'                 type.relation = c("tolerance", "eq.1"), type.aggregation =
#'                 c("t.tnorm", "lukasiewicz"))
#' red.2 <- FS.quickreduct.FRST(decision.table, type.method = "fuzzy.dependency",
#'                             type.QR = "fuzzy.QR", control = control)
#'
#' ## generate new decision table
#' new.decTable <- SF.applyDecTable(decision.table, red.2)
#'
#' ###################################################
#' ## Example 3: The Instance selection by IS.FRPS and
#' ## generate new decision table
#' ###################################################
#' dt.ex1 <- data.frame(c(0.5, 0.2, 0.3, 0.7, 0.2, 0.2),
#'                   c(0.1, 0.4, 0.2, 0.8, 0.4, 0.4), c(0, 0, 0, 1, 1, 1))
#' colnames(dt.ex1) <- c("a1", "a2", "d")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 3)
#'
#' ## evaluate and select instances
#' res.1 <- IS.FRPS.FRST(decision.table, type.alpha = "FRPS.3")
#'
#' ## generate new decision table
#' new.decTable <- SF.applyDecTable(decision.table, res.1)
#'
#' #################################################################
#' ## Example 4: Discretization by determining cut values and
#' ## then generate new decision table
#' #################################################################
#' dt.ex2 <- data.frame(c(1, 1.2, 1.3, 1.4, 1.4, 1.6, 1.3), c(2, 0.5, 3, 1, 2, 3, 1),
#'                              c(1, 0, 0, 1, 0, 1, 1))
#' colnames(dt.ex2) <- c("a", "b", "d")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex2, decision.attr = 3,
#'                   indx.nominal = 3)
#'
#' ## get cut values using the local strategy algorithm
#' cut.values <- D.discretization.RST(decision.table, type.method = "global.discernibility")
#'
#' ## generate new decision table
#' new.decTable <- SF.applyDecTable(decision.table, cut.values)
#'
#' #################################################################
#' ## Example 5: Missing value completion
#' #################################################################
#' dt.ex1 <- data.frame(
#'      c(100.2, 102.6, NA, 99.6, 99.8, 96.4, 96.6, NA),
#'      c(NA, "yes", "no", "yes", NA, "yes", "no", "yes"),
#'      c("no", "yes", "no", "yes", "yes", "no", "yes", NA),
#'      c("yes", "yes", "no", "yes", "no", "no", "no", "yes"))
#' colnames(dt.ex1) <- c("Temp", "Headache", "Nausea", "Flu")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4,
#'                                     indx.nominal = c(2:4))
#'
#' ## missing value completion
#' val.NA = MV.missingValueCompletion(decision.table, type.method = "globalClosestFit")
#'
#' ## generate new decision table
#' new.decTable <- SF.applyDecTable(decision.table, val.NA)
#' new.decTable
#' @export
SF.applyDecTable <- function(decision.table, object, control = list()) {

  if(!inherits(decision.table, "DecisionTable")) {
    stop("Provided data should inherit from the \'DecisionTable\' class.")
  }

  if(!(inherits(object, c("FeatureSubset", "ReductSet", "InstanceSelection",
                          "Discretization", "MissingValue")))) {
    stop("Class of the object was not recognized.")
  }

	control <- setDefaultParametersIfMissing(control, list(indx.reduct = 1))

  if (inherits(object, "FeatureSubset")) {
    tmpIdx = c(object$reduct, attr(decision.table, "decision.attr"))
    new.data <- decision.table[, tmpIdx, drop = FALSE]
    attr(new.data, "nominal.attrs") = attr(decision.table, "nominal.attrs")[tmpIdx]
    attr(new.data, "desc.attrs") = attr(decision.table, "desc.attrs")[tmpIdx]
    if	(!is.null(attr(decision.table, "decision.attr"))) {
      attr(new.data, "decision.attr") = ncol(new.data)
    }
    else {
      attr(new.data, "decision.attr") = NULL
    }
    new.data = ObjectFactory(new.data, "DecisionTable")

  }
  if (inherits(object, "ReductSet")) {
    indx.reduct <- control$indx.reduct
    reducts <- object$decision.reduct
    if (is.null(reducts[indx.reduct]) || is.na(reducts[indx.reduct])) {
      stop("there is no reducts at the given indx.reduct")
    }

    tmpIdx = c(reducts[[indx.reduct]], names(attr(decision.table, "desc.attrs"))[attr(decision.table, "decision.attr")])
    new.data <- decision.table[, tmpIdx, drop = FALSE]
    attr(new.data, "nominal.attrs") = attr(decision.table, "nominal.attrs")[tmpIdx]
    attr(new.data, "desc.attrs") = attr(decision.table, "desc.attrs")[tmpIdx]
    if (!is.null(attr(decision.table, "decision.attr"))) {
      attr(new.data, "decision.attr") = ncol(new.data)
    } else {
      attr(new.data, "decision.attr") = NULL
    }
    new.data = ObjectFactory(new.data, "DecisionTable")
  }
  if (inherits(object, "InstanceSelection")) {
    indx.objects <- object$indx.objects
    if (length(indx.objects) > 0) {
      new.data <- decision.table[c(indx.objects), , drop = FALSE]
      attr(new.data, "nominal.attrs") = attr(decision.table, "nominal.attrs")
      attr(new.data, "desc.attrs") = attr(decision.table, "desc.attrs")
      attr(new.data, "decision.attr") = attr(decision.table, "decision.attr")
      new.data = ObjectFactory(new.data, "DecisionTable")
    }
  }
  if (inherits(object, "Discretization")) {
    ## sort the cut values
    cut.values <- lapply(object$cut.values, sort)

    ## get discrete values according to the cut values
    if (!is.null(attr(decision.table, "decision.attr"))) {
      if (length(cut.values) != (ncol(decision.table) - 1))
        stop("The discretization is not conforming with the decision table.")

      decision.attr = factor(decision.table[[attr(decision.table, "decision.attr")]])
      new.data = mapply(applyDiscretization,
                        decision.table[-attr(decision.table, "decision.attr")], cut.values,
                        attr(decision.table, "nominal.attrs")[-attr(decision.table, "decision.attr")],
                        SIMPLIFY = FALSE)
      new.data[[length(new.data) + 1]] = decision.attr
    } else {
      if(length(cut.values) != ncol(decision.table)) {
        stop("The discretization is not conforming with the decision table.")
      } else {
        new.data = mapply(applyDiscretization,
                          decision.table, cut.values,
                          attr(decision.table, "nominal.attrs"),
                          SIMPLIFY = FALSE)
      }
    }
    new.data = data.frame(new.data, stringsAsFactors = TRUE)
    colnames(new.data) = colnames(decision.table)

    ## generate a decision table object
    new.data <- SF.asDecisionTable(dataset = new.data,
                                   decision.attr = attr(decision.table, "decision.attr"),
                                   indx.nominal = 1:ncol(new.data))
  }
  if (inherits(object, "MissingValue")){
    new.data <- decision.table
    nominal.indx <- attr(decision.table, "nominal.attrs")
    for (i in 1:nrow(object$val.NA)){
      if (nominal.indx[object$val.NA[i, 2]] == FALSE){
        new.data[object$val.NA[i, 1], object$val.NA[i, 2]] <-  as.numeric(object$val.NA[i, 3])
      } else {
        new.data[object$val.NA[i, 1], object$val.NA[i, 2]] <- object$val.NA[i, 3]
      }
    }
    new.data <- stats::na.omit(new.data)
    new.data <- SF.asDecisionTable(dataset = new.data,
                                   decision.attr = attr(decision.table, "decision.attr"),
                                   indx.nominal = which(nominal.indx == TRUE))
  }
  return(new.data)
}


## checking missing parameters
# @param control parameter values of each method
# @param defaults default parameter values of each method
setDefaultParametersIfMissing <- function(control, defaults) {
  for(i in names(defaults)) {
    if(is.null(control[[i]])) control[[i]] <- defaults[[i]]
  }

	return(control)
}

# It is used to convert rules into string
# @param rules rules in numeric
# @param type.task a type of task
# @param nominal.att a list of types of attributes
toStr.rules <- function(rules, type.task = "classification", nominal.att = NULL, type.model = "FRST"){
	options(stringsAsFactors = FALSE)
	Str.rules <- list()
	if (type.model == "FRST"){
		for (h in 1 : length(rules)){
			rule <- rules[[h]]
			if (ncol(rule) > 1){
				ante <- paste(colnames(rule[1]), rule[1], sep = ifelse(nominal.att[1] == TRUE, c(" is "), c(" is around ")))
				for (i in 2 : (ncol(rule) - 1)){
					temp <- paste(colnames(rule[i]), rule[i], sep = ifelse(nominal.att[i] == TRUE, c(" is "), c(" is around ")))
					ante <- paste(ante, temp, sep = " and ")
				}
			}
			else {
				ante <- paste(colnames(rule[1]), rule[1], sep = ifelse(nominal.att[1] == TRUE, c(" is "), c(" is around ")))
			}

			if (type.task == "classification"){
				cons <- paste(colnames(rule[ncol(rule)]), rule[[ncol(rule)]], sep = c(" is "))
			}
			else {
				cons <- paste(colnames(rule[ncol(rule)]), rule[ncol(rule)], sep = c(" is around "))
			}

			rule <- paste("IF", ante, "THEN", cons)
			Str.rules <- append(Str.rules, rule)
		}
	}
	else {
		colNames = attr(rules, "colnames")
    for (i in 1 : length(rules)){
			ante <- paste(colNames[rules[[i]]$idx[1]], rules[[i]]$values[1], sep = " is ")
      if(length(rules[[i]]$values) > 1) {
  			for (j in 2 : length(rules[[i]]$values)){
  				temp <- paste(colNames[rules[[i]]$idx[j]], rules[[i]]$values[j], sep = " is ")
  				ante <- paste(ante, temp, sep = " and ")
  			}
      }
			cons <- paste(attr(rules, "dec.attr"), paste(rules[[i]]$consequent, ";\n\t\t(supportSize=",
                                                   length(rules[[i]]$support), "; ", "laplace=",
                                                   rules[[i]]$laplace,")", sep=""), sep = c(" is "))
			rule <- paste("IF", ante, "THEN", cons)
			Str.rules <- append(Str.rules, rule)
		}
	}
	return(Str.rules)
}

#' The function can be used to change a custom set of attribute names from
#' a decision table into an object of the FeatureSubset class. It can be useful
#' for converting results of discernibility matrix-based attribute selection
#' methods (i.e. functions FS.all.reducts.computation and FS.one.reduct.computation).
#' @title Converting custom attribute name sets into a FeatureSubset object
#' @author Andrzej Janusz
#'
#' @param colNames a character vector containing names of attributes from a decision table
#' @param decisionTable a decision table which contains attributes from colNames
#' @param type.method an indicator of the method used for selecting the attributes
#' @param model an indicator of the model used for selecting the attributes
#' @return an object of a class FeatureSubset
#' #############################################################
#' ## Example 1:
#' #############################################################
#' data(RoughSetData)
#' wine.data <- RoughSetData$wine.dt
#'
#' ## discretization and generation of a reduct
#' cut.values <- D.discretization.RST(wine.data,
#'                                    type.method = "unsupervised.quantiles",
#'                                    nOfIntervals = 3)
#' decision.table <- SF.applyDecTable(wine.data, cut.values)
#' disc.matrix <- BC.discernibility.mat.RST(wine.data)
#' reduct <- FS.one.reduct.computation(disc.matrix, greedy=TRUE)
#' class(reduct)
#' is(reduct$decision.reduct[[1]])
#'
#' ## convertion into a FeatureSubset object
#' reduct <- SF.asFeatureSubset(reduct$decision.reduct[[1]], decision.table,
#'                              type.method = "greedy reduct from a discernibility matrix",
#'                              model = reduct$type.model)
#' class(reduct)
#' @export
SF.asFeatureSubset = function(colNames, decisionTable,
                            type.method = "custom subset",
                            model = "custom") {

  if(length(colNames) == 0 | !inherits(colNames, "character")) {
    stop("No correct attribute names were provided.")
  }

  if(!inherits(decisionTable, "DecisionTable")) {
    stop("Provided data should inherit from the \'DecisionTable\' class.")
  }

  fs = list()
  fs$reduct = which(colnames(decisionTable) %in% colNames)

  if(length(fs$reduct) == 0) {
    stop("No attribute name was recognized in the provided decision table.")
	}
	else {
    if(length(fs$reduct) < length(colNames)) {
      warning("Some of the attribute names were not recognized in the provided decision table.")
    }
  }

  names(fs$reduct) = colnames(decisionTable)[fs$reduct]

  fs$type.method = type.method
  fs$type.task = "feature selection"
  fs$model = model

  class(fs) = unique(c("FeatureSubset", class(fs)))
	return(fs)
}


#############################################################################
#
#  This file is a part of the R package "frbs".
#
#  Author: Lala Septem Riza
#  Co-author: Christoph Bergmeir
#  Supervisors: Francisco Herrera Triguero and Jose Manuel Benitez
#  Copyright (c) DiCITS Lab, Sci2s group, DECSAI, University of Granada.
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
#' It is the main function used for generating the frbsPMML format. In this package, we provide interfaces for writing and reading frbsPMML to/from a text file. 
#' See \code{\link{write.frbsPMML}} and \code{\link{read.frbsPMML}}.
#'
#' frbsPMML is a universal framework for representing FRBS models, which is a format adopted from the Predictive Model Markup Language (PMML). 
#' PMML is a format constructed by an XML-based language to provide a standard for describing models produced 
#' by data mining and machine learning algorithms. A main contribution of PMML is 
#' to provide interoperable schemata of predictive models. 
#' Using PMML, we can easily perform these tasks as our models are documented
#' in an XML-based language. Human experts can also update and modify the model on the files directly.
#' 
#' Since PMML is an XML-based language, the specification is defined by an XML Schema as 
#' recommended by the World Wide Web Consortium (W3C). The PMML format is specified by the main 
#' tag \emph{PMML} that contains some components. In the following, we describe the main components:
#' \itemize{
#' \item \emph{Header}: It contains general information about the PMML document, 
#'        such as copyright information for the model, its description, application, 
#'        and timestamp of generation. 
#' \item \emph{DataDictionary}: It contains information related to fields or variables, 
#'        such as number, names, types, and value ranges of variables.
#' \item \emph{MODEL-ELEMENT}: It is a main part of the PMML document that consists of models 
#'        supported by PMML. In each model, there are several components embedded in the element, 
#'        such as \emph{MiningSchema} and \emph{Output}. 
#'        \emph{MiningSchema} specifies outlier treatment, a missing value replacement policy, 
#'        and missing value treatment, whereas \emph{Output} shows a description of the output variable. 
#'        For example, in a clustering model, we define a schema representing the cluster centers 
#'        that are included in the \emph{ClusteringModel} element.       
#' }
#' Besides these components, there are some optional elements, such as \emph{MiningBuildTask}, 
#' \emph{TransformationDictionary}, and \emph{Extension}. 
#' More detailed information about PMML can be found in (Guazzelli et al., 2009).
#' 
#' Three models, which can be used for handling regression and classification tasks, 
#' are specified by the proposed representations: Mamdani, Takagi Sugeno Kang, and fuzzy rule-based classification systems. 
#' There are the following benefits offered by frbsPMML, as follows:
#' \itemize{
#' \item Interoperability: It is a standard format for representing many models without 
#'       depending on any programming languages (e.g., Java, Python, and C++) and platforms (e.g., Windows, Linux, and Mac). 
#' \item Tranparency: Since it is formed based on XML Schema containing formal definitions of the available elements, we can 
#'       understand FRBS models as written in frbsPMML.
#' \item Interpretability: frbsPMML expresses rulebase, database, and inference schema in simple ways. For example, 
#'       rulebase is constructed recursively, so that besides it meets to the mathematical logic (predicate), we can
#'       define different operators (i.e., \code{and} and \code{or}) in one rule.  
#' \item Flexibility: Since frbsPMML is based XML, human experts can easily modify and improve a model in the text file directly. 
#' \item Reproducibility: Sicen frbsPMML is a universal representation, it allows us to store, share, execute, and reproduce an FRBS model.	
#' }
#'
#' @title The frbsPMML generator
#'
#' @param model an frbs model. 
#' @param model.name a string representing the model name.
#' @param app.name a string representing an application name.
#' @param description a string representing the simulation description.
#' @param copyright a copyright of simulation.
#' @param algorithm.name a string representing the algorithm name.
#' @param ... other parameters
#' @return FRBS model in frbsPMML format
#' @examples
#' ## This example shows how to construct a frbsPMML file of the frbs model
#' ## Even though we are using MAMDANI model, other models have the same way
#' ## 
#' ## 1. Produce frbs model, for example: we perform Wang & Mendel's technique (WM)
#' ## Input data
#' \dontrun{data(frbsData)
#' data.train <- frbsData$GasFurnance.dt[1 : 204, ]
#' data.fit <- data.train[, 1 : 2]
#' data.tst <- frbsData$GasFurnance.dt[205 : 292, 1 : 2]
#' real.val <- matrix(frbsData$GasFurnance.dt[205 : 292, 3], ncol = 1)
#' range.data <- matrix(c(-2.716, 2.834, 45.6, 60.5, 45.6, 60.5), nrow = 2)
#' 
#' ## Set the method and its parameters
#' method.type <- "WM"
#' control <- list(num.labels = 3, type.mf = "GAUSSIAN", type.defuz = "WAM", 
#'                 type.tnorm = "MIN", type.implication.func = "ZADEH", 
#'                 name = "sim-0") 
#' 
#' ## Generate fuzzy model
#' object <- frbs.learn(data.train, range.data, method.type, control)
#' 
#' ## 2. Write frbsPMML file
#' ## by calling frbsPMML(), the frbsPMML format will be displayed in R console
#' frbsPMML(object)}
#'
#' @references
#' A. Guazzelli, M. Zeller, W.C. Lin, and G. Williams., 
#' "pmml: An open standard for sharing models", The R Journal, Vol. 1, No. 1, pp. 60-65 (2009).
#'
#' Data Mining Group, http://www.dmg.org/.
#'  
#' @export  
frbsPMML <- function(model,
                        model.name="frbs_model",
                        app.name="frbs",
                        description= NULL, 
                        copyright=NULL,
                        algorithm.name=model$method.type,
                        ...)
{
	requireNamespace("XML", quietly = TRUE)
	
	if (! inherits(model, "frbs")) stop("Not a legitimate frbs object")
	
	###################################################
	# Collect the required information.
	################################################### 
	range.data.ori <- model$range.data.ori  
	shadow.data <- as.data.frame(range.data.ori)
	
	## regression tasks
	if (any(model$type.model == c("MAMDANI", "TSK"))){ 	
		colnames(shadow.data) <- model$colnames.var
		funcName = "regression"
		optype="continuous"
		
		## Mamdani model
		if (model$type.model == "MAMDANI"){
			varinp.mf <- model$varinp.mf
			varout.mf <- model$varout.mf
			rule <- model$rule
			rule.data.num <- model$rule.data.num
			number.of.rules <- nrow(rule)
			num.fvalinput <- model$num.labels[, -ncol(model$num.labels), drop = FALSE]
			number.of.fields <- ncol(model$num.labels)
			num.varinput <- number.of.fields - 1
		}
		
		## Takagi Sugeno Kang model
		else {
			varinp.mf <- model$varinp.mf
			rule <- model$rule
			rule.data.num <- model$rule.data.num
			number.of.rules <- nrow(rule)
			func.tsk <- model$func.tsk
			num.fvalinput <- model$num.labels
			number.of.fields <- ncol(model$num.labels) + 1
			num.varinput <- ncol(model$num.labels)
		}
	}
	
	## classification tasks
	else if (model$type.model == "FRBCS") {		
		colnames(shadow.data) <- model$colnames.var[-length(model$colnames.var)]
		funcName = "classification"
		optype="string"
		varinp.mf <- model$varinp.mf
		rule <- model$rule
		rule.data.num <- model$rule.data.num
		number.of.rules <- nrow(rule)
		class <- model$class
		num.fvalinput <- model$num.labels[, -ncol(model$num.labels), drop = FALSE]
		number.of.fields <- ncol(model$num.labels)
		num.varinput <- number.of.fields - 1

	} else {
		stop("Currently, this model has not been supported yet")
	}
		
	if (algorithm.name == c("FS.HGD")){
			alpha.heuristic <- model$alpha.heuristic
	}
			
    field <- NULL
    field$name <-  model$colnames.var
  
    ## get number of attributes/variables
    if (funcName == "regression"){
		field$class <- rep("numeric", number.of.fields) # All fields are numeric
	}
	else {
		field$class <- c(rep("numeric", (number.of.fields - 1)), "factor")
	}
	names(field$class) <- field$name
  
	## get the name of output variable
	target <- model$colnames.var[length(model$colnames.var)]
  
	## it is used for classification tasks only
	for (i in 1:number.of.fields)
	{
		if (field$class[[field$name[i]]] == "factor") {
			field$levels[[field$name[i]]] <- c(as.character(seq(1, model$num.labels[1, ncol(model$num.labels)])))
		}
	}

	## Start constructing the frbsPMML format
	pmml <- .pmmlRootNode("1.0")
  
	## PMML -> Header
	pmml <- XML::append.XMLNode(pmml, .pmmlHeader(description, copyright, app.name))
   
	## PMML -> DataDictionary
	data.dictionary <- .pmmlDataDictionary(field = field, dataset = shadow.data)
	pmml <- XML::append.XMLNode(pmml, data.dictionary)
  
	## Construct a part of the frbs model 
	the.model <- XML::xmlNode("FrbsModel",
                      attrs=c(modelName=model$type.model,
                        functionName=funcName, 
                        algorithmName=algorithm.name,
                        targetFieldName=target
                        )) 

	## PMML -> frbs -> MiningSchema
	the.model <- XML::append.XMLNode(the.model, .pmmlMiningSchema(field = field, target = target))
  
	## PMML -> frbs -> Output
	if (funcName == "regression") {
		the.model <- XML::append.XMLNode(the.model, .pmmlOutput(field, target, optype))
	}
	else {
		the.model <- XML::append.XMLNode(the.model, .pmmlOutput(field, target))
	}
  
	## PMML -> frbs -> InferenceSchema&CLusterSchema
	the.model <- XML::append.XMLNode(the.model, .frbsInference(model))
  
	## For TSK Model
	if (model$type.model == "TSK") {	 
		## PMML -> FrbsModel
		the.model <- XML::append.XMLNode(the.model, .buildDatabaseTSK(model))	
		the.model <- XML::append.XMLNode(the.model, .buildRulesTSK(model))
	}
  
	## For Mamdani Model
	else if (model$type.model == "MAMDANI") {		
		# PMML -> FrbsModel
		the.model <- XML::append.XMLNode(the.model, .buildDatabaseMamdani(model))
		the.model <- XML::append.XMLNode(the.model, .buildRulesMamdani(model))
	}

	## For Classification Model
	else if (model$type.model == "FRBCS") { 
		# PMML -> FrbsModel
		the.model <- XML::append.XMLNode(the.model, .buildDatabaseFRBCS(model))
		the.model <- XML::append.XMLNode(the.model, .buildRulesFRBCS(model))
	}
  
	pmml <- XML::append.XMLNode(pmml, the.model)
  
	return(pmml)
}

# This is a internal function for constructing database based on Mamdani in frbsPMML
#  
# @title constructor of Mamdani Database 
#
# @param model an frbs object.
# @export
.buildDatabaseMamdani <- function(model)
{
	PMMLdata <- XML::xmlNode("Database")
	varinp.mf <- model$varinp.mf
	varout.mf <- model$varout.mf
	num.labels <- model$num.labels
	colnames.var <- model$colnames.var
	j <- 1
	
	## construct for varinp.mf into pmml
	for (h in 1 : (ncol(num.labels) - 1)){
		MF <- XML::xmlNode("MembershipFunction", attrs=c(name=colnames.var[h], numberOfLabels = num.labels[1, h]))	
		for (i in j:(j + num.labels[1, h] - 1)){
			shape.MF <- .nameMF(varinp.mf)					
			MF <- XML::append.XMLNode(MF,
						XML::xmlNode("FuzzyTerm",
						   attrs=c(name=colnames(varinp.mf)[[i]], 
							   type=shape.MF[i]),
							   .constructMF(varinp.mf[, i])))
		}
		j <- i + 1
		PMMLdata <- XML::append.XMLNode(PMMLdata, MF)
	}
	
	MF <- XML::xmlNode("MembershipFunction", attrs=c(name=colnames.var[length(colnames.var)], numberOfLabels = num.labels[1, ncol(num.labels)]))
	for (i in 1:ncol(varout.mf))
	{
		shape.MF <- .nameMF(varout.mf)	
		MF <- XML::append.XMLNode(MF,
						XML::xmlNode("FuzzyTerm",
						   attrs=c(name=colnames(varout.mf)[[i]], 
								   type=shape.MF[i]),	           
									.constructMF(varout.mf[, i])))
	}
	PMMLdata <- XML::append.XMLNode(PMMLdata, MF)
	return(PMMLdata)
}

# This is a internal function for constructing database based on Takagi Sugeno Kang in frbsPMML
#  
# @title Constructor of TSK Database
#
# @param model an frbs object.
# @export
.buildDatabaseTSK <- function(model)
{
	PMMLdata <- XML::xmlNode("Database")
	varinp.mf <- model$varinp.mf
	num.labels <- model$num.labels
	colnames.var <- model$colnames.var
	j <- 1
	for (h in 1 : ncol(num.labels)){
		MF <- XML::xmlNode("MembershipFunction", attrs=c(name=colnames.var[h], numberOfLabels = num.labels[1, h]))	
		for (i in j:(j + num.labels[1, h] - 1)){
			shape.MF <- .nameMF(varinp.mf)					
			MF <- XML::append.XMLNode(MF,
						XML::xmlNode("FuzzyTerm",
						   attrs=c(name=colnames(varinp.mf)[[i]], 
							   type=shape.MF[i]),				           
								.constructMF(varinp.mf[, i])))
		}
		j <- i + 1
		PMMLdata <- XML::append.XMLNode(PMMLdata, MF)
	}

	return(PMMLdata)
}

# This is a internal function for constructing database based on fuzzy rule-based classification system in frbsPMML
#  
# @title  constructor of FRBCS Database
#
# @param model an frbs object.
# @export
.buildDatabaseFRBCS <- function(model)
{
	PMMLdata <- XML::xmlNode("Database")
	varinp.mf <- model$varinp.mf
	num.labels <- model$num.labels
	colnames.var <- model$colnames.var
	j <- 1
	for (h in 1 : (ncol(num.labels) - 1)){
		MF <- XML::xmlNode("MembershipFunction", attrs=c(name=colnames.var[h], numberOfLabels = num.labels[1, h]))	
		for (i in j:(j + num.labels[1, h] - 1)){
			shape.MF <- .nameMF(varinp.mf)					
			MF <- XML::append.XMLNode(MF,
						XML::xmlNode("FuzzyTerm",
						   attrs=c(name=colnames(varinp.mf)[[i]], 
							   type=shape.MF[i]),				           
							   .constructMF(varinp.mf[, i])))
		}
		j <- i + 1
		PMMLdata <- XML::append.XMLNode(PMMLdata, MF)
	}
	
	return(PMMLdata)
}

# This is a internal function for constructing fuzzy rules based on Mamdani in frbsPMML
#  
# @title Rules constructor
#
# @param model an frbs object.
# @export
.buildRulesMamdani <- function(model)
{
	rules <- model$rule
	PMMLrules <- XML::xmlNode("Rulebase",
					  attrs=c(numberOfRules=nrow(rules))) 				  
	#seqq <- seq(2, (ncol(rules) - 4), by = 8) 
	for (i in 1:nrow(rules)){ 	
		rule <- XML::xmlNode("Rule", attrs=c(id=i))	
		
		dt.rule <- rules[i, 2:(ncol(rules)-4), drop = FALSE]
		rule <- XML::append.XMLNode(rule, .genIFPart(dt.rule))
	
		thenPart <- XML::xmlNode("Then")
		simplePredicate <- XML::xmlNode("SimplePredicate", attrs = c(field = rules[i, ncol(rules) - 2], value = rules[i, ncol(rules)]))
		thenPart <- XML::append.XMLNode(thenPart, simplePredicate)			
		rule <- XML::append.XMLNode(rule, thenPart)
		
		if (model$type.model == "FRBCS"){
			if (is.null(model$grade.cert[i])){
				gradePart <- XML::xmlNode("Grade", value = 1)
			}
			else {
				gradePart <- XML::xmlNode("Grade", value = model$grade.cert[i])			
			}
			rule <- XML::append.XMLNode(rule, gradePart)
		}
		PMMLrules <- XML::append.XMLNode(PMMLrules, rule)
	}	
	return(PMMLrules)
}

# This is a internal function for constructing fuzzy rules based on Takagi Sugeno Kang in frbsPMML
#  
# @title Rules constructor
#
# @param model an frbs object.
# @export
.buildRulesTSK <- function(model)
{
	rules <- model$rule
	colnames.var <- model$colnames.var
	func.tsk <- model$func.tsk
	PMMLrules <- XML::xmlNode("Rulebase",
					  attrs=c(numberOfRules=nrow(rules))) 				  
	#seqq <- seq(2, (ncol(rules) - 1), by = 4) 
	for (i in 1:nrow(rules)){ 	
		rule <- XML::xmlNode("Rule", attrs=c(id=i))	
		
		dt.rule <- rules[i, 2:(ncol(rules)-1),drop = FALSE]
		rule <- XML::append.XMLNode(rule, .genIFPart(dt.rule))

		## get THEN part
		thenPart <- XML::xmlNode("Then", attrs =c(type="LinearFunction"))		
		
		for (ii in 1:length(func.tsk[i, ])){
			if (ii == length(func.tsk[i, ])){
				funcLinear <- XML::xmlNode("Constant", attrs =c(value=func.tsk[i, ii]))
			}
			else {
				funcLinear <- XML::xmlNode("Coefficient", attrs=c(field=colnames.var[ii], value=func.tsk[i, ii]))
			}
			thenPart <- XML::append.XMLNode(thenPart, funcLinear)
		}
				
		rule <- XML::append.XMLNode(rule, thenPart)	
		PMMLrules <- XML::append.XMLNode(PMMLrules, rule)
	}	
		
	return(PMMLrules)
}

# This is a internal function for constructing fuzzy rules based on fuzzy rule-based classification systems (FRBCS) in frbsPMML
#  
# @title Rules constructor
#
# @param model an frbs object.
# @export
.buildRulesFRBCS <- function(model)
{
	return(.buildRulesMamdani(model))
}

# This is a internal function for writing InferenceSchema in frbsPMML
#  
# @title InferenceSchema writer
#
# @param model an frbs object.
# @export
.frbsInference <- function(model)
{ 
  InfSchema <- XML::xmlNode("InferenceSchema")
  InfSchema.fields <- list()
	
	if (!is.null(model$type.tnorm))								
		InfSchema.fields[[1]] <- XML::xmlNode("ConjunctionOperator", attrs=c(value=model$type.tnorm))
		
	if (!is.null(model$type.snorm))								
		InfSchema.fields[[2]] <- XML::xmlNode("DisjunctionOperator", attrs=c(value=model$type.snorm))								
	if (!is.null(model$type.implication.func))
		InfSchema.fields[[3]] <- XML::xmlNode("ImplicationOperator", attrs=c(value=model$type.implication.func))
	if (!is.null(model$type.defuz))
							
		InfSchema.fields[[4]] <- XML::xmlNode("AggregationOperator", attrs=c(value=model$type.defuz))								
  InfSchema$children <- InfSchema.fields
  return(InfSchema)
}

# This is a internal function for writing names of shape of membership functions in frbsPMML
#  
# @title InferenceSchema writer
#
# @param model an frbs object.
# @export
.nameMF <- function(var.mf)
{
	shape.MF <- list()
	shape.MF.str <- c("TRIANGLE", "TRAPEZOID", "TRAPEZOID", "TRAPEZOID",
					 "GAUSSIAN", "SIGMOID", "BELL")
	for (i in 1 : ncol(var.mf)){
		shape.MF[i] <- shape.MF.str[var.mf[1, i]]
	}	
	return(shape.MF)
}


# It is used to construct parameters in the frbsPMML format
# @param var.MF a vector of a parameter of membership functions
# @param shape.MF a string representing a name of shape of membership functions
.constructMF <- function(var.MF){
	mf <- XML::xmlNode("Parameters")
	if (var.MF[1] == 1){
		name.param <- c("Left", "Middle", "Right")
		for (i in 2 : 4){
			temp <- XML::xmlNode(name.param[(i - 1)], var.MF[i])
			mf <- XML::append.XMLNode(mf, temp)
		}
	}
	if (var.MF[1] == 4){
		name.param <- c("Left", "LeftMiddle", "RightMiddle", "Right")
		for (i in 2 : 5){
			temp <- XML::xmlNode(name.param[(i - 1)], var.MF[i])
			mf <- XML::append.XMLNode(mf, temp)
		}
	}
	if (var.MF[1] == 5){
		name.param <- c("Mean", "Variance")
		for (i in 2 : 3){
			temp <- XML::xmlNode(name.param[(i - 1)], var.MF[i])
			mf <- XML::append.XMLNode(mf, temp)
		}
	}
	if (var.MF[1] == 2){
		name.param <- c("Left", "LeftMiddle", "RightMiddle", "Right")
		for (i in 2 : 5){
			if (i == 2){
				temp <- XML::xmlNode(name.param[(i - 1)], var.MF[i])
			}
			else {
				temp <- XML::xmlNode(name.param[(i - 1)], var.MF[i - 1])
			}
			mf <- XML::append.XMLNode(mf, temp)
		}
	}
	if (var.MF[1] == 3){
		name.param <- c("Left", "LeftMiddle", "RightMiddle", "Right")
		for (i in 2 : 5){
			if (i == 5){
				temp <- XML::xmlNode(name.param[(i - 1)], var.MF[i - 1])
			}
			else {
				temp <- XML::xmlNode(name.param[(i - 1)], var.MF[i])
			}
			mf <- XML::append.XMLNode(mf, temp)
		}
	}
	return(mf)
}

# It is used to build the IF part of rules
# @param rule a matrix of rule
.genIFPart <- function(rule){
	## get IF part
	ifPart <- XML::xmlNode("If")	
	if (ncol(rule) == 3) {
		simplePredicate <- XML::xmlNode("SimplePredicate", attrs = c(field =  rule[1, 1], value = rule[1, 3]))
		ifPart <- XML::append.XMLNode(ifPart, simplePredicate)
	}
	else {
		ifPart <- XML::append.XMLNode(ifPart, .genCompoundPredicate(rule[1, ,drop =FALSE]))
	}
	return(ifPart)
}

# It is used to construct a compound of predicate
# @param rule a matrix of rule
.genCompoundPredicate <- function(rule){
	if (length(rule) <= 3) {
		return(XML::xmlNode("SimplePredicate", attrs = c(field =  rule[1, 1], value = rule[1, 3])))
	}
	else {
		compoundPredicate <- XML::xmlNode("CompoundPredicate", attrs = c(booleanOperator = rule[1, 4]))
		simplePredicate1 <- XML::xmlNode("SimplePredicate", attrs = c(field = rule[1, 1], value = rule[1, 3]))
		return(XML::append.XMLNode(compoundPredicate, simplePredicate1, .genCompoundPredicate(rule[1, -c(1:4), drop = FALSE])))
	}
}
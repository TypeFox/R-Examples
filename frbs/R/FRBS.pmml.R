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
#' It is a function used to save an FRBS model to the .frbsPMML file. Detailed information about frbsPMML can be seen in \code{\link{frbsPMML}}.
#'  
#' @title The frbsPMML writer 
#'
#' @param object a frbsPMML object which is an object produced by \code{\link{frbsPMML}}.
#' @param fileName a file name with extension \code{.frbsPMML}.
#' @seealso \code{\link{read.frbsPMML}} and \code{\link{frbsPMML}}.
#' @return a file containing an FRBS model in frbsPMML format 
#' @examples
#' ## This example shows how to construct frbsPMML file of frbs model
#' ## Even though we are using MAMDANI model, other models have the same way
#' ## 
#' ## 1. Produce frbs model, for example: we perform Wang & Mendel's technique (WM)
#' ##
#' ## Input data
#' data(frbsData)
#' data.train <- frbsData$GasFurnance.dt[1 : 204, ]
#' data.fit <- data.train[, 1 : 2]
#' data.tst <- frbsData$GasFurnance.dt[205 : 292, 1 : 2]
#' real.val <- matrix(frbsData$GasFurnance.dt[205 : 292, 3], ncol = 1)
#' range.data<-matrix(c(-2.716, 2.834, 45.6, 60.5, 45.6, 60.5), nrow = 2)
#' 
#' ## Set the method and its parameters
#' method.type <- "WM"
#' control <- list(num.labels = 15, type.mf = "GAUSSIAN", type.defuz = "WAM", 
#'                 type.tnorm = "MIN", type.implication.func = "ZADEH", 
#'                 name = "sim-0") 
#' 
#' ## Generate fuzzy model
#' \dontrun{object <- frbs.learn(data.train, range.data, method.type, control)}
#' 
#' ## 2. Write frbsPMML file
#' ## In this step, we provide two steps as follows:
#' ## a. by calling frbsPMML() function directly. 
#' ## b. by calling write.frbsPMML() function. 
#' 
#' ## 2a. by calling frbsPMML(), the frbsPMML format will be displayed in R console
#' \dontrun{pmml.obj <- frbsPMML(object)}
#' 
#' ## 2b. by calling write.frbsPMML(), the result will be saved as a file
#' ##     in the working directory.
#' \dontrun{write.frbsPMML(pmml.obj, file = "MAMDANI.GasFur")}
#'
#' @references
#' A. Guazzelli, M. Zeller, W.C. Lin, and G. Williams., 
#' "pmml: An open standard for sharing models", The R Journal, Vol. 1, No. 1, pp. 60-65 (2009).
#'
#' Data Mining Group, http://www.dmg.org/.
#'  
#' @export
write.frbsPMML <- function(object, fileName = NULL) {	
	requireNamespace("XML", quietly = TRUE)
			
	if (is.null(fileName))
		fileName = paste(object$name, "frbsPMML", sep =".")
	else fileName = paste(fileName, "frbsPMML", sep = ".")
    XML::saveXML(object, file=fileName, prefix = NULL)
}

#' It is used to read the frbsPMML format into an frbs model in R. Detailed information about frbsPMML can be seen in \code{\link{frbsPMML}}.
#'  
#' @title The frbsPMML reader
#'
#' @param fileName a file name with extension \code{.frbsPMML}.
#' @return an object representing the frbs model.
#' @seealso \code{\link{write.frbsPMML}} and \code{\link{frbsPMML}}.
#' @return an frbs object
#' @examples
#' ## This example shows how to construct and read frbsPMML file of frbs model
#' ## Even though we are using MAMDANI model, other models have the same way
#' ## 
#' ## 1. Produce frbs model, for example: we perform Wang & Mendel's technique (WM)
#' ##
#' ## Input data
#' data(frbsData)
#' data.train <- frbsData$GasFurnance.dt[1 : 204, ]
#' data.fit <- data.train[, 1 : 2]
#' data.tst <- frbsData$GasFurnance.dt[205 : 292, 1 : 2]
#' real.val <- matrix(frbsData$GasFurnance.dt[205 : 292, 3], ncol = 1)
#' range.data<-matrix(c(-2.716, 2.834, 45.6, 60.5, 45.6, 60.5), nrow = 2)
#' 
#' ## Set the method and its parameters
#' method.type <- "WM"
#' control <- list(num.labels = 15, type.mf = "GAUSSIAN", type.defuz = "WAM", 
#'                 type.tnorm = "MIN", type.implication.func = "ZADEH", 
#'                 name="sim-0") 
#' 
#' ## Generate fuzzy model
#' \dontrun{object <- frbs.learn(data.train, range.data, method.type, control)}
#' 
#' ## 2. Write frbsPMML file
#' ## In this step, we provide two ways as follows.
#' ## a. by calling frbsPMML() function directly. 
#' ## b. by calling write.frbsPMML() function. 
#' 
#' ## 2a. by calling frbsPMML(), the format will be displayed in R console
#' \dontrun{frbsPMML(object)}
#' 
#' ## 2b. by calling write.frbsPMML(), the result will be saved as a file
#' ##     in the working directory.
#' \dontrun{write.frbsPMML(object, file = "MAMDANI.GasFur")}
#' 
#' ## 3. Read frbsPMML file
#' \dontrun{object <- read.frbsPMML("MAMDANI.GasFur.frbsPMML")}
#' 
#' ## 4. Perform predicting step
#' \dontrun{res.test <- predict(object, data.tst)}
#'
#' @export
read.frbsPMML <- function(fileName) {
	requireNamespace("XML", quietly = TRUE)
	
    doc <- XML::xmlRoot(XML::xmlTreeParse(fileName))
    ## check model type
    if(is.element("FrbsModel", names(doc))){
		return(read.pmml.frbs(doc))
	} else {
		stop("File does not contain an frbs Model.")
	}
}

# This is a internal frbs frbsPMML reader
#  
# @title frbs frbsPMML reader
#
# @param model an object.
# @return frbs object
# @export
read.pmml.frbs <- function(model) {
	requireNamespace("XML", quietly = TRUE)
	
	## Reading frbs tag
    frbs.model <- model[["FrbsModel"]]
	if (is.null(frbs.model)){
		stop("It is not the frbs model, please make sure the FrbsModel tag is available.") 
	}
	
	attrs.model <- XML::xmlAttrs(frbs.model)
	type.model <- attrs.model[[1]]
		
	if (type.model == "MAMDANI") {
		object <- read.Mamdani.pmml(model)
	}
	else if (type.model == "TSK"){
		object <- read.TSK.pmml(model)
	}
	else if (type.model == "FRBCS"){
		object <- read.FRBCS.pmml(model)
	}
	else {
		stop("The model has not been supported yet.")
	}
		
	return(object)
}

# This is a internal function for reading Mamdani model in frbsPMML
#  
# @title Mamdani reader
#
# @param model an frbs object.
# @return frbs object
# @export
read.Mamdani.pmml <- function(model){
	type.model <- "MAMDANI"
	## Reading the header
    header <- read.HEADER(model)
	datadictionary <- read.DATADICTIONARY(model, type.model)
	mining <- read.MININGSCHEMA(model)
	inference <- read.INFERENCESCHEMA(model)
	database <- read.DATABASE(model)
	num.labels = database$num.labels
	varinp.mf <- database$var.MF[, 1:sum(num.labels[1, -c(ncol(num.labels))])]
	varout.mf <- database$var.MF[, (ncol(varinp.mf) + 1): ncol(database$var.MF)]
	rules <- read.RULEBASE(model, type.model = "MAMDANI")		
	
	pmml <- TRUE
	
	if (ncol(varinp.mf) != sum(num.labels[, -ncol(num.labels)]))
		stop("The number of labels is not the same as the number of membership functions.")
	
	if (ncol(varout.mf) != sum(num.labels[, ncol(num.labels)]))
		stop("The number of labels is not the same as the number of membership functions.")
	
	object <- list(name = header$name, method.type = toupper(header$method.type), range.data.ori = datadictionary$range.data.ori, num.labels = num.labels, 
				   varinp.mf = varinp.mf, varout.mf = varout.mf, rule = rules,  
				   type.defuz = toupper(inference$type.defuz), type.model = type.model,
				   type.tnorm = toupper(inference$type.tnorm), type.snorm = toupper(inference$type.snorm),
				   type.implication.func = toupper(inference$type.implication.func), colnames.var = mining$colnames.var,
				   pmml = pmml)
	
	mod <- frbsObjectFactory(object)
	return (mod)
}

# This is a internal function for reading Takagi Sugeno Kang model in frbsPMML
#  
# @title Takagi Sugeno Kang reader
#
# @param model an frbs object.
# @return frbs object
# @export
read.TSK.pmml <- function(model){
	type.model <- "TSK"
	## Reading the header
    header <- read.HEADER(model)
	datadictionary <- read.DATADICTIONARY(model, type.model)
	mining <- read.MININGSCHEMA(model)
	inference <- read.INFERENCESCHEMA(model)
	database <- read.DATABASE(model)
	num.labels = database$num.labels
	varinp.mf <- database$var.MF
	rules <- read.RULEBASE(model, type.model = "TSK")
	
	pmml <- TRUE
		
	if (ncol(varinp.mf) != sum(num.labels))
		stop("The number of labels is not the same as the number of membership functions.")

	object <- list(name = header$name, method.type = toupper(header$method.type), range.data.ori = datadictionary$range.data.ori, num.labels = num.labels, 
				   varinp.mf = varinp.mf, rule = rules$rules, func.tsk = rules$func.tsk, type.model = type.model, 
				   type.tnorm = toupper(inference$type.tnorm), type.snorm = toupper(inference$type.snorm),
				   type.implication.func = toupper(inference$type.implication.func), colnames.var = mining$colnames.var,
				   pmml = pmml)
	
	mod <- frbsObjectFactory(object)
	return (mod)
}


# This is a internal function for reading Takagi Sugeno Kang model in frbsPMML
#  
# @title Takagi Sugeno Kang reader
#
# @param model an frbs object.
# @return frbs object
# @export
read.FRBCS.pmml <- function(model){
	type.model <- "FRBCS"
	## Reading the header
    header <- read.HEADER(model)
	datadictionary <- read.DATADICTIONARY(model, type.model)
	mining <- read.MININGSCHEMA(model)
	inference <- read.INFERENCESCHEMA(model)
	database <- read.DATABASE(model)
	num.labels = database$num.labels
	varinp.mf <- database$var.MF
	rules <- read.RULEBASE(model, type.model = "FRBCS")	
	num.class <- max(datadictionary$num.class)
	num.labels <- cbind(num.labels, num.class)

	pmml <- TRUE
	
	if (ncol(varinp.mf) != sum(num.labels[, -ncol(num.labels)]))
		stop("The number of labels is not the same as the number of membership functions.")

	object <- list(name = header$name, method.type = toupper(header$method.type), range.data.ori = datadictionary$range.data.ori, 
				   num.labels = num.labels, varinp.mf = varinp.mf, rule = rules$rules, 
				   type.model = type.model, grade.cert = rules$grade.cert, class = rules$class,
				   type.tnorm = toupper(inference$type.tnorm), type.snorm = toupper(inference$type.snorm), 
				   type.implication.func = toupper(inference$type.implication.func), colnames.var = mining$colnames.var,
				   pmml = pmml)
	
	mod <- frbsObjectFactory(object)
	
	return(mod)
}

# This is a internal function for reading header of frbsPMML
#  
# @title header reader
#
# @param model an frbs object.
# @export
read.HEADER <- function(model){
	requireNamespace("XML", quietly = TRUE)
	
	## Reading the header
    frbs.model <- model[["FrbsModel"]]
	attrs.model <- XML::xmlAttrs(frbs.model)
	method.type <- attrs.model[[3]]	
	header <-  XML::xmlAttrs(model[["Header"]])
	name <- strsplit(header[2], " ")$description[2]
	if (is.null(name)) 
		name <- c("sim-0")
	
	if (is.null(method.type))
		method.type <- c("NA")
		
	return(list(method.type = method.type, name = name))
}

# This is a internal function for reading DataDictionary of frbsPMML
#  
# @title DataDictionary reader
#
# @param model an frbs object.
# @param type.model a type of frbs model
# @export
read.DATADICTIONARY <- function(model, type.model = "MAMDANI"){
	requireNamespace("XML", quietly = TRUE)
	
	## Reading the header
    frbs.model <- model[["FrbsModel"]]
	
	## Reading DataDictionary tag
	DataDict.model <-  model[names(model) == "DataDictionary"]
	if (is.null(DataDict.model))
		stop("DataDictionary is not available, please define it.")
	
	## Data Dictionary
	if (any(type.model == c("MAMDANI", "TSK"))) {
		interval.data <- lapply(DataDict.model, FUN = function(x) sapply(XML::xmlChildren(x), XML::xmlChildren))   
		num.var <- XML::xmlSize(interval.data$DataDictionary)	
		num.var.valid <- as.numeric(XML::xmlAttrs(DataDict.model$DataDictionary)[1])
		if (num.var != num.var.valid)
			stop("please check the number of variables.")
			
		range.data.ori <- matrix(nrow = 2, ncol = num.var)
		for (i in 1 : num.var) {
			interval.i <- XML::xmlAttrs(interval.data$DataDictionary[i]$DataField.Interval)		
			range.data.ori[1, i] <- as.numeric(interval.i[2])
			range.data.ori[2, i] <- as.numeric(interval.i[3])
		}		
		return(list(range.data.ori = range.data.ori))
	}
	else if (type.model == "FRBCS"){
		interval.data <- lapply(DataDict.model, FUN = function(x) sapply(XML::xmlChildren(x), XML::xmlChildren))   
		num.var <- XML::xmlSize(interval.data$DataDictionary)
		num.var.valid <- as.numeric(XML::xmlAttrs(DataDict.model$DataDictionary)[1])
		if (num.var != num.var.valid)
			stop("please check the number of variables.")
			
		range.data.ori <- matrix(nrow = 2, ncol = (num.var - 1))
		for (i in 1 : num.var) {
			optype <- XML::xmlAttrs(XML::xmlChildren(DataDict.model$DataDictionary)[i]$DataField)[[2]]
			if (optype == "categorical"){
				num.class <- XML::xmlSize(XML::xmlChildren(XML::xmlChildren(DataDict.model$DataDictionary)[i]$DataField))
			} else {
				interval.i <- XML::xmlAttrs(interval.data$DataDictionary[i]$DataField$Interval)		
				range.data.ori[1, i] <- as.numeric(interval.i[2])
				range.data.ori[2, i] <- as.numeric(interval.i[3])
			}
		}	
		return(list(range.data.ori = range.data.ori, num.class = num.class))
	} else {
		stop("The model has not been supported yet.")
	}
}

# This is a internal function for reading MiningSchema of frbsPMML
#  
# @title MiningSchema reader
#
# @param model an frbs object.
# @export
read.MININGSCHEMA <- function(model){
	requireNamespace("XML", quietly = TRUE)
	
	## Reading the header
    frbs.model <- model[["FrbsModel"]]
	
	## Reading MiningSchema tag
	MiningSchema.model <-  frbs.model[names(frbs.model) == "MiningSchema"]
	
	## Mining Schema/colnames.var
	MiningSchema <- XML::xmlChildren(MiningSchema.model$MiningSchema)
    names.var <- sapply(MiningSchema, FUN = function(x) XML::xmlAttrs(x))
	colnames.var <- names.var[1, ]
	
	return(list(colnames.var = colnames.var))
}

# This is a internal function for reading InferenceSchema of frbsPMML
#  
# @title InferenceSchema reader
#
# @param model an frbs object.
# @export
read.INFERENCESCHEMA <- function(model){
	requireNamespace("XML", quietly = TRUE)
	
	## Reading the header
    frbs.model <- model[["FrbsModel"]]
	
	## Reading InferenceSchema tag
	InferenceSchema.model <- frbs.model[names(frbs.model) == "InferenceSchema"]
	if (is.null(InferenceSchema.model))
		warning("InferenceSchema is not available.")
	
	## Get parameters and their values in Inference Schema
	InferenceSchema <- XML::xmlChildren(InferenceSchema.model$InferenceSchema)
	Inf.params <- sapply(InferenceSchema, FUN = function(x) XML::xmlAttrs(x))
	num.params <- length(Inf.params)
	
	type.implication.func <- NULL
	type.defuz <- NULL
	
	## Get type of model for checking
	for (i in 1 : num.params){
		if (names(Inf.params[i]) == "ConjunctionOperator.value")
			type.tnorm <- Inf.params[i]
		else if (names(Inf.params[i]) == "DisjunctionOperator.value")
			type.snorm <- Inf.params[i]
		else if (names(Inf.params[i]) == "ImplicationOperator.value")
			type.implication.func <- Inf.params[i]
		else if (names(Inf.params[i]) == "AggregationOperator.value")
			type.defuz <- Inf.params[i]
	}
	
	if (is.null(type.tnorm))
		type.tnorm <- c("MIN")
	if (is.null(type.snorm))
		type.snorm <- c("MAX")
	if (is.null(type.implication.func))
		type.implication.func <- c("ZADEH")
	if (is.null(type.defuz))
		type.defuz <- c("WAM")
		
	return(list(type.tnorm = type.tnorm, type.snorm = type.snorm,
	            type.implication.func = type.implication.func, type.defuz = type.defuz))
}

read.RULEBASE <- function(model, type.model = "MAMDANI"){
	requireNamespace("XML", quietly = TRUE)
	
	## Reading the header
    frbs.model <- model[["FrbsModel"]]

	Rulebase <- frbs.model[names(frbs.model) == "Rulebase"]$Rulebase
	
	## get Database tag
	Rules.model <- XML::xmlChildren(Rulebase)
	if (length(Rules.model) == 0)
		stop("There is no any rule.")
	
	## the function is used to generate one rule recursively.
	genRule <- function(CompoundPred){	
		operand <- XML::xmlAttrs(CompoundPred)[[1]]
		if (is.null(CompoundPred$children$CompoundPredicate)){
			rule <- c(XML::xmlAttrs(CompoundPred$children[[1]])[[1]], "is", XML::xmlAttrs(CompoundPred$children[[1]])[[2]], operand, XML::xmlAttrs(CompoundPred$children[[2]])[[1]], "is", XML::xmlAttrs(CompoundPred$children[[2]])[[2]])
			return(rule)
		}
		else {			
			rule <- c(XML::xmlAttrs(CompoundPred$children[[1]])[[1]], "is", XML::xmlAttrs(CompoundPred$children[[1]])[[2]], operand, genRule(CompoundPred$children$CompoundPredicate))
			return(rule)
		}
	}
	
	## check the number of rules
	if (length(Rules.model) != XML::xmlAttrs(Rulebase))
		stop("The number of rules is not correct.")
	
	rule.complete <- c()
	func.tsk <- c()
	grade.cert <- c()
	for (i in 1 : length(Rules.model)){
		## get tag for one rule
		ruleMod <- XML::xmlChildren(Rules.model[i]$Rule)		
		CompoundPred <- XML::xmlChildren(ruleMod$If)$CompoundPredicate
		
		## construct rules recursively
		rule.i <- genRule(CompoundPred)
		IfPart <- matrix(c("IF", rule.i), nrow = 1, byrow = TRUE)
		
		## Mamdani model
		if (type.model == "MAMDANI"){
			ThenPart <- matrix(c("THEN", XML::xmlAttrs(XML::xmlChildren(ruleMod$Then)$SimplePredicate)[[1]], "is", XML::xmlAttrs(XML::xmlChildren(ruleMod$Then)$SimplePredicate)[[2]]), nrow = 1)
			ruleComp.i <- cbind(IfPart, ThenPart)
			
			## accumulate the values
			rule.complete <- append(rule.complete, ruleComp.i)
		}
		
		## TSK model
		else if (type.model == "TSK"){
			## put THEN in the if part
			ruleComp.i <- cbind(IfPart, "THEN")
			
			## check the type of function
			if (XML::xmlAttrs(ruleMod$Then)[[1]]!="LinearFunction"){
				stop("The function is not supported")
			}
			
			## get list of the THEN part
			ThenValue <- XML::xmlChildren(ruleMod$Then)
			
			## inisialitation
			func.tsk.i <- c()
			
			## consider zero-order type (constant only)
			if (length(ThenValue) == 1){
				func.tsk.i <- c(func.tsk.i, ThenValue[1][[1]]$attributes[1])
			}
			else {
				## iterate as long as number of attributes + 1
				for (j in 1:length(ThenValue)){
					if (j == length(ThenValue)){
						func.tsk.i <- c(func.tsk.i, ThenValue[j][[1]]$attributes[1])
					}
					else{
						func.tsk.i <- c(func.tsk.i, ThenValue[j][[1]]$attributes[2])
					}
				}
			}
			## accumulate the values
			rule.complete <- append(rule.complete, ruleComp.i)
			func.tsk <- append(func.tsk, func.tsk.i)
		}
		else if (type.model == "FRBCS"){
			ThenPart <- matrix(c("THEN", XML::xmlAttrs(XML::xmlChildren(ruleMod$Then)$SimplePredicate)[[1]], "is", XML::xmlAttrs(XML::xmlChildren(ruleMod$Then)$SimplePredicate)[[2]]), nrow = 1)
			ruleComp.i <- cbind(IfPart, ThenPart)
			grade.cert.i <- as.numeric(XML::xmlChildren(ruleMod$Grade)[[1]]$value)
			
			## accumulate the values
			rule.complete <- append(rule.complete, ruleComp.i)
			grade.cert <- append(grade.cert, grade.cert.i)
		}
	}
	
	if (type.model == "MAMDANI"){
		## construct the matrix
		rule.complete <- matrix(rule.complete, nrow = length(Rules.model), byrow = TRUE)
		return(rule.complete)
	}
	else if (type.model == "TSK"){
		## construct the matrix
		rule.complete <- matrix(rule.complete, nrow = length(Rules.model), byrow = TRUE)
		func.tsk <- matrix(as.numeric(func.tsk), nrow = length(Rules.model), byrow = TRUE)
		if (is.null(func.tsk))
			stop("The functions of TSK must be defined.")
			
		return(list(rules = rule.complete, func.tsk = func.tsk))
	}
	else if (type.model == "FRBCS"){
		## construct the matrix
		rule.complete <- matrix(rule.complete, nrow = length(Rules.model), byrow = TRUE)
		grade.cert <- matrix(grade.cert, nrow = length(Rules.model), byrow = TRUE)
		class <- matrix(as.numeric(rule.complete[, ncol(rule.complete), drop = FALSE]), ncol = 1)
		return(list(rules = rule.complete, grade.cert = grade.cert, class = class))
	}
}

read.DATABASE <- function(model){
	requireNamespace("XML", quietly = TRUE)
	
	## Reading the header
    frbs.model <- model[["FrbsModel"]]
	
	## It is used to construct membership function for each fuzzy term
	construct.param <- function(type.mf, param){
		mf.var <- matrix(NA, nrow = 5, ncol = 1)
		if (type.mf == "TRIANGLE"){
			mf.var[1,1] <- 1
			mf.var[2,1] <- as.numeric(XML::xmlChildren(param)$Left[[1]]$value)
			mf.var[3,1] <- as.numeric(XML::xmlChildren(param)$Middle[[1]]$value)
			mf.var[4,1] <- as.numeric(XML::xmlChildren(param)$Right[[1]]$value)
		}
		else if (any(type.mf == c("TRAPEZOID"))){
			mf.var[1,1] <- 4
			mf.var[2,1] <- as.numeric(XML::xmlChildren(param)$Left[[1]]$value)
			mf.var[3,1] <- as.numeric(XML::xmlChildren(param)$LeftMiddle[[1]]$value)
			mf.var[4,1] <- as.numeric(XML::xmlChildren(param)$RightMiddle[[1]]$value)
			mf.var[5,1] <- as.numeric(XML::xmlChildren(param)$Right[[1]]$value)
		}
		else if (type.mf == "GAUSSIAN"){
			mf.var[1,1] <- 5
			mf.var[2,1] <- as.numeric(XML::xmlChildren(param)$Mean[[1]]$value)
			mf.var[3,1] <- as.numeric(XML::xmlChildren(param)$Variance[[1]]$value)
		}
		else {
			stop("The type of membership is not supported yet.")
		}
		return(mf.var)
	}
	
	## get Database tag
	MF.model <- XML::xmlChildren(frbs.model[names(frbs.model) == "Database"]$Database)
	num.var <- length(MF.model)
	num.labels <- matrix(nrow = 1, ncol = num.var)
	var.mf <- matrix(NA, nrow = 5, ncol = 1)
	name.fuzzyTerm <- c()
	for (i in 1:num.var){
		temp <- MF.model[i]$MembershipFunction
		num.labels[1, i] <- as.numeric(XML::xmlAttrs(temp)[2])
		mf.temp <- XML::xmlChildren(temp)
		num.mf <- length(mf.temp)
		if (length(mf.temp) != num.labels[1, i]){
			stop("The number of labels is not the same as membership functions.")
		}
		matrix.mf.temp <- matrix(NA, nrow = 5, ncol = num.mf)
		for (j in 1:num.mf){
			name.fuzzyTerm <- c(name.fuzzyTerm, XML::xmlAttrs(mf.temp[j]$FuzzyTerm)[[1]])
			type.mf <- XML::xmlAttrs(mf.temp[j]$FuzzyTerm)[[2]]
			param <- XML::xmlChildren(mf.temp[j]$FuzzyTerm)$Parameters
			matrix.mf.temp[, j] <- construct.param(type.mf, param)
		}
		var.mf <- cbind(var.mf, matrix.mf.temp)
	}
	var.mf <- var.mf[, -1, drop = FALSE]
	colnames(var.mf) <- name.fuzzyTerm
		
	return(list(var.MF = var.mf, num.labels = num.labels))
}


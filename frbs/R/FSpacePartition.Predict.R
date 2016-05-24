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
#' This function is one of the main internal functions of the package. 
#' It determines the values within the prediction phase.  
#'
#' This function involves four different processing steps on fuzzy rule-based systems. 
#' Firstly, the rulebase (see \code{\link{rulebase}}) validates 
#' the consistency of the fuzzy IF-THEN rules form. Then, the fuzzification 
#' (see \code{\link{fuzzifier}}) transforms crisp values 
#' into linguistic terms. Next, the inference calculates the degree of rule strengths using 
#' the t-norm and the s-norm. 
#' Finally, the defuzzification process calculates the results of the model using the Mamdani 
#' or the Takagi Sugeno Kang model.  
#'
#' @title The prediction phase
#' @param object the \code{\link{frbs-object}}.
#' @param newdata a matrix (\eqn{m \times n}) of data for the prediction process, 
#' where \eqn{m} is the number of instances and \eqn{n} is the number of input variables.
#' @seealso \code{\link{fuzzifier}}, \code{\link{rulebase}}, \code{\link{inference}} 
#' and \code{\link{defuzzifier}}.
#' @return A list with the following items:
#' \item{\code{rule}}{the fuzzy IF-THEN rules}
#' \item{\code{varinp.mf}}{a matrix to generate the shapes of the membership functions for 
#' the input variables}
#' \item{\code{MF}}{a matrix of the degrees of the membership functions}
#' \item{\code{miu.rule}}{a matrix of the degrees of the rules}
#' \item{\code{func.tsk}}{a matrix of the Takagi Sugeno Kang model for the consequent part of 
#' the fuzzy IF-THEN rules}
#' \item{\code{predicted.val}}{a matrix of the predicted values}
#' 
# @export
frbs.eng <- function(object, newdata){

	if (object$type.model == "TSK"){
		object$num.labels = cbind(object$num.labels, object$num.labels[1,1])
	}
	
	## get all of parameters
	range.output <- object$range.data.ori[, ncol(object$range.data.ori), drop = FALSE]
	num.varinput <- ncol(newdata)

	## change linguistic terms/labels to be unique
	temp <- ch.unique.fuzz(object$type.model, object$rule, object$varinp.mf, object$varout.mf, 
	        num.varinput, object$num.labels)

	rule <- temp$rule
	varinp.mf <- temp$varinp.mf
	varout.mf <- temp$varout.mf	
	names.fvalinput <- temp$names.fvalinput
	names.fvaloutput <- temp$names.fvaloutput	
	names.fvalues <- temp$names.fvalues	
	num.labels.input <- object$num.labels[, -ncol(object$num.labels), drop = FALSE]
	type.defuz <- object$type.defuz
	type.tnorm <- object$type.tnorm
	type.snorm <- object$type.snorm
	type.model <- object$type.model
	func.tsk <- object$func.tsk
	if (object$method.type != "MANUAL"){
		range.output[1, ] <- 0
		range.output[2, ] <- 1	
	}

	##################
	### I. Rule Based Module
	### In this function, Checking of the rule given by user will be done.
	### There are two kinds of model used which are Mamdani and TSK rule model.
	##################
	rule <- rulebase(type.model, rule, func.tsk)

	###################
	### II. Fuzzification Module
	### In this function, we convert crisp value into linguistic value based on the data and parameter of membership function.
	### There are several membership function can be used such as triangular, trapezoid, gaussian and logistic/sigmoid.
	###################
	
	MF <- fuzzifier(newdata, num.varinput, num.labels.input, varinp.mf)
	###################
	### III. Inference Module
	### In this function, we will calculate the confidence factor on antecedent for each rule. We use AND, OR, NOT operator. 
	###################
	ncol.MF <- ncol(MF)
	names.var <- names.fvalues[1 : ncol.MF]
	colnames(MF) <- c(names.var)
	
	miu.rule <- inference(MF, rule, names.fvalinput, type.tnorm, type.snorm)
	
	if(is.null(object$method.type) == FALSE && object$method.type == "FS.HGD"){
		if (!is.null(object$alpha.heuristic))
			miu.rule <- miu.rule ^ object$alpha.heuristic
		else 
			miu.rule <- miu.rule
	}
	
	if(is.null(object$method.type) == FALSE && object$method.type == "HYFIS"){
		if (!is.null(object$degree.rule))
			degree.rule <- object$degree.rule
		else 
			degree.rule <- 1
			
		for (i in 1 : nrow(miu.rule)){
			miu.rule[i, ] <- miu.rule[i, ] * t(degree.rule)
		}
	}

	###################
	### IV. Defuzzification Module
	### In this function, we calculate and convert linguistic value back into crisp value. 
	###################
	def <- defuzzifier(newdata, rule, range.output, names.fvaloutput, varout.mf, miu.rule, type.defuz, type.model, func.tsk)
	
	res <- list(rule = rule, varinp.mf = varinp.mf, varout.mf = varout.mf, MF = MF, miu.rule = miu.rule, func.tsk = func.tsk, predicted.val = def)
	return(res)
	
}

#' This function is the internal function of the fuzzy rule-based classification systems (FRBCS) to compute the predicted values.  
#'
#' @title FRBCS: prediction phase
#' @param object the \code{\link{frbs-object}}.
#' @param newdata a matrix (\eqn{m \times n}) of data for the prediction process, 
#' where \eqn{m} is the number of instances and \eqn{n} is the number of input variables.
#' @return A matrix of predicted values.
# @export
FRBCS.eng <- function(object, newdata){
	varinp.mf <- object$varinp.mf
	num.varinput <- ncol(object$num.labels) - 1
	num.labels.input <- object$num.labels[, -ncol(object$num.labels), drop = FALSE]
	num.fvaloutput <- object$num.labels[, ncol(object$num.labels), drop = FALSE]
	
	## Change the linguistic values
	seq.names <- rep(1:num.varinput, num.labels.input)
	names.fvalinput <- paste(seq.names, matrix(colnames(object$varinp.mf), nrow = 1), sep=".")
	names(varinp.mf) <- names.fvalinput
	
	rule <- object$rule
	if (rule[1, 1] == "IF"){
		k <- 1
		new.rule <- matrix(NA, nrow = nrow(rule), ncol = (2 * num.varinput + 1))
		for (i in 1 : num.varinput) {
			new.rule[, k] <- paste(i, rule[, (4 * i), drop = FALSE],sep=".")
			#new.rule[, k + 1] <- "and"
			new.rule[, k + 1] <- rule[, ((4 * i) + 1)] 
			k <- k + 2
		}
		new.rule[, (ncol(new.rule) - 1)] <- "->"
		new.rule[, ncol(new.rule)] <- rule[, ncol(rule), drop = FALSE]
		rule <- new.rule
	}

	classes <- as.matrix(as.numeric(rule[, ncol(rule), drop = FALSE]))
	
	if (!is.null(object$grade.cert)) {
		grade.certainty <- object$grade.cert
	} else { 
		grade.certainty <- cbind(as.numeric(rule[, ncol(rule), drop = FALSE]), 1)
	}
	
	type.tnorm <- object$type.tnorm
	type.snorm <- object$type.snorm
	type.model <- object$type.model
	
	##################
	### I. Rule Based Module
	### In this function, Checking of the rule given by user will be done.
	### There are two kinds of model used which are Mamdani and TSK rule model.
	##################
	rule <- rulebase(type.model, rule, classes)

	###################
	### II. Fuzzification Module
	### In this function, we convert crisp value into linguistic value based on the data and parameter of membership function.
	### There are several membership function can be used such as triangular, trapezoid, gaussian and logistic/sigmoid.
	###################
	MF <- fuzzifier(newdata, num.varinput, num.labels.input, varinp.mf)

	###################
	### III. Inference Module
	### In this function, we will calculate the confidence factor on antecedent for each rule. We use AND, OR, NOT operator. 
	################### 	
	miu.rule <- inference(MF, rule, names.fvalinput, type.tnorm, type.snorm)
	
	alpha.class.all <- miu.rule
	
	for (i in 1 : nrow(miu.rule)){
		alpha.class.all[i, ] <- t(miu.rule[i, ] * grade.certainty[, 1])
	}
	indx.max <- matrix()
	for (i in 1 : nrow(miu.rule)){
		indx.max[i] <- which.max(alpha.class.all[i, ])
	}	
	result <- matrix()
	
	for (i in 1 : length(indx.max)){
		result[i] <- object$class[indx.max[i], 1]
	}	
	res <- matrix(result)
	
	return(res)	
}
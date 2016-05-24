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
#' This is the internal function that implements the fuzzy inference rules by descent method (FIR.DM).
#' It is used to solve regression tasks. Users do not need to call it directly, 
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method was proposed by H. Nomura, I. Hayashi, and N. Wakami. FIR.DM uses simplified fuzzy 
#' reasoning where the consequent part is a real number (a particular case within the Takagi Sugeno Kang model),
#' while the membership function on the 
#' antecedent part is expressed by an isosceles triangle. So, in the learning phase, FIR.DM updates three parameters 
#' which are center and width of the triangular and a real number on the consequent part using a descent method.
#'
#' @title FIR.DM model building 
#'
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for training, where \eqn{m} is the number of instances and 
#'        \eqn{n} is the number of variables. The last column is the output variable.
#'        Note the data must be normalized between 0 and 1. 
#' @param num.labels a matrix (\eqn{1 \times n}) whose elements represent the number of labels (fuzzy terms),
#' where \eqn{n} is the number of variables.
#' @param max.iter the maximal number of iterations.
#' @param step.size the step size of the descent method, between 0 and 1.
#' @param type.tnorm the type of t-norm. For more detail, please have a look at \code{\link{inference}}. 
#' @param type.snorm the type of s-norm. For more detail, please have a look at \code{\link{inference}}. 
#' @param type.implication.func a value representing type of implication function. For more detail, please have a look at \code{\link{WM}}
#' @seealso \code{\link{DM.update}}, \code{\link{frbs.learn}}, and \code{\link{predict}}.
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking complete components. 
#' @references
#' H. Nomura, I. Hayashi and N. Wakami, "A learning method of fuzzy inference rules by descent method", 
#' IEEE International Conference on Fuzzy Systems, pp. 203 - 210 (1992).
# @export
FIR.DM <- function(data.train, num.labels, max.iter, step.size, type.tnorm = "MIN",
           type.snorm = "MAX", type.implication.func = "ZADEH"){
	
	## create progress bar
	progressbar <- txtProgressBar(min = 0, max = max.iter, style = 3)
	
	type.mf <- "TRIANGLE"
	range.data <- matrix(nrow = 2, ncol = ncol(data.train))
	range.data[1, ] <- 0
	range.data[2, ] <- 1
	
	## generate initial model
	mod <- generate.rule.GD(range.data, data.train, num.labels, type.mf, type.implication.func)
	
	## get data test from data training
	data.test <- as.matrix(data.train[, 1 : (ncol(data.train) - 1)])

	## get values of model	
	varinp.mf <- mod$varinp.mf
	names.varinp.mf <- colnames(varinp.mf)
	
	rule <- mod$rule[, -ncol(mod$rule), drop = FALSE]
	rule.data.num <- mod$rule.data.num
	mod$rule <- rule
	
	n.rowrule <- nrow(rule)	
	gal.iter <- matrix(nrow = max.iter)
	
	#generate func.tsk by random number
	num.ele <- n.rowrule 
	rand <- runif(num.ele, min = 0, max = 1)
	func.tsk <- matrix(rand, nrow = n.rowrule, ncol= 1, byrow = TRUE)
	
	mod$type.model <- "TSK"
	mod$func.tsk <- func.tsk
	mod$type.tnorm <- type.tnorm
	mod$type.snorm <- type.snorm
	mod$method.type <- "FIR.DM"
	mod$num.labels <- num.labels[, -ncol(num.labels), drop = FALSE]
	
	## calculate frbs.eng in order to obtain predicted values (def)
	init.res <- frbs.eng(mod, data.test)		
	init.def <- init.res$predicted.val
	
	init.miu.rule <- init.res$miu.rule
	
	## update initial parameters on consequent part
	for (h in 1 : nrow(data.test)) {
		init.gal <- init.def[h] - data.train[h, ncol(data.train)]
		for (i in 1 : nrow(func.tsk)){
			if (sum(init.miu.rule[h, ]) != 0)
				func.tsk[i, 1] <- func.tsk[i, 1] - step.size * init.miu.rule[h, i]/sum(init.miu.rule[h, ]) * init.gal 
		}
	}
	
	## calcualte error
	init.gal.l <- init.def - data.train[, ncol(data.train)]
	
	## update linear eq.
	mod$func.tsk <- func.tsk
	mod$num.labels <- num.labels[, -ncol(num.labels), drop = FALSE]
	galat <- matrix()
	error <- matrix()
	
	## update processes
	for (iter in 1 : max.iter){
		for (i in 1 : nrow(data.test)){
		
			data <- data.test[i, , drop = FALSE]
			dt.train.i <- data.train[i, ,drop = FALSE]
			res <- frbs.eng(mod, data)			
			def <- res$predicted.val
			
			error[i] <- def - dt.train.i[1, ncol(dt.train.i)]
			
			miu.rule <- res$miu.rule
			MF <- res$MF
			rule <- res$rule
			func.tsk <- res$func.tsk
			varinp.mf <- res$varinp.mf
			colnames(varinp.mf) <- names.varinp.mf
			
			## procedure for getting new parameters
			param.new <- DM.update(dt.train.i, rule.data.num, miu.rule, func.tsk, varinp.mf, step.size, def)

			func.tsk <- param.new$func.tsk.new
			varinp.mf <- param.new$varinp.mf.n
			
			## update with new params
			mod$varinp.mf <- varinp.mf
			colnames(mod$varinp.mf) <- names.varinp.mf			
			mod$func.tsk <- func.tsk
					
		}

		galat[iter] <- sqrt(mean(error^2))
		
		if (galat[iter] < 0.00001) {
			iter <- max.iter
			break
		}
		
		## progress bar
		setTxtProgressBar(progressbar, iter)
	}
	
	close(progressbar)
	rule <- mod$rule
	num.labels <- num.labels[, -ncol(num.labels), drop = FALSE]
	mod <- list(num.labels = num.labels, rule = rule, rule.data.num = rule.data.num, 
	          varinp.mf = varinp.mf, func.tsk = func.tsk, type.tnorm = type.tnorm, type.snorm = type.snorm,
			  type.mf = "TRIANGLE", type.model = "TSK", type.implication.func = type.implication.func)
	return (mod)
}

#' This is the internal function that implements the simplified TSK fuzzy rule generation method 
#' using heuristics and gradient descent method (FS.HGD). It is used to solve regression tasks. 
#' Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}.
#' 
#' This method was proposed by K. Nozaki, H. Ishibuchi, and H. Tanaka. It 
#' uses fuzzy IF-THEN rules with nonfuzzy singletons (i.e. real numbers) in
#' the consequent parts. The techniques of space partition are implemented to 
#' generate the antecedent part, while the initial consequent part of each 
#' rule is determined by the weighted mean value of the given training data. 
#' Then, the gradient descent method updates the value of the consequent part. 
#' Futhermore, the heuristic value given by the user affects the value of weight 
#' of each data. 
#' 
#' @title FS.HGD model building 
#'
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, where \eqn{m} is the number of instances and 
#'        \eqn{n} is the number of variables; the last column is the output variable.
#'        Note the data must be normalized between 0 and 1. 
#' @param num.labels a matrix (\eqn{1 \times n}), whose elements represent the number of labels (fuzzy terms); 
#' \eqn{n} is the number of variables.
#' @param max.iter maximal number of iterations.
#' @param step.size step size of the descent method. 
#' @param alpha.heuristic a positive real number which is the heuristic parameter.
#' @param type.tnorm the type of t-norm. For more detail, please have a look at \code{\link{inference}}. 
#' @param type.snorm the type of s-norm. For more detail, please have a look at \code{\link{inference}}. 
#' @param type.implication.func a value representing type of implication function. 
#' For more detail, please have a look at \code{\link{WM}}. 
#' @seealso \code{\link{frbs.learn}}, \code{\link{predict}}, and \code{\link{HGD.update}}
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking its complete components. 
#' @references
#' H. Ishibuchi, K. Nozaki, H. Tanaka, Y. Hosaka, and M. Matsuda, "Empirical study on learning in fuzzy systems by rice taste analysis",
#' Fuzzy Set and Systems, vol. 64, no. 2, pp. 129 - 144 (1994).
# @export
FS.HGD <- function(data.train, num.labels, max.iter = 100, step.size = 0.01, alpha.heuristic = 1, 
                 type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH"){
	
	## create progress bar
	progressbar <- txtProgressBar(min = 0, max = max.iter, style = 3)
		
	type.model <- "MAMDANI"
	func.tsk = NULL
	type.mf <- "TRIANGLE"
	range.data <- matrix(nrow = 2, ncol = ncol(data.train))
	range.data[1, ] <- 0
	range.data[2, ] <- 1
	
	## generate initial model
	mod <- generate.rule.GD(range.data, data.train, num.labels, type.mf, type.tnorm, type.implication.func)

	data.test <- as.matrix(data.train[, 1 : (ncol(data.train) - 1)])
	target.dt <- as.matrix(data.train[, ncol(data.train)], ncol = 1)
	
	##	get values from model
	num.varinput <- ncol(num.labels) - 1
	num.fvalinput <- num.labels[, -ncol(num.labels), drop = FALSE]
	
	#names.varinput <- mod$names.varinput
	names.fvalinput <- colnames(mod$varinp.mf)
	
	rule <- mod$rule
	varinp.mf <- mod$varinp.mf
	names.varinp.mf <- colnames(varinp.mf)
	rule.data.num <- mod$rule.data.num
	rule <- mod$rule[, 1 : (ncol(mod$rule) - 1)]
	mod$rule <- rule
	
	#number of rule
	n.rowrule <- nrow(rule)
	
	gal.iter <- matrix(nrow = max.iter)
	
	#generate func.tsk by heuristics
	#in order to initialize W
	data <- data.test
	rule <- rulebase(type.model, rule, func.tsk)
	MF <- fuzzifier(data, num.varinput, num.fvalinput, varinp.mf)
	
	miu.rule <- inference(MF, rule, names.fvalinput, type.tnorm, type.snorm)
	miu.rule <- miu.rule ^ alpha.heuristic
	
	func.tsk <- matrix(nrow = n.rowrule, ncol= 1, byrow=T)
	func.tsk <- (t(miu.rule) %*% target.dt) / colSums(miu.rule)		
	func.tsk[which(is.nan(func.tsk))] <- runif(1, 0, 1)
	mod$num.labels <- num.labels[, -ncol(num.labels), drop = FALSE]
	mod$type.model <- "TSK"
	mod$func.tsk <- func.tsk
	mod$type.tnorm <- type.tnorm
	mod$type.snorm <- type.snorm
	mod$method.type <- "FS.HGD"
	mod$alpha.heuristic <- alpha.heuristic
	
	galat <- matrix()
	error <- matrix()

	for (iter in 1 : max.iter){
		for (i in 1 : nrow(data.test)){
		
			dt.i <- data.test[i, ,drop = FALSE]
			dt.train.i <- data.train[i, ,drop = FALSE]
		
			data <- dt.i
		
			res <- frbs.eng(mod, data)
	
			def <- res$predicted.val
			
			error[i] <- def - dt.train.i[1, ncol(dt.train.i)]
			miu.rule <- res$miu.rule
			MF <- res$MF
			rule <- res$rule
			func.tsk <- res$func.tsk
			varinp.mf <- res$varinp.mf
			colnames(varinp.mf) <- names.varinp.mf
			
			## update procedure
			param.new <- HGD.update(dt.train.i, miu.rule, func.tsk, varinp.mf, step.size, def)
			
			## get new parameters
			func.tsk <- unlist(param.new[[1]])			
			mod$func.tsk <- func.tsk		
		}
	
		galat[iter] <- sqrt(mean(error^2))
	
		## progress bar
		setTxtProgressBar(progressbar, iter)
	}
	
	close(progressbar)
	rule <- mod$rule
	num.labels <- num.labels[, -ncol(num.labels), drop = FALSE]	
	mod <- list(num.labels = num.labels, rule = rule, rule.data.num = rule.data.num, varinp.mf = varinp.mf,
          	func.tsk = func.tsk, alpha.heuristic = alpha.heuristic, type.tnorm = type.tnorm, type.snorm = type.snorm,
			type.mf = "TRIANGLE", type.model = "TSK", type.implication.func = type.implication.func)
		
	return (mod)
}

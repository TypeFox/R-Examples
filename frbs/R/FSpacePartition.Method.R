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
#' This is the internal function that implements the model proposed by L. X. Wang and J. M. 
#' Mendel. It is used to solve regression task. Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}
#' 
#' The fuzzy rule-based system for learning from L. X. Wang and J. M. 
#' Mendel's paper is implemented in this function. For the learning process, there are four stages as follows: 
#' \itemize{
#' \item \code{Step 1:} Divide equally the input and output spaces of the given numerical data into 
#'      fuzzy regions as the database. In this case, fuzzy regions refers to intervals for each 
#'      linguistic term. Therefore, the length of fuzzy regions represents the number of 
#'      linguistic terms. For example, the linguistic term "hot" has the fuzzy region \eqn{[1, 3]}. 
#'      We can construct a triangular membership function having the corner points 
#'      \eqn{a = 1}, \eqn{b = 2}, and \eqn{c = 3} where \eqn{b} is a middle point 
#'      that its degree of the membership function equals one. 
#' \item \code{Step 2:} Generate fuzzy IF-THEN rules covering the training data, 
#'      using the database from Step 1. First, we calculate degrees of the membership function
#'      for all values in the training data. For each instance in the training data, 
#'      we determine a linguistic term having a maximum degree in each variable. 
#'      Then, we repeat the process for each instance in the training data to construct 
#'      fuzzy rules covering the training data.
#' \item \code{Step 3:} Determine a degree for each rule. 
#'      Degrees of each rule are determined by aggregating the degree of membership functions in 
#'      the antecedent and consequent parts. In this case, we are using the product aggregation operators. 
#' \item \code{Step 4:} Obtain a final rule base after deleting redundant rules. 
#'      Considering degrees of rules, we can delete the redundant rules having lower degrees. 
#' }
#' The outcome is a Mamdani model. In the prediction phase, 
#' there are four steps: fuzzification, checking the rules, inference, and defuzzification.
#' 
#' @title WM model building 
#'
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, where \eqn{m} is the number of instances and 
#'        \eqn{n} is the number of variables; the last column is the output variable. 
#'        Note the data must be normalized between 0 and 1. 
#' @param num.labels a matrix (\eqn{1 \times n}), whose elements represent the number of labels (linguistic terms); 
#' \eqn{n} is the number of variables.
#' @param type.mf the type of the membership function. See \code{\link{frbs.learn}}.
#' @param type.tnorm a value which represents the type of t-norm. See \code{\link{inference}}.
#' @param type.implication.func a value representing type of implication function. Let us consider a rule, \eqn{a \to b},  
#' \itemize{
#' \item \code{DIENES_RESHER} means \eqn{(b > 1 - a? b : 1 - a)}.
#' \item \code{LUKASIEWICZ} means \eqn{(b < a ? 1 - a + b : 1)}.
#' \item \code{ZADEH} means \eqn{(a < 0.5 || 1 - a > b ? 1 - a : (a < b ? a : b))}.
#' \item \code{GOGUEN} means \eqn{(a < b ? 1 : b / a)}.
#' \item \code{GODEL} means \eqn{(a <= b ? 1 : b)}.
#' \item \code{SHARP} means \eqn{(a <= b ? 1 : 0)}.
#' \item \code{MIZUMOTO} means \eqn{(1 - a + a * b)}.
#' \item \code{DUBOIS_PRADE} means \eqn{(b == 0 ? 1 - a : (a == 1 ? b : 1))}.
#' \item \code{MIN} means \eqn{(a < b ? a : b)}.
#' }
#' @param classification a boolean representing whether it is a classification problem or not.
#' @param range.data a matrix representing interval of data.
#' @seealso \code{\link{frbs.learn}}, \code{\link{predict}} and \code{\link{frbs.eng}}.
# @return an object of type \code{frbs}. Please have a look at An \code{\link{frbs-object}} for looking its complete components.
#' @references 
#' L.X. Wang and J.M. Mendel, "Generating fuzzy rule by learning from examples", IEEE Trans. Syst., Man, and Cybern.,
#' vol. 22, no. 6, pp. 1414 - 1427 (1992).
# @export
WM <- function(data.train, num.labels, type.mf = "GAUSSIAN", type.tnorm = "PRODUCT", type.implication.func = "ZADEH", classification = FALSE, range.data = NULL) {
	
	if (!is.null(range.data)) {
		 range.data = range.data
	}
	else if (classification == FALSE && is.null(range.data)) {
		range.data <- matrix(nrow = 2, ncol = ncol(data.train))
		range.data[1, ] <- 0
		range.data[2, ] <- 1
	} else {
		## make range of data according to class on data training
	    range.data.out <- matrix(c(0.5001, num.labels[1, ncol(num.labels)] + 0.4999), nrow = 2)
		range.data.inp <- matrix(nrow = 2, ncol = (ncol(data.train) - 1))
		range.data.inp[1, ] <- 0
		range.data.inp[2, ] <- 1
		range.data <- cbind(range.data.inp, range.data.out)
	}
	
	## build the names of linguistic values
	fuzzyTerm <- create.fuzzyTerm(classification, num.labels)
	
	names.fvalinput <- fuzzyTerm$names.fvalinput
	names.fvaloutput <- fuzzyTerm$names.fvaloutput
	names.fvalues <- c(names.fvalinput, names.fvaloutput)
	
	#get number of row and column
	nrow.data.train <- nrow(data.train)
	num.varinput <- ncol(data.train)
	
	##initialize
	seg <- matrix(ncol = 1)
	center.value <- matrix(ncol = 1)
	var.mf <- matrix(0, nrow = 5, ncol = sum(num.labels))
	
	## use data.train as data
	data <- data.train
	
	####### Wang and Mendel's Steps begin ####################
	###Step 1: Divide the Input and Output Spaces Into Fuzzy Regions
	
	## initialize
	jj <- 0
	
	## loop for all variable
	for (i in 1 : num.varinput){
		## initialize
		seg <- num.labels[1, i]
		kk <- 1
			
		## Make the depth on each linguistic value, assumed it's similar on all region. 
		delta.point <- (range.data[2, i] - range.data[1, i]) / (seg + 1)
		
			##contruct matrix of parameter of membership function (var.mf) on each variable (max(var) is num.varinput
			for (j in 1 : num.labels[1, i]){
			
				##counter to continue index of column var.mf					
				jj <- jj + 1
				
				## type.mf <- 1 means using trapezoid and triangular
				if (type.mf == 1 || type.mf == "TRIANGLE") {
				
					delta.tri <- (range.data[2, i] - range.data[1, i]) / (seg - 1)
					## on the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 1	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- var.mf[2, jj]
						var.mf[4, jj] <- var.mf[3, jj] + delta.tri
						var.mf[5, jj] <- NA						
					}
					
					## on the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 1	
						var.mf[4, jj] <- range.data[2, i]
						var.mf[3, jj] <- var.mf[4, jj] 
						var.mf[2, jj] <- var.mf[3, jj] - delta.tri
						var.mf[5, jj] <- NA						
					}
					
					## on the middle
					else{
						var.mf[1, jj] <- 1	
						var.mf[3, jj] <- range.data[1, i] + (j - 1) * delta.tri
						var.mf[2, jj] <-  var.mf[3, jj] - delta.tri
						var.mf[4, jj] <- var.mf[3, jj] + delta.tri
						var.mf[5, jj] <- NA						
					}	
				}
				## type.mf == 2 means we use trapezoid
				else if (type.mf == 2 || type.mf == "TRAPEZOID") {
					delta.tra <- (range.data[2, i] - range.data[1, i]) / (seg + 2)
					
					## on the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 2	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- var.mf[2, jj] + delta.tra
						var.mf[4, jj] <- var.mf[3, jj] + delta.tra
						var.mf[5, jj] <- NA												
					}
					
					## on the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 3	
						var.mf[4, jj] <- range.data[2, i]
						var.mf[3, jj] <- var.mf[4, jj] - delta.tra
						var.mf[2, jj] <- var.mf[3, jj] - delta.tra
						var.mf[5, jj] <- NA												
					}
					
					## on the middle
					else{
						var.mf[1, jj] <- 4	
						var.mf[2, jj] <- range.data[1, i] + (j - 1) * 1.15 * delta.tra
						var.mf[3, jj] <- var.mf[2, jj] + delta.tra
						var.mf[4, jj] <- var.mf[3, jj] + 0.5 * delta.tra					
						var.mf[5, jj] <- var.mf[4, jj] + delta.tra											
					}
				}
				##Type 5: Gaussian
				##parameter=(mean a, standard deviation b)
				##a=var.mf[2,]
				##b=var.mf[3,]
				else if (type.mf == 3 || type.mf == "GAUSSIAN") {
					delta.gau <- (range.data[2, i] - range.data[1, i]) / (seg - 1)
					##On the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 5	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- 0.35 * delta.gau
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					##On the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 5	
						var.mf[2, jj] <- range.data[2, i]
						var.mf[3, jj] <- 0.35 * delta.gau
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					## On the middle
					else {
					var.mf[1, jj] <- 5	
					var.mf[2, jj] <- range.data[1, i] + (j - 1) * delta.gau
					var.mf[3, jj] <- 0.35 * delta.gau
					var.mf[4, jj] <- NA
					var.mf[5, jj] <- NA
					}
				}
				##Type 6: Sigmoid/logistic
				##parameter=(gamma,c)
				##gamma=var.mf[2,]
				##c=var.mf[3,]
				else if (type.mf == 4 || type.mf == "SIGMOID") {
					delta.sig <- (range.data[2, i] - range.data[1, i]) / (seg + 1)
					##On the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 6	
						var.mf[2, jj] <- range.data[1, i]
						var.mf[3, jj] <- delta.sig
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					##On the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 6	
						var.mf[2, jj] <- range.data[2, i]
						var.mf[3, jj] <- delta.sig
						var.mf[4, jj] <- NA
						var.mf[5, jj] <- NA
					}
					## On the middle
					else {
					var.mf[1, jj] <- 6	
					var.mf[2, jj] <- range.data[1, i] + (j - 1) * delta.sig
					var.mf[3, jj] <- delta.sig
					var.mf[4, jj] <- NA
					var.mf[5, jj] <- NA
					}
				}
				##checking for type 7: Generalized Bell
				##parameter=(a,b,c)
				##a=var.mf[2,]
				##b=var.mf[3,]
				##c=var.mf[4,]
				else if (type.mf == 5 || type.mf == "BELL") {
					##On the left side
					if (kk %% num.labels[1, i] == 1){   
						var.mf[1, jj] <- 7	
						var.mf[2, jj] <- 0.6 * delta.point
						var.mf[3, jj] <- 0.6 * delta.point
						var.mf[4, jj] <- range.data[1, i]
						var.mf[5, jj] <- NA
					}
					##On the right side
					else if (kk %% num.labels[1, i] == 0){
						var.mf[1, jj] <- 7	
						var.mf[2, jj] <- 0.6 * delta.point
						var.mf[3, jj] <- 0.6 * delta.point
						var.mf[4, jj] <- range.data[2, i]
						var.mf[5, jj] <- NA
					}
					## On the middle
					else {
					
					var.mf[1, jj] <- 7	
					var.mf[2, jj] <- 0.6 * delta.point
					var.mf[3, jj] <- 0.6 * delta.point
					var.mf[4, jj] <- range.data[1, i] + delta.point * j
					var.mf[5, jj] <- NA
					}
				}						
				kk <- kk + 1
			}
	}

	if (classification == TRUE){
		range.data.out <- range.data[, ncol(range.data), drop = FALSE]
		num.fvaloutput <- num.labels[1, ncol(num.labels)]
		seg <- num.labels[1, ncol(num.labels)]
		delta.tra <- (range.data.out[2, 1] - range.data.out[1, 1]) / seg
		ncol.var.mf <- ncol(var.mf)
	
		k <- 1
		for (j in (ncol.var.mf - num.fvaloutput + 1) : ncol.var.mf){
			##On the left side	
			var.mf[1, j] <- 4	
			var.mf[2, j] <- range.data.out[1, 1] + (k - 1) * delta.tra 
			var.mf[3, j] <- var.mf[2, j] 
			var.mf[4, j] <- var.mf[3, j] + delta.tra
			var.mf[5, j] <- var.mf[4, j] 					
			k <- k + 1
		}
	}
	
	## Step 2: Generate Fuzzy Rules from Given Data Pairs.
	## Step 2a: Determine the degree of data pairs.
	## MF is matrix membership function. The dimension of MF is n x m, where n is number of data and m is num.labels * input variable (==ncol(var.mf))
		
	## get degree of membership by fuzzification
	MF <- fuzzifier(data, num.varinput, num.labels, var.mf)	
	colnames(MF) <- c(names.fvalues)
	
	####get max value of degree on each variable to get one rule.
	MF.max <- matrix(0, nrow = nrow(MF), ncol = ncol(MF))
	k <- 1
	for (i in 1 : length(num.labels)){
		start <- k
		end <- start + num.labels[1, i] - 1
		MF.temp <- MF[, start : end]
		
		for (m in 1 : nrow(MF)){
			max.MFTemp <- max(MF.temp[m, ], na.rm = TRUE)
			max.loc <- which.max(MF.temp[m, ])
			MF.max[m, k + max.loc - 1] <- max.MFTemp		
		}
		k <- k + num.labels[1, i]	
	}
	
	colnames(MF.max) <- c(names.fvalues)
	rule.matrix <- MF.max
	
	## Step III
	## determine the degree of the rule
	degree.rule <- matrix(nrow = nrow(rule.matrix), ncol =1)
	degree.ante <- matrix(nrow = nrow(rule.matrix), ncol =1)

	## calculate t-norm in antecedent part and implication function with consequent part
	for (i in 1:nrow(rule.matrix)){
		temp.ant.rule.degree <- c(1)
		max.ant.indx <- c(ncol(rule.matrix) - num.labels[1, ncol(num.labels)])
		for (j in 1 : max.ant.indx){
			if (rule.matrix[i, j] != 0){
				if (type.tnorm == "PRODUCT"){
					temp.ant.rule.degree <- temp.ant.rule.degree * rule.matrix[i, j]					
				}			
				else if (type.tnorm == "MIN"){
					temp.ant.rule.degree <- min(temp.ant.rule.degree, rule.matrix[i, j])
				}
				else if (type.tnorm == "HAMACHER"){
					temp.ant.rule.degree <- (temp.ant.rule.degree * rule.matrix[i, j]) / (temp.ant.rule.degree + rule.matrix[i, j] - temp.ant.rule.degree * rule.matrix[i, j])
				}
				else if (type.tnorm == "YAGER"){
					temp.ant.rule.degree <- 1 - min(1, ((1 - temp.ant.rule.degree) + (1 - rule.matrix[i, j])))
				}
				else if (type.tnorm == "BOUNDED"){
					temp.ant.rule.degree <- max(0, (temp.ant.rule.degree + rule.matrix[i, j] - 1))
				}
			}
		}
		degree.ante[i] <- temp.ant.rule.degree
		
		for (k in (max.ant.indx + 1) : ncol(rule.matrix)){
			if (rule.matrix[i, k] != 0){
				temp.rule.degree <- calc.implFunc(degree.ante[i], rule.matrix[i, k], type.implication.func)
			}
		}
		degree.rule[i] <- temp.rule.degree
	}

	rule.matrix[rule.matrix > 0] <- 1
	rule.matrix.bool <- rule.matrix

	## find the same elements on matrix rule considering degree of rules
	rule.red.mat <- rule.matrix
	temp <- cbind(degree.rule, degree.ante, rule.matrix.bool)
	temp <- temp[order(temp[,c(1)], decreasing = TRUE),]
	indx.nondup <- which(duplicated(temp[, -c(1,2)]) == FALSE, arr.ind = TRUE)
	rule.complete <- temp[indx.nondup, ,drop = FALSE] 
	degree.ante <- rule.complete[, 2, drop = FALSE] 
	degree.rule <- rule.complete[, 1, drop = FALSE]
	rule.red.mat <- rule.complete[, 3:ncol(rule.complete)]
	
	## delete incomplete rule
	ind.incomplete <- which(rowSums(rule.red.mat) != num.varinput)
	rule.red.mat[ind.incomplete, ] <- NA	
	degree.rule[ind.incomplete] <- NA
	degree.ante[ind.incomplete] <- NA
	
	## create rule in numeric
	seqq <- seq(1:ncol(rule.red.mat))
	rule.data.num <- t(apply(rule.red.mat, 1, function(x) x * seqq))
	
	rule.data.num <- na.omit(rule.data.num)
	degree.ante <- na.omit(degree.ante)
	degree.rule <- na.omit(degree.rule)
	
	## rule.data.num is the numbers representing the sequence of string of variable names 	
	rule.data.num[which(rule.data.num == 0)] <- NA
	rule.data.num <- t(apply(rule.data.num, 1, na.omit))
	
	## build the rule into list of string
	res <- generate.rule(rule.data.num, num.labels, names.fvalinput, names.fvaloutput)
	rule <- res$rule
	
	#############################################################
	### Collect the Input data for "frbs for testing"
	#############################################################
		
	## indx is index of output variable on num.labels
	indx <- length(num.labels)
	length.namesvar <- length(names.fvalues)
		
	## cut membership function of output variable on var.mf
	ll <- (ncol(var.mf) - num.labels[1, indx])
	varinp.mf <- var.mf[, 1 : ll, drop = FALSE]
	colnames(varinp.mf) <- names.fvalinput
	## get varout.mf
	varout.mf <- var.mf[, (ll + 1) : ncol(var.mf), drop = FALSE]
	colnames(varout.mf) <- names.fvaloutput
	
	## clean rule
	rule <- na.omit(rule)
	
	## range of original data
	range.data.ori <- range.data
	
	mod <- list(num.labels = num.labels, varout.mf = varout.mf, rule = rule, varinp.mf = varinp.mf, degree.ante = degree.ante, rule.data.num = rule.data.num, 
			   degree.rule = degree.rule, range.data.ori = range.data.ori, type.mf = type.mf, type.tnorm = type.tnorm, type.implication.func = type.implication.func)

	return (mod)
}  

#' This is the internal function that implements the fuzzy rule-based classification 
#' system using Chi's technique (FRBCS.CHI). It is used to solve classification tasks. 
#' Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}. This method is
#' suitable only for classification problems.
#' 
#' This method was proposed by Z. Chi, H. Yan, and T. Pham that extends
#' Wang and Mendel's method for tackling classification problems. 
#' Basically, the algorithm is quite similar as Wang and Mendel's technique. 
#' However, since it is based on the FRBCS model, Chi's method only takes class labels on each data
#' to be consequent parts of fuzzy IF-THEN rules. In other words, we generate rules as in 
#' Wang and Mendel's technique (\code{\link{WM}}) and then we replace consequent parts with their classes. 
#' Regarding calculating degress of each rule, they are determined by antecedent parts of the rules. 
#' Redudant rules can be deleted by considering their degrees. Lastly, we obtain fuzzy IF-THEN rules 
#' based on the FRBCS model.
#'
#' @title FRBCS.CHI model building 
#'
#' @param range.data a matrix (\eqn{2 \times n}) containing the range of the normalized data, where \eqn{n} is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, where \eqn{m} is the number of instances and 
#' \eqn{n} is the number of variables; the last column is the output variable. Note the data must be normalized between 0 and 1. 
#' @param num.labels a matrix (\eqn{1 \times n}), whose elements represent the number of labels (linguistic terms); 
#' \eqn{n} is the number of variables.
#' @param num.class an integer number representing the number of labels (linguistic terms).
#' @param type.mf the type of the shape of the membership functions. See \code{\link{fuzzifier}}.
#' @param type.tnorm the type of t-norm. See \code{\link{inference}}.
#' @param type.snorm the type of s-norm. See \code{\link{inference}}.
#' @param type.implication.func the type of implication function. See \code{\link{WM}}.
#' @seealso \code{\link{FRBCS.eng}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking its complete components. 
#' @references
#' Z. Chi, H. Yan, T. Pham, "Fuzzy algorithms with applications to image processing 
#' and pattern recognition", World Scientific, Singapore (1996).
# @export
FRBCS.CHI <- function(range.data, data.train, num.labels, num.class, type.mf = "TRIANGLE", type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH") {

	num.labels[1, ncol(num.labels)] <- num.class

	## generate initial model using WM
	mod <- WM(data.train, num.labels, type.mf, type.tnorm = type.tnorm, type.implication.func = type.implication.func, classification = TRUE)
	
	names.fvaloutput <- as.matrix(colnames(mod$varout.mf))
	rule <- mod$rule
	rule.constant <- mod$rule
	degree.ante <- mod$degree.ante
	
	rule.constant <- cbind(rule.constant[, -ncol(rule.constant), drop = FALSE], 
	                      as.numeric(factor(rule.constant[, ncol(rule.constant)])))
	
	classes <- matrix(strtoi(rule.constant[, ncol(rule.constant)])) 	
	grade.cert <- degree.ante
	
	mod$class <- classes
	mod$rule[, ncol(mod$rule)] <- classes
	mod$grade.cert <- grade.cert
	mod$varout.mf <- NULL
	mod$type.tnorm <- type.tnorm
	mod$type.snorm <- type.snorm
	mod$type.model <- "FRBCS"
	mod$type.mf <- type.mf
	mod$type.implication.func
	return (mod)	
}

#' This is the internal function that implements the fuzzy rule-based classification 
#' system with weight factor (FRBCS.W). It is used to solve classification tasks. 
#' Users do not need to call it directly,
#' but just use \code{\link{frbs.learn}} and \code{\link{predict}}. This method is
#' suitable only for classification problems.
#' 
#' This method is adopted from Ishibuchi and Nakashima's paper. 
#' Each fuzzy IF-THEN rule consists of antecedent linguistic values and a single consequent class with certainty grades 
#' (weights). The antecedent part is determined by a grid-type fuzzy partition from 
#' the training data. The consequent class is defined as the dominant class in 
#' the fuzzy subspace corresponding to the antecedent part of each fuzzy IF-THEN rule and 
#' the certainty grade is calculated from the ratio among the consequent class. 
#' A class of the new instance is determined by the consequent class of the rule with 
#' the maximal product of the compatibility grade and the certainty grade. 
#'
#' @title FRBCS.W model building 
#'
#' @param range.data a matrix (\eqn{2 \times n}) containing the range of the normalized data, where \eqn{n} is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, where \eqn{m} is the number of instances and 
#' \eqn{n} is the number of variables; the last column is the output variable. Note the data must be normalized between 0 and 1. 
#' @param num.labels a matrix (\eqn{1 \times n}), whose elements represent the number of labels (linguistic terms); 
#' \eqn{n} is the number of variables.
#' @param num.class an integer number representing the number of labels (linguistic terms).
#' @param type.mf the type of the shape of the membership functions.
#' @param type.tnorm the type of t-norm. See \code{\link{inference}}.
#' @param type.snorm the type of s-norm. See \code{\link{inference}}.
#' @param type.implication.func the type of implication function. See \code{\link{WM}}.
#' @seealso \code{\link{FRBCS.eng}}, \code{\link{frbs.learn}}, and \code{\link{predict}}
# @return a list of the model data. Please have a look at \code{\link{frbs.learn}} for looking its complete components. 
#' @references
#' H. Ishibuchi and T. Nakashima, "Effect of rule weights in fuzzy rule-based classification systems", 
#' IEEE Transactions on Fuzzy Systems, vol. 1, pp. 59 - 64 (2001).
# @export
FRBCS.W <- function(range.data, data.train, num.labels, num.class, type.mf, type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH") {

	num.labels[1, ncol(num.labels)] <- num.class
	
	## generate initial model using WM
	mod <- WM(data.train, num.labels, type.mf, type.tnorm = type.tnorm, type.implication.func = type.implication.func, 
	          classification = TRUE)
	
	names.fvaloutput <- as.matrix(colnames(mod$varout.mf))
	rule <- mod$rule
	rule.constant <- mod$rule
	degree.ante <- mod$degree.ante
	
	rule.constant <- cbind(rule.constant[, -ncol(rule.constant), drop = FALSE], 
	                      as.numeric(factor(rule.constant[, ncol(rule.constant)])))
	
	classes <- matrix(strtoi(rule.constant[, ncol(rule.constant)])) 
	
	grade.cert.temp <- cbind(classes, degree.ante)
	class.fact <- factor(grade.cert.temp[, 1], exclude = "")
	beta.class <- as.double(as.matrix(aggregate(grade.cert.temp[, 2], by = list(class.fact), sum)))
	beta.class <- matrix(beta.class, nrow = num.class)	
	
	CF <- matrix()
	 
	 for(i in 1 : nrow(beta.class)){
		CF[i] <- 1 + (beta.class[i, 2] - (sum(beta.class[, 2]) - (beta.class[i, 2] / (nrow(beta.class) - 1))))/sum(beta.class[, 2])
	 }
	grade.cert <- grade.cert.temp 
	for(i in 1 : nrow(grade.cert.temp)){
		grade.cert[i, 2] <- CF[grade.cert[i,1]] 
	}
	mod$class <- classes
	mod$rule[, ncol(mod$rule)] <- classes
	mod$grade.cert <- grade.cert[, 2, drop = FALSE]
	mod$varout.mf <- NULL
	mod$type.tnorm <- type.tnorm
	mod$type.snorm <- type.snorm
	mod$type.model <- "FRBCS"
	mod$type.mf <- type.mf
	mod$type.implication.func <- type.implication.func
	
	return (mod)
}  

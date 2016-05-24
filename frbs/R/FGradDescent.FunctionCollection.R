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
#' The role of this function is to update the parameters of the fuzzy inference rules by descent method (FIR.DM). 
#' This function is called by the main function of the FIR.DM method, \code{\link{FIR.DM}}.
#'
#' @title FIR.DM updating function
#'
#' @param data.train a matrix (\eqn{m \times n}) of normalized data, where \eqn{m} is the number of instances and 
#' \eqn{n} is the number of variables; the last column is the output variable.
#' @param rule.data.num a matrix containing the rulebase. Its elements are integers, see \code{\link{rulebase}}.
#' @param miu.rule a matrix with the degrees of rules which is a result of the \code{\link{inference}}.
#' @param func.tsk a matrix of parameters of the functions on the consequent part of the Takagi Sugeno Kang model.
#' @param varinp.mf a matrix of parameters of the membership functions of the input variables.
#' @param step.size the step size of the descent method, between 0 and 1. 
#' @param def a matrix which is obtained from the defuzzification. Please have a look at \code{\link{defuzzifier}}.
#' @seealso \code{\link{frbs.learn}}, \code{\link{predict}}, and \code{\link{FIR.DM}}.
# @return A list containing the following components:
# \item{\code{func.tsk.new}}{a new function of Takagi Sugeno Kang}
# \item{\code{varinp.mf.n}}{a new matrix of membership function}
# @export
DM.update <- function(data.train, rule.data.num, miu.rule, func.tsk, varinp.mf, step.size = 0.01, def){

data.m <- data.train[1, 1 : (ncol(data.train) - 1), drop = FALSE]
	
ff <- func.tsk 
miu.rule.t <- as.matrix(miu.rule)

#get margin between predicted and real value
gal <- (def - data.train[1, ncol(data.train)])
	
func.tsk.upd <- func.tsk

## procedure update w or contant on consequent part	
for (mm in 1 : nrow(func.tsk.upd)){
	sum.miu <- sum(miu.rule.t, na.rm = TRUE)
	if (sum.miu != 0)
		func.tsk.upd[mm, 1] <- func.tsk.upd[mm, 1] - step.size * gal * (miu.rule[mm]/sum.miu)
}

## update
func.tsk.new <- func.tsk.upd
num.varinput <- ncol(data.train) - 1


####################################################
#=== update input variable using gradient descent
###################################################
varinp.mf.n <- varinp.mf
for (ii in 1 : ncol(miu.rule.t)) {
	for (j in 1 : num.varinput){
		gap <- (varinp.mf[4, rule.data.num[ii, j]] - varinp.mf[2, rule.data.num[ii, j]])
			
		# in case we got bad shape.
		if (gap <= 0)
			gap <- 0.001
		
		##equation for getting a degree of membership function
		A <- 1 - (2 * abs(data.m[1, j] - varinp.mf[3, rule.data.num[ii, j]])/(gap))
		
		if (A == 0){
			A <- 0.001
		}
					
		if (sum.miu != 0) {
			varinp.mf.n[3, rule.data.num[ii, j]] <- varinp.mf[3, rule.data.num[ii, j]] - ((step.size * miu.rule[1, ii]) / sum.miu) * gal * (func.tsk[ii] - def) * 
							sign(data.m[1, j] - varinp.mf[3, rule.data.num[ii, j]]) * 2 / (gap * A)
							
			temp.b <- gap - step.size * miu.rule[1, ii] / sum.miu  * gal * (func.tsk[ii] - def) * ((1 - A) / A) * 1 / gap	
		}
		
		else {
				varinp.mf.n[3, rule.data.num[ii, j]] <- varinp.mf[3, rule.data.num[ii, j]]
				temp.b <- gap 
		}
			
		if (varinp.mf.n[3, rule.data.num[ii, j]] < 0)
			varinp.mf.n[3, rule.data.num[ii, j]] <- varinp.mf[3, rule.data.num[ii, j]]
		
		if (temp.b < 0 || temp.b == "NaN") {
			temp.b <- 0.001
		}
		
		varinp.mf.n[2, rule.data.num[ii, j]] <- varinp.mf[3, rule.data.num[ii, j]] - (0.5 * temp.b)
		varinp.mf.n[4, rule.data.num[ii, j]] <- varinp.mf[3, rule.data.num[ii, j]] + (0.5 * temp.b)
		
		if (varinp.mf.n[2, rule.data.num[ii, j]] < 0 ) {
			varinp.mf.n[2, rule.data.num[ii, j]] <- 0
		}
		
		if (varinp.mf.n[4, rule.data.num[ii, j]] < 0){
			varinp.mf.n[4, rule.data.num[ii, j]] <- 0
		}
		
		if (varinp.mf.n[2, rule.data.num[ii, j]] > 1 ) {
			varinp.mf.n[2, rule.data.num[ii, j]] <- 1
		}
		
		if (varinp.mf.n[4,rule.data.num[ii, j]] > 1) {
			varinp.mf.n[4, rule.data.num[ii, j]] <- 1
		}		
	}
}


param.new <- list(func.tsk.new = func.tsk.new, varinp.mf.n = varinp.mf.n)
return(param.new)
}

#' The role of this function is to update parameters within the simplified TSK fuzzy rule 
#' generation method using heuristics and the gradient descent method (FS.HGD).
#' This function is called by the main function 
#' of the FS.HGD method, see \code{\link{FS.HGD}}.
#'
#' @title FS.HGD updating function
#'
#' @param data.train a matrix (\eqn{m \times n}) of normalized data for the training process, where \eqn{m} is the number of instances and 
#' \eqn{n} is the number of variables; the last column is the output variable.
#' @param miu.rule a matrix with the degrees of rules which is the result of the \code{\link{inference}}.
#' @param func.tsk a matrix of parameters of the function on the consequent part using the Takagi Sugeno Kang model. See \code{\link{rulebase}}.
#' @param varinp.mf a matrix of parameters of membership functions of the input variables.
#' @param step.size a real number between 0 and 1 representing the step size of the gradient descent. 
#' @param def a matrix which is obtained by the \code{\link{defuzzifier}}.
#'  
#' @seealso \code{\link{FS.HGD}}
# @return a list of the new parameter
# @export
HGD.update <- function(data.train, miu.rule, func.tsk, varinp.mf, step.size = 0.01, def){

data <- data.train
alpha <- step.size

#TSK type

data.k <-(data[1, 1 : (ncol(data.train) - 1)])
data.m <- t(t(data.k))
	
ff <- func.tsk 
miu.rule.t <- as.matrix(miu.rule)

gal <- (def - data.train[1, ncol(data.train)])
	
func.tsk.upd <- func.tsk
for (mm in 1 : nrow(func.tsk.upd)){
	sum.miu <- sum(miu.rule.t, na.rm = TRUE)
	if (sum.miu != 0)
		func.tsk.upd[mm, 1] <- func.tsk.upd[mm, 1] - alpha * gal * (miu.rule[mm])
}


func.tsk.new <- func.tsk.upd

param.new <- list(func.tsk.new)
######################################
#=====end of code gradient descent 
#########################################

return(param.new)
}

# The role of this function is to generate initial rule for frbs based on gradient descent. 
#
# @title generate.rule.GD rule generating function
#
# @param data.train a matrix(m x n) of data for the training process, where m is the number of instances and 
# n is the number of variables; the last column is the output variable.
# @param range.data a matrix(2 x n) containing the range of the normalized data, where n is the number of variables, and
# first and second rows are the minimum and maximum value, respectively. 
# @param num.labels a matrix(1 x n), whose elements represent the number of labels (linguistic terms);
# n is the number of variables.
# @param type.mf a shape of membership functions.
# @param type.tnorm a value which represents the type of t-norm. See \code{\link{inference}}.
# @param type.implication.func a value representing type of implication function. 
generate.rule.GD <- function(range.data, data.train, num.labels, type.mf = "TRIANGLE", type.tnorm = "MIN", type.implication.func = "ZADEH"){

	## initialize mod
	mod <- NULL

	## generate initial model using WM
	mod <- WM(data.train, num.labels, type.mf, type.tnorm, type.implication.func)
	
	## get values of model		
	rule <- mod$rule
	rule.data.num <- mod$rule.data.num		

	## delete duplication on the antecedent part
	ind.nonDuplicate <- which(duplicated(rule[, -ncol(rule)]) == FALSE, arr.ind = TRUE)
	rule <- rule[ind.nonDuplicate, ,drop = FALSE]
	rule.data.num <- rule.data.num[ind.nonDuplicate, ,drop = FALSE]

	mod$rule <- rule
	mod$rule.data.num <- rule.data.num
	
	return(mod)
}
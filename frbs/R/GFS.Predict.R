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
#' This function is the internal function of the GFS.FR.MOGUL method to compute the predicted values.  
#'
#' @title GFS.FR.MOGUL: The prediction phase
#' @param object the \code{\link{frbs-object}}.
#' @param newdata a matrix (\eqn{m \times n}) of data for the prediction process, where \eqn{m} is the number of instances and 
#' \eqn{n} is the number of input variables.
#' @return A matrix of predicted values.
# @export
GFS.FR.MOGUL.test <- function(object, newdata){
	mod <- object
	data.test <- newdata
	rule.gen <- mod$rule
	res <- calc.pred.val(data.test = data.test, rule.gen = rule.gen, method.type = "GFS.FR.MOGUL")
	
	return(res)
}

#' This function is the internal function of the GFS.Thrift method to compute the predicted values.  
#'
#' @title GFS.Thrift: The prediction phase
#' @param object the \code{\link{frbs-object}}.
#' @param newdata a matrix (\eqn{m \times n}) of data for the prediction process, where \eqn{m} is the number of instances and 
#' \eqn{n} is the number of input variables.
#' @return A matrix of predicted values.
# @export
GFS.Thrift.test <- function(object, newdata){
	data.test <- newdata
	rule.gen <- object$rule.data.num
	var.mf <- cbind(object$varinp.mf, object$varout.mf)
	method.type <- object$method.type
	res <- calc.pred.val(data.test = data.test, rule.gen = rule.gen, method.type = method.type, var.mf = var.mf)
	
	return(res)
}


#' This function is the internal function of the GFS.LT.RS method to compute the predicted values.  
#'
#' @title GFS.LT.RS: The prediction phase
#' @param object the \code{\link{frbs-object}}.
#' @param newdata a matrix (\eqn{m \times n}) of data for the prediction process, where \eqn{m} is the number of instances and 
#' \eqn{n} is the number of input variables.
#' @return A matrix of predicted values.
# @export
GFS.LT.RS.test <- function(object, newdata){
	data.test <- newdata
	rule.gen <- object$rule.data.num
	var.mf <- object$var.mf
	method.type <- "GFS.LT.RS"
	mode.tuning <- object$mode.tuning
	var.mf.tune <- object$var.mf.tune
	num.labels <- object$num.labels
	params <- list(var.mf.tune = var.mf.tune, mode.tuning = mode.tuning, num.labels = num.labels)
	res <- calc.pred.val(data.test = data.test, rule.gen = rule.gen, method.type = method.type, var.mf = var.mf, params = params)
	
	return(res)
}

#' This function is the internal function of the GFS.GCCL and FH.GBML method to compute the predicted values.  
#'
#' @title GFS.GCCL.test: The prediction phase
#' @param object the \code{\link{frbs-object}}.
#' @param newdata a matrix (\eqn{m \times n}) of data for the prediction process, where \eqn{m} is the number of instances and 
#' \eqn{n} is the number of input variables.
#' @return A matrix of predicted values.
# @export
GFS.GCCL.eng <- function(object, newdata){
	classifier.rule <- list(rule = object$rule.data.num, grade.cert = object$grade.cert)
	varinp.mf <- object$varinp.mf	
	resClass.train <- GFS.GCCL.test(newdata, classifier.rule, varinp.mf)	
	return(resClass.train)
}

#' This function is the internal function of the SLAVE method to compute the predicted values.  
#'
#' @title SLAVE.test: The prediction phase
#' @param object the \code{\link{frbs-object}}.
#' @param newdata a matrix (\eqn{m \times n}) of data for the prediction process, where \eqn{m} is the number of instances and 
#' \eqn{n} is the number of input variables.
#' @return A matrix of predicted values.
# @export
SLAVE.test <- function(object, newdata){
	rule <- object$rule.data.num
	varinp.mf <- object$varinp.mf
	
	res <- matrix(nrow = nrow(newdata), ncol = 1)
	
	for (i in 1 : nrow(newdata)){
		data.test <- newdata[i, ,drop = FALSE]
		res.temp <- det.class(data.test, rule, grade.cert = NULL, varinp.mf, type = "NON-WEIGHTED")
		res[i, 1] <- res.temp$res.class
	}
	
	return(res)
}
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
#' This is a function that implements the fuzzy-rough instance selection (FRIS) proposed by  
#' (Jensen and Cornelis, 2010) which is used to perform instance selection.
#' 
#' FRIS is used to remove training instances that cause conflicts with other instances as 
#' determined by the fuzzy positive region. 
#'
#' This algorithm evaluates the degree of membership of each instance to the fuzzy positive region.
#' If there is a instance less than the threshold, then the instance can be removed. 
#' Additionally, it uses a fuzzy indiscernibility relation \eqn{R_a} to express the approximate equality between 
#' objects \eqn{x} and \eqn{y} on attribute \eqn{a} in the training set:
#'
#' \eqn{R_{a}^{\alpha}(x,y)=max(0, 1 - \alpha \frac{|a(x) - a(y)|}{l(a)})}
#' 
#' where parameter \eqn{\alpha} (\eqn{\alpha \ge 0}) determines the granularity of \eqn{R_{a}^{\alpha}}.
#' Then, the fuzzy \eqn{B}-indiscernibility relation, fuzzy lower approximation, positive region and degree of dependency are calculated based on 
#' the concept in \code{\link{B.Introduction-FuzzyRoughSets}}. 
#'  
#' It should be noted that this function does not give the new decision table directly. 
#' The other additional function called \code{\link{SF.applyDecTable}} is used to produce the new decision table based on 
#' the output of this function.
#'
#' @title The fuzzy rough instance selection algorithm
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param control a list of other parameters which are 
#'        \itemize{
#'        \item \code{threshold.tau}: a value determining whether an object can be removed or not. 
#'              The object can be removed if it is less than the threshold. The default value is 0.95.
#'        \item \code{alpha}: a parameter determining the granularity of the fuzzy similarity measure, which has positive values (>= 0). The default value is 1.
#'        \item \code{type.aggregation}: a list representing the type of aggregation and its value. 
#'              The default value is \code{type.aggregation = c("t.tnorm", "lukasiewicz")}. See \code{\link{BC.IND.relation.FRST}}.
#'        \item \code{t.implicator}: a string representing the value of implicator function. The default value is \code{"lukasiewicz"}. See \code{\link{BC.LU.approximation.FRST}}.
#'        }
#' @seealso \code{\link{IS.FRPS.FRST}}.
#' @return A class \code{"InstanceSelection"} that contains the following components:
#' \itemize{
#' \item \code{indx.objects}: it contains all indexes of objects that are selected. 
#' \item \code{type.method}: a string representing the type of method. In this case, it is \code{"IS.FRIS"}.
#' \item \code{type.task}: a string showing the type of task which is \code{"instance selection"}.
#' \item \code{type.model}: a string representing the type of model which is \code{"FRST"}. 
#' } 
#' @examples
#' #############################################
#' ## Example: Evaluate instances/objects and
#' ##          generate new decision table
#' #############################################
#' dt.ex1 <- data.frame(c(0.1, 0.5, 0.2, 0.3, 0.2, 0.2, 0.8), 
#'                   c(0.1, 0.1, 0.4, 0.2, 0.4, 0.4, 0.5), c(0, 0, 0, 0, 1, 1, 1))
#' colnames(dt.ex1) <- c("a1", "a2", "d")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 3, indx.nominal = c(3))
#'
#' ## evaluate index of objects
#' res.1 <- IS.FRIS.FRST(decision.table = decision.table, control = 
#'                         list(threshold.tau = 0.5, alpha = 0.8, 
#'                         type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                         t.implicator = "lukasiewicz"))
#' 
#' ## generate new decision table
#' new.decTable <- SF.applyDecTable(decision.table, res.1)
#'
#' @references
#' R. Jensen and C. Cornelis, "Fuzzy-rough Instance Selection",
#' Proceedings of the 19th International Conference on Fuzzy Systems (FUZZ-IEEE 2010), 
#' p. 1776 - 1782 (2010).
#' @export
IS.FRIS.FRST <- function(decision.table, control = list()){
	## set default values of all parameters
	control <- setDefaultParametersIfMissing(control, list(threshold.tau = 0.95, type.aggregation = c("t.tnorm", "lukasiewicz"),
                                     	t.implicator = "kleene_dienes", alpha = 1))

	## get parameters
	objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	num.att <- ncol(objects)

	## get the data
	threshold.tau <- control$threshold.tau
	num.instances <- nrow(objects)
	alpha <- control$alpha
	type.aggregation <- control$type.aggregation
	t.implicator <- control$t.implicator
	t.similarity = "eq.instance.selection"

	## Indiscernibility for all attributes
	attributes <- c(seq(1, num.att - 1))
	control.ind <- list(type.aggregation = type.aggregation, type.relation = c("tolerance", t.similarity), alpha = alpha)
	IND.all <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.ind)$IND.relation 

	## Indiscernibility for decision attribute
	attributes <- c(num.att)
	control.ind <- list(type.aggregation = type.aggregation, type.relation = c("tolerance", t.similarity), alpha = alpha)
	IND.dec.table <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.ind)$IND.relation 
	
	temp.res <- do.call(calc.implFunc, list(IND.all, IND.dec.table, t.implicator))
	pos <- matrix(apply(temp.res, 1, function(x) min(x)), ncol = 1)

	## create class of InstanceSelection
	indx.objects <- which(pos >= threshold.tau)
	mod <- list(indx.objects = indx.objects, type.method = "IS.FRIS", type.task = "instance selection", type.model = "FRST")					
	class.mod <- ObjectFactory(mod, classname = "InstanceSelection")	
		
	return(class.mod)
}


#' This is a function for implementing instance selection using prototype selection method (FRPS) proposed by (Verbiest et al, 2013).
#' 
#' This algorithm uses prototype selection (PS) to improve the accuracy of the k-nearest neighbour (kNN) method. 
#' It selects a subset of instances \eqn{S \subseteq X} and then classifies a new instance \eqn{t} using the kNN rule acting over
#' \eqn{S} instead of over \eqn{X}. Based on fuzzy rough set theory, \eqn{S} is built. Basically, two steps are performed in the algorithm.
#' First, the instances are ordered according to a measure called \eqn{\alpha} which is based on fuzzy rough set theory that evaluates the lack of 
#' predictive ability of the instances, and the instances for which the values exceeds a certain threshold are removed from training set.
#' To determine this threshold, we consider several equations which are applied on all instances. 
#' These equations are constructed by considering fuzzy positive region membership degree on several implication and t-norm operators.
#' The following is the equations of \eqn{\alpha}:
#' \itemize{
#' \item \code{"FRPS.1"}: \eqn{max_{y \notin [x]_{d}} \frac{1}{max_{i=1}^{m} \delta_{a_i}(x,y)}}
#' \item \code{"FRPS.2"}: \eqn{OW A_w \left(\frac{1}{OW A_w \delta_{a_i}(x,y)}\right)}
#' \item \code{"FRPS.3"}: \eqn{max_{y \notin [x]_{d}} \frac{1}{\displaystyle\sum\limits_{i=1}^m {\delta_{a_i}(x,y)}}}
#' \item \code{"FRPS.4"}: \eqn{OW A_w \left(\frac{1}{\displaystyle\sum\limits_{i=1}^m {\delta_{a_i}(x,y)}}\right)}
#' }
#' where \eqn{[x]_d} and \eqn{OW A_w} are equivalence class on the decision attribute and ordered weighted average operator, respectively.
#' And \eqn{\delta_{a_i}(x,y)} is a distance measure of the attribute \eqn{a_i} between \eqn{x} and \eqn{y}.
#' After that, \code{1knn} will be performed in order to select and get prototype \eqn{S} which produces the highest training accuracy.
#'
#' It should be noted that this function does not give the new decision table directly. 
#' The other additional function called \code{\link{SF.applyDecTable}} is used to produce new decision table based on 
#' the output of this function.
#'
#' @title The fuzzy rough prototype selection method
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param type.alpha type of FRPS expressing the equations of \eqn{\alpha}. The default value is \code{"FRPS.1"}.
#'        See Section \code{Details}.
#' @seealso \code{\link{IS.FRIS.FRST}}
#' @return A class \code{"InstanceSelection"} that contains the following components:
#' \itemize{
#' \item \code{indx.objects}: it contains all indexes of objects that are selected. 
#' \item \code{type.method}: a string representing a type of method. In this case, it is \code{"IS.FRPS"}.
#' \item \code{type.task}: a string showing the type of task which is \code{"instance selection"}.
#' \item \code{type.model}: a string representing the type of model which is \code{"FRST"}. It means fuzzy rough set theory.
#' }  
#' @examples
#' #############################################
#' ## Example: Evaluate instances/objects and
#' ## generate new decision table
#' #############################################
#' dt.ex1 <- data.frame(c(0.5, 0.2, 0.3, 0.7, 0.2, 0.2), c(0.1, 0.4, 0.2, 0.8, 0.4, 0.4), 
#'                       c(0, 0, 0, 1, 1, 1))
#' colnames(dt.ex1) <- c("a1", "a2", "d")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 3, indx.nominal = c(3))
#'
#' ## evaluate instances
#' res.1 <- IS.FRPS.FRST(decision.table, type.alpha = "FRPS.3")
#'
#' ## generate new decision table
#' new.decTable <- SF.applyDecTable(decision.table, res.1)
#' @references
#' N. Verbiest, C. Cornelis, and F. Herrera, "A Fuzzy Rough Prototype Selection Method", 
#' Pattern Recognition, Vol. 46, no. 10, p. 2770 - 2782 (2013).
#'
#' @export
IS.FRPS.FRST <- function(decision.table, type.alpha = "FRPS.1"){
	req.suc <- requireNamespace("class", quietly=TRUE)
	if(!req.suc) stop("In order to use this function, you need to install the package class.")

	## get parameters
	objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	num.att <- ncol(objects)

	## get range data from decision.table
	range.data.input <- matrix(unlist(desc.attrs[-length(desc.attrs)]), nrow = 2)
	
	## normalize data
	temp <- norm.data(objects[, -ncol(objects), drop = FALSE], range.data.input, min.scale = 0, max.scale = 1)
	objects <- cbind(temp, objects[, ncol(objects), drop = FALSE])
	
	## calculate alpha
	alpha <- matrix()
	for (i in 1 : nrow(objects)){
		## get one instance
		object.i <- objects[i, ,drop = FALSE]
		
		## get the condition attribute of object i
		cond.object.i <- object.i[1, -ncol(object.i), drop = FALSE]
		
		## get discernibility of decision attribute
		non.equivalence.cls <- which(objects[, ncol(objects)] != object.i[1, ncol(object.i)])
		
		all.temp <- matrix()
		for (j in non.equivalence.cls){
			## get condition attribute of object j
			cond.object.j <- objects[j, -ncol(objects), drop = FALSE]
			
			## combine object i and j
			val <- matrix()
			for (k in 1 : ncol(cond.object.i)){
					## check whether it is nominal attribute or not
					if (nominal.att[k] == FALSE){
						## non nominal attribtue
						val[k] <- (cond.object.i[1, k] - cond.object.j[1, k])^2
					} else {
						## for nominal attribute
						if (cond.object.i[1, k] != cond.object.j[1, k])
							val[k] <- 1
						else 
							val[k] <- 0
					}
			}
			
			##	check equation alpha
			if (type.alpha == "FRPS.1"){
				temp = c(1 / max(val))
				all.temp <- cbind(all.temp, temp)				
			}
			else if (type.alpha == "FRPS.2"){
				temp = calc.OWA(data = val, type.method = "weight")
				all.temp <- cbind(all.temp, temp)
			}
			else if (any(type.alpha == c("FRPS.3", "FRPS.4"))){
				temp = c(1 / sum(val))
				all.temp <- cbind(all.temp, temp)
			}			
		}
		## delete the first position
		all.temp <- all.temp[-1]
		
		## get max value
		if (any(type.alpha == c("FRPS.1", "FRPS.3"))){
			alpha[i] <- max(all.temp)
		}
		else if (any(type.alpha == c("FRPS.2", "FRPS.4"))){
			alpha[i] <- calc.OWA(data = all.temp, type.method = "weight")
		}
	}
	
	## combine alpha with sequential number in order to define the index of objects
	alpha <- cbind(seq(1, nrow(objects)), alpha)
	res.i <- c()
	rate <- 0
	for (i in 1 : nrow(objects)){
		## perform 1knn
		res <- class::knn1(train= objects[-i, -ncol(objects), drop = FALSE], 
		                 test = objects[i, -ncol(objects), drop = FALSE], 
						 cl= factor(objects[-i, ncol(objects)]))
		## check accuracy
		if (res == objects[i, ncol(objects)])
			rate <- rate + 1
	}
	
	## calculate total accuracy
	rate.init <- rate/nrow(objects)
	
	## collect current alpha 
	list.alpha <- list(alpha[which(alpha[, 2] == max(alpha[, 2])), 1])
	
	## initialization
	new.objects <- objects
	exit <- FALSE
	while (exit == FALSE){
		## delete some objects with max current alpha
		new.objects <- new.objects[-which(alpha[, 2] == max(alpha[, 2])), ,drop = FALSE]
	    
		## delete max current alpha
		alpha <- alpha[-which(alpha[, 2] == max(alpha[, 2])), ,drop = FALSE]
		
		## check if the number of objects over than 1
		if (nrow(new.objects) > 1){
			## perform 1knn
			res.knn <- class::knn1(train = new.objects[, -ncol(new.objects), drop = FALSE], 
			                test = objects[, -ncol(objects), drop = FALSE], 
							cl= factor(new.objects[, ncol(new.objects)]))
			## calculate accuracy
			rate <- sum(res.knn == objects[, ncol(objects)])/nrow(objects)
			
			## check if we got higher accuracy
			if (rate > rate.init){
				## replace alpha with this new one
				list.alpha <- c()
				list.alpha <- c(alpha[which(alpha[, 2] == max(alpha[, 2])), 1])
				rate.init <- rate
			}			
			## check if we got the same accuracy
			else if (rate == rate.init){
				## add new alpha
				list.alpha <- append(list.alpha, c(alpha[which(alpha[, 2] == max(alpha[, 2])), 1]))
				rate.init <- rate
			}
		}
		else if (nrow(new.objects) == 1){
			## add the last object
			list.alpha <- append(list.alpha, c(alpha[, 1]))
			exit <- TRUE
		} else {
			exit <- TRUE
		}
	}
	## sort the final alpha
	list.alpha <- sort(as.vector(unlist(list.alpha)))
	
	## construct class
	mod <- list(indx.objects = list.alpha, type.method = "IS.FRPS", type.task = "instance selection", type.model = "FRST")					
	class.mod <- ObjectFactory(mod, classname = "InstanceSelection")	
		
	return(class.mod)
}
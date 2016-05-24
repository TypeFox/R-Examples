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

# This function is to get k nearest neighbor
#
# @param decision.table data are used for training
# @param newdata data are used for testing
# @param k number of neighbor
# @param method.type a type of method used to calculate distance
get.knearest <- function(decision.table, newdata.i, k, method.type){
	nearest.dt <- matrix(nrow = 1, ncol = (ncol(decision.table) + 1))

	if (k == nrow(decision.table)){
		dist.dt <- apply(decision.table[, -ncol(decision.table)],1,function(x)sqrt(sum((x-newdata.i)^2)))
		nearest.dt <- rbind(nearest.dt, cbind(decision.table, dist.dt))	
		
	} else {
		for (j in 1 : nrow(decision.table)){
			## calculate distance based on method.type
			dist.dt <- stats::dist(rbind(decision.table[j, -ncol(decision.table), drop = FALSE], newdata.i), method = method.type)		
		
			## check if less than k (note. k + 1 because first row is NA
			if (nrow(nearest.dt) <= (k + 1)){
				nearest.dt <- rbind(nearest.dt, cbind(decision.table[j, ,drop = FALSE], dist.dt))
			} 
			
			## replace with lower value
			else if (max(nearest.dt[, ncol(nearest.dt)], na.rm = TRUE) > dist.dt){
				nearest.dt[which.max(nearest.dt[, ncol(nearest.dt)]), ] <- cbind(decision.table[j, ,drop = FALSE], dist.dt)
			}
		}
	}
	## delete first row
	nearest.dt <- stats::na.omit(nearest.dt)		
	
	return(nearest.dt)
}

# This function is to calculate degree of similarity
#
# @param nearest.dt a matrix of data of nearest neighbor
# @param control a list of parameters
calc.similarityDegree <- function(nearest.dt, control){	
	num.class <- control$num.class
	m <- control$m
	
	## calculate membership of newdata to each class
	## miu.class is a matrix (k, num.class) where k is number of neighbor
	miu.class <- calc.membershipClass(nearest.dt, num.class, type = "gradual")
	
	## calculate sum on denominator of the similarity eq.
	sum.denominator <- sum(1 / (nearest.dt[, ncol(nearest.dt)] ^ (2 / (m - 1))))
	
	## calculate membership function of class for determining class of newdata
	miu <- t(miu.class) %*% (1 / (nearest.dt[, ncol(nearest.dt), drop = FALSE] ^ (2 / (m - 1)))) / sum.denominator
	
	return(miu)
}

# This function is to calculate membership class of nearest neighbor 
#
# @param nearest.dt a matrix of data of nearest neighbor
# @param num.class the number of class
# @param control a list of parameters
calc.membershipClass <- function(nearest.dt, num.class, type = "gradual"){

	numClass <- matrix(nrow = 1, ncol = num.class)
	
	## get sum of class
	for (i in 1 : num.class){
		numClass[1, i] <- nrow(nearest.dt[which(nearest.dt[,(ncol(nearest.dt) - 1)] == i), ,drop = FALSE])
	}		
	
	miu.class <- matrix(nrow = nrow(nearest.dt), ncol = num.class)
	
	## calculate membership of each class based on Keller et. al.'s technique.
	for (i in 1 : nrow(nearest.dt)){
		for (j in 1 : num.class){
			if (nearest.dt[i, ncol(nearest.dt) - 1] == j) {
				if (type == "gradual") {
					miu.class[i, j] <- 0.51 + (numClass[1, nearest.dt[i, ncol(nearest.dt) - 1]]/nrow(nearest.dt) * 0.49)
				}
				else if (type == "crisp"){
					miu.class[i, j] <- 1
				}
			} else {
				if (type == "gradual"){
					miu.class[i, j] <- (numClass[1, j]/nrow(nearest.dt) * 0.49)
				}
				else if (type == "crisp"){
					miu.class[i, j] <- 0
				}
			}
		}
	}

	
	## Note: sum of miu.class must be 1
	return(miu.class)
}

# This is a function that implements the fuzzy-rough nearest neighbour algorithm.
# 
# @title The fuzzy-rough nearest neighbour algorithm
#
# @param decision.table a data frame representing the decision table. See \code{\link{BC.IND.relation.FRST}}. 
# @param newdata a data frame or matrix (p x n) of data for the test process, where p is the number of instances and 
#        n is the number of conditional attributes (input variables); 
# @param control a list of other parameters which depends on two approaches we are going to choose, as follows:
#        \itemize{
#        \item implicator/t-norm based model: we should set \code{type.LU = "implicator.tnorm"} and following parameters.
#                \itemize{
#                  \item \code{k}: the number of neighbours. The default value is 5.
#                  \item \code{type.aggregation}: the type of aggregation.
#                  \item \code{type.relation}: the type of relation. The default value is \code{c("tolerance", "eq.1")}.
#                  \item \code{type.implicator}: the type of implicator operator. The default value is "lukasiewicz".
#                 }
#        \item vaquely quantified rough set model: for using this type, we should set \code{type.LU = "vqrs"} and following parameters.
#                 \itemize{
#                  \item \code{k}:  the number of neighbor. The default value is 5.
#                  \item \code{q.some}: a numeric of alpha and beta parameter of VQRS. The default value is \code{c(0.1, 0.6)}.
#                  \item \code{q.most}: a numeric of alpha and beta parameter of VQRS. The default value is \code{c(0.2, 1)}.
#                  \item \code{type.relation}:  the type of relation. The default value is \code{c("tolerance", "eq.1")}.
#                  \item \code{type.aggregation}: the type of aggregation.
#                  } 
#        }
#        The description of those parameters could be seen at \code{\link{BC.LU.approximation.FRST}}.
# @seealso \code{\link{C.FRNN.O.FRST}}, 
# \code{\link{C.POSNN.FRST}} 
# @return a matrix of predicted classes of newdata. 
# @references
# Richard Jensen and Chris Cornelis, "Fuzzy-rough Nearest Neighbour Classificition and Prediction",
# Theoretical Computer Science, vol. 412, p. 5871 - 5884 (2011). 
FRNN.alg <- function(decision.table, newdata, type.method = "C.FRNN.FRST", control = list()){
	## get the data
	objects <- as.matrix(decision.table)
	newdata <- as.matrix(newdata)
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	num.att <- ncol(objects)
	classes.dt <- as.numeric(unlist(desc.attrs[length(desc.attrs)]))
	
	## check default parameters 
	control <- setDefaultParametersIfMissing(control, list(k = 5, type.aggregation = c("t.tnorm", "lukasiewicz"), 
	                 type.relation = c("tolerance", "eq.1"), t.implicator = "lukasiewicz", 
					 q.some = c(0.1, 0.6), q.most = c(0.2, 1), type.LU = "implicator.tnorm"))
	
	## get range data from decision.table
	range.data.input <- matrix(unlist(desc.attrs[-length(desc.attrs)]), nrow = 2)
	
	## get data
	k <- control$k
	method.type <- "euclidean"
	type.aggregation <- control$type.aggregation
	t.tnorm <- type.aggregation[2]
	type.relation <- control$type.relation
	t.implicator <- control$t.implicator
	type.LU <- control$type.LU
	q.some <- control$q.some
	q.most <- control$q.most
	num.inputvar <- num.att - 1
	P.attributes <- c(seq(1, num.inputvar))
	res.class <- matrix(nrow = nrow(newdata), ncol = 1)
	
	## iterate for each testing data
	for (ii in 1 : nrow (newdata)){
		## calculate and get K-nearest neighbor
		nearest.dt <- get.knearest(decision.table = objects, newdata.i = newdata[ii, ,drop = FALSE], k = k, method.type = method.type)
		objects.knn <- nearest.dt[, -ncol(nearest.dt), drop = FALSE]

		threshold.val <- 0
		
		## iterate over classes
		for (i in 1 : length(classes.dt)) {	
			## put testing data into decision table (training data)
			decision.table <- rbind(cbind(newdata[ii, ,drop = FALSE], classes.dt[i]), objects.knn)
			colnames(decision.table) <- names(desc.attrs)
			attr(decision.table, "nominal.attrs") = nominal.att
			attr(decision.table, "desc.attrs") = desc.attrs

			## calculate indiscernibility relation over conditional attributes and the decision attribute
			control.ind <- list(type.aggregation = type.aggregation, type.relation = type.relation)
			obj.IND <- BC.IND.relation.FRST(decision.table, attributes = P.attributes, control = control.ind)
			IND <- obj.IND$IND.relation
			obj.IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = c(num.att), control = control.ind)

			## Perform fuzzy lower and upper approximation using type.LU : "implicator.tnorm" and "vqrs"
			decision.attr = c(num.att)
			if (type.LU == "implicator.tnorm"){
				control <- list(t.implicator = t.implicator, t.tnorm = t.tnorm)				
			}
			else if (type.LU == "vqrs"){
				control <- list(q.some = q.some, q.most = q.most, t.tnorm = t.tnorm)
			}			
			LU <- BC.LU.approximation.FRST(decision.table, obj.IND, obj.IND.decAttr, 
					type.LU = type.LU, control = control)

			## check the method that is going to be used
			if (type.method == "C.FRNN.FRST"){						
				fuzzy.lower <- LU$fuzzy.lower[[1]][1]
				fuzzy.upper <- LU$fuzzy.upper[[1]][1]				
				new.value <- (fuzzy.lower + fuzzy.upper)/2				
			}
			else if (type.method == "C.POSNN.FRST"){
				#### Perform fuzzy indiscernibility relation ####
				control.ind <- list(type.aggregation = type.aggregation, type.relation = type.relation)
				IND.dec <- BC.IND.relation.FRST(decision.table, attributes = c(num.att), control = control.ind)$IND.relation				
				pos <- BC.positive.reg.FRST(decision.table, LU)$positive.freg
				new.value <- sum(IND.dec[1, ,drop = FALSE] * pos * IND[1, ,drop = FALSE])

			}
			
			## checking threshold value
			if (new.value >= threshold.val){
					res.class[ii] <- classes.dt[i]
					threshold.val <- new.value
			}
		}	
	}	
	
	return(res.class)
}


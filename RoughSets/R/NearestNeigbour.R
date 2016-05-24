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
#' It is used to predict classes of new datasets/patterns based on the fuzzy-rough ownership nearest neighbor algorithm (FRNN.O) 
#' proposed by (Sarkar, 2007). 
#' 
#' This method improves fuzzy \eqn{k}-nearest neighbors (FNN) by introducing rough sets into it. 
#' To avoid determining \eqn{k} by trial and error procedure, this method uses all training data. Uncertainties in data are accommodated by 
#' introducing the rough ownership function. It is the following equation \eqn{o_c} of each class expressing a squared weighted distance between a test pattern and all training data \eqn{d}
#' and constrained fuzzy membership \eqn{\mu_{C_c}}. 
#'
#' \eqn{o_c(y) = \frac{1}{|X|}\mu_{C_c}(x)\exp{(-d^{1/(q-1)})}}
#'
#' where \eqn{d = \sum_{j=1}^{N}K_j(y_j-x_{ij})^2}
#'
#' The predicted value of \eqn{y} is obtained by selecting class \eqn{c} where \eqn{o_c(y)} is maximum.
#'  
#' @title The fuzzy-rough ownership nearest neighbor algorithm
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}.
#'        It should be noted that the data must be numeric values instead of strings/characters.  
#' @param newdata a \code{"DecisionTable"} class representing data for the test process. 
#'
#'        See \code{\link{SF.asDecisionTable}}. 
#' @param control a list of other parameters.
#'        \itemize{
#'        \item \code{m}: the weight of distance. The default value is 2.
#'        }
#' @seealso \code{\link{C.FRNN.FRST}},  
#' \code{\link{C.POSNN.FRST}} 
#' @return A matrix of predicted classes of newdata. 
#' @examples
#' #############################################################
#' ## In this example, we are using Iris dataset.
#' ## It should be noted that since the values of the decision attribute are strings,
#' ## they should be transformed into numeric values using unclass()
#' #############################################################
#' data(iris)
#' ## shuffle the data
#' set.seed(2)
#' irisShuffled <- iris[sample(nrow(iris)),]
#'
#' ## transform values of the decision attribute into numerics
#' irisShuffled[,5] <- unclass(irisShuffled[,5])
#'
#' ## split the data into training and testing data
#' iris.training <- irisShuffled[1:105,]
#' iris.testing <- irisShuffled[106:nrow(irisShuffled),1:4]
#'
#' ## convert into the standard decision table
#' colnames(iris.training) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", 
#'                              "Petal.Width", "Species")
#' decision.table <- SF.asDecisionTable(dataset = iris.training, decision.attr = 5, 
#'                                     indx.nominal = c(5))
#' tst.iris <- SF.asDecisionTable(dataset = iris.testing)
#'
#' ## in this case, we are using "gradual" for type of membership					   
#' control <- list(m = 2)
#'
#' \dontrun{res.test.FRNN.O <- C.FRNN.O.FRST(decision.table = decision.table, newdata = tst.iris, 
#'                                  control = control)}
#'
#' @references
#' M. Sarkar, "Fuzzy-Rough Nearest-Neighbor Algorithm in Classification" 
#' Fuzzy Sets and Systems, vol. 158, no. 19, p. 2123 - 2152 (2007).
#'
#' @export
C.FRNN.O.FRST <- function(decision.table, newdata, control = list()){
	
	control <- setDefaultParametersIfMissing(control, list(m = 2, type.membership = "gradual"))
	objects <- as.matrix(decision.table)
	newdata <- as.matrix(newdata)
	
	## get the data
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	num.att <- ncol(objects)
	m <- control$m
	method.type <- "euclidean"
	type <- control$type.membership
	
	## get range data from decision.table
	range.data.input <- matrix(unlist(desc.attrs[-length(desc.attrs)]), nrow = 2)	
	num.class <- length(unlist(desc.attrs[length(desc.attrs)]))	
	objects.norm <- norm.data(objects[, -ncol(objects), drop = FALSE], range.data.input, min.scale = 0, max.scale = 1)
	objects <- cbind(objects.norm, objects[, ncol(objects), drop = FALSE])	
	newdata <- norm.data(newdata, range.data.input, min.scale = 0, max.scale = 1)

	num.inputvar <- ncol(objects) - 1
	num.instances <- nrow(objects)
	k.a <- matrix(nrow = 1, ncol = num.inputvar)
	data.inputvar <- objects[, -ncol(objects), drop = FALSE]
	res.class <- matrix(nrow = nrow(newdata), ncol = 1)
	
	for (ii in 1 : nrow(newdata)){

		distance <- apply(data.inputvar,1,function(x)sqrt(sum((x-newdata[ii, ,drop = FALSE])^2)))
		k.averDist <- 1 / ((1 / num.instances) * sum (distance ^ 2 /(m - 1))) 
		
		## calculate and get K-nearest neighbor (in this case, K = num.instances)
		nearest.dt <- get.knearest(decision.table = objects, 
		                           newdata.i = newdata[ii, ,drop = FALSE],
                                   k = num.instances, method.type = method.type)
		
		## change distance = 0 into 0.000001 to avoid devide by zero
		nearest.dt[nearest.dt[, ncol(nearest.dt)] == 0, ncol(nearest.dt)] <- 0.000001		
		
		## calculate membership of each class based on nearest.dt (in this case: all training data)
		miu <- calc.similarityDegree(nearest.dt, control = list(m = m, num.class = num.class))

		 ## calculate fuzzy-rough ownership function
		tau.C <- matrix(nrow = 1, ncol = num.class)		 
		 for (i in 1 : num.class){
			temp = 0
			for (j in 1 : num.instances){
				distance.i <- apply(data.inputvar[j, ,drop = FALSE],1,function(x)sqrt(sum((x-newdata[ii, ,drop = FALSE])^2)))
				temp = temp + miu[i, 1] / (1 + k.averDist * distance.i ^ (2 / (m - 1)))				
			}
				
			tau.C[1, i] <- (1 / num.instances) * temp
		 }
		
		## choose class which has max membership function		
		res.class[ii, 1] <- which.max(tau.C)
	}
	
	return(res.class)
}

#' It is used to predict new datasets/patterns based on the fuzzy-rough nearest neighbor algorithm (FRNN)
#' proposed by (Jensen and Cornelis, 2011).  
#' 
#' This method uses the fuzzy lower and upper approximations to improve the fuzzy nearest neighbor (FNN) algorithm.
#' This algorithm assigns a class to a target instance \eqn{t} as follows.
#' \itemize{
#' \item Determine \eqn{k} nearest neighbors considering their similarity to new patterns.
#' \item Assign new patterns to the class based on maximal value of fuzzy lower and upper approximations. 
#' If a value of fuzzy lower approximation is high, it shows that neighbors of newdata belong to a particular class, e.g. \code{C}. On the other hand, 
#' a high value of fuzzy upper approximation means that at least one neighbor belongs to that class. 
#' }  
#'
#' In this function, we provide two approaches based on types of fuzzy lower and upper approximations. The following is 
#' a list of the considered approximations:
#' \itemize{
#' \item \code{"implicator.tnorm"}: It refers to lower and upper approximations based on implicator/t-norm approach. 
#'       For more detail, it can be seen in \code{\link{BC.LU.approximation.FRST}}. When using this approach, 
#'       we need to assign the \code{control} parameter as follows:
#'
#'       \code{control <- list(type.LU = "implicator.tnorm", k,}
#' 
#'                        \code{type.aggregation, type.relation, t.implicator)}
#'
#'      The detailed description of the components in the \code{control} parameter can be seen in 
#'
#'      \code{\link{BC.LU.approximation.FRST}}.
#'
#' \item \code{"vqrs"}: It refers to lower and upper approximations based on vaguely quantified rough sets. 
#'       For more detail, it can be seen in \code{\link{BC.LU.approximation.FRST}}. When using this approach, 
#'       we need to assign the \code{control} parameter as follows:
#' 
#' 		\code{control <- list(type.LU = "vqrs", k, q.some, q.most,}
#'
#'                       \code{type.relation, type.aggregation)}
#'
#'      The detailed description of the components in the \code{control} parameter can be seen in 
#'
#'      \code{\link{BC.LU.approximation.FRST}}.
#' } 
#'
#' @title The fuzzy-rough nearest neighbor algorithm
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#'        It should be noted that the data must be numeric values instead of string/char. 
#' @param newdata a \code{"DecisionTable"} class representing data for the test process. 
#'
#'        See \code{\link{SF.asDecisionTable}}. 
#' @param control a list of other parameters as follows.
#'        \itemize{
#'        \item \code{type.LU}: a type of lower and upper approximations. See Section \code{Details}. The default value is \code{type.LU = "implicator.tnorm"}.
#'        \item \code{k}: the number of neighbors. It should be taken into account that 
#'                        this value could affect the accuracy. The default value is 5.
#'        \item \code{type.aggregation}: the type of the aggregation operator. See \code{\link{BC.IND.relation.FRST}}.
#'                   The default value is \code{type.aggregation = c("t.tnorm", "lukasiewicz")}.
#'        \item \code{type.relation}: the type of relation. See \code{\link{BC.LU.approximation.FRST}}.
#'
#'              The default value is \code{c("tolerance", "eq.1")}.
#'        \item \code{type.implicator}: the type of implicator operator. 
#'
#'             See \code{\link{BC.LU.approximation.FRST}}. The default value is \code{"lukasiewicz"}.
#'        \item \code{q.some}: a vector of values of alpha and beta parameters of VQRS. 
#'
#'             See \code{\link{BC.LU.approximation.FRST}}. The default value is \code{c(0.1, 0.6)}.
#'        \item \code{q.most}: a vector of values of alpha and beta parameter of VQRS. 
#'
#'             See \code{\link{BC.LU.approximation.FRST}}. The default value is \code{c(0.2, 1)}.
#'        }
#'
#' @seealso \code{\link{C.FRNN.O.FRST}}, 
#' \code{\link{C.POSNN.FRST}} 
#' @return A matrix of predicted classes of newdata. 
#' @examples
#' #############################################################
#' ## In this example, we are using Iris dataset.
#' ## It should be noted that since the values of the decision attribute are strings,
#' ## they should be transformed into numeric values using unclass()
#' #############################################################
#' data(iris)
#' ## shuffle the data
#' set.seed(2)
#' irisShuffled <- iris[sample(nrow(iris)),]
#'
#' ## transform values of the decision attribute into numerics
#' irisShuffled[,5] <- unclass(irisShuffled[,5])
#'
#' ## split the data into training and testing data
#' iris.training <- irisShuffled[1:105,]
#' iris.testing <- irisShuffled[106:nrow(irisShuffled),1:4]
#' 
#' colnames(iris.training) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", 
#'                        "Petal.Width", "Species")
#'
#' ## convert into a standard decision table
#' decision.table <- SF.asDecisionTable(dataset = iris.training, decision.attr = 5, 
#'                                      indx.nominal = c(5))
#' tst.iris <- SF.asDecisionTable(dataset = iris.testing)
#'
#' ###### FRNN algorithm using lower/upper approximation: 
#' ###### Implicator/tnorm based approach
#' control <- list(type.LU = "implicator.tnorm", k = 20, 
#'                 type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                 type.relation = c("tolerance", "eq.1"), t.implicator = "lukasiewicz") 									   
#' \dontrun{res.1 <- C.FRNN.FRST(decision.table = decision.table, newdata = tst.iris,
#'                              control = control)}
#'
#' ###### FRNN algorithm using VQRS
#' control <- list(type.LU = "vqrs", k = 20, q.some = c(0.1, 0.6), q.most = c(0.2, 1), 
#'                  type.relation = c("tolerance", "eq.1"), 
#'                  type.aggregation = c("t.tnorm","lukasiewicz"))
#' \dontrun{res.2 <- C.FRNN.FRST(decision.table = decision.table, newdata = tst.iris,
#'                              control = control)}
#' 
#' @references
#' R. Jensen and C. Cornelis, "Fuzzy-rough Nearest Neighbour Classification and Prediction",
#' Theoretical Computer Science, vol. 412, p. 5871 - 5884 (2011). 
#'
#' @export
C.FRNN.FRST <- function(decision.table, newdata, control = list()){
	res.class <- FRNN.alg(decision.table, type.method = "C.FRNN.FRST", newdata, control = control)
	return(res.class)
}

#' It is a function used to implement the positive region based fuzzy-rough nearest neighbor algorithm (POSNN)
#' which was proposed by (Verbiest et al, 2012) for predicting classes of new data. 
#' 
#' This method is aimed to improve the fuzzy-rough nearest neighbor algorithm (\code{\link{C.FRNN.FRST}}) algorithm by considering the fuzzy positive region. 
#' Basically the following steps are used to classify an instance \eqn{t}:
#' \itemize{
#'	\item determine the set of \eqn{k}-nearest neighbor of \eqn{t}, \eqn{NN}.
#'  \item assign \eqn{t} to the class \eqn{C} for which
#'
#' \eqn{\frac{\displaystyle\sum\limits_{x \in NN} R(x,t)C(x)POS(x)}{\displaystyle\sum\limits_{x \in NN} R(x,t)}}
#'
#' is maximal.
#' }
#' 
#' @title The positive region based fuzzy-rough nearest neighbor algorithm 
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#'        It should be noted that the data must be numeric values instead of string/char. 
#' @param newdata a \code{"DecisionTable"} class representing data for the test process. 
#'
#'        See \code{\link{SF.asDecisionTable}}. 
#' @param control a list of other parameters which is the same as \code{\link{C.FRNN.FRST}}.
#' @seealso \code{\link{C.FRNN.FRST}}, \code{\link{C.FRNN.O.FRST}}
#' @return A matrix of predicted classes of newdata. 
#' @examples
#' #############################################################
#' ## In this example, we are using Iris dataset.
#' ## It should be noted that since the values of the decision attribute are strings,
#' ## they should be transformed into numeric values using unclass()
#' #############################################################
#' data(iris)
#' ## shuffle the data
#' set.seed(2) 
#' irisShuffled <- iris[sample(nrow(iris)),]
#'
#' ## transform values of the decision attribute into numerics
#' irisShuffled[,5] <- unclass(irisShuffled[,5])
#'
#' ## split the data into training and testing data
#' iris.training <- irisShuffled[1:105,]
#' iris.testing <- irisShuffled[106:nrow(irisShuffled),1:4]
#' 
#' colnames(iris.training) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", 
#'                        "Petal.Width", "Species")
#'
#' ## convert into the standard decision table
#' decision.table <- SF.asDecisionTable(dataset = iris.training, decision.attr = 5, 
#'                                      indx.nominal = c(5))
#' tst.iris <- SF.asDecisionTable(dataset = iris.testing)
#'	   
#' ## FRNN algorithm using lower/upper approximation: Implicator/tnorm based approach
#' control <- list(type.LU = "implicator.tnorm", k = 20, t.tnorm = "lukasiewicz", 
#'                 type.relation = c("tolerance", "eq.1"), t.implicator = "lukasiewicz")
#'
#' \dontrun{res.test.POSNN <- C.POSNN.FRST(decision.table = decision.table, 
#'                               newdata = tst.iris, control = control)}
#'
#' @references
#' N. Verbiest, C. Cornelis and R. Jensen, "Fuzzy-rough Positive Region Based Nearest Neighbour Classification",
#' In Proceedings of the 20th International Conference on Fuzzy Systems (FUZZ-IEEE 2012), p. 1961 - 1967 (2012). 
#'
#' @export
C.POSNN.FRST <- function(decision.table, newdata, control = list()){
	res.class <- FRNN.alg(decision.table, newdata, type.method = "C.POSNN.FRST", control = control)	
	return(res.class)
}


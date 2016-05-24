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
#  This package is a free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' This function implements a fundamental part of RST: the indiscernibility relation.
#' This binary relation indicates whether it is possible to discriminate any given pair of objects from an information system. 
#'  
#' This function can be used as a basic building block for development of other RST-based methods.
#' A more detailed explanation of the notion of indiscernibility relation can be found in \code{\link{A.Introduction-RoughSets}}.
#'
#' @title Computation of indiscernibility classes based on the rough set theory
#' @author Andrzej Janusz
#'
#' @param decision.table an object inheriting from the \code{"DecisionTable"} class, which represents a decision system. 
#'        See \code{\link{SF.asDecisionTable}}.
#' @param feature.set an integer vector indicating indexes of attributes which should be used or an object inheriting from
#'                  the \code{FeatureSubset} class.
#'                  The computed indiscernibility classes will be relative to this attribute set. 
#'                  The default value is \code{NULL} which means that 
#'                  all conditional attributes should be considered. It is usually reasonable 
#'                  to discretize numeric attributes before the computation of indiscernibility classes.
#'
#' @return An object of a class \code{"IndiscernibilityRelation"} which is a list with the following components:
#'          \itemize{
#'          \item \code{IND.relation}: a list of indiscernibility classes in the data. Each class is represented by indices 
#'                of data instances which belong to that class
#'          \item \code{type.relation}: a character vector representing a type of relation used in computations. Currently, 
#'                only \code{"equivalence"} is provided. 
#'          \item \code{type.model}: a character vector identifying the type of model which is used. 
#'                In this case, it is \code{"RST"} which means the rough set theory.
#'          }
#'          
#' @seealso \code{\link{BC.LU.approximation.RST}}, \code{\link{FS.reduct.computation}}, \code{\link{FS.feature.subset.computation}}
#' 
#' @references
#' Z. Pawlak, "Rough Sets", International Journal of Computer and Information Sciences, 
#' vol. 11, no. 5, p. 341 - 356 (1982).
#'
#' @examples
#' #############################################
#' data(RoughSetData)
#' hiring.data <- RoughSetData$hiring.dt
#'
#' ## In this case, we only consider the second and third attribute:
#' A <- c(2,3)
#' ## We can also compute a decision reduct:
#' B <- FS.reduct.computation(hiring.data)
#' 
#' ## Compute the indiscernibility classes:
#' IND.A <- BC.IND.relation.RST(hiring.data, feature.set = A)
#' IND.A
#' 
#' IND.B <- BC.IND.relation.RST(hiring.data, feature.set = B)
#' IND.B
#'
#' @export
BC.IND.relation.RST <- function(decision.table, feature.set = NULL){
	
  if (!inherits(decision.table, "DecisionTable")) {
    stop("Provided data should inherit from the \'DecisionTable\' class.")
  }
  
	## get data
	objects <- decision.table
	nominal.att <- attr(decision.table, "nominal.attrs")
	decision.attr <- attr(decision.table, "decision.attr")
	
	## initialize
	if (is.null(feature.set)) {
		if(length(decision.attr) > 0) feature.set <- (1:ncol(objects))[-decision.attr]
    else feature.set <- 1:ncol(objects)
	} else {
    if (inherits(feature.set, "FeatureSubset")) feature.set <- feature.set$reduct
	}

	## check for non nominal attribute
	if (!all(nominal.att[c(feature.set)])){
		stop("please discretize attributes before computing an equivalence-based indiscernibility relation")
	}
	
	#compute the indiscernibility classes
	if (length(feature.set) == 1){
		IND = split(1:nrow(objects), do.call(paste, list(objects[ , feature.set])))
	} else {
		IND = split(1:nrow(objects), do.call(paste, objects[ , feature.set]))
	}
	
	## construct class
	mod <- list(IND.relation = IND, type.relation = "equivalence", type.model = "RST")	
	class.mod <- ObjectFactory(mod, classname = "IndiscernibilityRelation")
	
	return(class.mod)
}

#' This function implements a fundamental part of RST: computation of lower and upper approximations. 
#' The lower and upper approximations determine whether the objects can be certainty or possibly classified 
#' to a particular decision class on the basis of available knowledge.
#' 
#' This function can be used as a basic building block for development of other RST-based methods.
#' A more detailed explanation of this notion can be found in \code{\link{A.Introduction-RoughSets}}.
#' 
#' @title Computation of lower and upper approximations of decision classes
#' @author Andrzej Janusz
#'
#' @param decision.table an object inheriting from the \code{"DecisionTable"} class, which represents a decision system. 
#'        See \code{\link{SF.asDecisionTable}}.
#' @param IND an object inheriting from the \code{"IndiscernibilityRelation"} class, which represents indiscernibility clasees in the data. 
#' 
#' @return An object of a class \code{"LowerUpperApproximation"} which is a list with the following components:
#'         \itemize{
#'          \item \code{lower.approximation}: a list with indices of data instances included in lower approximations of decision classes.
#'          \item \code{upper.approximation}: a list with indices of data instances included in upper approximations of decision classes.
#'          \item \code{type.model}: a character vector identifying the type of model which was used. 
#'                In this case, it is \code{"RST"} which means the rough set theory.
#'          } 
#'          
#' @seealso \code{\link{BC.IND.relation.RST}}, \code{\link{BC.LU.approximation.FRST}}
#' 
#' @references
#' Z. Pawlak, "Rough Sets", International Journal of Computer and Information Sciences, 
#' vol. 11, no. 5, p. 341 - 356 (1982).
#'
#' @examples
#' #######################################
#' data(RoughSetData)
#' hiring.data <- RoughSetData$hiring.dt
#'
#' ## We select a single attribute for computation of indiscernibility classes:
#' A <- c(2)
#' 
#' ## Compute the indiscernibility classes:
#' IND.A <- BC.IND.relation.RST(hiring.data, feature.set = A)
#'
#' ## Compute the lower and upper approximations:
#' roughset <- BC.LU.approximation.RST(hiring.data, IND.A)
#' roughset
#' 
#' @export
BC.LU.approximation.RST <- function(decision.table, IND){
	
  if (!inherits(decision.table, "DecisionTable")) {
    stop("Provided data should inherit from the \'DecisionTable\' class.")
  }
  
	## get the data
	objects <- decision.table
	nominal.att <- attr(decision.table, "nominal.attrs")
  
  if (!inherits(IND, "IndiscernibilityRelation")) {
    stop("The list representing the indiscernibility relation should inherit from the \'IndiscernibilityRelation\' class.")
  }
  
	IND <- IND$IND.relation
	
	if (is.null(attr(decision.table, "decision.attr"))){
		decision.attr <- ncol(objects)
	}
	else {
		decision.attr <- attr(decision.table, "decision.attr")
		## check if index of decision attribute is not in last index, then change their position
		if (decision.attr != ncol(objects)){
			objects <- cbind(objects[, -decision.attr, drop = FALSE], objects[, decision.attr, drop = FALSE])
			nominal.att <- c(nominal.att[-decision.attr], nominal.att[decision.attr])
			decision.attr = ncol(objects)
		}
	}	
	num.att <- ncol(objects)

	## get indiscernibility of decision attribute
	if (!all(nominal.att[-decision.attr])){
		stop("please, discretize attributes before computing an equivalence-based indiscernibility relation")
	} else {
		## get unique decisions
		uniqueDecisions <- as.character(unique(objects[[decision.attr]]))
	}
	
	IND.decision.attr <- lapply(IND, function(x, splitVec) split(x, splitVec[x]), as.character(objects[[decision.attr]]))
	
	## initialization
	lower.appr <- list()
	upper.appr <- list()
	for(i in 1:length(uniqueDecisions)) {
		tmpIdx1 = which(sapply(IND.decision.attr, function(x) (uniqueDecisions[i] %in% names(x))))
		upper.appr[[i]] = unlist(IND[tmpIdx1])
		tmpIdx2 = which(sapply(IND.decision.attr[tmpIdx1], function(x) (length(x) == 1)))
		if(length(tmpIdx2) > 0) {
			lower.appr[[i]] = unlist(IND[tmpIdx1][tmpIdx2])
		}
		else lower.appr[[i]] = integer()
		colnames(lower.appr[[i]]) <- NULL
		colnames(upper.appr[[i]]) <- NULL
	}
	rm(tmpIdx1, tmpIdx2)
  	
	## give the names of list
	names(lower.appr) <- uniqueDecisions
	names(upper.appr) <- uniqueDecisions
	
	## build class
	res <- list(lower.approximation = lower.appr, upper.approximation = upper.appr, type.model = "RST")
	class.mod <- ObjectFactory(res, classname = "LowerUpperApproximation")
	
	return(class.mod)
}


#' This function implements a fundamental part of RST: computation of a positive region and the
#' degree of dependency. This function can be used as a basic building block for development 
#' of other RST-based methods. A more detailed explanation of this notion can be found 
#' in \code{\link{A.Introduction-RoughSets}}.
#' 
#' @title Computation of a positive region
#' @author Andrzej Janusz
#'
#' @param decision.table an object inheriting from the \code{"DecisionTable"} class, which represents a decision system. 
#'        See \code{\link{SF.asDecisionTable}}.
#' @param roughset an object inheriting from the \code{"LowerUpperApproximation"} class, which represents
#'        lower and upper approximations of decision classes in the data. Such objects are typically produced by calling 
#'        the \code{\link{BC.LU.approximation.RST}} function.

#' @return An object of a class \code{"PositiveRegion"} which is a list with the following components:
#'         \itemize{
#'           \item \code{positive.reg}: an integer vector containing indices of data instances belonging 
#'                 to the positive region,
#'           \item \code{degree.dependency}: a numeric value giving the degree of dependency,
#'           \item \code{type.model}: a varacter vector identifying the utilized model. In this case, 
#'                 it is \code{"RST"} which means the rough set theory.       
#'         } 
#'         
#' @seealso \code{\link{BC.IND.relation.RST}}, \code{\link{BC.LU.approximation.RST}}, \code{\link{BC.LU.approximation.FRST}}
#' 
#' @references
#' Z. Pawlak, "Rough Sets", International Journal of Computer and Information Sciences, 
#' vol. 11, no. 5, p. 341 - 356 (1982).
#'
#' @examples
#' ########################################################
#' data(RoughSetData)
#' hiring.data <- RoughSetData$hiring.dt
#'
#' ## We select a single attribute for computation of indiscernibility classes:
#' A <- c(2)
#' 
#' ## compute the indiscernibility classes:
#' IND.A <- BC.IND.relation.RST(hiring.data, feature.set = A)
#'
#' ## compute the lower and upper approximation:
#' roughset <- BC.LU.approximation.RST(hiring.data, IND.A)
#'
#' ## get the positive region:
#' pos.region = BC.positive.reg.RST(hiring.data, roughset)
#' pos.region
#' 
#' @export
BC.positive.reg.RST <- function(decision.table, roughset) {
	
	## get all objects from the lower approximations of decision classes
	positive.reg <- unlist(roughset$lower.approximation)
	names(positive.reg) <- NULL
	
	## get degree of dependecy 
	degree.depend <- length(positive.reg)/nrow(decision.table)
	
	res <- list(positive.reg = positive.reg[order(positive.reg)],
	            degree.dependency = degree.depend, type.model = "RST")
	
	## build class
	class.mod <- ObjectFactory(res, classname = "PositiveRegion")

	return(class.mod)
}

#' This function implements a fundamental part of RST: a decision-relative discernibility matrix. This notion
#' was proposed by (Skowron and Rauszer, 1992) as a middle-step in many RST algorithms for computaion of reducts, 
#' discretization and rule induction. A more detailed explanation of this notion can be found 
#' in \code{\link{A.Introduction-RoughSets}}.
#'
#' @title Computation of a decision-relative discernibility matrix based on the rough set theory
#' @author Lala Septem Riza and Andrzej Janusz
#'
#' @param decision.table an object inheriting from the \code{"DecisionTable"} class, which represents a decision system. 
#'        See \code{\link{SF.asDecisionTable}}.
#' @param range.object an integer vector indicating objects for construction of the \eqn{k}-relative discernibility matrix. 
#'                The default value is \code{NULL} which means that all objects in the decision table are used.
#' @param return.matrix a logical value determining whether the discernibility matrix should be retunred in the output. 
#'        If it is set to FALSE (the default) only a list containing unique clauses from the CNF representation 
#'        of the discernibility function is returned. 
#' 
#' @return An object of a class \code{"DiscernibilityMatrix"} which has the following components: 
#' \itemize{
#' \item \code{disc.mat}: the decision-relative discernibility matrix which for pairs of objects from different 
#'       decision classes stores names of attributes which can be used to discern between them. Only pairs of 
#'       objects from different decision classes are considered. For other pairs the \code{disc.mat} contains
#'       \code{NA} values. Moreover, since the classical discernibility matrix is symmetric only the pairs 
#'       from the lower triangular part are considered.
#' \item \code{disc.list}: a list containing unique clauses from the CNF representation of the discernibility function,
#' \item \code{discernibility.type}: a type of discernibility relation used in the computations,
#' \item \code{type.model}: a character vector identifying the type of model which was used. 
#'                In this case, it is \code{"RST"} which means the rough set theory.
#' }
#' 
#' @seealso \code{\link{BC.IND.relation.RST}}, \code{\link{BC.LU.approximation.RST}}, \code{\link{BC.LU.approximation.FRST}}
#'          and \code{\link{BC.discernibility.mat.FRST}}
#' 
#' @examples
#' #######################################################################
#' ## Example 1: Constructing the decision-relative discernibility matrix
#' #######################################################################
#' data(RoughSetData)
#' hiring.data <- RoughSetData$hiring.dt
#'
#' ## building the decision-relation discernibility matrix
#' disc.matrix <- BC.discernibility.mat.RST(hiring.data, return.matrix = TRUE)
#' disc.matrix
#'
#' @references
#' A. Skowron and C. Rauszer,  
#' "The Discernibility Matrices and Functions in Information Systems", 
#' in: R. Słowinski (Ed.), Intelligent Decision Support: Handbook of Applications and
#' Advances of Rough Sets Theory, Kluwer Academic Publishers, Dordrecht, Netherland,  
#' p. 331 - 362 (1992).
#' @export
BC.discernibility.mat.RST <- function(decision.table, range.object = NULL, return.matrix = FALSE){

  if(!inherits(decision.table, "DecisionTable")) {
    stop("Provided data should inherit from the \'DecisionTable\' class.")
  }

	## get the data
	objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	if (is.null(attr(decision.table, "decision.attr"))){
	  stop("A decision attribute is not indicated.")
	}	else {
		decision.attr <- attr(decision.table, "decision.attr")
		## check if index of decision attribute is not in last index, then change their position
		if (decision.attr != ncol(objects)){
			objects <- cbind(objects[, -decision.attr, drop = FALSE], objects[, decision.attr, drop = FALSE])
			nominal.att <- c(nominal.att[-decision.attr], nominal.att[decision.attr])
		}
	}
	num.object <- nrow(objects)
	names.attr <- t(colnames(objects)[-decision.attr])
	
	## replace if range.object = NULL
	if (is.null(range.object)){
		range.object <- matrix(c(1, nrow(objects)), nrow = 1)
		num.object <- range.object[1, 2] - range.object[1, 1] + 1
	}
	
	## initialize the discernibility matrix
  disc.mat <- array(list(NA), dim = c(num.object, num.object, 1))
	
	decVector = as.character(objects[, ncol(objects)])
	dataMatrix = as.matrix(objects[, -ncol(objects)])
  
   ## construct the discernibility matrix
	for (i in 1 : (num.object-1)){
		tmpIdx = which(decVector[(i+1) : num.object] != decVector[i])
		if(length(tmpIdx) > 0)  {
			tmpIdx = tmpIdx + i
			for (j in tmpIdx){
    			## select different values on decision attribute only
  				## construct one element has multi value/list which object is discerned.
  				disc.attr <- names.attr[which(dataMatrix[i,] != dataMatrix[j,])]
  				disc.mat[j, i, 1] <- list(disc.attr)
			}
		}
		rm(tmpIdx)
	}
	disc.mat = as.data.frame(disc.mat)
	disc.list = unique(do.call(c, disc.mat))[-1]
	
	## build class
	if (return.matrix){
		discernibilityMatrix = list(disc.mat = disc.mat, disc.list = disc.list, 
                              names.attr = colnames(decision.table), type.discernibility = "RST", type.model = "RST")
	}
	else {
		discernibilityMatrix = list(disc.list = disc.list, 
                              names.attr = colnames(decision.table), type.discernibility = "RST", type.model = "RST")
	}
	discernibilityMatrix = ObjectFactory(discernibilityMatrix, classname = "DiscernibilityMatrix")
	return(discernibilityMatrix)
}


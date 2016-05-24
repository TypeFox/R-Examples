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
#' It is used for handling missing values by deleting instances. It should be noted that the output of the function is \code{val.NA} 
#' which contains indices of missing values and their values (i.e., NA). Therefore, in order to generate a new decision table (dataset) 
#' the user need to execute \code{\link{SF.applyDecTable}}. 
#'
#' @title Missing value completion by deleting instances
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}. 
#'        Note: missing values are recognized as NA. 
#' @seealso \code{\link{MV.missingValueCompletion}}
#' @return A class \code{"MissingValue"}. See \code{\link{MV.missingValueCompletion}}.
#' @references
#' J. Grzymala-Busse and W. Grzymala-Busse, "Handling Missing Attribute Values," in Data Mining and Knowledge Discovery Handbook, 
#' O. Maimon and L. Rokach, Eds. New York : Springer, 2010, pp. 33-51
#'
#' @examples
#' #############################################
#' ## Example : Deletion Cases
#' #############################################
#' dt.ex1 <- data.frame(
#'      c("high", "very_high", NA, "high", "high", "normal", "normal", NA), 
#'      c(NA, "yes", "no", "yes", NA, "yes", "no", "yes"), 
#'      c("no", "yes", "no", "yes", "yes", "no", "yes", NA),
#'      c("yes", "yes", "no", "yes", "no", "no", "no", "yes"))
#' colnames(dt.ex1) <- c("Temp", "Headache", "Nausea", "Flu")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4, 
#'                                     indx.nominal = c(1:4))
#' indx = MV.deletionCases(decision.table)
#' @export
MV.deletionCases <- function(decision.table){
	## search NA and Inf position
	indx.na <- suppressWarnings(as.data.frame(which(is.na(decision.table), arr.ind=TRUE)))
	
	if (nrow(indx.na) != 0){
		indx.na <- cbind(indx.na, val = NA)
	}
	
	##new.decTable <- decision.table[-c(indx.na[, 1], indx.inf[, 1]), ,drop = FALSE]
	mod = list(val.NA = indx.na, type.method = "deletionCases")
	class.mod <- ObjectFactory(mod, classname = "MissingValue")
	
	return(class.mod)
}

#' It is used for handling missing values by assigning the most common value of an attribute restricted to a concept. 
#' If an attributes containing missing values is continuous/real, the method uses mean of the attribute instead of the most common value.
#' In order to generate a new decision table, we need to execute \code{\link{SF.applyDecTable}}. 
#'
#' @title The most common value or mean of an attribute restricted to a concept
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}. 
#'        Note: missing values are recognized as NA. 
#' @seealso \code{\link{MV.missingValueCompletion}}
#' @return A class \code{"MissingValue"}. See \code{\link{MV.missingValueCompletion}}.
#' @references
#' J. Grzymala-Busse and W. Grzymala-Busse, "Handling Missing Attribute Values," in Data Mining and Knowledge Discovery Handbook, 
#' O. Maimon and L. Rokach, Eds. New York : Springer, 2010, pp. 33-51
#'
#' @examples
#' #############################################
#' ## Example: The most common value
#' #############################################
#' dt.ex1 <- data.frame(
#'      c(100.2, 102.6, NA, 99.6, 99.8, 96.4, 96.6, NA), 
#'      c(NA, "yes", "no", "yes", NA, "yes", "no", "yes"), 
#'      c("no", "yes", "no", "yes", "yes", "no", "yes", NA),
#'      c("yes", "yes", "no", "yes", "no", "no", "no", "yes"))
#' colnames(dt.ex1) <- c("Temp", "Headache", "Nausea", "Flu")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4, 
#'                                     indx.nominal = c(2:4))
#' indx = MV.mostCommonValResConcept(decision.table)
#' @export
MV.mostCommonValResConcept <- function(decision.table){
	## search NA and Inf position
	indx.na <- suppressWarnings(data.frame(which(is.na(decision.table), arr.ind=TRUE)))
	indx.dec <- attr(decision.table, "decision.attr")
	i = 1
	while(i <= nrow(indx.na)){
		indx.consider <- which(decision.table[,  indx.dec] == decision.table[indx.na[i, 1], indx.dec]) 	
		if (attr(decision.table, "nominal.attr")[indx.na[i, 2]] == FALSE){
			temp <- summary(as.factor(stats::na.omit(decision.table[indx.consider, indx.na[i, 2]])))
			indx.na[i, 3] <- names(which.max(temp))[1]
		}
		else {
			data <- stats::na.omit(decision.table[indx.consider, indx.na[i, 2]])		
			temp <- summary(as.factor(data))
			indx.na[i, 3] <- names(which.max(temp))[1]
		}
		i = i + 1
	}
	
	if (nrow(indx.na)>=1){
		colnames(indx.na) <- c("Row", "Col", "Value")
	}
	
	mod = list(val.NA = indx.na, type.method = "mostCommonValResConcept")
	class.mod <- ObjectFactory(mod, classname = "MissingValue")
	
	return(class.mod)
}

#' It is used for handling missing values by replacing the attribute mean or common values. If an attributes containing missing values is continuous/real, the method uses mean of the attribute instead of the most common value.
#' In order to generate a new decision table, we need to execute \code{\link{SF.applyDecTable}}. 
#'
#' @title Replacing missing attribute values by the attribute mean or common values
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}. 
#'        Note: missing values are recognized as NA. 
#' @seealso \code{\link{MV.missingValueCompletion}}
#' @return A class \code{"MissingValue"}. See \code{\link{MV.missingValueCompletion}}.
#' @references
#' J. Grzymala-Busse and W. Grzymala-Busse, "Handling Missing Attribute Values," in Data Mining and Knowledge Discovery Handbook, 
#' O. Maimon and L. Rokach, Eds. New York : Springer, 2010, pp. 33-51
#'
#' @examples
#' #############################################
#' ## Example: Replacing missing attribute values
#' ##          by the attribute mean/common values
#' #############################################
#' dt.ex1 <- data.frame(
#'      c(100.2, 102.6, NA, 99.6, 99.8, 96.4, 96.6, NA), 
#'      c(NA, "yes", "no", "yes", NA, "yes", "no", "yes"), 
#'      c("no", "yes", "no", "yes", "yes", "no", "yes", NA),
#'      c("yes", "yes", "no", "yes", "no", "no", "no", "yes"))
#' colnames(dt.ex1) <- c("Temp", "Headache", "Nausea", "Flu")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4, 
#'                                     indx.nominal = c(2:4))
#' indx = MV.mostCommonVal(decision.table)
#' @export
MV.mostCommonVal <- function(decision.table){
	## search NA and Inf position
	indx.na <- suppressWarnings(data.frame(which(is.na(decision.table), arr.ind=TRUE)))
	
	i = 1
	while(i <= nrow(indx.na)){	 
		if (attr(decision.table, "nominal.attr")[indx.na[i, 2]] == FALSE){
			data <- stats::na.omit(decision.table[, indx.na[i, 2]])		
			indx.na[i, 3] <- mean(data)
		}
		else {
			data <- stats::na.omit(decision.table[, indx.na[i, 2]])		
			temp <- summary(as.factor(data))
			indx.na[i, 3] <- names(which.max(temp))[1]
		}
		
		i = i + 1
	}
	
	if (nrow(indx.na)>=1){
		colnames(indx.na) <- c("Row", "Col", "Value")
	}
	
	mod = list(val.NA = indx.na, type.method = "mostCommonVal")
	class.mod <- ObjectFactory(mod, classname = "MissingValue")
	
	return(class.mod)
}

#' It is used for handling missing values based on the global closest fit. 
#'
#' The global closes fit method is based on replacing
#' a missing attribute value by the known value in another case that resembles
#' as much as possible the case with the missing attribute value. In searching
#' for the closest fit case we compare two vectors of attribute values, one vector
#' corresponds to the case with a missing attribute value, the other vector is a candidate
#' for the closest fit. The search is conducted for all cases, hence the name
#' global closest fit. For each case a distance is computed, the case for which the
#' distance is the smallest is the closest fitting case that is used to determine the
#' missing attribute value.
#'
#' @title Global Closest Fit
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}. 
#'        Note: missing values are recognized as NA. 
#' @seealso \code{\link{MV.missingValueCompletion}}
#' @return A class \code{"MissingValue"}. See \code{\link{MV.missingValueCompletion}}.
#' @references
#' J. Grzymala-Busse and W. Grzymala-Busse, "Handling Missing Attribute Values," in Data Mining and Knowledge Discovery Handbook, 
#' O. Maimon and L. Rokach, Eds. New York : Springer, 2010, pp. 33-51
#'
#' @examples
#' #############################################
#' ## Example: Global Closest Fit
#' #############################################
#' dt.ex1 <- data.frame(
#'      c(100.2, 102.6, NA, 99.6, 99.8, 96.4, 96.6, NA), 
#'      c(NA, "yes", "no", "yes", NA, "yes", "no", "yes"), 
#'      c("no", "yes", "no", "yes", "yes", "no", "yes", NA),
#'      c("yes", "yes", "no", "yes", "no", "no", "no", "yes"))
#' colnames(dt.ex1) <- c("Temp", "Headache", "Nausea", "Flu")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4, 
#'                                     indx.nominal = c(2:4))
#' indx = MV.globalClosestFit(decision.table)
#' @export
MV.globalClosestFit <- function(decision.table){
	options(stringsAsFactors = FALSE)	
	
	## search NA and Inf position
	indx.na <- suppressWarnings(data.frame(which(is.na(decision.table), arr.ind=TRUE)))
	
	h = 1
	while(h <= nrow(indx.na)){
		seqq <- seq(1:nrow(decision.table))
		seqq.row <- seqq[-indx.na[h, 1]]
		dt.ref <- decision.table[indx.na[h, 1], ,drop=FALSE]
		distance <- matrix()
		for (i in seqq.row){
			temp <- 0
			for (j in 1 : (ncol(decision.table)-1)){
				if ((attr(decision.table, "nominal.attr")[j] && dt.ref[1, j] != decision.table[i, j]) || is.na(decision.table[i, j]) || is.na(dt.ref[1, j])){
					dis <- 1
				}
				else if (dt.ref[1, j] == decision.table[i, j]){
					dis <- 0
				}				
				else {
					dis <- abs(dt.ref[1, j] - decision.table[i, j])/(max(decision.table[, j], na.rm = TRUE) - min(decision.table[, j], na.rm = TRUE))
				}
				temp <- temp + dis
			}
			distance[i] <- temp
		}
		indx.min <- which.min(distance)[1]
		if (attr(decision.table, "nominal.attr")[indx.na[h, 2]] == FALSE){
			indx.na[h, 3] <- decision.table[seqq.row[indx.min], indx.na[h, 2]]
		}
		else {
			indx.na[h, 3] <- as.character(decision.table[seqq.row[indx.min], indx.na[h, 2]])
		}
		
		h = h + 1
	}
		
	if (nrow(indx.na)>=1){
		colnames(indx.na) <- c("Row", "Col", "Value")
	}
	
	mod = list(val.NA = indx.na, type.method = "globalClosestFit")
	class.mod <- ObjectFactory(mod, classname = "MissingValue")
	
	return(class.mod)
}

#' It is used for handling missing values based on the concept closest fit.
#'
#' This method is similar to the global closest fit method. The difference is
#' that the original data set, containing missing attribute values, is first split into
#' smaller data sets, each smaller data set corresponds to a concept from the original
#' data set. More precisely, every smaller data set is constructed from one of
#' the original concepts, by restricting cases to the concept.
#'
#' @title Concept Closest Fit
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}. 
#'        Note: missing values are recognized as NA. 
#' @seealso \code{\link{MV.missingValueCompletion}}
#' @return A class \code{"MissingValue"}. See \code{\link{MV.missingValueCompletion}}.
#' @references
#' J. Grzymala-Busse and W. Grzymala-Busse, "Handling Missing Attribute Values," in Data Mining and Knowledge Discovery Handbook, 
#' O. Maimon and L. Rokach, Eds. New York : Springer, 2010, pp. 33-51
#'
#' @examples
#' #############################################
#' ## Example: Concept Closest Fit
#' #############################################
#' dt.ex1 <- data.frame(
#'      c(100.2, 102.6, NA, 99.6, 99.8, 96.4, 96.6, NA), 
#'      c(NA, "yes", "no", "yes", NA, "yes", "no", "yes"), 
#'      c("no", "yes", "no", "yes", "yes", "no", "yes", NA),
#'      c("yes", "yes", "no", "yes", "no", "no", "no", "yes"))
#' colnames(dt.ex1) <- c("Temp", "Headache", "Nausea", "Flu")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4, 
#'                                     indx.nominal = c(2:4))
#' indx = MV.conceptClosestFit(decision.table)
#' @export
MV.conceptClosestFit <- function(decision.table){
	options(stringsAsFactors = FALSE)	
	
	## search NA and Inf position
	indx.na <- suppressWarnings(data.frame(which(is.na(decision.table), arr.ind=TRUE)))
	
	h = 1
	while(h <= nrow(indx.na)){
		dt.ref <- decision.table[indx.na[h, 1], ,drop=FALSE]
		seqq <- seq(1:nrow(decision.table))
		indx.del <- which(decision.table[, attr(decision.table, "decision.attr")] == dt.ref[1, attr(decision.table, "decision.attr")])
		seqq.row <- seqq[-c(indx.na[h, 1], indx.del)]		
		distance <- matrix()
		for (i in seqq.row){
			temp <- 0
			for (j in 1 : (ncol(decision.table)-1)){
				if ((attr(decision.table, "nominal.attr")[j] && dt.ref[1, j] != decision.table[i, j]) || is.na(decision.table[i, j]) || is.na(dt.ref[1, j])){
					dis <- 1
				}
				else if (dt.ref[1, j] == decision.table[i, j]){
					dis <- 0
				}				
				else {
					dis <- abs(dt.ref[1, j] - decision.table[i, j])/(max(decision.table[, j], na.rm = TRUE) - min(decision.table[, j], na.rm = TRUE))
				}
				temp <- temp + dis
			}
			distance[i] <- temp
		}
		indx.min <- which.min(distance)[1]
		if (attr(decision.table, "nominal.attr")[indx.na[h, 2]] == FALSE){
			indx.na[h, 3] <- decision.table[seqq.row[indx.min], indx.na[h, 2]]
		}
		else {
			indx.na[h, 3] <- as.character(decision.table[seqq.row[indx.min], indx.na[h, 2]])
		}
		
		h = h + 1
	}
	
	if (nrow(indx.na)>=1){
		colnames(indx.na) <- c("Row", "Col", "Value")
	}
	
	mod = list(val.NA = indx.na, type.method = "conceptClosestFit")
	class.mod <- ObjectFactory(mod, classname = "MissingValue")
	
	return(class.mod)
}


#' It is a wrapper function for missing value completion. 
#'
#' @title Wrapper function of missing value completion
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}. 
#'        Note: missing values are recognized as NA. 
#' @param type.method one of the following methods:
#'         \itemize{
#'         \item \code{"deletionCases"}: See \code{\link{MV.deletionCases}}.
#'         \item \code{"mostCommonValResConcept"}: See \code{\link{MV.mostCommonValResConcept}}.
#'         \item \code{"mostCommonVal"}: See \code{\link{MV.mostCommonVal}}.
#'         \item \code{"globalClosestFit"}: See \code{\link{MV.globalClosestFit}}.
#'         \item \code{"conceptClosestFit"}: See \code{\link{MV.conceptClosestFit}}.
#'         }
#' @seealso \code{\link{MV.deletionCases}}, \code{\link{MV.mostCommonValResConcept}}, \code{\link{MV.mostCommonVal}},
#'       \code{\link{MV.globalClosestFit}}, and \code{\link{MV.conceptClosestFit}}.
#' @return A class \code{"MissingValue"} which contains
#'          \itemize{
#'          \item \code{val.NA}: a matrix containing indices of missing value (i.e., unknown values) positions and their values.
#'          \item \code{type.method}: a string showing the type of used method. In this case, it is \code{"deleteCases"}. 
#'          }
#' @references
#' J. Grzymala-Busse and W. Grzymala-Busse, "Handling Missing Attribute Values," in Data Mining and Knowledge Discovery Handbook, 
#' O. Maimon and L. Rokach, Eds. New York : Springer, 2010, pp. 33-51
#'
#' @examples
#' #############################################
#' ## Example : 
#' #############################################
#' dt.ex1 <- data.frame(
#'      c(100.2, 102.6, NA, 99.6, 99.8, 96.4, 96.6, NA), 
#'      c(NA, "yes", "no", "yes", NA, "yes", "no", "yes"), 
#'      c("no", "yes", "no", "yes", "yes", "no", "yes", NA),
#'      c("yes", "yes", "no", "yes", "no", "no", "no", "yes"))
#' colnames(dt.ex1) <- c("Temp", "Headache", "Nausea", "Flu")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4, 
#'                                     indx.nominal = c(2:4))
#' indx = MV.missingValueCompletion(decision.table, type.method = "deletionCases")
#'
#' ## generate new decision table
#' new.decTable <- SF.applyDecTable(decision.table, indx)
#' @export
MV.missingValueCompletion <- function(decision.table, type.method = "deletionCases"){
  
  if (!(type.method %in% c("deletionCases", "mostCommonValResConcept", "mostCommonVal",
                           "globalClosestFit", "conceptClosestFit"))) {
    stop("Unrecognized method.")
  }
  
  if (!inherits(decision.table, "DecisionTable")) {
    stop("Provided data should inherit from the \'DecisionTable\' class.")
  }
  
	val.NA <- switch(type.method,  
    	             deletionCases = MV.deletionCases(decision.table),
    	             mostCommonValResConcept = MV.mostCommonValResConcept(decision.table),
                   mostCommonVal = MV.mostCommonVal(decision.table),
    	             globalClosestFit = MV.globalClosestFit(decision.table),
    	             conceptClosestFit = MV.conceptClosestFit(decision.table) )
  	
	return(val.NA)
}

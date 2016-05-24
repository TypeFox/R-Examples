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
#' This is a function used to implement a fundamental concept of FRST which is fuzzy indiscernibility relations. 
#' It is used for any fuzzy relations that determine the degree to which two objects are indiscernibility. 
#' The detailed description about basic concepts of FRST 
#' can be found in \code{\link{B.Introduction-FuzzyRoughSets}}.
#' 
#' Briefly, the indiscernibility relation is a relation that shows a degree of similarity among the objects.
#' For example, \eqn{R(x_i, x_j) = 0} means the object \eqn{x_i} is completely different from \eqn{x_j}, 
#' and otherwise if \eqn{R(x_i, x_j) = 1}, while between those values we consider a degree of similarity. 
#' To calculate this relation, several methods have been implemented
#' in this function which are approaches based on fuzzy tolerance, equivalence and \eqn{T}-equivalence relations. 
#' The fuzzy tolerance relations proposed by (Jensen and Shen, 2009) include three equations while
#' (Hu, 2004) proposed five \eqn{T_{cos}}-transitive kernel functions as fuzzy \eqn{T}-equivalence relations. The simple
#' algorithm of fuzzy equivalence relation is implemented as well. Furthermore, we facilitate users to define their own equation for similarity relation.
#' 
#' To calculate a particular relation, we should pay attention to several components in 
#' the parameter \code{control}. The main component in the \code{control} parameter is \code{type.relation} that defines 
#' what kind of approach we are going to use. The detailed explanation about the parameters and their equations 
#' is as follows:
#' \itemize{
#' \item \code{"tolerance"}: It refers to fuzzy tolerance relations proposed by (Jensen and Shen, 2009). 
#' 		In order to represent the \code{"tolerance"} relation, we must set \code{type.relation} as follows:
#'
#' 		\code{type.relation = c("tolerance", <chosen equation>)} 
#' 
#'      where the chosen equation called as \code{t.similarity} is one of the
#' 		\code{"eq.1"}, \code{"eq.2"}, and \code{"eq.3"} equations which have been explained in \code{\link{B.Introduction-FuzzyRoughSets}}.
#'      
#' \item \code{"transitive.kernel"}: It refers to the relations employing kernel functions (Genton, 2001).
#' 		In order to represent the relation, we must set the \code{type.relation} parameter as follows.
#'
#' 		\code{type.relation = c("transitive.kernel", <chosen equation>, <delta>)}      
#'
#'      where the chosen equation is one of five following equations (called \code{t.similarity}):
#'       \itemize{
#'       \item \code{"gaussian"}: It means Gaussian kernel which is \eqn{R_G(x,y) = \exp (-\frac{|x - y|^2}{\delta})}
#'       \item \code{"exponential"}: It means exponential kernel which is \eqn{R_E(x,y) = \exp(-\frac{|x - y|}{\delta})}
#'       \item \code{"rational"}: It means rational quadratic kernel which is \eqn{R_R(x,y) = 1 - \frac{|x - y|^2}{|x - y|^2 + \delta}}
#'       \item \code{"circular"}: It means circular kernel which is if \eqn{|x - y| < \delta}, 
#'              \eqn{R_C(x,y) = \frac{2}{\pi}\arccos(\frac{|x - y|}{\delta}) - \frac{2}{\pi}\frac{|x - y|}{\delta}\sqrt{1 - (\frac{|x - y|}{\delta})^2}}
#'       \item \code{"spherical"}: It means spherical kernel which is if \eqn{|x - y| < \delta},
#'              \eqn{R_S(x,y) = 1 - \frac{3}{2}\frac{|x - y|}{\delta} + \frac{1}{2}(\frac{|x - y|}{\delta})^3}
#'       }
#'       and \code{delta} is a specified value. 
#'       For example: let us assume we are going to use \code{"transitive.kernel"} as the fuzzy relation,
#'       \code{"gaussian"} as its equation and the delta is 0.2. So, we assign the \code{type.relation} parameter as follows:
#' 
#'       \code{type.relation = c("transitive.kernel", "gaussian", 0.2)}
#'
#'       If we omit the \code{delta} parameter then we are using \code{"gaussian"} defined as \eqn{R_E(x,y) = \exp(-\frac{|x - y|}{2\sigma^2})}, where \eqn{\sigma} is the variance.
#'       Furthermore, when using this relation, usually we set 
#'
#'       \code{type.aggregation = c("t.tnorm", "t.cos")}.  
#'
#' \item \code{"kernel.frst"}: It refers to \eqn{T}-equivalence relation proposed by (Hu, 2004).
#' 		In order to represent the relation, we must set \code{type.relation} parameter as follows.
#'
#' 		\code{type.relation = c("kernel.frst", <chosen equation>, <delta>)}      
#'
#'      where the chosen equation is one of the kernel functions, but they have different names corresponding to previous ones: 
#'      \code{"gaussian.kernel"}, \code{"exponential.kernel"}, \code{"rational.kernel"}, \code{"circular.kernel"}, and \code{"spherical.kernel"}. And
#'      \code{delta} is a specified value. 
#'       For example: let us assume we are going to use \code{"gaussian.kernel"} as its equation and the delta is 0.2. 
#'       So, we assign the \code{type.relation} parameter as follows:
#' 
#'       \code{type.relation = c("kernel.frst", "gaussian.kernel", 0.2)}
#'
#'       In this case, we do not need to define type of aggregation. Furthemore, regarding the distance used in the equations if objects \eqn{x} and \eqn{y} contains mixed values (nominal and continuous)
#'       then we use the Gower distance and we use the euclidean distance for continuous only. 
#'
#'  \item \code{"transitive.closure"}: It refers to similarity relation (also called fuzzy equivalence relation). 
#'       We consider a simple algorithm to calculate this relation as follows:
#'
#'       Input: a fuzzy relation R
#'
#'       Output: a min-transitive fuzzy relation \eqn{R^m}
#'
#'       Algorithm:
#'
#'       1. For every x, y: compute
#'
#'       \eqn{R'(x,y) = max(R(x,y), max_{z \in U}min(R(x,z), R(z,y)))}
#'
#'       2. If \eqn{R' \not= R}, then \eqn{R \gets R'} and go to 1, else \eqn{R^m \gets R'}
#'
#'      For interested readers, other algorithms can be seen in (Naessens et al, 2002). Let \code{"eq.1"} be the \eqn{R} fuzzy relations, to define it as parameter is
#'      
#'      \code{type.relation = c("transitive.closure", "eq.1")}. We can also use other equations that have been explained in \code{"tolerance"} and \code{"transitive.kernel"}.
#'
#' \item \code{"crisp"}: It uses crisp equality for all attributes and 
#'       we set the parameter \code{type.relation = "crisp"}. In this case, we only have \eqn{R(x_i, x_j) = 0} 
#'       which means the object \eqn{x_i} is completely different from \eqn{x_j}, 
#'       and otherwise if \eqn{R(x_i, x_j) = 1}. 
#'       
#'  \item \code{"custom"}: this value means that we define our own equation for the indiscernibility relation.
#'        The equation should be defined in parameter \code{FUN.relation}. 
#'
#'       \code{type.relation = c("custom", <FUN.relation>)}
#'
#'        The function \code{FUN.relation} should consist of three arguments which are \code{decision.table}, 
#'        \code{x}, and \code{y}, where \code{x} and \code{y} represent two objects which we want to compare. 
#'        It should be noted that the user must ensure that the values of this equation are always between 0 and 1. 
#'        An example can be seen in Section \code{Examples}.
#' }
#'
#' Beside the above \code{type.relation}, we provide several options of values for the \code{type.aggregation} parameter. 
#' The following is a description about it.
#' \itemize{
#'    \item \code{type.aggregation = c("crisp")}: It uses crisp equality for all attributes.
#'    \item \code{type.aggregation = c("t.tnorm", <t.tnorm operator>)}: It means we are using \code{"t.tnorm"} aggregation 
#'               which is a triangular norm operator with 
#'               a specified operator \code{t.tnorm} as follows: 
#'     			\itemize{
#'     				\item \code{"min"}: standard t-norm i.e., \eqn{min(x_1, x_2)}.
#'     				\item \code{"hamacher"}: hamacher product i.e., \eqn{(x_1 * x_2)/(x_1 + x_2 - x_1 * x_2)}.
#'     				\item \code{"yager"}: yager class i.e., \eqn{1 - min(1, ((1 - x_1) + (1 - x_2)))}.
#'     				\item \code{"product"}: product t-norm i.e., \eqn{(x_1 * x_2)}.
#'     				\item \code{"lukasiewicz"}: lukasiewicz's t-norm (default) i.e., \eqn{max(x_2 + x_1 - 1, 0)}. 
#'     				\item \code{"t.cos"}: \eqn{T_{cos}}t-norm i.e., \eqn{max(x_1 * x_2 - \sqrt{1 - x_1^2}\sqrt{1 - x_2^2, 0})}.
#'                  \item \code{FUN.tnorm}: It is a user-defined function for \code{"t.tnorm"}. It has to have two arguments, for example:
#'                              
#' 								\code{FUN.tnorm <- function(left.side, right.side)}
#'
#'               					 \code{if ((left.side + right.side) > 1)}
#'
#'                  					  \code{return(min(left.side, right.side))}
#'
#'                					 \code{else return(0)}
#'               	}
#' 		  			 The default value is \code{type.aggregation = c("t.tnorm", "lukasiewicz")}.
#'    \item \code{type.aggregation = c("custom", <FUN.agg>)}: It is used to define our own function for a type of aggregations. \code{<FUN.agg>} is 
#'          a function having one argument representing data that is produced by fuzzy similarity equation calculation.
#'          The data is a list of one or many matrices which depend on the number of considered attributes and has dimension: 
#'          the number of object \eqn{\times} the number of object. For example:
#'
#'         \code{FUN.agg <- function(data) return(Reduce("+", data)/length(data))}
#'
#'         which is a function to calculate average along all attributes. Then,
#'         we can set \code{type.aggregation} as follows:
#'
#'        \code{type.aggregation = c("general.custom", <FUN.agg>)}. An example can be seen in Section \code{Examples}.        
#'    }
#'
#' Furthermore, the use of this function has been illustrated in Section \code{Examples}. 
#' Finally, this function is important since it is a basic function needed by other functions, such as \code{\link{BC.LU.approximation.FRST}} and 
#' \code{\link{BC.positive.reg.FRST}} for calculating lower and upper approximation and determining positive regions.
#'
#' @title The indiscernibility relation based on fuzzy rough set theory
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing a decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param attributes a numerical vector expressing indexes of subset of attributes to be considered. 
#'                 The default value is \code{NULL} which means that 
#'                 all conditional attributes will be considered.
#' @param control a list of other parameters consisting of the following parameters:
#'        \itemize{ 
#'        \item \code{type.relation}: a list containing string values that express the 
#'              type of the fuzzy relation and its equation. The default value is \code{type.relation = c("tolerance", "eq.1")}. See in the Section \code{Details}.
#'        \item \code{type.aggregation}: a list expressing type of aggregation. The default value is \code{type.aggregation = c("t.tnorm", "lukasiewicz")}. 
#'              See in the Section \code{Details}.
#'        }
#' @return A class \code{"IndiscernibilityRelation"} which contains
#'          \itemize{
#'          \item \code{IND.relation}: a matrix representing the indiscernibility relation over all objects. 
#'          \item \code{type.relation}: a vector representing the type of relation. 
#'          \item \code{type.aggregation}: a vector representing the type of aggregation operator.
#'          \item \code{type.model}: a string showing the type of model which is used. In this case it is \code{"FRST"} which means fuzzy rough set theory.
#'          }
#' @seealso \code{\link{BC.LU.approximation.FRST}}, \code{\link{BC.IND.relation.RST}}, \code{\link{BC.LU.approximation.RST}}, 
#'
#' and \code{\link{BC.positive.reg.FRST}}
#' @references
#' M. Genton, "Classes of Kernels for Machine Learning: a Statistics Perspective", 
#' J. Machine Learning Research, vol. 2, p. 299 - 312 (2001). 
#'
#' H. Naessens, H. De Meyer, and B. De Baets, 
#' "Algorithms for the Computation of T-Transitive Closures",
#' IEEE Trans. on Fuzzy Systems, vol. 10, No. 4, p. 541 - 551 (2002).
#'
#' R. Jensen and Q. Shen,  
#' "New Approaches to Fuzzy-Rough Feature Selection", 
#' IEEE Trans. on Fuzzy Systems, vol. 19, no. 4,
#' p. 824 - 838 (2009).
#'
#' Q. Hu, D. Yu, W. Pedrycz, and D. Chen, "Kernelized Fuzzy Rough Sets and Their Applications",
#' IEEE Trans. Knowledge Data Eng., vol. 23, no. 11, p. 1649 - 1471 (2011).
#' 
#' @examples
#' ###########################################################
#' ## Example 1: Dataset containing nominal values for 
#' ## all attributes.
#' ###########################################################
#' ## Decision table is represented as data frame
#' dt.ex1 <- data.frame(c(1,0,2,1,1,2,2,0), c(0, 1,0, 1,0,2,1,1), 
#'                         c(2,1,0,0,2,0,1,1), c(2,1,1,2,0,1,1,0), c(0,2,1,2,1,1,2,1))
#' colnames(dt.ex1) <- c("aa", "bb", "cc", "dd", "ee")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 5, 
#'       indx.nominal = c(1:5))
#'
#' ## In this case, we only consider the second and third attributes.
#' attributes <- c(2, 3)
#' 
#' ## calculate fuzzy indiscernibility relation ##
#' ## in this case, we are using "crisp" as a type of relation and type of aggregation
#' control.ind <- list(type.relation = c("crisp"), type.aggregation = c("crisp"))
#' IND <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.ind)
#' 
#' ###########################################################
#' ## Example 2: Dataset containing real-valued attributes
#' ###########################################################
#' dt.ex2 <- data.frame(c(-0.4, -0.4, -0.3, 0.3, 0.2, 0.2), 
#'                      c(-0.3, 0.2, -0.4, -0.3, -0.3, 0),
#'				        c(-0.5, -0.1, -0.3, 0, 0, 0),
#'				        c("no", "yes", "no", "yes", "yes", "no"))
#' colnames(dt.ex2) <- c("a", "b", "c", "d")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex2, decision.attr = 4)
#' 
#' ## in this case, we only consider the first and second attributes
#' attributes <- c(1, 2)
#'
#' ## Calculate fuzzy indiscernibility relation 
#' ## in this case, we choose "tolerance" relation and "eq.1" as similarity equation
#' ## and "lukasiewicz" as t-norm of type of aggregation
#' control.1 <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                     type.relation = c("tolerance", "eq.1"))
#' IND.1 <- BC.IND.relation.FRST(decision.table, attributes = attributes, 
#'                               control = control.1) 
#'
#' ## Calculate fuzzy indiscernibility relation: transitive.kernel
#' control.2 <- list(type.aggregation = c("t.tnorm", "t.cos"), 
#'                     type.relation = c("transitive.kernel", "gaussian", 0.2))
#' IND.2 <- BC.IND.relation.FRST(decision.table, attributes = attributes, 
#'                               control = control.2) 
#'
#' ## Calculate fuzzy indiscernibility relation: kernel.frst 
#' control.3 <- list(type.relation = c("kernel.frst", "gaussian.kernel", 0.2))
#' IND.3 <- BC.IND.relation.FRST(decision.table, attributes = attributes, 
#'                               control = control.3) 
#'
#' ## calculate fuzzy transitive closure
#' control.4 <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                     type.relation = c("transitive.closure", "eq.1"))
#' IND.4 <- BC.IND.relation.FRST(decision.table, attributes = attributes, 
#'                               control = control.4) 
#'
#' ## Calculate fuzzy indiscernibility relation: using user-defined relation
#' ## The customized function should have three arguments which are : decision.table 
#' ## and object x, and y.
#' ## This following example shows that we have an equation for similarity equation: 
#' ## 1 - abs(x - y) where x and y are two objects that will be compared.
#' ## In this case, we do not consider decision.table in calculation.
#' FUN.relation <- function(decision.table, x, y) {
#'            return(1 - (abs(x - y)))
#'        }
#' control.5 <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                      type.relation = c("custom", FUN.relation))
#' IND.5 <- BC.IND.relation.FRST(decision.table, attributes = attributes, 
#'                               control = control.5) 
#'
#' ## In this case, we calculate aggregation as average of all objects 
#' ## by executing the Reduce function
#' FUN.average <- function(data){
#'   	 return(Reduce("+", data)/length(data))
#' }
#' control.6 <- list(type.aggregation = c("custom", FUN.average), 
#'                      type.relation = c("tolerance", "eq.1"))
#' IND.6 <- BC.IND.relation.FRST(decision.table, attributes = attributes, 
#'                               control = control.6)
#'	
#' ## using user-defined tnorms 
#' FUN.tnorm <- function(left.side, right.side) {
#'                if ((left.side + right.side) > 1)
#'                    return(min(left.side, right.side))
#'                else return(0)}
#' control.7 <- list(type.aggregation = c("t.tnorm", FUN.tnorm), 
#'                     type.relation = c("tolerance", "eq.1"))
#' IND.7 <- BC.IND.relation.FRST(decision.table, attributes = attributes, 
#'                               control = control.7) 
#'
#' ## Calculate fuzzy indiscernibility relation: kernel fuzzy rough set 
#' control.8 <- list(type.relation = c("kernel.frst", "gaussian.kernel", 0.2))
#' IND.8 <- BC.IND.relation.FRST(decision.table, attributes = attributes, 
#'                               control = control.8) 						   
#' ##################################################################
#' ## Example 3: Dataset containing continuous and nominal attributes
#' ## Note. we only consider type.relation = c("tolerance", "eq.1")
#' ## but other approaches have the same way.
#' ##################################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$housing7.dt 
#' 
#' ## in this case, we only consider the attribute: 1, 2, 3, 4 
#' attributes <- c(1,2,3,4)
#'
#' ## Calculate fuzzy indiscernibility relation
#' control.9 <- list(type.aggregation = c("t.tnorm", "lukasiewicz"),
#'                     type.relation = c("tolerance", "eq.1"))
#' IND.9 <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.9) 
#' 
#' @export
BC.IND.relation.FRST <- function(decision.table, attributes = NULL, control = list()){
	## set default values of all parameters
	control <- setDefaultParametersIfMissing(control, list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
	           type.relation = c("tolerance", "eq.1"), disc.mat = FALSE, type.MVCompletion = FALSE, alpha = 1))
	
	## check missing values
	if (any(is.na(decision.table))){
		type.MVCompletion <- TRUE
	}
	else {
		type.MVCompletion <- FALSE
	}
	
	## get the data
	type.relation <- control$type.relation
	disc.mat <- control$disc.mat

	## get type of aggregation
	type.aggregation <- ch.typeAggregation(control$type.aggregation)

	## check whether attributes is exist or not
	if (is.null(attributes)){
		attributes <- c(seq(1, (ncol(decision.table) - 1)))
	}
	
	## calculate fuzzy similarity relation based on the chosen type.relation
	if (type.relation[1] == "crisp"){
		miu.Ra <- calc.fsimilarity(decision.table, attributes, t.tnorm = type.aggregation[[2]], disc.mat = disc.mat, 
		                           t.aggregation = type.aggregation[[1]], type.MVCompletion = type.MVCompletion, alpha = control$alpha)	
	}
	else if (type.relation[1] == "tolerance"){
		## set default value
		if (is.na(type.relation[2])) t.similarity <- c("eq.1")
		else t.similarity <- type.relation[2]
		type.relation <- list(type = type.relation[1], equation = t.similarity)
		
		## calculate relation
		miu.Ra <- calc.fsimilarity(decision.table, attributes, t.similarity, t.tnorm = type.aggregation[[2]], FUN = NULL, 
		                           disc.mat = disc.mat, t.aggregation = type.aggregation[[1]], type.MVCompletion = type.MVCompletion, alpha = control$alpha)	
	}
	else if (type.relation[1] == "transitive.kernel"){
		## set the default values
		if (is.na(type.relation[2])) t.similarity <- c("gaussian")
		else t.similarity <- type.relation[2]		
		if (is.na(type.relation[3])) delta <- NULL
		else delta <- as.numeric(type.relation[3])

		## calculate relation
		type.relation <- list(type = type.relation[1], equation = t.similarity, delta = delta)
		miu.Ra <- calc.fsimilarity(decision.table, attributes, t.similarity, t.tnorm = type.aggregation[[2]], FUN = NULL, disc.mat = disc.mat, 
		                           delta = delta, t.aggregation = type.aggregation[[1]], type.MVCompletion = type.MVCompletion, alpha = control$alpha)	
	}
	else if (type.relation[1] == "kernel.frst"){
		## set the default values
		if (is.na(type.relation[2])) t.similarity <- c("gaussian.kernel")
		else t.similarity <- type.relation[2]		
		if (is.na(type.relation[3])) delta <- NULL
		else delta <- as.numeric(type.relation[3])		
		if (t.similarity == "gaussian.kernel"){
			t.tnorm <- func.gaussian.kernel
		}
		else if (t.similarity == "exponential.kernel"){
			t.tnorm <- func.exponential.kernel
		}
		else if (t.similarity == "rational.kernel"){
			t.tnorm <- func.rational.kernel
		}
		else if (t.similarity == "circular.kernel"){
			t.tnorm <- func.circular.kernel
		}
		else if (t.similarity == "spherical.kernel"){
			t.tnorm <- func.spherical.kernel
		}
		
		## calculate relation
		type.relation <- list(type = type.relation[1], equation = t.similarity, delta = delta)
		miu.Ra <- calc.fsimilarity(decision.table, attributes, t.similarity, t.tnorm = t.tnorm, FUN = NULL, disc.mat = disc.mat, 
		                           delta = delta, t.aggregation = "kernel.frst", type.MVCompletion = type.MVCompletion, alpha = control$alpha)	
	}
	else if (type.relation[1] == "custom"){
		FUN <- type.relation[[2]]
		type.relation <- list(type = type.relation[1], equation = FUN)
		miu.Ra <- calc.fsimilarity(decision.table = decision.table, attributes = attributes, FUN = FUN, t.tnorm = type.aggregation[[2]], 
		                          disc.mat = disc.mat, t.aggregation = type.aggregation[[1]], type.MVCompletion = type.MVCompletion, alpha = control$alpha) 
	}
	else if (type.relation[1] == "transitive.closure"){
		if (is.na(type.relation[2])) t.similarity <- c("eq.1")
		else t.similarity <- type.relation[2]
			
		type.relation <- list(type = type.relation[1], equation = t.similarity)
		IND.relation.tolerance <- calc.fsimilarity(decision.table, attributes, t.similarity, t.tnorm = type.aggregation[[2]], 
		                                          FUN = NULL, disc.mat = disc.mat, t.aggregation = type.aggregation[[1]], 
												  type.MVCompletion = type.MVCompletion, alpha = control$alpha)	
		miu.Ra <- calc.transitive.closure(IND.relation.tolerance, type.MVCompletion = type.MVCompletion)
	}
	else {
		stop("Please define the type of relation which is going to be used.")
	}
	
	## build class
	if (type.relation[1] == "crisp"){
		mod <- list(IND.relation = miu.Ra, type.relation = "crisp", type.aggregation = "crisp", missing.value = type.MVCompletion, type.model = "FRST")
	}
	else if (type.relation[1] == "kernel.frst"){
		mod <- list(IND.relation = miu.Ra, type.relation = type.relation, missing.value = type.MVCompletion, type.model = "FRST")
	}
	else {
		mod <- list(IND.relation = miu.Ra, type.relation = type.relation, type.aggregation = type.aggregation, missing.value = type.MVCompletion, type.model = "FRST")
	}	
	class.mod <- ObjectFactory(mod, classname = "IndiscernibilityRelation")
	
	return(class.mod)

}


#' This is a function implementing a fundamental concept of FRST: fuzzy lower and upper approximations. 
#' Many options have been considered for determining lower and upper approximations, 
#' such as techniques based on implicator and t-norm functions proposed by 
#' (Radzikowska and Kerre, 2002).
#' 
#' Fuzzy lower and upper approximations as explained in \code{\link{B.Introduction-FuzzyRoughSets}} are used
#' to define to what extent the set of elements can be classified into a certain class strongly or weakly. We can perform various methods by choosing the parameter \code{type.LU}. 
#' The following is a list of all \code{type.LU} values: 
#'    \itemize{
#'    \item \code{"implicator.tnorm"}: It means implicator/t-norm based model proposed by (Radzikowska and Kerre, 2002). 
#'          The explanation has been given in \code{\link{B.Introduction-FuzzyRoughSets}}.
#'          Other parameters in \code{control} related with this approach are \code{t.tnorm} and \code{t.implicator}.
#'          In other words, when we are using \code{"implicator.tnorm"} as \code{type.LU}, 
#'          we should consider parameters \code{t.tnorm} and \code{t.implicator}.
#'          The possible values of these parameters can be seen in the description of parameters. 
#'    \item \code{"vqrs"}: It means vaguely quantified rough sets proposed by 
#'          (Cornelis et al, 2007). Basically, this concept proposed to replace fuzzy lower and upper approximations 
#'          based on Radzikowska and Kerre's technique (see \code{\link{B.Introduction-FuzzyRoughSets}})
#'          with the following equations, respectively. 
#' 
#'          \eqn{(R_{Q_u} \downarrow A)(y) = Q_u(\frac{|R_y \cap A|}{|R_y|})} 
#'
#'          \eqn{(R_{Q_l} \uparrow A)(y) = Q_l(\frac{|R_y \cap A|}{|R_y|})}
#'
#'          where the quantifier \eqn{Q_u} and \eqn{Q_l} represent the terms \code{most} and \code{some}. 
#'
#'    \item \code{"owa"}: It refers to ordered weighted average based fuzzy rough sets. 
#'            This method was introduced by (Cornelis et al, 2010) and computes the approximations
#'            by an aggregation process proposed by (Yager, 1988). The OWA-based lower and upper approximations of 
#'            \eqn{A} under \eqn{R} with weight vectors \eqn{W_l} and \eqn{W_u} are defined as
#'
#'           \eqn{(R \downarrow W_l A)(y) = OW A_{W_l}\langle I(R(x, y), A(y))\rangle}
#'
#'           \eqn{(R \uparrow W_u A)(y) = OW A_{W_u}\langle T(R(x, y), A(y))\rangle}
#'
#'           We provide two ways to define the weight vectors as follows:
#'
#'           \itemize{
#'           \item \code{m.owa}: Let \eqn{m.owa} be \eqn{m} and \eqn{m \le n}, this model is defined by
#'
#'           \eqn{W_l = <w_i^l> = w_{n+1-i}^l = \frac{2^{m-i}}{2^{m}-1}} for \eqn{i = 1,\ldots, m} and \eqn{0} for \eqn{i = m + 1, \ldots, n}
#'
#'           \eqn{W_u = <w_i^u> = w_i^u = \frac{2^{m - i}}{2^{m} - 1}} for \eqn{i = 1, \ldots, m} and \eqn{0} for \eqn{i = m + 1, \ldots, n}
#'
#'           where \eqn{n} is the number of data.
#'
#'           \item \code{custom}: In this case, users define the own weight vector.
#'           It should be noted that the weight vectors \eqn{<w_i>} should satisfy \eqn{w_i \in [0, 1]} and 
#'           their summation is 1.
#'           }
#'
#    \item \code{"gaussian.kernel"}: It means by using the gaussian kernel function based fuzzy rough sets
#            to calculate the approximations. It was introduced by (D. G. Chen et al, 2011). 
#            The following are fuzzy lower and upper approximations of a fuzzy set \eqn{A} in \eqn{U}
# 
# 			 \eqn{(R \downarrow A)(y) = inf_{x \in U} \mathcal{V}_{cos}(R(x,y), A(y))}
#
#           \eqn{(R \uparrow A)(y) = sup_{x \in U} \mathcal{T}_{cos}(R(x,y), A(y))}
# 
#            where \eqn{\mathcal{V}_{cos}(a,b)} is \eqn{ab + \sqrt{(1 - a^2)(1 - b^2)}} if \eqn{a > b} and
#            it is 1 when \eqn{a \le b}. And, 
#            \eqn{T_{cos}(a,b) = max(ab - \sqrt{1 - a^2}\sqrt{1 - b^2}, 0)}.  
#
#'    \item \code{"fvprs"}: It refers to fuzzy variable precision rough sets (FVPRS) introduced by 
#'            (Zhao et al, 2009). It is a combination between variable precision rough sets (VPRS)
#'             and FRST. This function implements the construction of lower and upper approximations as follows.
#'       
#'			  \eqn{(R_{\alpha} \downarrow A)(y) = inf_{A(y) \le \alpha} \mathcal{I}(R(x,y), \alpha) \wedge inf_{A(y) > \alpha} \mathcal{I}(R(x,y), A(y))}
#'
#'           \eqn{(R_{\alpha} \uparrow A)(y) = sup_{A(y) \ge N(\alpha)} \mathcal{T}(R(x,y), N(\alpha)) \vee sup_{A(y) < N(\alpha)} \mathcal{T}(R(x,y), A(y))}
#'             
#'           where \eqn{\alpha}, \eqn{\mathcal{I}} and \eqn{\mathcal{T}} are the variable precision parameter, implicator and t-norm operators, respectively.
#	  \item \code{"sfrs"}: It refers to soft fuzzy rough sets (SFRS) proposed by (Hu et al, 2010). It is inspired by idea of
#           soft-margin support vector machine (SVM) which can reduce the influence of noise. So, this concept introduced 
#           soft distance between an object \eqn{x} and a set of objects \eqn{Y}, which is defined as 
# 
#           \eqn{SD(x, Y) = arg_{d(x,y)} sup_{y \in Y}\{d(x,y) - \beta m_Y\}} 
#
#           where \eqn{d} is a distance function, \eqn{\beta} is a penalty factor and \eqn{m_Y = |\{y_i|d(x, y_i) < d(x, y)\}|} 
#
#           Therefore, let \eqn{U} be a nonempty universe, R be a fuzzy similarity relation on \eqn{U} and \eqn{F(U)} be
#           the fuzzy power set of \eqn{U}. The soft fuzzy lower and upper approximations of \eqn{A \in F(U)} are defined as
#
#           \eqn{(R_{\alpha} \downarrow A)(x) = 1 - R(x, arg_y sup_{A(y) \le A(y_L)}\{1 - R(x, y) - \beta m_{Y_L}\})},
#
#           \eqn{(R_{\alpha} \uparrow A)(x) = R(x, arg_y inf_{A(y) \ge A(y_U)}\{R(x,y) + \beta n_{Y_U}\})}
#
#           where
#
#           \eqn{Y_L = \{y|A(y) \le A(y_L), y \in U\}, y_L = arg_y inf_{y \in U} max\{1 - R(x,y), A(y)\}},
#
#           \eqn{Y_U = \{y|A(y) \ge A(y_U), y \in U\}, y_U = arg_y sup_{y \in U} min\{1 - R(x,y), A(y)\}}.
# 
#'    \item \code{"rfrs"}: It refers to robust fuzzy rough sets (RFRS) proposed by (Hu et al, 2012). 
#'            This package provides six types of RFRS which are k-trimmed minimum, k-mean minimum, k-median minimum, 
#'            k-trimmed maximum, k-mean maximum, and k-median maximum. 
#'            Basically, these methods are a special case of ordered weighted average (OWA) where they consider 
#'            the weight vectors as follows.
#'            \itemize{
#'            \item \code{"k.trimmed.min"}: \eqn{w_i^l = 1} for \eqn{i = n - k} and \eqn{w_i^l = 0} otherwise.
#'            \item \code{"k.mean.min"}: \eqn{w_i^l = 1/k} for \eqn{i > n - k} and \eqn{w_i^l = 0} otherwise.
#'            \item \code{"k.median.min"}: \eqn{w_i^l = 1} if k odd, \eqn{i = n - (k-1)/2} and \eqn{w_i^l = 1/2} if k even, \eqn{i = n - k/2} 
#'                                 and \eqn{w_i^l = 0} otherwise.
#'            \item \code{"k.trimmed.max"}: \eqn{w_i^u = 1} for \eqn{i = k + 1} and \eqn{w_i^u = 0} otherwise. 
#'            \item \code{"k.mean.max"}: \eqn{w_i^u = 1/k} for \eqn{i < k + 1} and \eqn{w_i^u = 0} otherwise.   
#'            \item \code{"k.median.max"}: \eqn{w_i^u = 1} if k odd, \eqn{i = (k + 1)/2} and \eqn{w_i^u = 1/2} if k even, \eqn{i = k/2 + 1} 
#'                                 or \eqn{w_i^u = 0} otherwise.
#'            }
#'
#'    \item \code{"beta.pfrs"}: It refers to \eqn{\beta}-precision fuzzy rough sets (\eqn{\beta}-PFRS) proposed by 
#'            (Salido and Murakami, 2003). This algorithm uses \eqn{\beta}-precision quasi-\eqn{\mathcal{T}}-norm and 
#'            \eqn{\beta}-precision quasi-\eqn{\mathcal{T}}-conorm. The following are the \eqn{\beta}-precision versions of fuzzy lower and upper approximations of a fuzzy set \eqn{A} in \eqn{U}
#' 
#' 				\eqn{(R_B \downarrow A)(y) = T_{\beta_{x \in U}} \mathcal{I}(R_B(x,y), A(x))}
#'
#' 				\eqn{(R_B \uparrow A)(y) = S_{\beta_{x \in U}} \mathcal{T}(R_B(x,y), A(x))} 
#'
#'            where \eqn{T_{\beta}} and \eqn{S_{\beta}} are \eqn{\beta}-precision quasi-\eqn{\mathcal{T}}-norm and \eqn{\beta}-precision quasi-\eqn{\mathcal{T}}-conorm.
#'            Given a t-norm \eqn{\mathcal{T}}, a t-conorm \eqn{S}, \eqn{\beta \in [0,1]} and \eqn{n \in N \setminus \{0, 1\}}, the corresponding 
#'            \eqn{\beta}-precision quasi-t-norm \eqn{T_{\beta}} and \eqn{\beta}-precision-\eqn{\mathcal{T}}-conorm \eqn{S_{\beta}} of order \eqn{n} are
#'            \eqn{[0,1]^n \to [0,1]} mappings such that for all \eqn{x = (x_1,...,x_n)} in \eqn{[0,1]^n},
#'
#'             \eqn{T_{\beta}(x) = \mathcal{T}(y_1,...,y_{n-m})},
#'
#'             \eqn{S_{\beta}(x) = \mathcal{T}(z_1,...,z_{n-p})},
#'
#'             where \eqn{y_i} is the \eqn{i^{th}} greatest element of \eqn{x} and \eqn{z_i} is the \eqn{i^{th}} smallest element of \eqn{x}, and
#'
#'             \eqn{m = max(i \in \{0,...,n\}|i \le (1-\beta)\sum_{j=1}^{n}x_j)},
#'
#'             \eqn{p = max(i \in \{0,...,n\}|i \le (1-\beta)\sum_{j=1}^{n}(a - x_j))}.
#' 
#'            In this package we use \code{min} and \code{max} for \eqn{\mathcal{T}}-norm and \eqn{\mathcal{T}}-conorm, respectively. 
#' 
#'    \item \code{"custom"}: It refers to user-defined lower and upper approximations. An example can be seen in Section \code{Examples}.
#'    }
#' The parameter \code{type.LU}, which is explained above, is related with parameter \code{control}. 
#' In other words, when choosing a specific value of \code{type.LU}, we should take into account to set values of related components in \code{control}.
#' The components that are considered depend on what kind of lower and upper approximations are used. 
#' So, we do not need to assign all components for a particular approach but only components related with \code{type.LU}.
#' The following is a list showing the components of each approaches.
#'                \itemize{
#'                \item \code{type.LU = "implicator.tnorm"}: 
#'
#'                 \code{control <- list(t.implicator, t.tnorm)}
#' 
#'                 \item \code{type.LU = "vqrs"}:
#'
#'                 \code{control <- list(q.some, q.most, type.aggregation, t.tnorm)}
#'
#'                 \item \code{type.LU = "owa"}:
#'
#'                 \code{control <- list(t.implicator, t.tnorm, m.owa)} 
#'                 
#'                 or
#'
#'                \code{control <- list(t.implicator, t.tnorm, w.owa)} 
#'
#                 \item \code{type.LU = "gaussian.kernel"}:
#
#                 \code{control <- list(type.relation)}
# 
#'                 \item \code{type.LU = "fvprs"}:
#'
#'                 \code{control <- list(t.implicator, t.tnorm, alpha)}
#'
#                 \item \code{type.LU = "sfrs"}:
#
#                 \code{control <- list(penalty.fact)}
#
#'                \item \code{type.LU = "beta.pfrs"}: 
#'
#'                 \code{control <- list(t.implicator, t.tnorm, beta.quasi)}
#'
#'                 \item \code{type.LU = "rfrs"}:
#'
#'                 \code{control <- list(t.implicator, t.tnorm, type.rfrs, k.rfrs)}
#'
#'                 \item \code{type.LU = "custom"}:
#'
#'                 \code{control <- list(t.implicator, t.tnorm, FUN.lower, FUN.upper)}
#'                }
#' The description of the components can be seen in the \code{control} parameter.
#' In Section \code{Examples}, we provide two examples showing different cases which are
#' when we have to handle a nominal decision attribute and a continuous one. 
#'
#' It should be noted that this function depends on \code{\link{BC.IND.relation.FRST}}
#' which is a function used to calculate the fuzzy indiscernibility relation as input data. 
#' So, it is obvious that before performing this function, users must execute \code{\link{BC.IND.relation.FRST}} first. 
#'
#' @title The fuzzy lower and upper approximations based on fuzzy rough set theory
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param IND.condAttr a \code{"IndiscernibilityRelation"} class of the conditional attributes which is produced by \code{\link{BC.IND.relation.FRST}}.
#' @param IND.decAttr a \code{"IndiscernibilityRelation"} class of the decision attribute which is produced by \code{\link{BC.IND.relation.FRST}}.
#' @param type.LU a string representing a chosen method to calculate lower and upper approximations. See the explanation in Section \code{Details}. 
#' @param control a list of other parameters. In order to understand how to express the \code{control} parameter, 
#'         see the explanation in Section \code{Details}. 
#'         The descriptions of those components and their values is as follows. 
#'         \itemize{ 
#'              \item \code{t.tnorm}: a type of triangular functions which have been explained 
#'
#'                   in \code{\link{BC.IND.relation.FRST}}.
#'        		\item \code{t.implicator}: a type of implicator functions.  
#'               	 The following are values of this parameter:               
#' 						\itemize{
#' 						\item \code{"kleene_dienes"} means \eqn{max(1 - x_1, x_2)}.
#' 						\item \code{"lukasiewicz"} means \eqn{min(1 - x_1 + x_2, 1)}. It is the default value. 
#' 						\item \code{"zadeh"} means \eqn{max(1 - x_1, min(x_1, x_2))}.
#' 						\item \code{"gaines"} means \eqn{(x_1 <= x_2 ? 1 : x_2 / x_1)}.
#' 						\item \code{"godel"} means \eqn{(x_1 <= x_2 ? 1 : x_2)}.
#' 						\item \code{"kleene_dienes_lukasiewicz"} means \eqn{1 - x_1 + x_1 * x_2}.
#' 						\item \code{"mizumoto"} means \eqn{(1 - x_1 + x_1 * x_2)}.
#' 						\item \code{"dubois_prade"} means \eqn{(x_2 == 0 ? 1 - x_1 : (x_1 == 1 ? x_2 : 1))}.
#' 						}
#'              	Where we consider the following rule: \eqn{x_1 -> x_2}. 
#'        		\item \code{q.some}: a vector of alpha and beta parameters of vaguely quantified rough set  
#'              	for quantifier \code{some}. The default value is \code{q.some = c(0.1, 0.6)}.
#'        		\item \code{q.most}: a vector of alpha and beta parameters of vaguely quantified rough set 
#'              	for quantifier \code{most}. The default value is \code{q.most = c(0.2, 1)}.
#'        		\item \code{alpha}: a numeric between 0 and 1 representing the threshold parameter of the fuzzy variable precision rough sets 
#'              	 (FVPRS) (see Section \code{Details}). The default value is 0.05.
#              \item \code{penalty.fact}: a numeric representing penalty factor which is used in soft fuzzy rough sets (SFRS) (see Section \code{Details}).
#                    The default value is 0.8.
#'              \item \code{m.owa}: an integer number (\eqn{m}) which is used in the OWA fuzzy rough sets (see Section \code{Details}). 
#'                    
#'                   The default value is \code{m.owa = round(0.5 * ncol(decision.table))}.
#'              \item \code{w.owa}: a vector representing the weight vector in the OWA fuzzy rough sets (see Section \code{Details}).
#'                    The default value is \code{NULL}, which means we use the \code{m.owa} type.
#'              \item \code{type.rfrs}: a type of robust fuzzy rough sets which is one of the following methods:
#'                    \code{"k.trimmed.min"}, \code{"k.mean.min"}, \code{"k.median.min"}, \code{"k.trimmed.max"},
#'                    \code{"k.mean.max"}, and \code{"k.median.max"} (see Section \code{Details}). The default value is \code{"k.trimmed.min"}.
#'              \item \code{k.rfrs}: a number between 0 and the number of data which is used to define considered data on 
#'                    robust fuzzy rough sets (RFRS) (see Section \code{Details}). The default value is 
#'                    \code{k.rfrs = round(0.5*nrow(decision.table))}.
#'              \item \code{beta.quasi}: a number between 0 and 1 representing \eqn{\beta}-precision t-norms and t-conorms in \eqn{\beta}-PFRS.
#'                    The default value is 0.05.
#'        }
#' @return A class \code{"LowerUpperApproximation"} representing fuzzy rough set (fuzzy lower and upper approximations). It contains the following components:
#'         \itemize{
#'          \item \code{fuzzy.lower}: a list showing the lower approximation classified 
#'                based on decision concepts for each index of objects. The value refers to
#'                the degree of objects included in the lower approximation.  
#'                In case the decision attribute is continuous, the result is in a data frame 
#'                with dimension (number of objects x number of objects) and the value on position \eqn{(i,j)} 
#'                shows the membership of object \eqn{i} to the lower approximation of the similarity class of object \eqn{j}.
#'          \item \code{fuzzy.upper}: a list showing the upper approximation classified 
#'                based on decision concepts for each index of objects. The value refers to
#'                the degree of objects included in the upper approximation. 
#'                In case the decision attribute is continuous values, the result is in data frame 
#'                with dimension (number of objects x number of objects) and the value on position \eqn{(i,j)} 
#'                shows the membership of object \eqn{i} to the upper approximation of the similarity class of object \eqn{j}.
#'          \item \code{type.LU}: a string representing the type of lower and upper approximation approaches.
#'          \item \code{type.model}: a string showing the type of model which is used. In this case, it is \code{"FRST"} which means fuzzy rough set theory.
#'          } 
#'         
#' @seealso \code{\link{BC.IND.relation.RST}}, \code{\link{BC.LU.approximation.RST}}, 
#'           and \code{\link{BC.positive.reg.FRST}}
#' @references
#' A. M. Radzikowska and E. E. Kerre, "A Comparative Study of Fuzzy Rough Sets", 
#' Fuzzy Sets and Systems, vol. 126, p. 137 - 156 (2002). 
#' 
#' C. Cornelis, M. De Cock, and A. Radzikowska, "Vaguely Quantified Rough Sets",
#' Proceedings of 11th International Conference on Rough Sets, Fuzzy Sets,
#' Data Mining and Granular Computing (RSFDGrC2007), Lecture Notes in
#' Artificial Intelligence 4482, p. 87 - 94 (2007).
#'
#' C. Cornelis, N. Verbiest, and R. Jensen, "Ordered Weighted Average Based Fuzzy
#' Rough Sets", Proceedings of the 5th International Conference on Rough Sets
#' and Knowledge Technology (RSKT 2010), p. 78 - 85 (2010).
#' 
# D. G. Chen, Q. H. Hu, Y. P. Yang, "Parameterized Attribute Reduction with 
# Gaussian Kernel Based Fuzzy Rough Sets", Information Sciences, vol. 181, no. 23,
# p. 5169 - 5179 (2011). 
#'
#' J. M. F. Salido and S. Murakami, "Rough Set Analysis of a General Type of Fuzzy Data
#' Using Transitive Aggregations of Fuzzy Similarity Relations", 
#' Fuzzy Sets Syst., vol. 139, p. 635 - 660 (2003).
#'
# Q. Hu, S. An and D. Yu, "Soft Fuzzy Rough Sets for Robust Feature Evaluation and Selection",
# Information Sciences, vol. 180, p. 4384 - 4400 (2010).
#'
#' Q. Hu, L. Zhang, S. An, D. Zhang, and D. Yu, "On Robust Fuzzy Rough Set Models",
#' IEEE Trans. on Fuzzy Systems, vol. 20, no. 4, p. 636 - 651 (2012).
#'
#' R. Jensen and Q. Shen,  
#' "New Approaches to Fuzzy-Rough Feature Selection", 
#' IEEE Trans. on Fuzzy Systems, vol. 19, no. 4,
#' p. 824 - 838 (2009).
#' 
#' R. R. Yager, "On Ordered Weighted Averaging Aggregation Operators in Multicriteria
#' Decision Making", IEEE Transactions on Systems, Man, and Cybernetics, vol. 18, p. 183 - 190 (1988).
#'
#' S. Y. Zhao, E. C. C. Tsang, and D. G. Chen, 
#' "The Model of Fuzzy Variable Precision Rough Sets",
#' IEEE Trans. Fuzzy Systems, vol. 17, no. 2,
#' p. 451 - 467 (2009).
#'
#' @examples
#' ###########################################################
#' ## 1. Example: Decision table contains nominal decision attribute
#' ## we are using the same dataset and indiscernibility 
#' ## relation along this example.
#' ###########################################################
#' dt.ex1 <- data.frame(c(-0.4, -0.4, -0.3, 0.3, 0.2, 0.2), 
#'                      c(-0.3, 0.2, -0.4, -0.3, -0.3, 0),
#'				        c(-0.5, -0.1, -0.3, 0, 0, 0),
#'				        c("no", "yes", "no", "yes", "yes", "no"))
#' colnames(dt.ex1) <- c("a", "b", "c", "d")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4)
#'
#' ## let us consider the first and second attributes 
#' ## only as conditional attributes
#' condAttr <- c(1, 2)
#' 
#' ## let us consider the fourth attribute as decision attribute
#' decAttr <- c(4)
#' 
#' #### calculate fuzzy indiscernibility relation ####
#' control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                     type.relation = c("tolerance", "eq.1"))
#' control.dec <- list(type.aggregation = c("crisp"), type.relation = "crisp")
#'
#' ## fuzzy indiscernibility relation of conditional attribute
#' IND.condAttr <- BC.IND.relation.FRST(decision.table, attributes = condAttr, 
#'                                      control = control.ind)
#'
#' ## fuzzy indiscernibility relation of decision attribute 
#' IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = decAttr, 
#'                                      control = control.dec)
#'
#' #### Calculate fuzzy lower and upper approximation using type.LU : "implicator.tnorm" ####	 
#' control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz")
#' FRST.LU <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "implicator.tnorm", control = control)
#'
#' #### Calculate fuzzy lower and upper approximation using type.LU : "vqrs" ####	 
#' control <- list(q.some = c(0.1, 0.6), q.most = c(0.2, 1), t.tnorm = "lukasiewicz")
#' FRST.VQRS <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "vqrs", control = control)
#'
#' #### Calculate fuzzy lower and upper approximation using type.LU : "owa" ####	 
#' control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", m.owa = 3) 
#' FRST.OWA.1 <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "owa", control = control)
#'
#' #### Calculate fuzzy lower and upper approximation using type.LU : 
#' #### "owa" with customized function 
#' #### In this case, we are using the same weight vector as
#' #### previous one with m.owa = 3
#' control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", 
#'                w.owa =  c(0, 0, 0, 0.14, 0.29, 0.57)) 
#' FRST.OWA.2 <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "owa", control = control)
#'
## #### Calculate fuzzy lower and upper approximation using type.LU : "gaussian.kernel" ####	 
## control <- list(type.relation = c("transitive.kernel", "gaussian"), delta = 0.2) 
## FRST.G <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
##               type.LU = "gaussian.kernel", control = control)
##
#' #### Calculate fuzzy lower and upper approximation using type.LU : "fvprs" ####	 
#' control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", alpha = 0.05)
#' FRST.fvprs <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "fvprs", control = control)
#'
# ### Calculate fuzzy lower and upper approximation using type.LU : "sfrs" ####	 
# control <- list(penalty.fact = 1)
# FRST.sfrs <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr,
#               type.LU = "sfrs", control = control)
#'
#' #### Calculate fuzzy lower and upper approximation using type.LU : "rfrs" ####	 
#' control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", 
#'                 type.rfrs = "k.trimmed.min", k.rfrs = 0)
#' FRST.rfrs <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "rfrs", control = control)
#'
#' #### Calculate fuzzy lower and upper approximation using type.LU : "beta.pfrs" ####	 
#' control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", beta.quasi = 0.05)
#' FRST.beta.pfrs <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "beta.pfrs", control = control)
#'
#' #### Calculate fuzzy lower and upper approximation using type.LU : "custom" ####	
#' ## In this case, we calculate approximations randomly. 
#' f.lower <- function(x){
#'         return(min(runif(1, min = 0, max = 1) * x))	
#' }
#' f.upper <- function(x){
#'         return(max(runif(1, min = 0, max = 1) * x))
#' }
#' control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", FUN.lower = f.lower, 
#'                 FUN.upper = f.upper)
#' FRST.custom <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "custom", control = control)
#' 
#' #### In this case, we use custom function for triangular norm and implicator operator
#' ## For example, let us define our implicator and t-norm operator as follows.
#' imp.lower <- function(antecedent, consequent){
#'	                 return(max(1 - antecedent, consequent))
#'               }
#' tnorm.upper <- function(x, y){
#'                 return (x * y)
#'              } 
#' control.custom <- list(t.implicator = imp.lower, t.tnorm = tnorm.upper)
#' FRST.LU.custom <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "implicator.tnorm", control = control.custom)
#'
#' ###########################################################
#' ## 2. Example: Decision table contains a continuous decision attribute.
#' ## It should be noted that in this example, we are using
#' ## the same dataset and indiscernibility relation.
#' ## We only show one method but for other approaches 
#' ## the procedure is analogous to the previous example
#' ###########################################################
#' ## In this case, we are using housing dataset containing 7 objects
#' data(RoughSetData)
#' decision.table <- RoughSetData$housing7.dt
#'
#' ## let us consider the first and second conditional attributes only,
#' ## and the decision attribute at 14.
#' cond.attributes <- c(1, 2)
#' dec.attributes <- c(14)
#' control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                type.relation = c("tolerance", "eq.1"))
#' IND.condAttr <- BC.IND.relation.FRST(decision.table, attributes = cond.attributes, 
#'                                      control = control.ind) 
#' IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = dec.attributes, 
#'                                     control = control.ind) 
#'
#' #### Calculate fuzzy lower and upper approximation using type.LU : "implicator.tnorm" ####	 
#' control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz")
#' FRST.LU <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "implicator.tnorm", control = control)
#' 
#' @export
BC.LU.approximation.FRST <- function(decision.table, IND.condAttr, IND.decAttr, type.LU = "implicator.tnorm", control = list()){
	############################################################
	### NOTE: The procedure in this function is based on 
	### R. Jensen and Q. Shen,  
	### "New Approaches to Fuzzy-Rough Feature Selection", 
	### IEEE Trans. on Fuzzy Systems, vol. 19, no. 4,
	### p. 824 - 838 (2009).
    ##########################################################	
	## set default values of all parameters
	control <- setDefaultParametersIfMissing(control, list(t.tnorm = "lukasiewicz", t.implicator = "lukasiewicz", q.some = c(0.1, 0.6), q.most = c(0.2, 1), 
										m.owa = 0.5 * ncol(decision.table), k.rfrs = round(0.5*nrow(decision.table)), 
										penalty.fact = 0.8, delta = 0.2, alpha = 0.05, beta.quasi = 0.05, type.rfrs = "k.trimmed.min", 
										FUN.lower = NULL, FUN.upper = NULL, w.owa = NULL))
		
	## get the data
	objects <- decision.table
	nominal.att <- attr(decision.table, "nominal.attrs")
	desc.attrs <- attr(decision.table, "desc.attrs")
	if (is.null(attr(decision.table, "decision.attr"))){
		decision.attr <- ncol(objects)
	}
	else {
		decision.attr <- attr(decision.table, "decision.attr")
		## check if index of decision attribute is not in last index, then change their position
		if (decision.attr != ncol(objects)){
			objects <- cbind(objects[, -decision.attr, drop = FALSE], objects[, decision.attr, drop = FALSE])
			nominal.att <- c(nominal.att[-decision.attr], nominal.att[decision.attr])
		}
	}
		
	#nominal.att <- decision.table$nominal.attrs
	t.tnorm <- control$t.tnorm	
	t.implicator <- control$t.implicator
	q.some <- control$q.some
	q.most <- control$q.most
	m.owa <- control$m.owa
	FUN.lower <- control$FUN.lower
	FUN.upper <- control$FUN.upper
	delta <- control$delta
	alpha <- control$alpha
	penalty.fact <- control$penalty.fact
	k.rfrs <- control$k.rfrs
	type.rfrs <- control$type.rfrs
	beta.quasi <- control$beta.quasi
	w.owa <- control$w.owa

	## check value of m
	if (m.owa > nrow(objects) && type.LU == "owa" && is.null(w.owa)){
		stop("please insert m.owa less than the number of instances")
	}

	## check for customized functions
	if (type.LU == "custom" && (is.null(FUN.lower) || is.null(FUN.upper))){
		stop("please insert your defined function for approximation equation")
	}
	
	## get indiscernibility relation
	if(!inherits(IND.condAttr, "IndiscernibilityRelation")) stop("not a legitimate IndiscernibilityRelation object")
	if(!inherits(IND.decAttr, "IndiscernibilityRelation")) stop("not a legitimate IndiscernibilityRelation object")
	miu.Ra <- IND.condAttr$IND.relation
	
	## if decision attribute is nominal values
	if (nominal.att[decision.attr] == TRUE){
		if (any(type.LU == c("implicator.tnorm", "owa", "rfrs", "custom", "beta.pfrs"))){
			## get the classes
			class.dec.attrs <- unique(objects[, decision.attr])	
	
			## initialization
			fuzzy.lower <- c()
			fuzzy.upper <- c()
			
			## looping based on decision concepts
			for (i in 1 : length(class.dec.attrs)){				
				## get indiscernibility relation for decision attribute toward to decision concepts
				## 1 means values are the same				
				temp.Indx <- which(objects[, decision.attr] == class.dec.attrs[i])
				new.val.dec.attrs <- as.numeric(IND.decAttr$IND.relation[temp.Indx[1], ,drop = FALSE])
				
				imp.val <- matrix(do.call(calc.implFunc, list(miu.Ra, new.val.dec.attrs, t.implicator)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
				tnorm.val <- matrix(do.call(func.tnorm, list(new.val.dec.attrs, miu.Ra, t.tnorm)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
					
				if (type.LU == "implicator.tnorm"){			
					## compute and save them into lower and upper
					## here, min == inf and max = sup					
					fuzzy.lower.Ra <- apply(imp.val, 2, min) 
					fuzzy.upper.Ra <- apply(tnorm.val, 2, max)
				}
				else if (type.LU == "owa"){
					## compute and save them into lower and upper
					fuzzy.lower.Ra <- apply(imp.val, 2, min) 
					fuzzy.upper.Ra <- apply(tnorm.val, 2, max)
					
					if (is.null(w.owa)){
						fuzzy.lower.Ra <- apply(imp.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "owa", type.apprx = "lower"))
						fuzzy.upper.Ra <- apply(tnorm.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "owa", type.apprx = "upper"))
					}
					else {
						fuzzy.lower.Ra <- apply(imp.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "owa", type.apprx = "lower", w.owa = w.owa))
						fuzzy.upper.Ra <- apply(tnorm.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "owa", type.apprx = "upper", w.owa = w.owa))
					}
				}
				else if (type.LU == "rfrs"){
					## compute and save them into lower and upper				
					fuzzy.lower.Ra <- apply(imp.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "rfrs", type.apprx = "lower"))
					fuzzy.upper.Ra <- apply(tnorm.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "rfrs", type.apprx = "upper"))
				}
				else if (type.LU == "custom"){	
					## compute and save them into lower and upper								
					fuzzy.lower.Ra <- apply(imp.val, 2, function(x) FUN.lower(x))
					fuzzy.upper.Ra <- apply(tnorm.val, 2, function(x) FUN.upper(x))
				}
				else if (type.LU == "beta.pfrs"){
					res.Temp <- calc.LU.betaPFRS(imp.val, tnorm.val, beta.quasi)
					fuzzy.lower.Ra <- res.Temp$fuzzy.lower.Ra
					fuzzy.upper.Ra <- res.Temp$fuzzy.upper.Ra
				}
				
				## give names
				names(fuzzy.lower.Ra) <- c(seq(1, nrow(miu.Ra)))
				names(fuzzy.upper.Ra) <- c(seq(1, nrow(miu.Ra)))
				
				## collect the values of approximation
				fuzzy.lower <- append(fuzzy.lower, list(fuzzy.lower.Ra))
				fuzzy.upper <- append(fuzzy.upper, list(fuzzy.upper.Ra))
				
			}
		}
		
		else if (type.LU == "sfrs"){
			## get the classes
			class.dec.attrs <- unique(objects[, decision.attr])	
	
			## initialization
			fuzzy.lower <- c()
			fuzzy.upper <- c()
			
			## looping based on decision concepts
			for (i in 1 : length(class.dec.attrs)){			
				
				miu.Ra.lower <- miu.Ra
				miu.Ra.upper <- miu.Ra
				miu.Ra.lower[, which(objects[, decision.attr] == class.dec.attrs[i])] <- NA
				miu.Ra.upper[, which(objects[, decision.attr] != class.dec.attrs[i])] <- NA
				
				## calculate soft distance
				indx.obj.lower <- apply(miu.Ra.lower, 1, function(x) calc.sd(x, penalty.fact = penalty.fact, type = "lower"))
				indx.obj.upper <- apply(miu.Ra.upper, 1, function(x) calc.sd(x, penalty.fact = penalty.fact, type = "upper"))
				
				## save lower and upper appr
				fuzzy.lower.Ra <- sapply(1:length(indx.obj.lower), function(x) 1 - miu.Ra[x, indx.obj.lower[x]])
				fuzzy.upper.Ra <- sapply(1:length(indx.obj.upper), function(x) miu.Ra[x, indx.obj.upper[x]])				
								
				## give names
				names(fuzzy.lower.Ra) <- c(seq(1, nrow(miu.Ra)))
				names(fuzzy.upper.Ra) <- c(seq(1, nrow(miu.Ra)))
				
				# collect the values of approximation
				fuzzy.lower <- append(fuzzy.lower, list(fuzzy.lower.Ra))
				fuzzy.upper <- append(fuzzy.upper, list(fuzzy.upper.Ra))
			}
		}
		
		else if (type.LU == "fvprs"){
			## get the classes
			class.dec.attrs <- unique(objects[, decision.attr])	
	
			## initialization
			fuzzy.lower <- c()
			fuzzy.upper <- c()
			
			## looping based on decision concepts
			for (i in 1 : length(class.dec.attrs)){	
				## get indiscernibility relation for decision attribute toward to decision concepts
				## 1 means values are the same				
				temp.Indx <- which(objects[, decision.attr] == class.dec.attrs[i])
				new.val.dec.attrs <- as.numeric(IND.decAttr$IND.relation[temp.Indx[1], ,drop = FALSE])

				lower.val.dec.attrs <- new.val.dec.attrs
				lower.val.dec.attrs[lower.val.dec.attrs <= alpha] <- alpha
			
				imp.val <- matrix(do.call(calc.implFunc, list(miu.Ra, lower.val.dec.attrs, t.implicator)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
				upper.val.dec.attrs <- new.val.dec.attrs
				upper.val.dec.attrs[upper.val.dec.attrs >= (1 - alpha)] <- (1 - alpha)
				
				tnorm.val <- matrix(do.call(func.tnorm, list(upper.val.dec.attrs, miu.Ra, t.tnorm)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
				
				## compute and save them into lower and upper	
				fuzzy.lower.Ra <- apply(imp.val, 2, min) 
				fuzzy.upper.Ra <- apply(tnorm.val, 2, max)	
				
				## give names
				names(fuzzy.lower.Ra) <- c(seq(1, nrow(miu.Ra)))
				names(fuzzy.upper.Ra) <- c(seq(1, nrow(miu.Ra)))
				
				## collect the values of approximation
				fuzzy.lower <- append(fuzzy.lower, list(fuzzy.lower.Ra))
				fuzzy.upper <- append(fuzzy.upper, list(fuzzy.upper.Ra))
			}
		}
		
		else if (type.LU == "gaussian.kernel"){
			## get the classes
			class.dec.attrs <- unique(objects[, decision.attr])	
	
			## initialization
			fuzzy.lower <- c()
			fuzzy.upper <- c()
			
			## looping based on decision concepts
			for (i in 1 : length(class.dec.attrs)){	
				## get indiscernibility relation for decision attribute toward to decision concepts
				## 1 means values are the same				
				temp.Indx <- which(objects[, decision.attr] == class.dec.attrs[i])
				new.val.dec.attrs <- as.numeric(IND.decAttr$IND.relation[temp.Indx[1], ,drop = FALSE])
				
				## calculate implication function for fuzzy lower appr
				imp.func <- function(x, y, M, N){
					if (M[x, y] <= N[y])
						imp.val <- 1
					else 
						imp.val <- M[x, y] * N[y] + sqrt((1 - M[x, y]^2) * (1 -  N[y]^2))
				}
				
				tnorm.val.func <- function(x, y, M, N){
					return (pmax(M[x, y] * N[y] - sqrt(1 - M[x, y]^2)*sqrt(1 - N[y]^2), 0))
				}
				
				Vec.imp.func <- Vectorize(imp.func,vectorize.args = c('x','y'))
				Vec.tnorm.val.func <- Vectorize(tnorm.val.func,vectorize.args = c('x','y'))
				
				imp.val <- outer(1:nrow(objects), 1:nrow(objects), Vec.imp.func, miu.Ra, new.val.dec.attrs)					
				tnorm.val <- outer(1:nrow(objects), 1:nrow(objects), Vec.tnorm.val.func, miu.Ra, new.val.dec.attrs)
				
				## compute and save them into lower and upper
				## here, min == inf and max = sup					
				fuzzy.lower.Ra <- apply(imp.val, 1, min) 
				fuzzy.upper.Ra <- apply(tnorm.val, 1, max)
				
				## give names
				names(fuzzy.lower.Ra) <- c(seq(1, nrow(miu.Ra)))
				names(fuzzy.upper.Ra) <- c(seq(1, nrow(miu.Ra)))
				
				## collect the values of approximation
				fuzzy.lower <- append(fuzzy.lower, list(fuzzy.lower.Ra))
				fuzzy.upper <- append(fuzzy.upper, list(fuzzy.upper.Ra))	
			}
		}
		
		else if (type.LU == "vqrs"){
			## get the classes
			class.dec.attrs <- unique(objects[, decision.attr])	
	
			## initialization
			fuzzy.lower <- c()
			fuzzy.upper <- c()
			
			## looping based on decision concepts
			for (i in 1 : length(class.dec.attrs)){	
			
				temp.Indx <- which(objects[, decision.attr] == class.dec.attrs[i])
				new.val.dec.attrs <- as.numeric(IND.decAttr$IND.relation[temp.Indx[1], ,drop = FALSE])
				
				## looping for each equivalence classes
				tnorm.val <- t(do.call(func.tnorm, list(new.val.dec.attrs, miu.Ra, t.tnorm)))
				
				sum.tnorm.val <- rowSums(tnorm.val)
				sum.miu.Ra <- rowSums(miu.Ra)
				dt <- sum.tnorm.val/sum.miu.Ra
				
				fuzzy.lower.Ra <- sapply(dt, function(x) if (is.infinite(x)) return(0) else return(calc.vqrs.relation(x, q.most)))
				fuzzy.upper.Ra <- sapply(dt, function(x) if (is.infinite(x)) return(0) else return(calc.vqrs.relation(x, q.some)))
				
				## give names
				names(fuzzy.lower.Ra) <- c(seq(1, nrow(miu.Ra)))
				names(fuzzy.upper.Ra) <- c(seq(1, nrow(miu.Ra)))
				
				## collect the values of approximation
				fuzzy.lower <- append(fuzzy.lower, list(fuzzy.lower.Ra))
				fuzzy.upper <- append(fuzzy.upper, list(fuzzy.upper.Ra))	
			}
		}
		
		names(fuzzy.lower) <- c(as.character(unique(objects[, decision.attr])))
		names(fuzzy.upper) <- c(as.character(unique(objects[, decision.attr])))				
	}
	else {
		IND.dec.att <- IND.decAttr$IND.relation
		
		if (any(type.LU == c("implicator.tnorm", "owa", "rfrs", "custom", "fvprs", "gaussian.kernel", "vqrs", "beta.pfrs"))){
			## get the classes
			## initialization
			fuzzy.lower <- c()
			fuzzy.upper <- c()
			
			## looping based on decision concepts
			for (i in 1 : nrow(IND.dec.att)){				
				## initialization
				new.val.dec.attrs <- IND.dec.att[i, , drop = FALSE]
				
				if (type.LU == "implicator.tnorm"){
					imp.val <- matrix(do.call(calc.implFunc, list(miu.Ra, c(new.val.dec.attrs), t.implicator)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
					tnorm.val <- matrix(do.call(func.tnorm, list(c(new.val.dec.attrs), miu.Ra, t.tnorm)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
					
					## compute and save them into lower and upper
					## here, min == inf and max = sup					
					fuzzy.lower.Ra <- apply(imp.val, 2, min) 
					fuzzy.upper.Ra <- apply(tnorm.val, 2, max)
				} 
				else if (type.LU == "owa"){
					imp.val <- matrix(do.call(calc.implFunc, list(miu.Ra, c(new.val.dec.attrs), t.implicator)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))				
					tnorm.val <- matrix(do.call(func.tnorm, list(c(new.val.dec.attrs), miu.Ra, t.tnorm)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
					
					## compute and save them into lower and upper
					fuzzy.lower.Ra <- apply(imp.val, 2, min) 
					fuzzy.upper.Ra <- apply(tnorm.val, 2, max)
					
					if (is.null(w.owa)){
						fuzzy.lower.Ra <- apply(imp.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "owa", type.apprx = "lower"))
						fuzzy.upper.Ra <- apply(tnorm.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "owa", type.apprx = "upper"))
					}
					else {
						fuzzy.lower.Ra <- apply(imp.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "owa", type.apprx = "lower", w.owa = w.owa))
						fuzzy.upper.Ra <- apply(tnorm.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "owa", type.apprx = "upper", w.owa = w.owa))
					}
				}
				else if (type.LU == "rfrs"){
					imp.val <- matrix(do.call(calc.implFunc, list(miu.Ra, c(new.val.dec.attrs), t.implicator)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))				
					tnorm.val <- matrix(do.call(func.tnorm, list(c(new.val.dec.attrs), miu.Ra, t.tnorm)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
					
					## compute and save them into lower and upper				
					fuzzy.lower.Ra <- apply(imp.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "rfrs", type.apprx = "lower"))
					fuzzy.upper.Ra <- apply(tnorm.val, 2, function (x) calc.OWA(x, m.owa = m.owa, type.method = "rfrs", type.apprx = "upper"))
				}
				else if (type.LU == "custom"){
					imp.val <- matrix(do.call(calc.implFunc, list(miu.Ra, c(new.val.dec.attrs), t.implicator)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))				
					tnorm.val <- matrix(do.call(func.tnorm, list(c(new.val.dec.attrs), miu.Ra, t.tnorm)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
					
					## compute and save them into lower and upper								
					fuzzy.lower.Ra <- apply(imp.val, 2, function(x) FUN.lower(x))
					fuzzy.upper.Ra <- apply(tnorm.val, 2, function(x) FUN.upper(x))
				}
				else if (type.LU == "fvprs"){
					lower.val.dec.attrs <- new.val.dec.attrs
					lower.val.dec.attrs[lower.val.dec.attrs <= alpha] <- alpha
				
					imp.val <- matrix(do.call(calc.implFunc, list(miu.Ra, c(lower.val.dec.attrs), t.implicator)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
					upper.val.dec.attrs <- new.val.dec.attrs
					upper.val.dec.attrs[upper.val.dec.attrs >= (1 - alpha)] <- (1 - alpha)
					
					tnorm.val <- matrix(do.call(func.tnorm, list(c(upper.val.dec.attrs), miu.Ra, t.tnorm)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
					
					## compute and save them into lower and upper	
					fuzzy.lower.Ra <- apply(imp.val, 2, min) 
					fuzzy.upper.Ra <- apply(tnorm.val, 2, max)		
				}
				else if (type.LU == "gaussian.kernel"){
					## calculate implication function for fuzzy lower appr
					imp.func <- function(x, y, M, N){
						if (M[x, y] <= N[y])
							imp.val <- 1
						else 
							imp.val <- M[x, y] * N[y] + sqrt((1 - M[x, y]^2) * (1 -  N[y]^2))
					}
					
					tnorm.val.func <- function(x, y, M, N){
						return (pmax(M[x, y] * N[y] - sqrt(1 - M[x, y]^2)*sqrt(1 - N[y]^2), 0))
					}
					
					Vec.imp.func <- Vectorize(imp.func,vectorize.args = c('x','y'))
					Vec.tnorm.val.func <- Vectorize(tnorm.val.func,vectorize.args = c('x','y'))
					
					imp.val <- outer(1:nrow(objects), 1:nrow(objects), Vec.imp.func, miu.Ra, c(new.val.dec.attrs))
					
					tnorm.val <- outer(1:nrow(objects), 1:nrow(objects), Vec.tnorm.val.func, miu.Ra, c(new.val.dec.attrs))
					
					## compute and save them into lower and upper
					## here, min == inf and max = sup					
					fuzzy.lower.Ra <- apply(imp.val, 1, min) 
					fuzzy.upper.Ra <- apply(tnorm.val, 1, max)
				}
				else if (type.LU == "vqrs"){
					## looping for each equivalence classes
					tnorm.val <- t(do.call(func.tnorm, list(c(new.val.dec.attrs), miu.Ra, t.tnorm)))
					
					sum.tnorm.val <- rowSums(tnorm.val)
					sum.miu.Ra <- rowSums(miu.Ra)
					dt <- sum.tnorm.val/sum.miu.Ra
					
					fuzzy.lower.Ra <- sapply(dt, function(x) if (is.infinite(x)) return(0) else return(calc.vqrs.relation(x, q.most)))
					fuzzy.upper.Ra <- sapply(dt, function(x) if (is.infinite(x)) return(0) else return(calc.vqrs.relation(x, q.some)))
				}
				else if (type.LU == "beta.pfrs"){
					imp.val <- matrix(do.call(calc.implFunc, list(miu.Ra, c(new.val.dec.attrs), t.implicator)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
					tnorm.val <- matrix(do.call(func.tnorm, list(c(new.val.dec.attrs), miu.Ra, t.tnorm)), nrow = nrow(miu.Ra), ncol = ncol(miu.Ra))
					
					res.Temp <- calc.LU.betaPFRS(imp.val, tnorm.val, beta.quasi)
					fuzzy.lower.Ra <- res.Temp$fuzzy.lower.Ra
					fuzzy.upper.Ra <- res.Temp$fuzzy.upper.Ra
				}
				
				## collect the values of approximation
				fuzzy.lower <- rbind(fuzzy.lower, fuzzy.lower.Ra)
				fuzzy.upper <- rbind(fuzzy.upper, fuzzy.upper.Ra)
			}
		}
		
		else if (type.LU == c("sfrs")){		
			stop("Currently, the package doesn't support for continuous decision attribute")
		}
		
		rownames(fuzzy.lower) <- NULL
		rownames(fuzzy.upper) <- NULL
	}	
	
	## build class
	res <- list(fuzzy.lower = fuzzy.lower, fuzzy.upper = fuzzy.upper, type.LU = type.LU, type.model = "FRST")
	class.mod <- ObjectFactory(res, classname = "LowerUpperApproximation")
	
	return(class.mod)
}


#' This is a function that implements a fundamental concept of fuzzy rough set theory which is
#' the positive region and the corresponding degree of dependency. The explanation about this concept can be seen 
#' in \code{\link{B.Introduction-FuzzyRoughSets}}.
#' 
#' In order to compute the function, we need to calculate the indiscernibility relation by executing \code{\link{BC.IND.relation.FRST}} 
#' and the lower and upper approximations by calling \code{\link{BC.LU.approximation.FRST}}.
#'
#' @title Positive region based on fuzzy rough set
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#' @param fuzzyroughset a \code{"LowerUpperApproximation"} class representing a fuzzy rough set that is produced by \code{\link{BC.LU.approximation.FRST}}.
#' @seealso \code{\link{BC.LU.approximation.FRST}}, \code{\link{BC.IND.relation.FRST}}, \code{\link{BC.IND.relation.RST}}, 
#'
#' \code{\link{BC.LU.approximation.RST}}, and \code{\link{BC.positive.reg.FRST}}.
#' @return A class \code{"PositiveRegion"} containing the following components:
#'         \itemize{
#'         \item \code{positive.freg}: a vector representing membership degrees to the fuzzy positive region for each index of objects.
#'         \item \code{degree.dependency}: a value expressing the degree of dependency. 
#'         \item \code{type.model}: a string representing type of models. In this case, it is \code{"FRST"} which means fuzzy rough set theory.
#'         }
#' @references
#' R. Jensen and Q. Shen,  
#' "New Approaches to Fuzzy-Rough Feature Selection", 
#' IEEE Trans. on Fuzzy Systems, vol. 19, no. 4,
#' p. 824 - 838 (2009).
#
#' @examples
#' ###########################################################
#' ##### 1. Example: Using a simple decision table containing 
#' #####             nominal values for the decision attribute
#' ###########################################################
#' dt.ex1 <- data.frame(c(-0.4, -0.4, -0.3, 0.3, 0.2, 0.2), 
#'                      c(-0.3, 0.2, -0.4, -0.3, -0.3, 0),
#'				        c(-0.5, -0.1, -0.3, 0, 0, 0),
#'				        c("no", "yes", "no", "yes", "yes", "no"))
#' colnames(dt.ex1) <- c("a", "b", "c", "d")
#' decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4)
#'
#' ## let us consider the first and second attributes only as conditional attribute
#' condAttr <- c(1, 2)
#' 
#' ## let us consider the fourth attribute as decision attribute
#' decAttr <- c(4)
#' 
#' #### Calculate fuzzy indiscernibility relation ####
#' control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                     type.relation = c("tolerance", "eq.1"))
#' control.dec <- list(type.aggregation = c("crisp"), type.relation = "crisp")
#'
#' IND.condAttr <- BC.IND.relation.FRST(decision.table, attributes = condAttr, 
#'                                      control = control.ind) 
#' IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = decAttr, 
#'                                      control = control.dec) 
#' 
#' #### Calculate fuzzy lower and upper approximation using type.LU : 
#' #### "implicator.tnorm" 
#' control <- list(t.implicator = "lukasiewicz")
#' FRST.LU <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "implicator.tnorm", control = control)
#'
#' #### Determine positive regions ####
#' res.1 <- BC.positive.reg.FRST(decision.table, FRST.LU)
#' 
#' ###########################################################
#' ##### 2. Example: Using the housing decision table containing 
#' #####             continuous values for the decision attribute
#' ###########################################################
#'
#' ## In this case, we are using the housing dataset containing 7 objects
#' data(RoughSetData)
#' decision.table <- RoughSetData$housing7.dt
#' 
#' conditional.attr <- c(1, 2)
#' decision.attr = c(14)
#' control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                      type.relation = c("tolerance", "eq.1"))
#'
#' #### Calculate fuzzy indiscernibility relation ####
#' IND.condAttr <- BC.IND.relation.FRST(decision.table, attributes = conditional.attr, 
#'                                      control = control.ind) 
#' IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = decision.attr, 
#'                                      control = control.ind) 
#'
#' #### Calculate fuzzy lower and upper approximation using type.LU : 
#' #### "implicator.tnorm" 
#' control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz")
#'
#' FRST.LU <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "implicator.tnorm", control = control)
#'
#' #### Determine fuzzy regions ####
#' res.2 <- BC.positive.reg.FRST(decision.table, FRST.LU)
#'
#' @export
BC.positive.reg.FRST <- function(decision.table, fuzzyroughset){
	
	res <- BC.def.region.FRST(decision.table, fuzzyroughset)
	positive.freg <- res$positive.freg
	degree.dependency <- res$degree.dependency

	res <- list(positive.freg = positive.freg, degree.dependency = degree.dependency, type.model = "FRST")
	
	## build class
	class.mod <- ObjectFactory(res, classname = "PositiveRegion")

	return(class.mod)
}

#' This is a function that is used to build the decision-relative discernibility matrix based on FRST.
#' It is a matrix whose elements contain discernible attributes among pairs of objects. 
#' By means of this matrix, we are able to produce all decision reducts of the given decision system.
#' 
#' In this function, we provide several approaches in order to generate the decision-relative discernibility matrix. 
#' Theoretically, all reducts are found by constructing 
#' the matrix that contains elements showing discernible attributes among objects. 
#' The discernible attributes are determined by a specific condition which depends on the selected algorithm. A particular approach can be executed by selecting
#' a value of the parameter \code{type.discernibility}. The following shows the different 
#' values of the parameter \code{type.discernibility} corresponding approaches considered in this function.
#' \itemize{
#' \item \code{"standard.red"}: It is adopted from (Tsang et al, 2008)'s approach. 
#' The concept has been explained briefly in \code{\link{B.Introduction-FuzzyRoughSets}}. 
#' In order to use this algorithm, we assign the \code{control} parameter
#'        with the following components:
#' 
#'        \code{control = list(type.aggregation, type.relation, type.LU, t.implicator)}
#'
#'        The detailed description of the components can be seen in \code{\link{BC.IND.relation.FRST}} and 
#'
#'        \code{\link{BC.LU.approximation.FRST}}. Furthermore, in this case the authors suggest to use the "min" t-norm  
#'       (i.e., \code{type.aggregation = c("t.tnorm", "min")}) and the implicator operator "kleene_dienes" (i.e., \code{t.implicator = "kleene_dienes"}).
#'   
#' \item \code{"alpha.red"}: It is based on (Zhao et al, 2009)'s approach where all reductions will 
#'       be found by building an \eqn{\alpha}-discernibility matrix. This matrix contains elements which are defined by
#'
#'       1) if \eqn{x_i} and \eqn{x_j} belong to different decision concept,
#'
#'       \eqn{c_{ij} = \{R : \mathcal{T}(R(x_i, x_j), \lambda) \le \alpha \}},
#'
#'       where \eqn{\lambda = (R_{\alpha} \downarrow A)(u)} which is lower approximation 
#'       of FVPRS (See \code{\link{BC.LU.approximation.FRST}}). 
#'
#'       2) \eqn{c_{ij}={\oslash}}, otherwise.
#'
#'       To generate the discernibility matrix based on this approach, we use the \code{control} parameter
#'       with the following components:
#' 
#'        \code{control = list(type.aggregation, type.relation, t.implicator, alpha.precision)} 
#'
#'        where the lower approximation \eqn{\lambda} is fixed to \code{type.LU = "fvprs"}. The detailed
#'        description of the components can be seen in \code{\link{BC.IND.relation.FRST}} and \code{\link{BC.LU.approximation.FRST}}.
#'        Furthermore, in this case the authors suggest to use \eqn{\mathcal{T}}-similarity relation 
#'
#'        e.g., \code{type.relation = c("transitive.closure", "eq.3")},
#'
#'        the "lukasiewicz" t-norm  (i.e., \code{type.aggregation = c("t.tnorm", "lukasiewicz")}), and \code{alpha.precision} from 0 to 0.5.
#'    
# \item \code{"consistence.degree"}: It is based on (E. C. C. Tsang and S. Y. Zhao, 2010)'s approach. This algorithm defines the 
#        discernibility matrix of \eqn{M_D(U, A)} where \eqn{M} is a \eqn{n \times n} matrix that constitutes the following \eqn{(c_{ij})}.
#        Suppose, \eqn{(U, A \cup D)} is a decision table/system, 
#        
#        1) if \eqn{x_i} and \eqn{x_j} belong to different decision concept,
#
#         \eqn{c_{ij}= \{R : \mathcal{T}(R(x_i, x_j), \lambda)=0\}},
# 
#        where \eqn{\lambda = inf_{u \in U}\mathcal{I}(R(x_i, u), A(u))}. 
#
#        2) \eqn{c_{ij}={\oslash}}, otherwise.
#
#        To generate fuzzy relation \eqn{R} and lower approximation \eqn{\lambda}, we use the \code{control} parameter
#        with the following components:
# 
#        \code{control = list(type.aggregation, type.relation, type.LU, t.implicator)}. 
#
#        The detailed description of the components can be seen in \code{\link{BC.IND.relation.FRST}} and 
#
#        \code{\link{BC.LU.approximation.FRST}}.
#
#' \item \code{"gaussian.red"}: It is based on (Chen et al, 2011)'s approach. The discernibility matrix contains elements which are defined by: 
#'        
#'        1) if \eqn{x_i} and \eqn{x_j} belong to different decision concept, 
#'
#'         \eqn{c_{ij}= \{R : R(x_i, x_j) \le \sqrt{1 - \lambda^2(x_i)}\}},
#' 
#'        where \eqn{\lambda = inf_{u \in U}\mathcal{I}_{cos}(R(x_i, u), A(u)) - \epsilon}. To generate fuzzy relation \eqn{R} , we use the fixed parameters as follows:
#' 
#'        \code{t.tnorm = "t.cos"} and \code{type.relation = c("transitive.kernel", "gaussian")}. 
#'        
#'        2) \eqn{c_{ij}={\oslash}}, otherwise.
#' 
#'        In this case, we need to define \code{control} parameter as follows.
#'
#'        \code{control <- list(epsilon)}
#'
#'       It should be noted that when having nominal values on all attributes then \code{epsilon} (\eqn{\epsilon}) should be 0. 
#'
#' \item \code{"min.element"}: It is based on (Chen et al, 2012)'s approach where we only consider finding 
#'       the minimal element of the discernibility matrix by introducing the binary relation \eqn{DIS(R)} the relative discernibility relation 
#'      of conditional attribute \eqn{R} with respect to decision attribute \eqn{d}, which is computed as
#'      
#'      \eqn{DIS(R) = \{(x_i, x_j) \in U \times U: 1 - R(x_i, x_j) > \lambda_i, x_j \notin [x_i]_d\}},
#'
#'      where \eqn{\lambda_i = (Sim(R) \downarrow [x_i]_d)(x_i)} with \eqn{Sim(R)} a fuzzy equivalence relation. 
#'      In other words, this algorithm does not need to build the discernibility matrix. 
#'      To generate the fuzzy relation \eqn{R} and lower approximation \eqn{\lambda}, we use the \code{control} parameter
#'        with the following components:
#' 
#'        \code{control = list(type.aggregation, type.relation, type.LU, t.implicator)}. 
#'
#'        The detailed description of the components can be seen in \code{\link{BC.IND.relation.FRST}} and 
#'
#'        \code{\link{BC.LU.approximation.FRST}}.
#' }
#'
#' @title The decision-relative discernibility matrix based on fuzzy rough set theory
#' @author Lala Septem Riza
#'
#' @param decision.table a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
#'        It should be noted that this case only supports the nominal/symbolic decision attribute.
#' @param type.discernibility a string representing a type of discernibility. See in Section \code{Details}. 
#' @param control a list of other parameters.
#'        \itemize{
#'        \item \code{type.relation}: a type of fuzzy indiscernibility relation. The default value is \code{type.relation = c("tolerance", "eq.1")}.
#'
#'        		 See \code{\link{BC.IND.relation.FRST}}.
#'        \item \code{type.aggregation}: a type of aggregation operator. The default value is \code{type.aggregation = c("t.tnorm", "lukasiewicz")}.
#'
#'       		 See \code{\link{BC.IND.relation.FRST}}.
#'        \item \code{t.implicator}: a type of implicator operator. The default value is \code{"lukasiewicz"}.
#'
#'        		See \code{\link{BC.LU.approximation.FRST}}.
#'        \item \code{type.LU}: a type of method of lower and upper approximations. The default value is \code{"implicator.tnorm"}.
#'        
#'        		See \code{\link{BC.LU.approximation.FRST}}.
#'        \item \code{alpha.precision}: a numeric value representing a precision variable. It is used when using \code{"alpha.red"} as \code{type.discernibility}.
#'              The default value is 0.05.
#'
#'         		See \code{\link{BC.LU.approximation.FRST}}.  
#'        \item \code{show.discernibilityMatrix}: a boolean value determining whether the discernibility matrix will be shown or not (NULL). The default value is \code{FALSE}.
#'        \item \code{epsilon}: a numeric between 0 and 1 representing the \eqn{\epsilon} value on 
#'
#'             \code{type.discernibility = "gaussian.red"}. It should be noted that when having nominal values on all attributes then \eqn{\epsilon} should be 0. 
#'             The default value is 0.
#'        \item \code{delta}: a numeric representing the \eqn{\delta} value on \code{"gaussian"} equations 
#'
#'             (see \code{\link{BC.IND.relation.FRST}}). The default value is 2.
#' 		  \item \code{range.object}: a vector representing considered objects to construct the \code{k}-relative discernibility matrix. 
#'                The default value is \code{NULL} which means that we are using all objects in the decision table.
#'        }
#' @return A class \code{"DiscernibilityMatrix"} containing the following components: 
#' \itemize{
#' \item \code{disc.mat}: a matrix showing the decision-relative discernibility matrix \eqn{M(\mathcal{A})} 
#'        which contains \eqn{n \times n} where \eqn{n} is the number of objects. It will be printed when choosing 
#'
#'        \code{show.discernibilityMatrix = TRUE}.
#' \item \code{disc.list}: the decision-relative discernibility represented in a list.
#' \item \code{discernibility.type}: a string showing the chosen type of discernibility methods.
#' \item \code{type.model}: in this case, it is \code{"FRST"}.
#' }
#' @seealso \code{\link{BC.discernibility.mat.RST}}, \code{\link{BC.LU.approximation.RST}}, and \code{\link{BC.LU.approximation.FRST}}
#' @references
#' D. Chen, L. Zhang, S. Zhao, Q. Hu, and P. Zhu, "A Novel Algorithm for Finding Reducts 
#' with Fuzzy Rough Sets", IEEE Trans. on Fuzzy Systems, vol. 20, no. 2, p. 385 - 389 (2012). 
#'
#' D. G. Chen, Q. H. Hu, and Y. P. Yang, "Parameterized Attribute Reduction with
#' Gaussian Kernel Based Fuzzy Rough Sets", Information Sciences, vol. 181, no. 23, 
#' p. 5169 - 5179 (2011).
#'
#' E. C. C. Tsang, D. G. Chen, D. S. Yeung, X. Z. Wang, and J. W. T. Lee, 
#' "Attributes Reduction Using Fuzzy Rough Sets", IEEE Trans. Fuzzy Syst., 
#' vol. 16, no. 5, p. 1130 - 1141 (2008).
#'
# E. C. C. Tsang and S. Y. Zhao, "Decision Table Reduction in KDD: Fuzzy Rough Based Approach",
# Transaction on Rough Sets, Lecture Notes in Computer Sciences, vol. 5946, p. 177 - 188 (2010).
#'
#' S. Zhao, E. C. C. Tsang, and D. Chen, "The Model of Fuzzy Variable Precision Rough Sets",
#' IEEE Trans. on Fuzzy Systems, vol. 17, no. 2, p. 451 - 467 (2009).
#'
#' @examples
#' #######################################################################
#' ## Example 1: Constructing the decision-relative discernibility matrix
#' ## In this case, we are using The simple Pima dataset containing 7 rows. 
#' #######################################################################
#' data(RoughSetData)
#' decision.table <- RoughSetData$pima7.dt 
#' 
#' ## using "standard.red"
#' control.1 <- list(type.relation = c("tolerance", "eq.1"), 
#'                 type.aggregation = c("t.tnorm", "min"), 
#'                 t.implicator = "kleene_dienes", type.LU = "implicator.tnorm")
#' res.1 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "standard.red", 
#'                                     control = control.1)
#' 
#' ## using "gaussian.red"
#' control.2 <- list(epsilon = 0)
#' res.2 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "gaussian.red",
#'                                     control = control.2)
#'
# control.3 <- list(type.relation = c("tolerance", "eq.1"), 
#                 type.aggregation = c("t.tnorm", "lukasiewicz"),
#                 t.implicator = "lukasiewicz", type.LU = "implicator.tnorm")
# res.3 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "consistence.degree",
#                                     control = control.3) 
#' ## using "alpha.red"
#' control.3 <- list(type.relation = c("tolerance", "eq.1"), 
#'                 type.aggregation = c("t.tnorm", "min"),
#'                 t.implicator = "lukasiewicz", alpha.precision = 0.05)
#' res.3 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "alpha.red", 
#'                                     control = control.3)
#'
#' ## using "min.element"
#' control.4 <- list(type.relation = c("tolerance", "eq.1"), 
#'                 type.aggregation = c("t.tnorm", "lukasiewicz"),
#'                 t.implicator = "lukasiewicz", type.LU = "implicator.tnorm")
#' res.4 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "min.element", 
#'                                     control = control.4)
#' 
#' #######################################################################
#' ## Example 2: Constructing the decision-relative discernibility matrix
#' ## In this case, we are using the Hiring dataset containing nominal values
#' #######################################################################
#' \dontrun{data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt 
#'
#' control.1 <- list(type.relation = c("crisp"), 
#'                 type.aggregation = c("crisp"),
#'                 t.implicator = "lukasiewicz", type.LU = "implicator.tnorm")
#' res.1 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "standard.red", 
#'                                     control = control.1)
#' 
#' control.2 <- list(epsilon = 0)
#' res.2 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "gaussian.red",
#'                                     control = control.2)}
#' @export
BC.discernibility.mat.FRST <- function(decision.table, type.discernibility = "standard.red", control = list()){

	if(!inherits(decision.table, "DecisionTable")) {
		stop("Provided data should inherit from the \'DecisionTable\' class.")
	}
  
	if (is.null(attr(decision.table, "decision.attr"))){
		stop("A decision attribute is not indicated.")
	}
	
	if (attr(decision.table, "nominal.attrs")[attr(decision.table, "decision.attr")] == FALSE){
		stop("The decision attribute must be nominal values")
	}
  
	control <- setDefaultParametersIfMissing(control, list(type.relation = c("tolerance", "eq.1"), type.aggregation = c("t.tnorm", "lukasiewicz"),
                                     	t.implicator = "lukasiewicz", type.LU = "implicator.tnorm", alpha.precision = 0.05, 
										show.discernibilityMatrix = FALSE, epsilon = 0, delta = 2, range.object = NULL))
	
	type.relation <- control$type.relation
	type.aggregation <- ch.typeAggregation(control$type.aggregation)
	t.implicator <- control$t.implicator
	type.LU <- control$type.LU
	alpha.precision <- control$alpha.precision
	show.discernibilityMatrix <- control$show.discernibilityMatrix
	epsilon <- control$epsilon
	delta <- control$delta
	range.object <- control$range.object
	if (!is.null(range.object)){
		decision.table <- decision.table[range.object, ,drop=FALSE]
	}
	
	if (type.discernibility == "min.element"){
		return(min.disc.mat.FRST(decision.table, t.tnorm = type.aggregation[2], type.relation, t.implicator, type.LU))
	}
	else{
		## build decision-relative discernibility matrix
		res.temp <- build.discMatrix.FRST(decision.table, type.discernibility, num.red = "all", alpha.precision, type.relation, 
		                               t.implicator, type.LU, show.discernibilityMatrix, epsilon = epsilon, delta = delta, type.aggregation = type.aggregation)
		disc.mat = res.temp$disc.mat
		disc.list = res.temp$disc.list
		
		disc.list = unique(disc.list)
		discernibilityMatrix = list(disc.mat = disc.mat, disc.list = disc.list, 
                                names.attr = colnames(decision.table), type.discernibility = type.discernibility, type.model = "FRST")
		discernibilityMatrix = ObjectFactory(discernibilityMatrix, classname = "DiscernibilityMatrix")
		return(discernibilityMatrix)
	}	
}
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
#' This part introduces briefly fuzzy rough set theory (FRST) and its application to data analysis. 
#' Since recently there are a lot of FRST variants that have been
#' proposed by researchers, in this introduction, we only provide some basic concepts of FRST based on (Radzikowska and Kerre, 2002). 
#'
#' Just like in RST (see \code{\link{A.Introduction-RoughSets}}),
#' a data set is represented as a table called an information system \eqn{\mathcal{A} = (U, A)}, where
#' \eqn{U} is a non-empty set of finite objects as the universe of discourse (note: it refers to all instances/experiments/rows 
#' in datasets) and \eqn{A} is a non-empty finite set of attributes, such that \eqn{a : U \to V_{a}} for every \eqn{a \in A}. 
#' The set \eqn{V_{a}} is the set of values that attribute \eqn{a} may take. Information systems that involve a decision attribute, 
#' containing classes or decision values of each objects, are called decision systems (or said as decision tables). More formally, it is a pair \eqn{\mathcal{A} = (U, A \cup \{d\})},
#' where \eqn{d \notin A} is the decision attribute. The elements of \eqn{A} are called conditional attributes. However, different from RST, FRST has several ways
#' to express indiscernibility. 
#' 
#' Fuzzy indiscernibility relation (FIR) is used for any fuzzy relation that determines the degree to which two objects are indiscernible. 
#' We consider some special cases of FIR.
#' \itemize{
#' \item fuzzy tolerance relation: this relation has properties which are reflexive and symmetric where
#' 
#'       reflexive: \eqn{R(x,x) = 1}
#'
#'       symmetric: \eqn{R(x,y) = R(y,x)}
#'
#' \item similarity relation (also called fuzzy equivalence relation): this relation has properties not only reflexive and symmetric but also
#'       transitive defined as
#'
#'       \eqn{min(R(x,y), R(y,z)) \le R(x,z)}
#'
#' \item \eqn{\mathcal{T}}-similarity relation (also called fuzzy \eqn{\mathcal{T}}-equivalence relation): this relation is a fuzzy tolerance relation that is also \eqn{\mathcal{T}}-transitive.
#'
#'       \eqn{\mathcal{T}(R(x,y), R(y,z)) \le R(x,z)}, for a given triangular norm \eqn{\mathcal{T}}. 
#' }
#'
#' The following equations are the tolerance relations on a quantitative attribute \eqn{a}, \eqn{R_a}, proposed by (Jensen and Shen, 2009).
#' \itemize{
#' \item \code{eq.1}: \eqn{R_a(x,y) = 1 - \frac{|a(x) - a(y)|}{|a_{max} - a_{min}|}}
#' \item \code{eq.2}: \eqn{R_a(x,y) = exp(-\frac{(a(x) - a(y))^2}{2 \sigma_a^2})}
#' \item \code{eq.3}: \eqn{R_a(x,y) = max(min(\frac{a(y) - a(x) + \sigma_a}{\sigma_a}, \frac{a(x) - a(y) + \sigma_a}{\sigma_a}), 0)}
#' }
#' where \eqn{\sigma_{a}^2} is the variance of feature \eqn{a} and \eqn{a_{min}} and \eqn{a_{max}} are the minimal and maximal values of data supplied by user. 
#' Additionally, other relations have been implemented in \code{\link{BC.IND.relation.FRST}}
#'
#' For a qualitative (i.e., nominal) attribute \eqn{a}, the classical manner of discerning objects is used, i.e., \eqn{R_a(x,y) = 1}
#' if \eqn{a(x) = a(y)} and \eqn{R_a(x,y) = 0}, otherwise. We can then define, for any subset \eqn{B} of \eqn{A}, the fuzzy \eqn{B}-indiscernibility relation by
#'
#' \eqn{R_B(x,y) = \mathcal{T}(R_a(x,y))},
#'
#' where \eqn{\mathcal{T}} is a t-norm operator, for instance minimum, product and Lukasiewicz t-norm. 
#' In general, \eqn{\mathcal{T}} can be replaced by any aggregation operator, like e.g., the average.
#'
#' In the context of FRST, according to (Radzikowska and Kerre, 2002) lower and upper approximation 
#' are generalized by means of an implicator \eqn{\mathcal{I}} and a t-norm \eqn{\mathcal{T}}. 
#' The following are the fuzzy \eqn{B}-lower and \eqn{B}-upper approximations of a fuzzy set \eqn{A} in \eqn{U}
#' 
#' \eqn{(R_B \downarrow A)(y) = inf_{x \in U} \mathcal{I}(R_B(x,y), A(x))}
#'
#' \eqn{(R_B \uparrow A)(y) = sup_{x \in U} \mathcal{T}(R_B(x,y), A(x))}
#' 
#' The underlying meaning is that \eqn{R_B \downarrow A} is the set of elements \emph{necessarily} satisfying
#' the concept (strong membership), while \eqn{R_B \uparrow A} is the set of elements \emph{possibly} belonging
#' to the concept (weak membership). Many other ways to define the approximations can be found in \code{\link{BC.LU.approximation.FRST}}.
#' Mainly, these were designed to deal with noise in the data and to make the approximations more robust.
#' 
#' Based on fuzzy \eqn{B}-indiscernibility relations, we define the fuzzy \eqn{B}-positive region by, for \eqn{y \in X},
#'
#' \eqn{POS_B(y) = (\cup_{x \in U} R_B \downarrow R_dx)(y)}
#' 
#' We can define the degree of dependency of \eqn{d} on \eqn{B}, \eqn{\gamma_{B}} by
#'
#' \eqn{\gamma_{B} = \frac{|POS_{B}|}{|U|} = \frac{\sum_{x \in U} POS_{B}(x)}{|U|}}
#'
#' A decision reduct is a set \eqn{B \subseteq A} such that \eqn{\gamma_{B} = \gamma_{A}} and \eqn{\gamma_{B'} = \gamma_{B}} for every \eqn{B' \subset B}. 
#'
#' As we know from rough set concepts (See \code{\link{A.Introduction-RoughSets}}), we are able to calculate the 
#' decision reducts by constructing the decision-relative discernibility matrix. Based on (Tsang et al, 2008), the discernibility matrix
#' can be defined as follows. 
#' The discernibility matrix is an \eqn{n \times n} matrix \eqn{(c_{ij})} where 
#' for \eqn{i,j = 1, \ldots, n}
#'        
#' 1) \eqn{c_{ij}= \{a \in A : 1 - R_{a}(x_i, x_j) \ge	\lambda_i\}} if \eqn{\lambda_j < \lambda_i}.
#'
#' 2) \eqn{c_{ij}={\oslash}}, otherwise. 
#'
#' with \eqn{\lambda_i = (R_A \downarrow R_{d}x_{i})(x_i)} and \eqn{\lambda_j = (R_A \downarrow R_{d}x_{j})(x_{j})}
#'
#' Other approaches of discernibility matrix can be read at \code{\link{BC.discernibility.mat.FRST}}.
#'
#' The other implementations of the FRST concepts can be seen at \code{\link{BC.IND.relation.FRST}}, 
#'
#' \code{\link{BC.LU.approximation.FRST}}, and \code{\link{BC.positive.reg.FRST}}. 
#' 
#' @name B.Introduction-FuzzyRoughSets
#' @aliases FuzzyRoughSets-intro
#' @docType package
#' @title Introduction to Fuzzy Rough Set Theory
#' @references 
#' A. M. Radzikowska and E. E. Kerre, "A Comparative Study of Fuzzy Rough Sets", 
#' Fuzzy Sets and Systems, vol. 126, p. 137 - 156 (2002). 
#'
#' D. Dubois and H. Prade, "Rough Fuzzy Sets and Fuzzy Rough Sets",
#' International Journal of General Systems, vol. 17, p. 91 - 209 (1990).
#'
#' E. C. C. Tsang, D. G. Chen, D. S. Yeung, X. Z. Wang, and J. W. T. Lee, 
#' "Attributes Reduction Using Fuzzy Rough Sets", IEEE Trans. Fuzzy Syst., 
#' vol. 16, no. 5, p. 1130 - 1141 (2008).
#'
#' L. A. Zadeh, "Fuzzy Sets",
#' Information and Control, vol. 8, p. 338 - 353 (1965).
#'
#' R. Jensen and Q. Shen,  
#' "New Approaches to Fuzzy-Rough Feature Selection", 
#' IEEE Trans. on Fuzzy Systems, vol. 19, no. 4,
#' p. 824 - 838 (2009).
#'
#' Z. Pawlak, "Rough Sets", International Journal of Computer and Information Sciences, 
#' vol. 11, no. 5, p. 341 - 356 (1982).
NULL



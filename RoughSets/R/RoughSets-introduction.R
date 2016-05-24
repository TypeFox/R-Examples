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
#' This part attempts to introduce rough set theory (RST) and its application to data analysis. 
#' While the classical RST proposed by Pawlak in 1982 is explained in detail in this section, 
#' some recent advancements will be treated in the documentation of the related functions.  
#' 
#' In RST, a data set is represented as a table called an information system \eqn{\mathcal{A} = (U, A)}, where
#' \eqn{U} is a non-empty set of finite objects known as the universe of discourse (note: it refers to all instances/rows 
#' in datasets) and \eqn{A} is a non-empty finite set of attributes, such that \eqn{a : U \to V_{a}} for every \eqn{a \in A}. 
#' The set \eqn{V_{a}} is the set of values that attribute \eqn{a} may take. Information systems that involve a decision attribute, 
#' containing classes for each object, are called decision systems or decision tables. More formally, it is a pair \eqn{\mathcal{A} = (U, A \cup \{d\})},
#' where \eqn{d \notin A} is the decision attribute. The elements of \eqn{A} are called conditional attributes. The information system
#' representing all data in a particular system may contain redundant parts. It could happen because there are the same 
#' or indiscernible objects or some superfluous attributes. The indiscernibility relation is a binary relation showing the relation between two objects.
#' This relation is an equivalence relation.  
#' Let \eqn{\mathcal{A} = (U, A)} be an information system, then for any \eqn{B \subseteq A} there is an equivalence 
#' relation \eqn{R_B(x,y)}:
#' 
#' \eqn{R_B(x,y)= \{(x,y) \in U^2 | \forall a \in B, a(x) = a(y)\}}
#' 
#' If \eqn{(x,y) \in R_B(x,y)}, then \eqn{x} and \eqn{y} are indiscernible by attributes from \eqn{B}. The equivalence
#' classes of the \eqn{B}-indiscernibility relation are denoted \eqn{[x]_{B}}. The indiscernibility relation will be further used to define basic concepts of rough
#' set theory which are lower and upper approximations.
#' 
#' Let \eqn{B \subseteq A} and \eqn{X \subseteq U},
#' \eqn{X} can be approximated using the information contained within \eqn{B} by constructing 
#' the \eqn{B}-lower and \eqn{B}-upper approximations of \eqn{X}:
#' 
#' \eqn{R_B \downarrow X = \{ x \in U | [x]_{B} \subseteq X \}}
#' 
#' \eqn{R_B \uparrow X = \{ x \in U | [x]_{B} \cap X \not= \emptyset \}}
#'
#' The tuple \eqn{\langle R_B \downarrow X, R_B \uparrow X \rangle} is called a rough set.
#' The objects in \eqn{R_B \downarrow X} mean that they can be with certainty classified as members of \eqn{X} on the basis of knowledge in \eqn{B}, while
#' the objects in \eqn{R_B \uparrow X} can be only classified as possible members of \eqn{X} on the basis of knowledge in \eqn{B}. 
#' 
#' In a decision system, for \eqn{X} we use decision concepts (equivalence classes of decision attribute) \eqn{[x]_d}. 
#' We can define \eqn{B}-lower and \eqn{B}-upper approximations as follows.
#'
#' \eqn{R_B \downarrow [x]_d = \{ x \in U | [x]_{B} \subseteq [x]_d \}}
#' 
#' \eqn{R_B \uparrow [x]_d = \{ x \in U | [x]_{B} \cap [x]_d \not= \emptyset \}}
#' 
#' The positive, negative and boundary of \eqn{B} regions can be defined as:
#' 
#' \eqn{POS_{B} = \bigcup_{x \in U } R_B \downarrow [x]_d}
#'
#' The boundary region, \eqn{BND_{B}}, is the set of objects that can possibly, but not certainly, be classified.
#' 
#' \eqn{BND_{B} = \bigcup_{x \in U} R_B \uparrow [x]_d - \bigcup_{x \in U} R_B \downarrow [x]_d}
#' 
#' Furthermore, we can calculate the degree of dependency of the decision on a set of attributes. The decision attribute \eqn{d}
#' depends totally on a set of attributes \eqn{B}, denoted \eqn{B \Rightarrow d},
#' if all  attribute values from \eqn{d} are uniquely determined by values of attributes from \eqn{B}. It can be defined as follows.
#' For \eqn{B \subseteq A}, it is said that \eqn{d} depends on \eqn{B} in a degree of dependency \eqn{\gamma_{B} = \frac{|POS_{B}|}{|U|}}. 
#' 
#' A decision reduct is a set \eqn{B \subseteq A} such that \eqn{\gamma_{B} = \gamma_{A}} and \eqn{\gamma_{B'} < \gamma_{B}} for every \eqn{B' \subset B}.
#' One algorithm to determine all reducts is by constructing the decision-relative discernibility matrix. 
#' The discernibility matrix \eqn{M(\mathcal{A})} is an \eqn{n \times n} matrix \eqn{(c_{ij})} where
#'
#' \eqn{c_{ij} = \{a \in A: a(x_i) \neq a(x_j) \}} if \eqn{d(x_i) \neq d(x_j)} and
#'
#' \eqn{c_{ij} = \oslash} otherwise
#'
#' The discernibility function \eqn{f_{\mathcal{A}}} for a decision system \eqn{\mathcal{A}} is a boolean function of \eqn{m} boolean variables \eqn{\bar{a}_1, \ldots, \bar{a}_m}
#' corresponding to the attributes \eqn{a_1, \ldots, a_m} respectively, and defined by
#'
#' \eqn{f_{\mathcal{A}}(\bar{a_1}, \ldots, \bar{a_m}) = \wedge \{\vee \bar{c}_{ij}: 1 \le j < i \le n, c_{ij} \neq \oslash \}}
#'
#' where \eqn{\bar{c}_{ij}= \{ \bar{a}: a \in c_{ij}\}}. The decision reducts of \eqn{A} are then the prime implicants of the function \eqn{f_{\mathcal{A}}}. 
#' The complete explanation of the algorithm can be seen in (Skowron and Rauszer, 1992).
#' 
#' The implementations of the RST concepts can be seen in \code{\link{BC.IND.relation.RST}}, 
#'
#' \code{\link{BC.LU.approximation.RST}}, \code{\link{BC.positive.reg.RST}}, and 
#'
#' \code{\link{BC.discernibility.mat.RST}}. 
#' 
#' @name A.Introduction-RoughSets
#' @aliases RoughSets-intro
#' @docType package
#' @title Introduction to Rough Set Theory
#' @references 
#' A. Skowron and C. Rauszer,  
#' "The Discernibility Matrices and Functions in Information Systems", 
#' in: R. Slowinski (Ed.), Intelligent Decision Support: Handbook of Applications and
#' Advances of Rough Sets Theory, Kluwer Academic Publishers, Dordrecht, Netherland,  
#' p. 331 - 362 (1992).
#'
#' Z. Pawlak, "Rough Sets", 
#' International Journal of Computer and Information System, 
#' vol. 11, no.5, p. 341 - 356 (1982).
#' 
#' Z. Pawlak, "Rough Sets: Theoretical Aspects of Reasoning about Data, System Theory, Knowledge Engineering and Problem Solving",
#' vol. 9, Kluwer Academic Publishers, Dordrecht, Netherlands (1991). 
#'
NULL



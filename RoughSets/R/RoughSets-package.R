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
#' This part contains global explanations about the implementation and use of the \code{RoughSets} package.
#' The package \code{RoughSets} attempts to provide a complete tool to model and analyze
#' information systems based on rough set theory (RST) and fuzzy rough set theory (FRST).
#' From fundamental point of view, this package allows to construct rough sets by defining lower and upper approximations.
#' Furthermore, recent methods for tackling common tasks in data mining, such as data preprocessing (e.g., discretization, feature selection, missing value completion,
#' and instance selection), rule induction, and prediction classes or decision values of new datasets
#' are available as well.
#'
#' There are two main parts considered in this package which are RST and FRST.
#' RST was introduced by (Pawlak, 1982; Pawlak, 1991) which provides sophisticated mathematical tools to model
#' and analyze information systems that involve uncertainty and imprecision. By employing indiscernibility relation among objects, RST does not require
#' additional parameters to extract information.
#' The detailed explanation about fundamental concepts of RST can be read in Section \code{\link{A.Introduction-RoughSets}}. Secondly, FRST, an extension of
#' RST, was introduced by (Dubois and Prade, 1990) as a combination between
#' fuzzy sets proposed by (Zadeh, 1965) and RST. This concept allows to analyze continuous attributes without
#' performing discretization on data first. Basic concepts of FRST
#' can be seen in \code{\link{B.Introduction-FuzzyRoughSets}}.
#'
#' Based on the above concepts, many methods have been proposed and applied for dealing with several different domains.
#' In order to solve the problems, methods employ the indiscernibility relation and lower and upper approximation concepts.
#' All methods that have been implemented in this package will be explained by grouping based on their domains. The following is
#' a list of domains considered in this package:
#' \itemize{
#' \item Basic concepts: This part, we can divide into four different tasks which are
#'       indiscernibility relation, lower and upper approximations, positive region and discernibility matrix.
#'       All of those tasks have been explained briefly in Section \code{\link{A.Introduction-RoughSets}} and
#'
#'      \code{\link{B.Introduction-FuzzyRoughSets}}.
#' \item Discretization: It is used to convert real-valued attributes into nominal/symbolic ones in an information system.
#'       In RST point of view, this task attempts to maintain the discernibility between objects.
#' \item Feature selection: It is a process for finding a subset of features which have the same quality as the complete feature set.
#'      In other words, its purpose is to select the significant features and eliminate the dispensible ones.
#'      It is a useful and necessary process when we are facing datasets containing large numbers of features. From RST and FRST perspective,
#'      feature selection refers to searching superreducts and reducts. The detailed information about reducts can be read in
#'      \code{\link{A.Introduction-RoughSets}} and \code{\link{B.Introduction-FuzzyRoughSets}}.
#' \item Instance selection: This process is aimed to remove noisy, superfluous, or inconsistent instances from training datasets but retain consistent ones.
#'      In other words, good accuracy of classification is achieved by removing instances which do not give positive contributions.
#' \item Prediction/classification: This task is used to predict decision values of a new dataset (test data).
#'      We consider implementing some methods to perform this task, such as fuzzy-rough nearest neighbor approaches, etc.
#' \item Rule induction: This task refers to generate IF - THEN rules. The rule represents knowledge which is contained in a dataset.
#'      One advantage of building rules is that naturally the model is easy to interpret. Then, predicted values over new datasets can be determined by
#'      considering the rules.
#' }
#'
#' As we mentioned before, we have embedded many well-known algorithms or techniques for handling the above domains. The algorithms were considered
#' since experimentally it has been proven that they were able to tackle complex tasks. They are implemented as functions that were organized
#' to work with the same data structures. So, users can perform various approaches for a particular task easily and then compare their results.
#' In order to be recognized quickly, generally we have chosen the names of the functions with some conventions. The names contain three parts
#' which are \code{prefix}, \code{suffix}, and \code{middle} that are separated by a point. The following is a description of each
#' part.
#' \itemize{
#' \item \code{prefix}: There are some different prefixes for names of functions expressing a kind of task to be performed.
#'                      The function names with prefix \code{BC} refer to \emph{basic concepts} which means that the functions are created for
#'                      implementing the basic concepts of RST and FRST.
#'                      While prefix \code{D} refers to \emph{discretization}, \code{FS}, \code{IS}, \code{RI}, \code{MV}, and \code{C} refer to \emph{feature selection},
#'                      \emph{instance selection}, \emph{rule induction}, \emph{missing value completion}, and \emph{classifier based on nearest neighbor} domains. Furthermore, \code{SF} and \code{X} mean that
#'                      functions are used as \emph{supporting functions} which are not related directly with RST and FRST and \emph{auxiliary} functions which are called as a parameter.
#' \item \code{suffix}: It is located at the end of names. There are two types available: \code{RST} and \code{FRST}. \code{RST} represents \emph{rough set theory}
#'                      while \code{FRST} shows that the function is applied to \emph{fuzzy rough set theory}. Additionally, some functions that do not have
#'                      \code{RST} or \code{FRST} suffix are used for both theories.
#' \item \code{middle}: All other words in the middle of the names are used to express the actual name of a particular method/algorithm or functionality.
#'                      In this case, it could consist of more than one word separated by points.
#' }
#' For instance, the function \code{\link{BC.IND.relation.RST}} is used to calculate the indiscernibility relation which is one of the basic concepts of RST.
#' Other functions that have names not based on the above rules are S3 functions e.g., \code{summary} and \code{predict} which are
#' used to summarize objects and predict new data, respectively.
#'
#' The following description explains domains and their algorithms implemented in the package:
#' \enumerate{
#' \item \bold{The implementations of RST}: This part outlines some considered algorihtms/methods based on RST.
#'              The approaches can be classified based on their tasks as follows:
#'             \enumerate{
#'             \item The basic concepts: The following is a list showing tasks and their implementations as functions.
#'                       \itemize{
#'                       \item Indiscernibility relation: It is a relation determining whether two objects are indiscernible by some attributes.
#'                             It is implemented in \code{\link{BC.IND.relation.RST}}.
#'                       \item Lower and upper approximations: These approximations show whether objects can be classified with certainty or not.
#'                             It is implemented in \code{\link{BC.LU.approximation.RST}}.
#'                       \item Positive region: It is used to determine objects that are included in positive region and the degree of dependency.
#'                             It is implemented in \code{\link{BC.positive.reg.RST}}.
#'                       \item Discernibility matrix: It is used to create a discernibility matrix showing attributes that discern each pair of objects.
#'                             It is implemented in \code{\link{BC.discernibility.mat.RST}}.
#'                       }
#'             \item Discretization: There are a few methods included in the package:
#'                      \itemize{
#'                       \item \code{\link{D.global.discernibility.heuristic.RST}}: It implements the global discernibility algorithm
#'                              which is computing globally semi-optimal cuts using the maximum discernibility heuristic.
#'						           \item \code{\link{D.discretize.quantiles.RST}}: It is a function used for computing cuts of the "quantile-based" discretization into \eqn{n} intervals.
#'                       \item \code{\link{D.discretize.equal.intervals.RST}}: It is a function used for computing cuts of the "equal interval size" discretization into \eqn{n} intervals.
#'                       }
#'      		The output of these functions is a list of cut values which are the values for converting real to nominal values.
#'              So, in order to generate a new decision table according to the cut values, we need to call \code{\link{SF.applyDecTable}}.
#'              Additionally, we have implemented \code{\link{D.discretization.RST}} as a wrapper function collecting all methods considered to perform discretization tasks.
#'             \item Feature selection: According to its output, it can be classified into the following groups:
#'                      \itemize{
#'                       \item Feature subset: It refers to a superreduct which is not necessarily minimal. In other words, the methods in this group
#'                              might generate just a subset of attributes.
#'                             \itemize{
#'                                  \item QuickReduct algorithm: It has been implemented in \code{\link{FS.quickreduct.RST}}.
#'                                  \item Superreduct generation: It is based on some criteria:
#'                                                    entropy, gini index, discernibility measure, size of positive region.
#'
#'                                                    It is implemented in \code{\link{FS.greedy.heuristic.superreduct.RST}}.
#'                             }
#'                             Furthermore, we provide a wrapper function \code{\link{FS.feature.subset.computation}} in order to give a user interface for many methods of RST and FRST that are included in this group.
#'                       \item Reduct: The following are methods that produce a single decision reduct:
#'                             \itemize{
#'                                  \item Reduct generation based on criteria: It is based on different criteria which are
#'                                                    entropy, gini index, discernibility measure, size of positive region.
#'                                                    It has been implemented in \code{\link{FS.greedy.heuristic.reduct.RST}}.
#'                                  \item Permutation reduct: It is based on a permutation schema over all attributes.
#'                                                     It has been implemented in \code{\link{FS.permutation.heuristic.reduct.RST}}.
#'                              }
#'                             Furthermore, we provide a wrapper function \code{\link{FS.reduct.computation}} in order to give a user interface toward many methods of RST and FRST that are included in this group.
#'                       \item All reducts: In order to generate all reducts, we execute \code{\link{FS.all.reducts.computation}}. However,
#'                             before doing that, we need to call \code{\link{BC.discernibility.mat.RST}} for
#'                             constructing a decision-relative discernibility matrix
#'                       }
#'                       It should be noted that the outputs of the functions are decision reducts. So, for generating a new decision table according to the decision reduct,
#'                       we need to call \code{\link{SF.applyDecTable}}.
#'            \item Rule induction: We provide several functions used to generate rules, as follows:
#'                       \itemize{
#'                           \item The function \code{\link{RI.indiscernibilityBasedRules.RST}}: This function requires the output of the feature selection functions.
#'                           \item The function \code{\link{RI.CN2Rules.RST}}: It is a rule induction method based on the CN2 algorithm.
#'                           \item The function \code{\link{RI.LEM2Rules.RST}}: It implements a rule induction method based on the LEM2 algorithm.
#'                           \item The function \code{\link{RI.AQRules.RST}}: It is a rule induction based on the AQ-style algorithm.
#'                       }
#'                        After obtaining the rules, we execute \code{\link{predict.RuleSetRST}} considering our rules and given newdata/testing data to obtain predicted values/classes.
#'              }
#' \item \bold{The implementations of FRST}: As in the \code{RST} part, this part contains several algorithms that can be classified into several groups based on their purpose.
#'           The following is a description of all methods that have been implemented in functions:
#'  \enumerate{
#'	  \item Basic concepts: The following is a list showing tasks and their implementations:
#'    \itemize{
#'       \item Indiscernibility relations: they are fuzzy relations determining to which degree two objects are similar.
#'             This package provides several types of relations which are implemented in a single function
#'             called \code{\link{BC.IND.relation.FRST}}. We consider several types of relations e.g.,
#'             fuzzy equivalence, tolerance, and \eqn{T}-similarity relations. These relations can be chosen by
#'             assigning \code{type.relation}. Additionally, in this function, we provide several options to
#'             calculate aggregation e.g., triangular norm operators (e.g., \code{"lukasiewicz"}, \code{"min"}, etc)
#'             and user-defined operators.
#' 	   	 \item Lower and upper approximations: These approximations show to what extent objects can be classified with certainty or not.
#'           This task has been implemented in
#'
#'           \code{\link{BC.LU.approximation.FRST}}. There are many approaches available in this package that can be selected by assigning the parameter \code{type.LU}.
#'           The considered methods are
#'           implication/t-norm, \eqn{\beta}-precision fuzzy rough sets (\eqn{\beta}-PFRS), vaguely quantified rough sets (VQRS), fuzzy variable precision rough sets (FVPRS), ordered weighted average (OWA),
#'           soft fuzzy rough sets (SFRS), and robust fuzzy rough sets (RFRS). Furthermore, we provide a facility, which is \code{"custom"}, where users can create their own approximations by
#'           defining functions to calculate lower and upper approximations. Many options to calculate implicator and triangular norm are also available.
#'      \item Positive region: It is used to determine the membership degree of each object to the positive region and the degree of dependency.
#'                   It is implemented in \code{\link{BC.positive.reg.FRST}}.
#'      \item Discernibility matrix: It is used to construct the decision-relative discernibility matrix. There are some approaches to construct the matrix,
#'                   e.g., based on standard approach, Gaussian reduction, alpha reduction, and minimal element in discernibility matrix. They have been implemented
#'                   in \code{\link{BC.discernibility.mat.FRST}}.
#'    }
#'    \item Feature selection: According to the output of functions,
#'       we may divide them into three groups: those that produce a superreduct, a set of reducts, or a single reduct. The following is a description of functions based on their types:
#'       \itemize{
#'           \item Feature subset: It refers to methods which produce a superreduct which is not necessarily a reduct. In other words methods in this group
#'                               might generate just a subset of attributes.
#'                 The following is a complete list of methods considered in this package:
#'                 \itemize{
#'                  \item positive region based algorithms: It refers to
#'                        positive regions, as a way to evaluate attributes to be selected. They are implemented in \code{\link{FS.quickreduct.FRST}}.
#'                        Furthermore, we provide several different measures based on the positive region in this function.
#'                        All methods included in this part employ the QuickReduct algorithm to obtain selected features.
#'                        In order to choose a particular algorithm, we need to assign parameter \code{type.method} in \code{\link{FS.quickreduct.FRST}}.
#'                \item boundary region based algorithm: This algorithm is based on the membership degree to the fuzzy boundary region.
#'                      This algorithm has been implemented in \code{\link{FS.quickreduct.FRST}}.
#'               }
#'               Furthermore, we provide a wrapper function \code{\link{FS.feature.subset.computation}} in order to give a user interface for many methods of RST and FRST.
#'           \item Reduct: It refers to a method that produces a single decision reduct. We provide one algorithm which is the near-optimal reduction proposed by Zhao et al.
#'                       It is implemented in \code{\link{FS.nearOpt.fvprs.FRST}}.
#'                Furthermore, we provide a wrapper function \code{\link{FS.reduct.computation}} in order to provide a user interface toward many methods of RST and FRST.
#'           \item All reducts:  In order to get all decision reducts, we execute \code{\link{FS.all.reducts.computation}}. However, before doing that, we firstly execute the \code{\link{BC.discernibility.mat.FRST}} function for
#'                             constructing a decision-relative discernibility matrix.
#'     }
#'     The output of the above methods is a class containing a decision reduct/feature subset and other descriptions.
#'     For generating a new decision table according to the decision reduct, we provide the function \code{\link{SF.applyDecTable}}.
#'
#'     \item Rule induction: It is a task used to generate
#'           rules representing knowledge of a decision table. Commonly, this process is called learning phase in machine learning.
#'           The following methods are considered to generate rules:
#'           \itemize{
#' 	    		 \item \code{\link{RI.hybridFS.FRST}}: It combines fuzzy-rough rule induction
#'                      and feature selection.
#'               \item \code{\link{RI.GFRS.FRST}}: It refers to rule induction based on generalized fuzzy rough sets (GFRS).
#'           }
#'           After generating rules, we can use them to predict decision values/classes of new data
#'           by executing the S3 function \code{\link{predict.RuleSetFRST}}.
#'     \item Instance selection: The following functions select instances to improve accuracy by
#'           removing noisy, superfluous or inconsistent ones from training datasets.
#'           \itemize{
#'             \item \code{\link{IS.FRIS.FRST}}: It refers to the fuzzy rough instance selection (FRIS). It evaluates the degree of membership to the positive region of each instance.
#'             		  If an instance's membership degree is less than the threshold, then the instance can be removed.
#'             \item \code{\link{IS.FRPS.FRST}}: It refers to the fuzzy-rough prototype selection (FRPS). It employs prototype selection (PS) to improve the accuracy of the $k$-nearest neighbor (kNN) method.
#'           }
#'       We provide the function \code{\link{SF.applyDecTable}} that is used to generate a new decision table according to the output of instance selection functions.
#'    \item Fuzzy-rough nearest neighbors: This part provides methods based on nearest neighbors for
#'          predicting decision values/classes of new datasets. In other words, by supplying a decision table as training data
#'          we can predict decision values of new data at the same time.
#'          We have considered the following methods:
#'          \itemize{
#'            \item \code{\link{C.FRNN.FRST}}: It refers to the fuzzy-rough nearest neighbors based on Jensen and Cornelis' technique.
#'            \item \code{\link{C.FRNN.O.FRST}}: It refers to the fuzzy-rough ownership nearest neighbor algorithm based on Sarkar's method.
#'            \item \code{\link{C.POSNN.FRST}}: The positive region based fuzzy-rough nearest neighbor algorithm based on Verbiest et al's technique.
#'          }
#' }
#' }
#'
#' Furthermore, we provide an additional feature which is missing value completion. Even though algorithms, included in this feature, are not based on RST and FRST, they will be usefull to do data analysis.
#' The following is a list of functions implemented for handling missing values in the data preprocessing step:
#' \itemize{
#' \item \code{\link{MV.deletionCases}}: it refers to the approach deleting instances.
#' \item \code{\link{MV.mostCommonValResConcept}}: it refers to the approach based on the most common value or mean of an attribute restricted to a concept.
#' \item \code{\link{MV.mostCommonVal}}: it refers to the approach replacing missing attribute values by the attribute mean or common values.
#' \item \code{\link{MV.globalClosestFit}}: it refers to the approach based on the global closest fit approach.
#' \item \code{\link{MV.conceptClosestFit}}: it refers to the approach based on the concept closest fit approach.
#' }
#' Additionally, we provide a wrapper function which is \code{\link{MV.missingValueCompletion}}
#' in order to give a user interface for the methods.
#'
#' To get started with the package, the user can have a look at the examples included in
#' the documentation on each function. Additionally, to show general usage of the package briefly,
#' we also provide some examples showing general usage in this section.
#'
#' If you have problems using the package, find a bug, or have suggestions,
#' please contact the package maintainer by email, instead of writing to the general R lists
#' or to other internet forums and mailing lists.
#'
#' There are many demos that ship with the package. To get a list of them, type:
#'
#' \code{demo()}
#'
#' Then, to start a demo, type \code{demo(<demo_name_here>)}. All the demos are presented as
#' R scripts in the package sources in the "demo" subdirectory.
#'
#' Currently, there are the following demos available:
#'
#' \itemize{
#' \item Basic concepts of RST and FRST:
#'
#' \code{demo(BasicConcept.RST)},
#' \code{demo(BasicConcept.FRST)},
#'
#' \code{demo(DiscernibilityMatrix.RST)},
#' \code{demo(DiscernibilityMatrix.FRST)}.
#'
#' \item Discretization based on RST:
#'
#' \code{demo(D.local.discernibility.matrix.RST)},
#' \code{demo(D.max.discernibility.matrix.RST)},
#'
#' \code{demo(D.global.discernibility.heuristic.RST)},
#' \code{demo(D.discretize.quantiles.RST)},
#'
#' \code{demo(D.discretize.equal.intervals.RST)}.
#'
#' \item Feature selection based on RST:
#'
#' \code{demo(FS.permutation.heuristic.reduct.RST)},
#' \code{demo(FS.quickreduct.RST)},
#'
#' \code{demo(FS.greedy.heuristic.reduct.RST)},
#' \code{demo(FS.greedy.heuristic.reduct.RST)}.
#'
#' \item Feature selection based on FRST:
#'
#' \code{demo(FS.QuickReduct.FRST.Ex1)},
#' \code{demo(FS.QuickReduct.FRST.Ex2)},
#'
#' \code{demo(FS.QuickReduct.FRST.Ex3)},
#' \code{demo(FS.QuickReduct.FRST.Ex4)},
#'
#' \code{demo(FS.QuickReduct.FRST.Ex5)},
#' \code{demo(FS.nearOpt.fvprs.FRST)}.
#'
#' \item Instance selection based on FRST:
#'
#' \code{demo(IS.FRIS.FRST)},
#' \code{demo(IS.FRPS.FRST)}
#'
#' \item Classification using the Iris dataset:
#'
#' \code{demo(FRNN.O.iris)},
#' \code{demo(POSNN.iris)},
#' \code{demo(FRNN.iris)}.
#'
#' \item Rule induction based on RST:
#'
#' \code{demo(RI.indiscernibilityBasedRules.RST)}.
#'
#' \item Rule induction based on FRST:
#'
#' \code{demo(RI.classification.FRST)},
#' \code{demo(RI.regression.FRST)}.
#'
#' \item Missing value completion:
#' \code{demo(MV.simpleData)}.
#' }
#'
#' Some decision tables have been embedded in this package which can be seen in
#' \code{\link{RoughSetData}}.
#'
#' Finally, you may visit the package webpage \url{http://sci2s.ugr.es/dicits/software/RoughSets},
#' where we provide a more extensive introduction as well as additional explanations of
#' the procedures.
#'
#' @name RoughSets-package
#' @aliases RoughSets
#' @docType package
#' @title Getting started with the RoughSets package
#' @references
#' D. Dubois and H. Prade, "Rough Fuzzy Sets and Fuzzy Rough Sets",
#' International Journal of General Systems, vol. 17, p. 91 - 209 (1990).
#'
#' L.A. Zadeh, "Fuzzy Sets",
#' Information and Control, vol. 8, p. 338 - 353 (1965).
#'
#' Z. Pawlak, "Rough Sets",
#' International Journal of Computer and Information System,
#' vol. 11, no. 5, p. 341 - 356 (1982).
#'
#' Z. Pawlak, "Rough Sets: Theoretical Aspects of Reasoning About Data, System Theory, Knowledge Engineering and Problem Solving",
#' vol. 9, Kluwer Academic Publishers, Dordrecht, Netherlands (1991).
#'
# @keywords rough sets fuzzy rough sets instance selection feature selection classification prediction
#' @author Lala Septem Riza \email{lala.s.riza@@decsai.ugr.es},
#'
#' Andrzej Janusz \email{andrzejanusz@@gmail.com},
#'
#' Chris Cornelis \email{chriscornelis@@decsai.ugr.es},
#'
#' Francisco Herrera \email{herrera@@decsai.ugr.es},
#'
#' Dominik Slezak \email{slezak@@mimuw.edu.pl},
#'
#' and Jose Manuel Benitez \email{j.m.benitez@@decsai.ugr.es}
#'
#' DiCITS Lab, SCI2S group, CITIC-UGR, DECSAI, University of Granada,
#'
#' \url{http://dicits.ugr.es}, \url{http://sci2s.ugr.es}
#'
#' Institute of Mathematics, University of Warsaw.
#'
#' @useDynLib RoughSets
#' @examples
#' ##############################################################
#' ## A.1 Example: Basic concepts of rough set theory
#' ##############################################################
#' ## Using hiring data set, see RoughSetData
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt
#'
#' ## define considered attributes which are first, second, and
#' ## third attributes
#' attr.P <- c(1,2,3)
#'
#' ## compute indiscernibility relation
#' IND <- BC.IND.relation.RST(decision.table, feature.set = attr.P)
#'
#' ## compute lower and upper approximations
#' roughset <- BC.LU.approximation.RST(decision.table, IND)
#'
#' ## Determine regions
#' region.RST <- BC.positive.reg.RST(decision.table, roughset)
#'
#' ## The decision-relative discernibility matrix and reduct
#' disc.mat <- BC.discernibility.mat.RST(decision.table, range.object = NULL)
#'
#' ##############################################################
#' ## A.2 Example: Basic concepts of fuzzy rough set theory
#' ##############################################################
#' ## Using pima7 data set, see RoughSetData
#' data(RoughSetData)
#' decision.table <- RoughSetData$pima7.dt
#'
#' ## In this case, let us consider the first and second attributes
#' conditional.attr <- c(1, 2)
#'
#' ## We are using the "lukasiewicz" t-norm and the "tolerance" relation
#' ## with "eq.1" as fuzzy similarity equation
#' control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"),
#'                     type.relation = c("tolerance", "eq.1"))
#'
#' ## Compute fuzzy indiscernibility relation
#' IND.condAttr <- BC.IND.relation.FRST(decision.table, attributes = conditional.attr,
#'                             control = control.ind)
#'
#' ## Compute fuzzy lower and upper approximation using type.LU : "implicator.tnorm"
#' ## Define index of decision attribute
#' decision.attr = c(9)
#'
#' ## Compute fuzzy indiscernibility relation of decision attribute
#' ## We are using "crisp" for type of aggregation and type of relation
#' control.dec <- list(type.aggregation = c("crisp"), type.relation = "crisp")
#'
#' IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = decision.attr,
#'                             control = control.dec)
#'
#' ## Define control parameter containing type of implicator and t-norm
#' control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz")
#'
#' ## Compute fuzzy lower and upper approximation
#' FRST.LU <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr,
#'               type.LU = "implicator.tnorm", control = control)
#'
#' ## Determine fuzzy positive region and its degree of dependency
#' fuzzy.region <- BC.positive.reg.FRST(decision.table, FRST.LU)
#'
#' ###############################################################
#' ## B Example : Data analysis based on RST and FRST
#' ## In this example, we are using wine dataset for both RST and FRST
#' ###############################################################
#' ## Load the data
#' \dontrun{data(RoughSetData)
#' dataset <- RoughSetData$wine.dt
#'
#' ## Shuffle the data with set.seed
#' set.seed(5)
#' dt.Shuffled <- dataset[sample(nrow(dataset)),]
#'
#' ## Split the data into training and testing
#' idx <- round(0.8 * nrow(dt.Shuffled))
#'   wine.tra <-SF.asDecisionTable(dt.Shuffled[1:idx,],
#' decision.attr = 14, indx.nominal = 14)
#'   wine.tst <- SF.asDecisionTable(dt.Shuffled[
#'  (idx+1):nrow(dt.Shuffled), -ncol(dt.Shuffled)])
#'
#' ## DISCRETIZATION
#' cut.values <- D.discretization.RST(wine.tra,
#' type.method = "global.discernibility")
#' d.tra <- SF.applyDecTable(wine.tra, cut.values)
#' d.tst <- SF.applyDecTable(wine.tst, cut.values)
#'
#' ## FEATURE SELECTION
#' red.rst <- FS.feature.subset.computation(d.tra,
#'   method="quickreduct.rst")
#' fs.tra <- SF.applyDecTable(d.tra, red.rst)
#'
#' ## RULE INDUCTION
#' rules <- RI.indiscernibilityBasedRules.RST(d.tra,
#'   red.rst)
#'
#' ## predicting newdata
#' pred.vals <- predict(rules, d.tst)
#'
#' #################################################
#' ## Examples: Data analysis using the wine dataset
#' ## 2. Learning and prediction using FRST
#' #################################################
#'
#' ## FEATURE SELECTION
#' reduct <- FS.feature.subset.computation(wine.tra,
#'  method = "quickreduct.frst")
#'
#' ## generate new decision tables
#' wine.tra.fs <- SF.applyDecTable(wine.tra, reduct)
#' wine.tst.fs <- SF.applyDecTable(wine.tst, reduct)
#'
#' ## INSTANCE SELECTION
#' indx <- IS.FRIS.FRST(wine.tra.fs,
#'  control = list(threshold.tau = 0.2, alpha = 1))
#'
#' ## generate a new decision table
#' wine.tra.is <- SF.applyDecTable(wine.tra.fs, indx)
#'
#' ## RULE INDUCTION (Rule-based classifiers)
#' control.ri <- list(
#'  type.aggregation = c("t.tnorm", "lukasiewicz"),
#'  type.relation = c("tolerance", "eq.3"),
#'  t.implicator = "kleene_dienes")
#'
#' decRules.hybrid <- RI.hybridFS.FRST(wine.tra.is,
#'   control.ri)
#'
#' ## predicting newdata
#' predValues.hybrid <- predict(decRules.hybrid,
#'   wine.tst.fs)
#'
#' #################################################
#' ## Examples: Data analysis using the wine dataset
#' ## 3. Prediction using fuzzy nearest neighbor classifiers
#' #################################################
#'
#' ## using FRNN.O
#' control.frnn.o <- list(m = 2,
#'   type.membership = "gradual")
#'
#' predValues.frnn.o <- C.FRNN.O.FRST(wine.tra.is,
#'   newdata = wine.tst.fs, control = control.frnn.o)
#'
#' ## Using FRNN
#' control.frnn <- list(type.LU = "implicator.tnorm",k=20,
#'   type.aggregation = c("t.tnorm", "lukasiewicz"),
#'   type.relation = c("tolerance", "eq.1"),
#'   t.implicator = "lukasiewicz")
#'
#' predValues.frnn <- C.FRNN.FRST(wine.tra.is,
#'   newdata = wine.tst.fs, control = control.frnn)
#'
#' ## calculating error
#' real.val <- dt.Shuffled[(idx+1):nrow(dt.Shuffled),
#'   ncol(dt.Shuffled), drop = FALSE]
#'
#' err.1 <- 100*sum(pred.vals!=real.val)/nrow(pred.vals)
#' err.2 <- 100*sum(predValues.hybrid!=real.val)/
#'   nrow(predValues.hybrid)
#' err.3 <- 100*sum(predValues.frnn.o!=real.val)/
#'   nrow(predValues.frnn.o)
#' err.4 <- 100*sum(predValues.frnn!=real.val)/
#'   nrow(predValues.frnn)
#'
#' cat("The percentage error = ", err.1, "\n")
#' cat("The percentage error = ", err.2, "\n")
#' cat("The percentage error = ", err.3, "\n")
#' cat("The percentage error = ", err.4, "\n")}
NULL

#' @import Rcpp
NULL


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
#' This is one of the central functions of the package. This function is used to 
#' generate/learn the model from numerical data using fuzzy rule-based systems.
#' 
#' This function makes accessible all learning methods that are implemented 
#' in this package. All of the methods use this function as interface for the learning 
#' stage, so users do not need to call other functions in the learning phase. 
#' In order to obtain good results, users need to adjust some parameters such as the 
#' number of labels, the type of the shape of the membership function, the maximal number of iterations, 
#' the step size of the gradient descent, or other method-dependent parameters which are collected in the \code{control}
#' parameter. After creating the model using this function, it can be used to predict new data with \code{\link{predict}}.
#'
#' @title The frbs model building function
#'
#' @param data.train a data frame or matrix (\eqn{m \times n}) of data for the training process, 
#'        where \eqn{m} is the number of instances and 
#'        \eqn{n} is the number of variables; the last column is the output variable. It should be noted that
#'        the training data must be expressed in numbers (numerical data). And,
#'        especially for classification tasks, the last column representing class names/symbols isn't allowed 
#'        to have values 0 (zero). In the other words, the categorical values 0 should be replaced with other values.
#' @param range.data a matrix (\eqn{2 \times n}) containing the range of the data, where \eqn{n} is the number of variables, and
#'        first and second rows are the minimum and maximum values, respectively. It should be noted that
#'        for \code{"FRBCS.W"}, \code{"FRBCS.CHI"}, \code{"GFS.GCCL"}, \code{"FH.GBML"}, and \code{"SLAVE"}, \eqn{n} represents the number of input variables only
#'        (without the output variable). It will be assigned as min/max of training data if it is omitted. 
#' @param method.type this parameter determines the learning algorithm to be used. 
#'        The following methods are implemented: 
#' \itemize{
#' \item \code{"WM"}: Wang and Mendel's technique to handle regression tasks. See \code{\link{WM}};
#' \item \code{"SBC"}: subtractive clustering method to handle regression tasks. See \code{\link{SBC}};
#' \item \code{"HYFIS"}: hybrid neural fuzzy inference systems to handle regression tasks. See \code{\link{HyFIS}};
#' \item \code{"ANFIS"}: adaptive neuro-fuzzy inference systems to handle regression tasks. See \code{\link{ANFIS}};
#' \item \code{"FRBCS.W"}: fuzzy rule-based classification systems with weight factor based on Ishibuchi's method 
#'                  to handle classification tasks. See \code{\link{FRBCS.W}}; 
#' \item \code{"FRBCS.CHI"}: fuzzy rule-based classification systems based on Chi's method to handle
#'                  classification tasks. See \code{\link{FRBCS.CHI}}; 
#' \item \code{"DENFIS"}: dynamic evolving neuro-fuzzy inference systems to handle regression tasks. See \code{\link{DENFIS}};
#' \item \code{"FS.HGD"}: fuzzy system using heuristic and gradient descent method to handle regression tasks. See \code{\link{FS.HGD}}; 
#' \item \code{"FIR.DM"}: fuzzy inference rules by descent method to handle regression tasks. See \code{\link{FIR.DM}}; 
#' \item \code{"GFS.FR.MOGUL"}: genetic fuzzy systems for fuzzy rule learning based on the MOGUL methodology 
#'                    to handle regression tasks. See \code{\link{GFS.FR.MOGUL}};
#' \item \code{"GFS.THRIFT"}: Thrift's technique based on genetic algorithms to handle regression tasks. See \code{\link{GFS.Thrift}};
#' \item \code{"GFS.GCCL"}: Ishibuchi's method based on genetic cooperative-competitive learning
#'                   to handle classification tasks. See \code{\link{GFS.GCCL}};
#' \item \code{"FH.GBML"}: Ishibuchi's method based on hybridization of genetic cooperative-competitive learning and Pittsburgh to handle
#'                   classification tasks. See \code{\link{FH.GBML}};
#' \item \code{"SLAVE"}: structural learning algorithm on vague environment to handle classification tasks. See \code{\link{SLAVE}};
#' \item \code{"GFS.LT.RS"}: genetic algorithm for lateral tuning and rule selection. See \code{\link{GFS.LT.RS}} 
#' }
#' @param control a list containing all arguments, depending on the learning algorithm to use. The following list are 
#'                  parameters required for each methods, whereas their descriptions will be explained later on.
#' \itemize{
#' \item \code{WM}: 
#' 
#'     \code{list(num.labels, type.mf, type.tnorm, type.defuz,}
#'
#'     \code{type.implication.func, name)}
#'
#' \item \code{HYFIS}: 
#'
#'     \code{list(num.labels, max.iter, step.size, type.tnorm,}
#'
#'     \code{type.defuz, type.implication.func, name)}
#'
#' \item \code{ANFIS} and \code{FIR.DM}: 
#'
#'     \code{list(num.labels, max.iter, step.size,}
#'
#'     \code{type.tnorm, type.implication.func , name)}
#'
#' \item \code{SBC}: 
#'
#'     \code{list(r.a, eps.high, eps.low, name)}
#'
#' \item \code{FS.HGD}: 
#'
#'     \code{list(num.labels, max.iter, step.size, alpha.heuristic,}
#'
#'     \code{type.tnorm, type.implication.func, name)}
#'
#' \item \code{FRBCS.W} and \code{FRBCS.CHI}: 
#'
#'     \code{list(num.labels, type.mf, type.tnorm,}
#'
#'     \code{type.implication.func, name)}
#'
#' \item \code{DENFIS} method: 
#'
#'     \code{list(Dthr, max.iter, step.size, d, name)}
#'
#' \item \code{GFS.FR.MOGUL}: 
#'
#'     \code{list(persen_cross, max.iter, max.gen, max.tune,}
#'
#'     \code{persen_mutant, epsilon, name)}
#'
#' \item \code{GFS.THRIFT} method: 
#'
#'     \code{list(popu.size, num.labels, persen_cross,}
#'
#'     \code{max.gen, persen_mutant, type.tnorm, type.defuz,}
#'
#'     \code{type.implication.func, name)}
#'
#' \item \code{GFS.GCCL}: 
#'
#'     \code{list(popu.size, num.class, num.labels, persen_cross,}
#'
#'     \code{max.gen, persen_mutant, name)}
#'
#' \item \code{FH.GBML}: 
#'
#'     \code{list(popu.size, max.num.rule, num.class, persen_cross,}
#'
#'     \code{max.gen, persen_mutant, p.dcare, p.gccl, name)}
#'
#' \item \code{SLAVE}: 
#'
#'     \code{list(num.class, num.labels, persen_cross, max.iter,}
#'
#'     \code{max.gen, persen_mutant, k.lower, k.upper, epsilon, name)}
#'
#' \item \code{GFS.LT.RS}: 
#'
#'     \code{list(popu.size, num.labels, persen_mutant, max.gen,}
#'
#'     \code{mode.tuning, type.tnorm, type.implication.func,} 
#'
#'     \code{type.defuz, rule.selection, name)}
#'
#' }
#' 
#' \bold{Description of the \code{control} Parameters}
#' \itemize{
#' \item \code{num.labels}: a positive integer to determine the number of labels (linguistic terms). 
#'       The default value is 7.
#' \item \code{type.mf}: the following type of the membership function. The default value is \code{GAUSSIAN}. For more detail, see \code{\link{fuzzifier}}.
#'        \itemize{
#'        \item \code{TRIANGLE}: it refers triangular shape.
#'        \item \code{TRAPEZOID}: it refers trapezoid shape.
#'        \item \code{GAUSSIAN}: it refers gaussian shape.
#'        \item \code{SIGMOID}: it refers sigmoid.
#'        \item \code{BELL}: it refers generalized bell.
#'        }
#' \item \code{type.defuz}: the type of the defuzzification method as follows. The default value is \code{WAM}. For more detail, see \code{\link{defuzzifier}}.
#'       \itemize{
#'       \item \code{WAM}: the weighted average method.
#'       \item \code{FIRST.MAX}: the first maxima.
#'       \item \code{LAST.MAX}: the last maxima.
#'       \item \code{MEAN.MAX}: the mean maxima.
#'       \item \code{COG}: the modified center of gravity (COG).
#'       }
#' \item \code{type.tnorm}: the type of conjunction operator (t-norm). The following are options of t-norm available. For more detail, please have a look at \code{\link{inference}}.
#'       The default value is \code{MIN}. 
#'       \itemize{
#'       \item \code{MIN} means standard type (minimum).
#'       \item \code{HAMACHER} means Hamacher product.
#'       \item \code{YAGER} means Yager class (with tao = 1).
#'       \item \code{PRODUCT} means product.
#'       \item \code{BOUNDED} mean bounded product.
#'       }
#' \item \code{type.snorm}: the type of disjunction operator (s-norm). The following are options of s-norm available. For more detail, please have a look at \code{\link{inference}}.
#'       The default value is \code{MAX}. 
#'       \itemize{ 
#'       \item \code{MAX} means standard type (maximum). 
#'       \item \code{HAMACHER} means Hamacher sum.
#'       \item \code{YAGER} means Yager class (with tao = 1).
#'       \item \code{SUM} means sum.
#'       \item \code{BOUNDED} mean bounded sum. 
#'       }
#' \item \code{type.implication.func}: the type of implication function. The following are options of implication function available:
#'                              \code{DIENES_RESHER}, \code{LUKASIEWICZ}, \code{ZADEH},
#'                              \code{GOGUEN}, \code{GODEL}, \code{SHARP}, \code{MIZUMOTO},
#'                              \code{DUBOIS_PRADE}, and \code{MIN}.
#'                              For more detail, please have a look at \code{\link{WM}}. The default value is \code{ZADEH}. 
#' \item \code{name}: a name for the model. The default value is \code{"sim-0"}.
#' \item \code{max.iter}: a positive integer to determine the maximal number of iterations. 
#'       The default value is 10.
#' \item \code{step.size}: the step size of the gradient descent, a real number between 0 and 1. 
#'       The default value is 0.01.
#' \item \code{r.a}: a positive constant which is effectively the radius defining a neighborhood. 
#'       The default value is 0.5.
#' \item \code{eps.high}: an upper threshold value. The default value is 0.5.
#' \item \code{eps.low}: a lower threshold value. The default value is 0.15.
#' \item \code{alpha.heuristic}: a positive real number representing a heuristic value. 
#'       The default value is 1.
#' \item \code{Dthr}: the threshold value for the envolving clustering method (ECM), between 0 and 1. 
#'       The default value is 0.1.
#' \item \code{d}: a parameter for the width of the triangular membership function. 
#'       The default value is 2.
#' \item \code{persen_cross}: a probability of crossover. The default value is 0.6.
#' \item \code{max.gen}: a positive integer to determine the maximal number of generations of the genetic algorithm. 
#'       The default value is 10.
#' \item \code{max.tune}: a positive integer to determine the maximal number of tuning iterations.
#'       The default value is 10.
#' \item \code{persen_mutant}: a probability of mutation. The default value is 0.3.
#' \item \code{epsilon}: a real number between 0 and 1 representing the level of generalization.
#'       A high epsilon can lead to overfitting. The default value is 0.9. 
#' \item \code{popu.size}: the size of the population which is generated in each generation. The default value is 10.
#' \item \code{max.num.rule}: the maximum size of the rules. The default value is 5.
#' \item \code{num.class}: the number of classes.
#' \item \code{p.dcare}: a probability of "don't care" attributes. The default value is 0.5.
#' \item \code{p.gccl}: a probability of the GCCL process. The default value is 0.5.
#' \item \code{k.lower}: a lower bound of the noise threshold with interval between 0 and 1. The default value is 0.
#' \item \code{k.upper}: an upper bound of the noise threshold with interval between 0 and 1. The default value is 1.
#' \item \code{mode.tuning}: a type of lateral tuning which are \code{"LOCAL"} or \code{"GLOBAL"}. The default value is \code{"GLOBAL"}.
#' \item \code{rule.selection}:a boolean value representing whether performs rule selection or not. 
#'       The default value is \code{"TRUE"}.
#' }
#'
#' @seealso \code{\link{predict}} for the prediction phase, and 
#' the following main functions of each of the methods for theoretical background and references:  \code{\link{WM}}, \code{\link{SBC}},
#' \code{\link{HyFIS}}, \code{\link{ANFIS}}, \code{\link{FIR.DM}}, \code{\link{DENFIS}}, 
#' \code{\link{FS.HGD}}, \code{\link{FRBCS.W}}, \code{\link{FRBCS.CHI}}, \code{\link{GFS.FR.MOGUL}},
#' \code{\link{GFS.Thrift}}, \code{\link{GFS.GCCL}}, \code{\link{FH.GBML}}, \code{\link{GFS.LT.RS}}, and \code{\link{SLAVE}}.
#' @return The \code{\link{frbs-object}}. 
#' @examples
#' ##################################
#' ## I. Regression Problem
#' ## Suppose data have two input variables and one output variable.  
#' ## We separate them into training, fitting, and testing data.
#' ## data.train, data.fit, data.test, and range.data are inputs 
#' ## for all regression methods.
#' ###################################
#' ## Take into account that the simulation might take a long time 
#' ## depending on the hardware you are using. The chosen parameters 
#' ## may not be optimal.
#' ## Data must be in data.frame or matrix form and the last column 
#' ## is the output variable/attribute.
#' ## The training data must be expressed in numbers (numerical data).
#' data.train <- matrix(c(5.2, -8.1, 4.8, 8.8, -16.1, 4.1, 10.6, -7.8, 5.5, 10.4, -29.0, 
#'                       5.0, 1.8, -19.2, 3.4, 12.7, -18.9, 3.4, 15.6, -10.6, 4.9, 1.9, 
#'                       -25.0, 3.7, 2.2, -3.1, 3.9, 4.8, -7.8, 4.5, 7.9, -13.9, 4.8, 
#'                       5.2, -4.5, 4.9, 0.9, -11.6, 3.0, 11.8, -2.1, 4.6, 7.9, -2.0, 
#'                       4.8, 11.5, -9.0, 5.5, 10.6, -11.2, 4.5, 11.1, -6.1, 4.7, 12.8, 
#'                       -1.0, 6.6, 11.3, -3.6, 5.1, 1.0, -8.2, 3.9, 14.5, -0.5, 5.7, 
#'                       11.9, -2.0, 5.1, 8.1, -1.6, 5.2, 15.5, -0.7, 4.9, 12.4, -0.8, 
#'                       5.2, 11.1, -16.8, 5.1, 5.1, -5.1, 4.6, 4.8, -9.5, 3.9, 13.2, 
#'                       -0.7, 6.0, 9.9, -3.3, 4.9, 12.5, -13.6, 4.1, 8.9, -10.0, 
#'                       4.9, 10.8, -13.5, 5.1), ncol = 3, byrow = TRUE)
#' colnames(data.train) <- c("inp.1", "inp.2", "out.1") 
#'
#' data.fit <- data.train[, -ncol(data.train)]
#'
#' data.test <- matrix(c(10.5, -0.9, 5.8, -2.8, 8.5, -0.6, 13.8, -11.9, 9.8, -1.2, 11.0,
#'                      -14.3, 4.2, -17.0, 6.9, -3.3, 13.2, -1.9), ncol = 2, byrow = TRUE)
#'
#' range.data <- matrix(apply(data.train, 2, range), nrow = 2)
#'
#' #############################################################
#' ## I.1 Example: Constructing an FRBS model using Wang & Mendel
#' #############################################################
#' method.type <- "WM" 
#' 
#' ## collect control parameters into a list
#' ## num.labels = 3 means we define 3 as the number of linguistic terms
#' control.WM <- list(num.labels = 3, type.mf = "GAUSSIAN", type.tnorm = "MIN",  
#' type.defuz = "WAM", type.implication.func = "ZADEH", name = "Sim-0") 
#' 
#' ## generate the model and save it as object.WM
#' object.WM <- frbs.learn(data.train, range.data, method.type, control.WM)
#'
#' #############################################################
#' ## I.2 Example: Constructing an FRBS model using SBC
#' #############################################################
#' \dontrun{method.type <- "SBC" 
#' control.SBC <- list(r.a = 0.5, eps.high = 0.5, eps.low = 0.15, name = "Sim-0")
#'
#' object.SBC <- frbs.learn(data.train, range.data, method.type, control.SBC)}
#'
#' #############################################################
#' ## I.3 Example: Constructing an FRBS model using HYFIS
#' #############################################################
#' \dontrun{method.type <- "HYFIS"
#' 
#' control.HYFIS <- list(num.labels = 5, max.iter = 50, step.size = 0.01, type.tnorm = "MIN", 
#'                       type.defuz = "COG", type.implication.func = "ZADEH", name = "Sim-0")
#' 
#' object.HYFIS <- frbs.learn(data.train, range.data, method.type, control.HYFIS)}
#'
#' #############################################################
#' ## I.4 Example: Constructing an FRBS model using ANFIS
#' #############################################################
#' \dontrun{method.type <- "ANFIS" 
#'
#' control.ANFIS <- list(num.labels = 5, max.iter = 10, step.size = 0.01, type.tnorm = "MIN", 
#'                       type.implication.func = "ZADEH", name = "Sim-0") 
#'
#' object.ANFIS <- frbs.learn(data.train, range.data, method.type, control.ANFIS)}
#'
#' #############################################################
#' ## I.5 Example: Constructing an FRBS model using DENFIS
#' #############################################################
#' 
#' \dontrun{control.DENFIS <- list(Dthr = 0.1, max.iter = 10, step.size = 0.001, d = 2, 
#'                        name = "Sim-0")
#' method.type <- "DENFIS"
#' 
#' object.DENFIS <- frbs.learn(data.train, range.data, method.type, control.DENFIS)}
#'
#' #############################################################
#' ## I.6 Example: Constructing an FRBS model using FIR.DM
#' #############################################################
#' \dontrun{method.type <- "FIR.DM"
#'  
#' control.DM <- list(num.labels = 5, max.iter = 10, step.size = 0.01, type.tnorm = "MIN", 
#'                      type.implication.func = "ZADEH", name = "Sim-0") 
#' object.DM <- frbs.learn(data.train, range.data, method.type, control.DM)}
#'
#' #############################################################
#' ## I.7 Example: Constructing an FRBS model using FS.HGD
#' #############################################################
#' \dontrun{method.type <- "FS.HGD" 
#'  
#' control.HGD <- list(num.labels = 5, max.iter = 10, step.size = 0.01, 
#'                alpha.heuristic = 1, type.tnorm = "MIN",  
#'                type.implication.func = "ZADEH", name = "Sim-0") 
#' object.HGD <- frbs.learn(data.train, range.data, method.type, control.HGD)}
#'
#' #############################################################
#' ## I.8 Example: Constructing an FRBS model using GFS.FR.MOGUL
#' #############################################################
#' \dontrun{method.type <- "GFS.FR.MOGUL" 
#'  
#' control.GFS.FR.MOGUL <- list(persen_cross = 0.6, 
#'                     max.iter = 5, max.gen = 2, max.tune = 2, persen_mutant = 0.3, 
#'                     epsilon = 0.8, name="sim-0") 
#' object.GFS.FR.MOGUL <- frbs.learn(data.train, range.data, 
#'                        method.type, control.GFS.FR.MOGUL)}
#'
#' #############################################################
#' ## I.9 Example: Constructing an FRBS model using Thrift's method (GFS.THRIFT)
#' #############################################################
#' \dontrun{method.type <- "GFS.THRIFT" 
#'  
#' control.Thrift <- list(popu.size = 6, num.labels = 3, persen_cross = 1, 
#'                       max.gen = 5, persen_mutant = 1, type.tnorm = "MIN", 
#'                       type.defuz = "COG", type.implication.func = "ZADEH", 
#'                       name="sim-0") 
#' object.Thrift <- frbs.learn(data.train, range.data, method.type, control.Thrift)}
#' 
#' ##############################################################
#' ## I.10 Example: Constructing an FRBS model using
#' ##      genetic for lateral tuning and rule selection (GFS.LT.RS)
#' #############################################################
#' ## Set the method and its parameters
#' \dontrun{method.type <- "GFS.LT.RS" 
#'   
#' control.lt.rs <- list(popu.size = 5, num.labels = 5, persen_mutant = 0.3,
#'	               max.gen = 10, mode.tuning = "LOCAL", type.tnorm = "MIN",  
#'                 type.implication.func = "ZADEH", type.defuz = "WAM", 
#'                 rule.selection = TRUE, name="sim-0")
#'
#' ## Generate fuzzy model
#' object.lt.rs <- frbs.learn(data.train, range.data, method.type, control.lt.rs)}
#'
#' #############################################################
#' ## II. Classification Problems 
#' #############################################################
#' ## The iris dataset is shuffled and divided into training and 
#' ## testing data. Bad results in the predicted values may result
#' ## from casual imbalanced classes in the training data.
#' ## Take into account that the simulation may take a long time 
#' ## depending on the hardware you use. 
#' ## One may get better results with other parameters. 
#' ## Data are in data.frame or matrix form and the last column is 
#' ## the output variable/attribute
#' ## The data must be expressed in numbers (numerical data).
#' 
#' data(iris)
#' irisShuffled <- iris[sample(nrow(iris)),]
#' irisShuffled[,5] <- unclass(irisShuffled[,5])
#' tra.iris <- irisShuffled[1:105,]
#' tst.iris <- irisShuffled[106:nrow(irisShuffled),1:4]
#' real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)
#'
#' ## Please take into account that the interval needed is the range of input data only.
#' range.data.input <- matrix(apply(iris[, -ncol(iris)], 2, range), nrow = 2)
#' 
#' ######################################################### 
#' ## II.1 Example: Constructing an FRBS model using 
#' ##      FRBCS with weighted factor based on Ishibuchi's method
#' ###############################################################
#' ## generate the model
#' \dontrun{method.type <- "FRBCS.W"
#' control <- list(num.labels = 3, type.mf = "TRIANGLE", type.tnorm = "MIN", 
#'                type.implication.func = "ZADEH", name = "sim-0") 
#' 
#' object <- frbs.learn(tra.iris, range.data.input, method.type, control)
#' 
#' ## conduct the prediction process
#' res.test <- predict(object, tst.iris)}
#'
#' ######################################################### 
#' ## II.2 Example: Constructing an FRBS model using 
#' ##      FRBCS based on Chi's method
#' ###############################################################
#' ## generate the model
#' \dontrun{method.type <- "FRBCS.CHI"
#' control <- list(num.labels = 7, type.mf = "TRIANGLE", type.tnorm = "MIN", 
#'                type.implication.func = "ZADEH", name = "sim-0") 
#' 
#' object <- frbs.learn(tra.iris, range.data.input, method.type, control)
#' 
#' ## conduct the prediction process
#' res.test <- predict(object, tst.iris)}
#'
#' ######################################################### 
#' ## II.3 The example: Constructing an FRBS model using GFS.GCCL
#' ###############################################################
#' \dontrun{method.type <- "GFS.GCCL" 
#' 
#' control <- list(popu.size = 5, num.class = 3, num.labels = 5, persen_cross = 0.9, 
#'                     max.gen = 2, persen_mutant = 0.3,
#'                     name="sim-0") 
#' ## Training process
#' ## The main result of the training is a rule database which is used later for prediction.
#' object <- frbs.learn(tra.iris, range.data.input, method.type, control)
#'
#' ## Prediction process
#' res.test <- predict(object, tst.iris)}
#'
#' ######################################################### 
#' ## II.4 Example: Constructing an FRBS model using FH.GBML
#' ###############################################################
#' \dontrun{method.type <- "FH.GBML" 
#'	 
#'	control <- list(popu.size = 5, max.num.rule = 5, num.class = 3, 
#'				persen_cross = 0.9, max.gen = 2, persen_mutant = 0.3, p.dcare = 0.5, 
#'              p.gccl = 1, name="sim-0") 
#'	 
#'	## Training process
#'	## The main result of the training is a rule database which is used later for prediction.
#'	object <- frbs.learn(tra.iris, range.data.input, method.type, control)
#'
#'	## Prediction process
#'	res.test <- predict(object, tst.iris)}
#'
#' ######################################################### 
#' ## II.5 The example: Constructing an FRBS model using SLAVE
#' ###############################################################
#' \dontrun{method.type <- "SLAVE" 
#'	 
#'	control <- list(num.class = 3, num.labels = 5,
#'				persen_cross = 0.9, max.iter = 5, max.gen = 3, persen_mutant = 0.3, 
#'              k.lower = 0.25, k.upper = 0.75, epsilon = 0.1, name="sim-0") 
#'	 
#'	## Training process
#'	## The main result of the training is a rule database which is used later for prediction.
#'	object <- frbs.learn(tra.iris, range.data.input, method.type, control)
#'
#'	## Prediction process
#'	res.test <- predict(object, tst.iris)}
#'
#' @export
frbs.learn <- function(data.train, range.data = NULL, method.type = c("WM"), control=list()){
	options(verbose=FALSE) 
	
	## get type of method 
	method.type <- toupper(method.type)

	## get names of variables
	colnames.var <- colnames(data.train)

	## initialize mod
	mod <- NULL

	## condition if data.train is in data frame type
	if (class(data.train) != "matrix"){
		data.train <- as.matrix(data.train)
	}

	## if user did not give range of data, calculate from data
	if (is.null(range.data)){
		## for classification methods, we only need range of input data
		if (any(method.type == c("FRBCS.W", "FRBCS.CHI", "GFS.GCCL", "FH.GBML", "SLAVE"))){
			dt.min <- matrix(do.call(pmin, lapply(1:nrow(data.train[, -ncol(data.train), drop = FALSE]), function(i)data.train[i, -ncol(data.train), drop = FALSE])), nrow = 1)
			dt.max <- matrix(do.call(pmax, lapply(1:nrow(data.train[, -ncol(data.train), drop = FALSE]), function(i)data.train[i, -ncol(data.train), drop = FALSE])), nrow = 1)
		} else {
			dt.min <- matrix(do.call(pmin, lapply(1:nrow(data.train), function(i)data.train[i,])), nrow = 1)
			dt.max <- matrix(do.call(pmax, lapply(1:nrow(data.train), function(i)data.train[i,])), nrow = 1)		
		}
		range.data <- rbind(dt.min, dt.max)
	}    
	
	## Wang & Mendel's technique and HYFIS
	if(any(method.type == c("WM", "HYFIS"))){

		## getting all of parameters
		control <- setDefaultParametersIfMissing(control, list(num.labels = 7, max.iter = 10,
					step.size = 0.01, type.mf = "GAUSSIAN", type.defuz = "WAM", type.tnorm = "MIN",
						type.snorm = "MAX", type.implication.func = "ZADEH", name="sim-0"))
		
		## get parameters
		range.data.ori <- range.data
		data.train.ori <- data.train
		num.labels <- control$num.labels
		type.mf <- control$type.mf
		name <- control$name
		type.tnorm <- control$type.tnorm
		type.snorm <- control$type.snorm
		type.defuz <- control$type.defuz
		type.implication.func <- control$type.implication.func
		## replicate num of labels for each attributes
		num.labels <- matrix(rep(num.labels, ncol(range.data)), nrow=1)

		## normalize data training
		data.tra.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
		
		## generate FRBS model
		if (method.type == "WM"){
			modelSpecific <- WM(data.tra.norm, num.labels, type.mf, type.tnorm, type.implication.func)
			
			## collect results as model
			mod <- modelSpecific
			mod$type.model <- "MAMDANI"
			mod$func.tsk <- NULL
			mod$type.defuz <- type.defuz
			mod$type.snorm <- type.snorm
			mod$range.data.ori <- range.data.ori
		}
		else if(method.type == "HYFIS"){
			
			max.iter <- control$max.iter
			step.size <- control$step.size
			modelSpecific <- HyFIS(data.tra.norm, num.labels, max.iter, step.size, type.tnorm, type.snorm, type.defuz, type.implication.func)
			mod <- modelSpecific
			mod$range.data.ori <- range.data.ori
		}
	}
	
	## Substractive clustering approach
	else if(method.type == "SBC"){
		
		## get all of parameters
		control <- setDefaultParametersIfMissing(control, list(r.a = 0.5, eps.high = 0.5, eps.low = 0.15, name ="sim-0"))
		r.a <- control$r.a
		eps.high <- control$eps.high
		eps.low <- control$eps.low
		name <- control$name
		range.data.ori <- range.data
		
		## generate FRBS model 
		modelSpecific <- SBC(data.train, range.data.ori, r.a, eps.high, eps.low)
		mod <- modelSpecific
	}

	## Takagi Sugeno Kang: ANFIS, FIR.DM, FS.HGD
	else if (any(method.type == c("ANFIS", "FIR.DM", "FS.HGD"))){

		## get all of parameters
		control <- setDefaultParametersIfMissing(control, list(num.labels = 7, max.iter = 10, 
				   step.size = 0.01, type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH", alpha.heuristic = 1, name="sim-0"))
						   
		range.data.ori <- range.data
		data.train.ori <- data.train
		n.labels <- control$num.labels
		max.iter <- control$max.iter
		step.size <- control$step.size
		type.tnorm <- control$type.tnorm
		type.snorm <- control$type.snorm
		type.implication.func <- control$type.implication.func
		name <- control$name
		
		## normalize data training
		data.tra.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
		
		## generate labels of each variables
		num.labels <- matrix(rep(n.labels, ncol(range.data.ori)), nrow=1)

		## generate FRBS model
		if (method.type == "ANFIS") {
			modelSpecific <- ANFIS(data.tra.norm, num.labels, max.iter, step.size, type.tnorm, type.snorm, type.implication.func)	
		}
		else if (method.type == "FIR.DM"){
			modelSpecific <- FIR.DM(data.tra.norm, num.labels, max.iter, step.size, type.tnorm, type.snorm, type.implication.func)
		}
		else if (method.type == "FS.HGD"){
			alpha.heuristic <- control$alpha.heuristic
			modelSpecific <- FS.HGD(data.tra.norm, num.labels, max.iter, step.size, alpha.heuristic, type.tnorm, type.snorm, type.implication.func)
		}
		mod <- modelSpecific
		mod$range.data.ori <- range.data.ori
	}

	## Fuzzy rule-based classification systems: FRBCS.W, FRBCS.CHI
	else if (any(method.type == c("FRBCS.W", "FRBCS.CHI"))){

		## get all of parameters 
		control <- setDefaultParametersIfMissing(control, list(num.labels = 7, type.mf = "GAUSSIAN", 
					  type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH", name="sim-0"))			
		range.data.input <- range.data
		n.labels <- control$num.labels
		## generate labels of each variables
		num.labels <- matrix(rep(n.labels, ncol(data.train)), nrow=1)
		
		type.mf <- control$type.mf
		type.tnorm <- control$type.tnorm
		type.snorm <- control$type.snorm
		type.implication.func <- control$type.implication.func
		name <- control$name

		## make range of data according to class on data training
		range.data.out <- matrix(c(min(data.train[, ncol(data.train)], na.rm = TRUE) - 0.4999, max(data.train[, ncol(data.train)], na.rm = TRUE) + 0.4999), nrow = 2)
		num.class <- floor(max(range.data.out))
		
		## normalize range of data and data training
		range.data.norm <- range.data.input
		range.data.norm[1, ] <- 0
		range.data.norm[2, ] <- 1	
		range.data.ori <- range.data.input
		range.data.inout <- cbind(range.data.norm, range.data.out)	
		data.tra.norm <- norm.data(data.train[, 1 : (ncol(data.train) - 1)], range.data.ori, min.scale = 0, max.scale = 1)
		data.train <- cbind(data.tra.norm, matrix(data.train[, ncol(data.train)], ncol = 1))
		
		## generate FRBS model
		## FRBCS.W
		if (method.type == "FRBCS.W"){
			modelSpecific <- FRBCS.W(range.data.inout, data.train, num.labels, num.class, type.mf, type.tnorm, type.snorm, type.implication.func)
		}
		## FRBCS.CHI
		else if (method.type == "FRBCS.CHI"){
			modelSpecific <- FRBCS.CHI(range.data.inout, data.train, num.labels, num.class, type.mf, type.tnorm, type.snorm, type.implication.func)
		}
		mod <- modelSpecific
		mod$range.data.ori <- range.data.ori
	}

	## Clustering approach: DENFIS
	else if (method.type == "DENFIS"){	
		## get all of parameters
		control <- setDefaultParametersIfMissing(control, list(Dthr = 0.1, max.iter = 100, step.size = 0.01, d = 2, name="sim-0"))
		Dthr <- control$Dthr
		max.iter <- control$max.iter
		step.size <- control$step.size
		d <- control$d
		name <- control$name
		range.data.ori <- range.data
		
		## generate FRBS model
		modelSpecific <- DENFIS(data.train, range.data.ori, Dthr, max.iter, step.size, d)
		mod <- modelSpecific
	}
	
	## Approximate model: GFS.FR.MOGUL
	else if (method.type == "GFS.FR.MOGUL"){
		## get all of parameters
		control <- setDefaultParametersIfMissing(control, list(persen_cross = 0.6, max.iter = 10, 
						 max.gen = 10, max.tune = 10, persen_mutant = 0.3, epsilon = 0.8, name="sim-0"))

		 ## getting all of parameters
		range.data.ori <- range.data
		data.train.ori <- data.train
		persen_cross <- control$persen_cross
		persen_mutant <- control$persen_mutant
		max.iter <- control$max.iter
		max.gen <- control$max.gen
		max.tune <- control$max.tune
		epsilon <- control$epsilon
		name <- control$name
		
		## normalize data training
		data.tra.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
			
		modelSpecific <- GFS.FR.MOGUL(data.tra.norm, persen_cross, persen_mutant, 
									 max.iter, max.gen, max.tune, range.data.ori, epsilon)
		mod <- modelSpecific
	}

	## Thrift's technique using genetic algorithms
	else if (method.type == "GFS.THRIFT"){
		## get all of parameters
		control <- setDefaultParametersIfMissing(control, list(popu.size = 10, num.labels = 3, 
						 persen_cross = 0.6, max.gen = 10, persen_mutant = 0.3, type.defuz = "WAM", 
						 type.tnorm = "MIN", type.snorm = "MAX", type.mf = "TRIANGLE", type.implication.func = "ZADEH", name="sim-0"))
		
		## getting all of parameters
		range.data.ori <- range.data
		data.train.ori <- data.train
		popu.size <- control$popu.size
		persen_cross <- control$persen_cross
		persen_mutant <- control$persen_mutant
		max.gen <- control$max.gen
		name <- control$name
		n.labels <- control$num.labels
		type.defuz <- control$type.defuz
		type.tnorm <- control$type.tnorm
		type.snorm <- control$type.snorm
		type.mf <- control$type.mf
		type.implication.func <- control$type.implication.func
		
		## normalize data training
		data.tra.norm <- norm.data(data.train.ori, range.data.ori, min.scale = 0, max.scale = 1)
		num.labels <- matrix(rep(n.labels, ncol(range.data.ori)), nrow=1)
		
		## generate FRBS model
		modelSpecific <- GFS.Thrift(data.tra.norm, popu.size, num.labels, persen_cross,
			  persen_mutant, max.gen, range.data.ori, type.defuz, type.tnorm, type.snorm, type.mf, type.implication.func)
		mod <- modelSpecific
	}
	
	else if (method.type == "GFS.GCCL"){
		## get all of parameters
		control <- setDefaultParametersIfMissing(control, list(popu.size = 10, num.class = 2, num.labels = 3, persen_cross = 0.6, 
						   max.gen = 10, persen_mutant = 0.3, name="sim-0"))
		
		## getting all of parameters
		range.data.input <- range.data
		data.train.ori <- data.train
		popu.size <- control$popu.size
		persen_cross <- control$persen_cross
		persen_mutant <- control$persen_mutant
		max.gen <- control$max.gen
		name <- control$name
		n.labels <- control$num.labels
		n.class <- control$num.class
			
		num.labels <- matrix(rep(n.labels, ncol(range.data)), nrow = 1)
		num.labels <- cbind(num.labels, n.class)
		
		## normalize range of data and data training
		range.data.norm <- range.data.input
		range.data.norm[1, ] <- 0
		range.data.norm[2, ] <- 1	
		range.data.input.ori <- range.data.input
		data.tra.norm <- norm.data(data.train[, 1 : ncol(data.train) - 1], range.data.input, min.scale = 0, max.scale = 1)
		data.train <- cbind(data.tra.norm, matrix(data.train[, ncol(data.train)], ncol = 1))
		
		## generate FRBS model
		modelSpecific <- GFS.GCCL(data.train, popu.size, range.data.norm, num.labels, persen_cross, persen_mutant, max.gen, range.data.input.ori)
		
		mod <- modelSpecific
	}

	else if (method.type == "FH.GBML"){
		## get all of parameters
		control <- setDefaultParametersIfMissing(control, list(popu.size = 10, max.num.rule = 5, num.class = 3, persen_cross = 0.6, 
									   max.gen = 10, persen_mutant = 0.3, p.dcare = 0.5, p.gccl = 0.5, name="sim-0"))
		
		## getting all of parameters
		range.data.input <- range.data
		data.train.ori <- data.train
		popu.size <- control$popu.size
		persen_cross <- control$persen_cross
		persen_mutant <- control$persen_mutant
		max.gen <- control$max.gen
		name <- control$name
		num.class <- control$num.class
		max.num.rule <- control$max.num.rule
		p.dcare <- control$p.dcare
		p.gccl <- control$p.gccl
			
		## normalize data.train excluded output attribute
		data.tra.norm <- norm.data(data.train.ori[, -ncol(data.train.ori), drop = FALSE], range.data.input, min.scale = 0, max.scale = 1)
		data.tra.norm <- cbind(data.tra.norm, data.train.ori[, ncol(data.train.ori), drop = FALSE])
		
		## generate FRBS model
		modelSpecific <- FH.GBML(data.tra.norm, popu.size, max.num.rule, persen_cross, persen_mutant, max.gen, 
								 num.class, range.data.input, p.dcare, p.gccl)	
		mod <- modelSpecific
	}

	else if (method.type == "SLAVE"){
		## get all of parameters
		control <- setDefaultParametersIfMissing(control, list(num.class = 3, num.labels = 3, persen_cross = 0.6, 
					max.iter = 10, max.gen = 10, persen_mutant = 0.3, k.lower = 0, k.upper = 1, epsilon = 0.5, name="sim-0"))
		
		## getting all of parameters
		range.data.input <- range.data
		data.train.ori <- data.train
		persen_cross <- control$persen_cross
		persen_mutant <- control$persen_mutant
		max.iter <- control$max.iter
		max.gen <- control$max.gen
		name <- control$name
		num.class <- control$num.class
		num.labels <- control$num.labels
		k.lower <- control$k.lower
		k.upper <- control$k.upper
		epsilon <- control$epsilon
		
		num.labels <- matrix(rep(num.labels, ncol(range.data)), nrow = 1)
		num.labels <- cbind(num.labels, num.class)
			
		## normalize data.train excluded output attribute
		data.tra.norm <- norm.data(data.train.ori[, -ncol(data.train.ori), drop = FALSE], 
								   range.data.input, min.scale = 0, max.scale = 1)
		data.tra.norm <- cbind(data.tra.norm, data.train.ori[, ncol(data.train.ori), drop = FALSE])
		
		## generate FRBS model
		modelSpecific <- SLAVE(data.tra.norm, persen_cross, persen_mutant, max.iter, max.gen, num.labels, range.data.input, k.lower, k.upper, epsilon)	
		mod <- modelSpecific
	}

	else if (method.type == "GFS.LT.RS"){

		## get all of parameters
		control <- setDefaultParametersIfMissing(control, list(popu.size = 10, num.labels = 3, persen_mutant = 0.3,
						   max.gen = 10, mode.tuning = "GLOBAL", type.tnorm = "MIN", type.snorm = "MAX", type.defuz = "WAM", 
						   type.implication.func = "ZADEH", rule.selection = TRUE, name="sim-0"))
		
		## getting all of parameters
		range.data.ori <- range.data
		data.train.ori <- data.train
		popu.size <- control$popu.size
		persen_mutant <- control$persen_mutant
		max.gen <- control$max.gen
		name <- control$name
		n.labels <- control$num.labels
		mode.tuning <- control$mode.tuning
		rule.selection <- control$rule.selection
		type.tnorm <- control$type.tnorm
		type.snorm <- control$type.snorm
		type.implication.func <- control$type.implication.func
		type.defuz <- control$type.defuz
			
		num.labels <- matrix(rep(n.labels, ncol(range.data)), nrow = 1)
		
		## normalize range of data and data training
		range.data.norm <- range.data
		range.data.norm[1, ] <- 0
		range.data.norm[2, ] <- n.labels - 1	
		data.tra.norm <- norm.data(data.train, range.data.ori, min.scale = 0, max.scale = (n.labels - 1))
		
		## generate FRBS model
		modelSpecific <- GFS.LT.RS(data.tra.norm, popu.size, range.data.norm, num.labels, persen_mutant, max.gen, mode.tuning, 
		                   type.tnorm, type.snorm, type.implication.func, type.defuz, rule.selection, range.data.ori)
		
		mod <- modelSpecific
	}

	mod$method.type <- method.type
	mod$name <- name

	## convert numeric values into string values of parameters
	mod <- convert.params(mod)
	
	## keep colnames of training data into mod
	if (!is.null(colnames.var)) {
		mod$colnames.var <- colnames.var 
	} else {
		mod$colnames.var <- paste("var", seq(1, ncol(data.train)), sep = ".")
	}

	mod <- frbsObjectFactory(mod)

	## change rule format into IF ... THEN ...
	if (!is.null(mod$rule) && !is.null(mod$num.labels))
		mod$rule <- rep.rule(mod)

	return(mod)
}

## checking missing parameters
# @param control parameter values of each method
# @param defaults default parameter values of each method
setDefaultParametersIfMissing <- function(control, defaults) {
  for(i in names(defaults)) {
    if(is.null(control[[i]])) control[[i]] <- defaults[[i]]
  }
  control
}

#' This function creates objects of type \code{frbs}. Currently, its 
#' implementation is very basic and does no argument checking, as 
#' it is only used internally.
#' 
#' The members of the \code{frbs} object depend on the used learning method. The following list describes all of the members that can be present. 
#' \describe{
#' \item{\code{num.labels}}{the number of linguistic terms for the variables}
#' \item{\code{varout.mf}}{a matrix to generate the shapes of the membership functions for the output variable. 
#'          The first row represents the shape of the membership functions, the other rows contain the parameters that have been generated. 
#'          Whether the values of parameters within the matrix are normalized to lie between 0 and 1 or not depends on the selected method.}
#' \item{\code{rule}}{the fuzzy IF-THEN rules; In the \code{GFS.FR.MOGUL} case, a rule refers to the parameter values of the membership function 
#'          which represents the rule.}
#' \item{\code{rule.data.num}}{the fuzzy IF-THEN rules in integer format.}
#' \item{\code{varinp.mf}}{a matrix to generate the shapes of the membership functions for the input variables. 
#'           The first row represents the shape of the membership functions, 
#'           the other rows contain the non \code{NA} values representing the parameters related with their type of membership function. 
#'           For example, \code{TRAPEZOID}, \code{TRIANGLE}, and \code{GAUSSIAN} have four, three, and two values as their parameters, respectively. 
#'           Whether the values of parameters within the matrix are normalized to lie between 0 and 1 or not depends on the selected method.}
#' \item{\code{type.model}}{the type of model. Here, \code{MAMDANI} refers to the Mamdani model, and \code{TSK} refers to the Takagi Sugeno Kang model on the consequence part.}
#' \item{\code{func.tsk}}{a matrix of the Takagi Sugeno Kang model consequent part of the fuzzy IF-THEN rules.}
#' \item{\code{class}}{a matrix representing classes of \code{FRBCS} model}
#' \item{\code{num.labels}}{a number of linguistic terms on each variables/attributes.}
#' \item{\code{type.defuz}}{the type of the defuzzification method.}
#' \item{\code{type.tnorm}}{the type of the t-norm method.}
#' \item{\code{type.snorm}}{the type of the s-norm method.}
#' \item{\code{type.mf}}{the type of shapes of membership functions.}
#' \item{\code{type.implication.func}}{the type of the implication function.}
#' \item{\code{method.type}}{the type of the selected method.}
#' \item{\code{name}}{the name given to the model.}
#' \item{\code{range.data.ori}}{range of the original data (before normalization).}
#' \item{\code{cls}}{cluster centers.}
#' \item{\code{Dthr}}{the boundary parameter of the \code{DENFIS} method.}
#' \item{\code{d}}{the multiplier parameters of the \code{DENFIS} method.}
#' \item{\code{r.a}}{the neighborhood factor of \code{SBC}.}
#' \item{\code{degree.rule}}{certainty degree of rules.}
#' \item{\code{rule.data.num}}{a matrix representing the rules in integer form.}
#' \item{\code{grade.cert}}{grade of certainty for classification problems.}
#' \item{\code{alpha.heuristic}}{a parameter for the heuristic of the \code{FS.HGD} method.}
#' \item{\code{var.mf.tune}}{a matrix of parameters of membership function for lateral tuning.}
#' \item{\code{mode.tuning}}{a type of lateral tuning.}
#' \item{\code{rule.selection}}{a boolean of rule selection.}
#' \item{\code{colnames.var}}{the names of variables.}
#' }
#' 
#' @title The object factory for frbs objects
#' @param mod a list containing all the attributes for the object
#' @return an object of type \code{frbs}
#' @aliases frbs-object
frbsObjectFactory <- function(mod){
	class(mod) <- "frbs"
	return(mod)
}

#' This is the main function to obtain a final result as predicted values for all methods in this package. 
#' In order to get predicted values, this function is run using an \code{\link{frbs-object}}, which is typically generated using \code{\link{frbs.learn}}.
#' 
#' @title The frbs prediction stage
#'
#' @param object an \code{\link{frbs-object}}.
#' @param newdata a data frame or matrix (\eqn{m \times n}) of data for the prediction process, where \eqn{m} is the number of instances and 
#' \eqn{n} is the number of input variables. It should be noted that the testing data must be expressed in numbers (numerical data).
#' @param ... the other parameters (not used)
#' @seealso \code{\link{frbs.learn}} and \code{\link{frbs.gen}} for learning and model generation, 
#' and the internal main functions of each method for the theory:  
#' \code{\link{WM}}, \code{\link{SBC}}, \code{\link{HyFIS}}, \code{\link{ANFIS}}, 
#' \code{\link{FIR.DM}}, \code{\link{DENFIS}}, \code{\link{FS.HGD}}, \code{\link{FRBCS.W}}, 
#' \code{\link{GFS.FR.MOGUL}}, \code{\link{GFS.Thrift}}, \code{\link{GFS.GCCL}}, \code{\link{FRBCS.CHI}}, 
#' \code{\link{FH.GBML}}, \code{\link{GFS.LT.RS}}, and \code{\link{SLAVE}}.
#' @return The predicted values. 
#' @aliases predict
#' @examples
#' ##################################
#' ## I. Regression Problem
#' ###################################
#' ## In this example, we just show how to predict using Wang and Mendel's technique but
#' ## users can do it in the same way for other methods.
#' data.train <- matrix(c(5.2, -8.1, 4.8, 8.8, -16.1, 4.1, 10.6, -7.8, 5.5, 10.4, -29.0, 
#'                        5.0, 1.8, -19.2, 3.4, 12.7, -18.9, 3.4, 15.6, -10.6, 4.9, 1.9, 
#'                        -25.0, 3.7, 2.2, -3.1, 3.9, 4.8, -7.8, 4.5, 7.9, -13.9, 4.8, 
#'                        5.2, -4.5, 4.9, 0.9, -11.6, 3.0, 11.8, -2.1, 4.6, 7.9, -2.0, 
#'                        4.8, 11.5, -9.0, 5.5, 10.6, -11.2, 4.5, 11.1, -6.1, 4.7, 12.8, 
#'                        -1.0, 6.6, 11.3, -3.6, 5.1, 1.0, -8.2, 3.9, 14.5, -0.5, 5.7, 
#'                        11.9, -2.0, 5.1, 8.1, -1.6, 5.2, 15.5, -0.7, 4.9, 12.4, -0.8, 
#'                        5.2, 11.1, -16.8, 5.1, 5.1, -5.1, 4.6, 4.8, -9.5, 3.9, 13.2, 
#'                        -0.7, 6.0, 9.9, -3.3, 4.9, 12.5, -13.6, 4.1, 8.9, -10.0, 
#'                        4.9, 10.8, -13.5, 5.1), ncol = 3, byrow = TRUE)
#' 
#' data.fit <- matrix(c(10.5, -0.9, 5.2, 5.8, -2.8, 5.6, 8.5, -0.2, 5.3, 13.8, -11.9,
#'                      3.7, 9.8, -1.2, 4.8, 11.0, -14.3, 4.4, 4.2, -17.0, 5.1, 6.9, 
#'                      -3.3, 5.1, 13.2, -1.9, 4.6), ncol = 3, byrow = TRUE)
#'
#' newdata <- matrix(c(10.5, -0.9, 5.8, -2.8, 8.5, -0.2, 13.8, -11.9, 9.8, -1.2, 11.0,
#'                       -14.3, 4.2, -17.0, 6.9, -3.3, 13.2, -1.9), ncol = 2, byrow = TRUE)
#'
#' range.data<-matrix(c(0.9, 15.6, -29, -0.2, 3, 6.6), ncol=3, byrow = FALSE)
#' #############################################################
#' ## I.1 Example: Implementation of Wang & Mendel
#' #############################################################
#' method.type <- "WM" 
#' 
#' ## collect control parameters into a list
#' ## num.labels = 3 means we define 3 as the number of linguistic terms
#' control.WM <- list(num.labels = 3, type.mf = "GAUSSIAN", type.tnorm = "MIN", 
#'                type.snorm = "MAX", type.defuz = "WAM", 
#'                type.implication.func = "ZADEH", name = "Sim-0") 
#' 
#' ## generate the model and save it as object.WM
#' object.WM <- frbs.learn(data.train, range.data, method.type, control.WM)
#'
#' ## the prediction process
#' ## The following code can be used for all methods
#' res <- predict(object.WM, newdata)
#' 
#' @export  
#' @method predict frbs
#' @S3method predict frbs
predict.frbs <- function(object, newdata, ...) {
	options(verbose=FALSE)
	
	## Initial checking
	init.check <- validate.params(object, newdata)
	object <- init.check$object
	newdata <- init.check$newdata
		
	## get the type of method
	m.type <- object$method.type
	
	## condition for PMML using FRBCS model
	if (object$type.model == "FRBCS" && !is.null(object$pmml)){
		m.type <- "FRBCS.CHI"
	}

	## 0. MANUAL
	if (m.type == "MANUAL"){
		if (any(object$type.model == c("MAMDANI", "TSK"))){
			res <- frbs.eng(object, newdata)
		} else {
			res <- FRBCS.eng(object, newdata)
		}
	}

	## 1. WM, ANFIS, HYFIS, FIR.DM, FS.HGD approach
	else if (any(m.type == c("WM", "ANFIS", "HYFIS", "FIR.DM", "FS.HGD")) 
	                   || any(object$type.model == c("MAMDANI", "TSK"))) {
		range.data.ori <- object$range.data.ori
		range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1)]
		range.output.ori <- range.data.ori[, ncol(range.data.ori), drop = FALSE]
		data.tst.norm <- norm.data(newdata, range.input.ori, min.scale = 0, max.scale = 1)
		res.comp <- frbs.eng(object, data.tst.norm)
		res.denorm <- denorm.data(res.comp$predicted.val, range.output.ori, min.scale = 0, max.scale = 1)
		res <- res.denorm
	}

	## 2.SBC approach
	else if(m.type == "SBC"){	
		res <- SBC.test(object, newdata)
	}

	## 3. FRBCS.W and FRBCS.CHI approach
	else if(any(m.type == c("FRBCS.W", "FRBCS.CHI"))){
		data.tst.norm <- norm.data(newdata, object$range.data.ori, min.scale = 0, max.scale = 1)
		res <- FRBCS.eng(object, data.tst.norm)
	}

	## 4. DENFIS approach
	else if(m.type == "DENFIS"){
		res <- DENFIS.eng(object, newdata)
	}

	## 5. GFS.FR.MOGUL approach
	else if(m.type == "GFS.FR.MOGUL" || object$type.model == "APPROXIMATE"){	
		range.data.ori <- object$range.data.ori
		range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1), drop = FALSE]
		range.output.ori <- range.data.ori[, ncol(range.data.ori), drop = FALSE]
		data.tst.norm <- norm.data(newdata, range.input.ori, min.scale = 0, max.scale = 1)
		res.comp <- GFS.FR.MOGUL.test(object, data.tst.norm)
		res.denorm <- denorm.data(res.comp, range.output.ori, min.scale = 0, max.scale = 1)
		res <- res.denorm
	}

	## 6. GFS.THRIFT
	else if(m.type == "GFS.THRIFT"){	
		range.data.ori <- object$range.data.ori
		range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1), drop = FALSE]
		range.output.ori <- range.data.ori[, ncol(range.data.ori), drop = FALSE]
		data.tst.norm <- norm.data(newdata, range.input.ori, min.scale = 0, max.scale = 1)
		res.comp <-	GFS.Thrift.test(object, data.tst.norm)
		res.denorm <- denorm.data(res.comp, range.output.ori, min.scale = 0, max.scale = 1)
		res <- res.denorm
	}

	## 7. GFS.GCCL, FH.GBML
	else if(any(m.type == c("GFS.GCCL", "FH.GBML"))){
		newdata <- norm.data(newdata, object$range.data.ori, min.scale = 0, max.scale = 1)
		res <- GFS.GCCL.eng(object, newdata)
	}

	## 8. SLAVE
	else if (m.type == "SLAVE"){
		data.tst.norm <- norm.data(newdata, object$range.data.ori, min.scale = 0, max.scale = 1)
		res <- SLAVE.test(object, data.tst.norm)
	}

	## 9. GFS.LT.RS
	else if (m.type == "GFS.LT.RS"){
		range.data.ori <- object$range.data.ori
		num.labels <- object$num.labels
		range.input.ori <- range.data.ori[, 1:(ncol(range.data.ori) - 1), drop = FALSE]
		range.output.ori <- range.data.ori[, ncol(range.data.ori), drop = FALSE]

		data.tst.norm <- norm.data(newdata, range.input.ori, min.scale = 0, max.scale = (num.labels[1,1] - 1))

		res.comp <-	GFS.LT.RS.test(object, data.tst.norm)
		res.denorm <- denorm.data(res.comp, range.output.ori, min.scale = 0, max.scale = (num.labels[1,1] - 1))
		res <- res.denorm
	}

	return(res)
}

#' This function enables the output of a summary of the \code{\link{frbs-object}}. 
#'
#' This function displays several components of the object. The components of 
#' one particular method can be different from components of other methods.
#' The following is a description of all components which might be printed.
#' \itemize{
#' \item The name of the model: A name given by the user representing the name of the simulation
#'                or data or model.
#' \item Model was trained using: It shows which method we have been used.
#' \item The names of attributes: a list of names of training data.
#' \item The interval of training data: It is a matrix representing the original 
#'            interval of data where the first and second rows are minimum and maximum of data, 
#'            respectively. The number of columns represents the number of variables.
#' \item Type of FRBS model: a description expresses one of the following FRBS model available such as \code{"MAMDANI"}, \code{"TSK"},
#'                           \code{"FRBCS"}, \code{"CLUSTERING"}, \code{"APPROXIMATE"}, and \code{"2TUPPLE"}.
#' \item Type of membership function: a description expresses one of the following shapes of membership functions: 
#'                           \code{"GAUSSIAN"}, code{"TRIANGLE"}, \code{"TRAPEZOID"}, \code{"SIGMOID"}, and \code{"BELL"}.
#' \item Type of t-norm method: a description expresses one of the following type of t-norm: \code{"MIN"}, \code{"PRODUCT"}, \code{"HAMACHER"},
#'                           \code{"YAGER"}, and \code{"BOUNDED"}.
#' \item Type of s-norm method: a description expresses one of the following type of s-norm: \code{"MAX"}, \code{"SUM"}, \code{"HAMACHER"},
#'                           \code{"YAGER"}, and \code{"BOUNDED"}.
#' \item Type of defuzzification technique: a description expresses one of the following types: \code{"WAM"}, \code{"FIRST_MAX"}, 
#'                          \code{"LAST_MAX"}, \code{"MEAN_MAX"}, and \code{"COG"}. 
#' \item Type of implication function: a description expresses one of the following types: 
#'                          \code{"DIENES_RESHER"}, \code{"LUKASIEWICZ"},
#'                          \code{"ZADEH"}, \code{"GOGUEN"}, \code{"GODEL"}, \code{"SHARP"}, \code{"MIZUMOTO"}, 
#'                          \code{"DUBOIS_PRADE"}, and \code{"MIN"}.
#' \item The names of linguistic terms of the input variables: These names are generated
#'           automatically by frbs expressing all linguistic terms considered. Generally,
#'           these names are built by two parts which are the name of variables expressed 
#'           by \code{"v"} and the name of linguistic terms of each variables represented by \code{"a"}. 
#'           For example, \code{"v.1_a.1"} means the linguistic value \code{"a.1"} of the first variable (v.1).
#'           However, we provide different format if we set the number of linguistic terms (\code{num.labels}) to 3, 5, 7. 
#'           For example, for the number of label 3, it will be \code{"small"}, \code{"medium"}, and \code{"large"}.
#' \item The names of linguistic terms of the output variable: For the Mamdani model, since the frbs package only considers
#'           single output, the names of the linguistic terms for the output variable 
#'           are simple and clear and start with \code{"c"}. However, for the Takagi Sugeno Kang model and
#'           fuzzy rule-based classification systems, this component is always \code{NULL}.
#' \item The parameter values of membership functions of the input variables (normalized):
#'          It is represented by a matrix (\eqn{5 \times n}) where n depends on the number of 
#'          linguistic terms on the input variables and the first row of the matrix describes 
#'          a type of membership function, and the rest of rows are 
#'          their parameter values. 
#'          For example, label \code{"v.1_a.2"} has value 
#'          {4.0, 0.23, 0.43, 0.53, 0.73} on its column. It means that the label a.2 of variable v.1 
#'           has a parameter as follows. 
#'          4.0 on the first row shows \code{TRAPEZOID} shape in the middle position, 
#'          while 0.23, 0.43, 0.53, and 0.73 are corner points of a \code{TRAPEZOID}. 
#'          Furthermore, the following is the complete list of shapes of membership functions:
#'          \itemize{
#'          \item \code{TRIANGLE}: 1 on the first row and rows 2, 3, and 4 represent corner points. 
#'          \item \code{TRAPEZOID}: 2, 3, or 4 on the first row means they are \code{TRAPEZOID} in left, right and middle side, respectively,
#'                     and rows 2, 3, 4, and 5 represent corner points. But for \code{TRAPEZOID} at left or right side the fifth row is \code{NA}. 
#'          \item \code{GAUSSIAN}: 5 on the first row means it uses \code{GAUSSIAN} and second and third row represent mean and variance.
#'          \item \code{SIGMOID}: 6 on the first row and two parameters (gamma and c) on second and third rows.
#'          \item \code{BELL}: 7 on the first row and three parameters (a, b, c) on second, third, and fourth rows.
#'          }
#' \item The fuzzy IF-THEN rules: In this package, there are several models for representing
#'          fuzzy IF-THEN rules based on the method used. 
#'          \itemize{
#'          \item the Mamdani model: they are represented as a knowledge base containing two parts: 
#'          antecedent and consequent parts which are separated by a sign "THEN", as for example in the
#'          following rule:
#' 
#'          \code{IF var.1 is v.1_a.1 and var.2 is v.2_a.2 THEN var.3 is c.2}
#'          
#'          \item the Takagi Sugeno Kang model: In this model, this component only represents the antecedent
#'          of rules while the consequent part will be represented by linear equations. 
#'          \item fuzzy rule-based classification systems (FRBCS): This model is quite similar to the Takagi Sugeno Kang model,
#'          but the consequent part expresses pre-defined classes instead of a simplify of linear equations.
#'          \item approximate approach: Especially for \code{GFS.FR.MOGUL}, a matrix of parameters
#'          of membership functions is used to represent the fuzzy IF-THEN rules as well.  
#'          The representation of rules and membership functions is a matrix (\eqn{n \times (p \times m)}) where
#'          n is the number of rules and m is the number of variables while p is the number of corner points 
#'          of the membership function, if we are using \code{TRIANGLE} or \code{TRAPEZOID} then p = 3 or 4, respectively. 
#'          For example, let us consider the triangular membership function and a number of variables of 3. 
#'          The representation of rules and membership functions is as follows:
#' 
#'          \code{<<a11 a12 a13>> <<b11 b12 b13>> <<c11 c12 c13>>}.
#'          
#'          }
#' \item The linear equations on consequent parts of fuzzy IF-THEN rules: It is used in
#'         the Takagi Sugeno Kang model.
#' \item The weight of the rules or the certainty factor: For the \code{FRBCS.W} method, this shows the weight related to the rules 
#'         representing the ratio of dominance among the rules.
#' \item The cluster centers: This component is used in clustering methods representing cluster centers.
#' }
#' @title The summary function for frbs objects
#' 
#' @param object the \code{\link{frbs-object}}
#' @param ... the other parameters (not used)
#' @export  
#' @method summary frbs
#' @S3method summary frbs
summary.frbs <- function(object, ...){
	
  if(!inherits(object, "frbs")) stop("not a legitimate frbs model")
	cat("The name of model: ", object$name, "\n")
	cat("Model was trained using: ", object$method.type, "\n") 
	cat("The names of attributes: ", object$colnames.var, "\n")
	
	## display for range.data
	if (any(object$method.type == c("FRBCS.W", "FRBCS.CHI", "GFS.GCCL", "FH.GBML", "SLAVE"))){
		colnames(object$range.data.ori) <- object$colnames.var[-length(object$colnames.var)]
		rownames(object$range.data.ori) <- c("min", "max")
		cat("The interval of input data: ", "\n")
		print(object$range.data.ori)
	}
	else if (object$method.type == "MANUAL") {
		if (object$type.model == "FRBCS"){
			colnames(object$range.data.ori) <- object$colnames.var[-length(object$colnames.var)]
			rownames(object$range.data.ori) <- c("min", "max")
			cat("The interval of input data: ", "\n")
			print(object$range.data.ori)
		}
		else if (any(object$type.model == c("MAMDANI", "TSK"))){
			range.data.ori <- object$range.data.ori
			colnames(range.data.ori) <- object$colnames.var
			rownames(range.data.ori) <- c("min", "max")
			cat("The interval of training data: ", "\n")
			print(range.data.ori)
		} else {
			stop("The model is not supported")
		}
	} else {
		range.data.ori <- object$range.data.ori
		colnames(range.data.ori) <- object$colnames.var
		rownames(range.data.ori) <- c("min", "max")
		cat("The interval of training data: ", "\n")
		print(range.data.ori)
	}
	
	## display for schema of Inference	
	if (any(object$method.type == c("WM", "HYFIS", "ANFIS", "FIR.DM", "FS.HGD", "GFS.THRIFT", "FRBCS.W", "FRBCS.CHI", 
	                                "SLAVE", "GFS.GCCL", "FH.GBML", "GFS.FR.MOGUL", "GFS.LT.RS", "MANUAL"))){
		cat("Type of FRBS model:", "\n")
		print(object$type.model)
		cat("Type of membership functions:", "\n")
		print(object$type.mf)
		cat("Type of t-norm method:", "\n")
		if (object$type.tnorm == 1 || object$type.tnorm == "MIN") print("Standard t-norm (min)")
		else if (object$type.tnorm == 2 || object$type.tnorm == "HAMACHER") print("Hamacher product")
		else if (object$type.tnorm == 3 || object$type.tnorm == "YAGER") print("Yager class (with tao = 1)")
		else if (object$type.tnorm == 4 || object$type.tnorm == "PRODUCT") print("Product")
		else print("Bounded product")
		cat("Type of s-norm method:", "\n")
		if (object$type.snorm == 1 || object$type.snorm == "MAX") print("Standard s-norm")
		else if (object$type.snorm == 2 || object$type.snorm == "HAMACHER") print("Hamacher sum")
		else if (object$type.snorm == 3 || object$type.snorm == "YAGER") print("Yager class (with tao = 1)")
		else if (object$type.snorm == 4 || object$type.snorm == "SUM") print("Sum")
		else print("Bounded sum")
		
		if (any(object$method.type == c("WM", "HYFIS", "GFS.THRIFT", "GFS.LT.RS"))){
			cat("Type of defuzzification technique:", "\n")
			if (object$type.defuz == 1 || object$type.defuz == "WAM") print("Weighted average method")
			else if (object$type.defuz == 2 || object$type.defuz == "FIRST.MAX") print("first of maxima")
			else if (object$type.defuz == 3 || object$type.defuz == "LAST.MAX") print("last of maxima")
			else if (object$type.defuz == 4 || object$type.defuz == "MEAN.MAX") print("mean of maxima")
			else print("modified COG")
		}
		cat("Type of implication function:", "\n")
		print(object$type.implication.func)
	}
	
	## display for parameters of membership function and rule
	if (any(object$method.type == c("WM", "HYFIS", "GFS.THRIFT", "ANFIS", "FIR.DM", "FS.HGD", "FRBCS.W", 
	                              "FRBCS.CHI", "SLAVE", "GFS.GCCL", "FH.GBML", "GFS.LT.RS", "MANUAL"))){
		num.labels <- object$num.labels
		rule <- as.data.frame(object$rule)
		cat("The names of linguistic terms on the input variables: ", "\n")
		print(colnames(object$varinp.mf))	
		cat("The parameter values of membership function on the input variable (normalized): ", "\n")
		print(object$varinp.mf)
		
		if (any(object$method.type == c("ANFIS", "FIR.DM", "FS.HGD"))){
			colnames(num.labels) <- object$colnames.var[-length(object$colnames.var)]
			cat("The number of linguistic terms on input variables", "\n")
			print(num.labels)
			cat("The fuzzy IF-THEN rules: ", "\n")
			print(rule)
			if (ncol(object$func.tsk) > 1){	
				seq.deg <- seq(from = 1, to = (ncol(object$func.tsk) - 1), by = 1)
				coef.var <- paste("var", seq.deg, sep = ".")
				names.func <- c(coef.var, "const")	
			} else {
				names.func <- c("const")	
			}
			colnames(object$func.tsk) <- names.func
			cat("The linear equations on consequent parts of fuzzy IF-THEN rules: ", "\n")
			print(object$func.tsk)		
		}
		# else if (any(object$method.type == c("FRBCS.W", "FRBCS.CHI"))){
			# colnames(num.labels) <- object$colnames.var
			# cat("The number of linguistic terms on each variables", "\n")
			# print(num.labels)
			# cat("The fuzzy IF-THEN rules: ", "\n")
			# print(rule)
			# if (any(object$method.type == c("FRBCS.W"))){
				# cat("The weight of the rules", "\n")
				# print(object$grade.cert[, 2, drop = FALSE])
			# }
		# }
		else if (any(object$method.type == c("WM", "HYFIS", "GFS.THRIFT", "GFS.LT.RS"))){
			cat("The names of linguistic terms on the output variable: ", "\n")
			print(colnames(object$varout.mf))
			cat("The parameter values of membership function on the output variable (normalized): ", "\n")
			print(object$varout.mf)		
			colnames(num.labels) <- object$colnames.var
			cat("The number of linguistic terms on each variables", "\n")
			print(num.labels)
			cat("The fuzzy IF-THEN rules: ", "\n")
			print(rule)
		}
		
		else if (any(object$method.type == c("FRBCS.W", "FRBCS.CHI", "GFS.GCCL", "FH.GBML"))){
			colnames(num.labels) <- object$colnames.var
			cat("The number of linguistic terms on each variables", "\n")
			print(num.labels)
			cat("The fuzzy IF-THEN rules: ", "\n")
			print(rule)
			cat("The certainty factor:", "\n")
			print(object$grade.cert)
		}
		
		else if (object$method.type == c("SLAVE")){
			colnames(num.labels) <- object$colnames.var
			cat("The number of linguistic terms on each variables", "\n")
			print(num.labels)
			cat("The fuzzy IF-THEN rules: ", "\n")
			print(rule)
		}
		
		else if (object$method.type == "MANUAL"){
			if (object$type.model == "MAMDANI"){
				cat("The names of linguistic terms on the output variable: ", "\n")
				print(colnames(object$varout.mf))
				cat("The parameter values of membership function on the output variable (normalized): ", "\n")
				print(object$varout.mf)		
				colnames(num.labels) <- object$colnames.var
				cat("The number of linguistic terms on each variables", "\n")
				print(num.labels)
				cat("The fuzzy IF-THEN rules: ", "\n")
				print(rule)
			}
			else if (object$type.model == "TSK"){
				colnames(num.labels) <- object$colnames.var[-length(object$colnames.var)]
				cat("The number of linguistic terms on input variables", "\n")
				print(num.labels)
				cat("The fuzzy IF-THEN rules: ", "\n")
				print(rule)
				if (ncol(object$func.tsk) > 1){	
					seq.deg <- seq(from = 1, to = (ncol(object$func.tsk) - 1), by = 1)
					coef.var <- paste("var", seq.deg, sep = ".")
					names.func <- c(coef.var, "const")	
				} else {
					names.func <- c("const")	
				}
				colnames(object$func.tsk) <- names.func
				cat("The linear equations on consequent parts of fuzzy IF-THEN rules: ", "\n")
				print(object$func.tsk)
			} else {
				colnames(num.labels) <- object$colnames.var
				cat("The number of linguistic terms on each variables", "\n")
				print(num.labels)
				cat("The fuzzy IF-THEN rules: ", "\n")
				print(rule)
			}
		}
		if (object$method.type == "GFS.LT.RS"){
			cat("The mode of lateral tuning:", "\n")
			print(object$mode.tuning)
			if (object$mode.tuning == "LOCAL"){
				cat("The values of lateral tuning of membership function parameters:", "\n")
				print(object$var.mf.tune)
			}
		}
	}
	
	## display for clustering approach
	if (any(object$method.type == c("DENFIS", "SBC"))){		
		if (any(object$method.type == c("DENFIS", "SBC"))){
			colnames(object$cls) <- object$colnames.var
			cat("The cluster centers: ", "\n")
			print(object$cls)
		}
	}
	
	## display for GFS.FR.MOGUL
	if (object$method.type == c("GFS.FR.MOGUL")){
		cat("The fuzzy IF-THEN rules: ", "\n")
		print(object$rule)	
	}
	
  invisible(object)	
}

#' This function can be used to plot the shapes of the membership functions.
#'
#' @title The plotting function
#' 
#' @param object an \code{\link{frbs-object}} or a list of parameters to plot membership functions when we build the frbs model without learning. 
#'        For plotting using the list, there are several parameters that must be inserted in params as follows.
#'        \itemize{
#'        \item \code{var.mf}: a matrix of membership function of input and output variables. Please see \code{\link{fuzzifier}}.
#'        \item \code{range.data.ori}: a matrix (\eqn{2 \times n}) containing the range of the data, where \eqn{n} is the number of variables, and
#'        first and second rows are the minimum and maximum values, respectively. 
#'        \item \code{num.labels}: the number of linguistic terms of the input and output variables. 
#' 
#'              For example: \code{num.labels <- matrix(c(3, 3, 3), nrow = 1)}
#'
#'              It means we have 3 linguistic values/labels for two input variables and one output variable.
#'        \item \code{names.variables}: a list of names of variables. 
#'
#'              For example: \code{names.variables <- c("input1", "input2", "output1")}
#'        }
#' @examples
#' ## The following examples contain two different cases which are
#' ## using an frbs-object and the manual way. 
#' ## 
#' ## 1. Plotting using frbs object.
#' data(iris)
#' irisShuffled <- iris[sample(nrow(iris)),]
#' irisShuffled[,5] <- unclass(irisShuffled[,5])
#' tra.iris <- irisShuffled[1:105,]
#' tst.iris <- irisShuffled[106:nrow(irisShuffled),1:4]
#' real.iris <- matrix(irisShuffled[106:nrow(irisShuffled),5], ncol = 1)
#'
#' ## Please take into account that the interval needed is the range of input data only.
#' range.data.input <- matrix(c(4.3, 7.9, 2.0, 4.4, 1.0, 6.9, 0.1, 2.5), nrow=2)
#' 
#' ## generate the model
#' method.type <- "FRBCS.W"
#' control <- list(num.labels = 7, type.mf = 1) 
#' \dontrun{object <- frbs.learn(tra.iris, range.data.input, method.type, control)} 
#' 
#' ## plot the frbs object
#' \dontrun{plotMF(object)}
#'
#' ## 2. Plotting using params.
#' ## Define shape and parameters of membership functions of input variables.
#' ## Please see the fuzzifier function of how to contruct the matrix.
#' varinp.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
#'                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
#'                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
#'                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
#'                       nrow = 5, byrow = FALSE)
#'
#' ## Define the shapes and parameters of the membership functions of the output variables.
#' varout.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
#'                       nrow = 5, byrow = FALSE)
#' var.mf <- cbind(varinp.mf, varout.mf)
#' range.data <- matrix(c(0,100, 0, 100, 0, 100, 0, 100, 0, 100), nrow=2)
#' num.labels <- matrix(c(3,3,3,3,3), nrow = 1)
#' names.variables <- c("input1", "input2", "input3", "input4", "output1")
#'
#' ## plot the membership function.
#' \dontrun{plotMF(object = list(var.mf = var.mf, range.data.ori = range.data, 
#'           num.labels = num.labels, names.variables = names.variables))}
#' @export
plotMF <- function(object) {
  
  ## define whether object as a list or frbs object
  if (inherits(object, "frbs")) { 
	method.type <- object$method.type
  }
  else if (class(object) == "list") {
	method.type <- c("MANUAL")
  } else {
	stop("please input a frbs object or a list of parameters")
  }
  
  ## check whether sort of method can plot or not and get some parameters
  if (any(method.type == c("WM", "HYFIS", "ANFIS", "FS.HGD", "GFS.THRIFT", "FIR.DM", "FRBCS.W", 
                           "FRBCS.CHI", "GFS.GCCL", "FH.GBML", "SLAVE", "GFS.LT.RS", "MANUAL"))){
	  if (any(method.type == c("WM", "HYFIS", "GFS.THRIFT"))){
		  range.data <- matrix(nrow = 2, ncol = ncol(object$num.labels))
		  range.data[1, ] <- 0
		  range.data[2, ] <- 1
		  num.varinput <- ncol(object$num.labels)
		  num.fvalinput <- object$num.labels 
		  ## get parameters of MF for all variables
		  var.mf <- cbind(object$varinp.mf, object$varout.mf) 
		  
	  }
	  else if (method.type == "GFS.LT.RS"){
		  range.data <- matrix(nrow = 2, ncol = ncol(object$num.labels))
		  range.data[1, ] <- 0
		  range.data[2, ] <- object$num.labels[1,1] - 1  
		  num.varinput <- ncol(object$num.labels)
		  num.fvalinput <- object$num.labels 
		  ## get parameters of MF for all variables
		  var.mf <- cbind(object$varinp.mf, object$varout.mf) 
		  if (object$mode.tuning == "LOCAL")
			  warning("The plot shows a original shape of membership functions only")
	  }
	 else if (any(method.type == c("ANFIS", "FS.HGD", "FIR.DM", "FRBCS.W", "FRBCS.CHI", "GFS.GCCL", "FH.GBML", "SLAVE"))){
		  var.mf <- object$varinp.mf
		  ## for TSK and FRBCS model, plotting can be done only for input variables		  
		  if (any(method.type == c("ANFIS", "FIR.DM", "FS.HGD"))){
			range.data <- matrix(nrow = 2, ncol = (ncol(object$range.data.ori) - 1))
			range.data[1, ] <- 0
			range.data[2, ] <- 1
			num.fvalinput <- object$num.labels
			num.varinput <- ncol(object$range.data.ori) - 1
		  } else { 
			range.data <- matrix(nrow = 2, ncol = ncol(object$range.data.ori))
			range.data[1, ] <- 0
			range.data[2, ] <- 1
			num.fvalinput <- object$num.labels[, -ncol(object$num.labels), drop = FALSE]  
			num.varinput <- ncol(object$range.data.ori)
		  }		  
	  } 
	  ## get parameters when plotting from manual condition
	  else {
		  if (is.null(object$num.labels) || is.null(object$range.data.ori)) {
			  stop("please input the matrix of num.labels and range.data")
		  } else {
			  if (!is.null(object$type.model) && object$type.model == "TSK") num.varinput <- ncol(object$range.data.ori) - 1
			  else  num.varinput <- ncol(object$range.data.ori)
			  if (!is.null(object$varout.mf)) var.mf <- cbind(object$varinp.mf, object$varout.mf)
			  else if (!is.null(object$var.mf)) var.mf <- object$var.mf
			  else var.mf <- object$varinp.mf
			  range.data <- object$range.data.ori
			  num.fvalinput <- object$num.labels
			  if (!is.null(object$names.variables)){ 
				names.variables <- object$names.variables
			  } else {
				names.variables <- paste("var", seq(1, ncol(object$range.data.ori)), sep = ".")
			  }
		  }
	  }
	  
	  ##get number of column of var.mf
	  col.var.mf <- ncol(var.mf)
	  
	  ## counter is used to make continue index j
	  counter <- 1
	  
	  ## k is used to index number linguistic value on each variable on varinput.
	  k <- 1
	  
	  ## make row plot 
	  op <- par(mfrow = c(ceiling(num.varinput/2), 2))
	  
	  ## loop as many as number of input variable
	  for (i in 1 : num.varinput){
		j <- counter
		
		## Initialize plot
		MF.degree <- function(x){
		  y <- x - x
		  return (y)
		}
		if (!is.null(object$colnames.var)) { 
			names <- object$colnames.var[i]
		} else {
			names <- names.variables[i]
		}
				
		curve(MF.degree, range.data[1, i], range.data[2, i], ylim = c(0, 1), col = "white")
		title(main = names)

		## loop as many as number of column of var.mf
		for (j in counter : col.var.mf){
		  
		  ## make boundary point
		  oo <- range.data[1, i]
		  aa <- var.mf[2, j]
		  bb <- var.mf[3, j]
		  cc <- var.mf[4, j]
		  dd <- var.mf[5, j]
		  mm <- range.data[2, i]
		  
		  ##condition for triangular type
		  if (var.mf[1, j] == 1){				
			
			## make a function for plotting, args (x) is sequence data
			f.y1 <- function(x){
			  
			  ## range of input data
			  p.0 <- x[x >= oo & x <= aa]
			  p.1 <- x[x > aa & x <= bb]
			  p.2 <- x[x > bb & x <= cc]
			  p.3 <- x[x > cc & x <= mm]
			  
			  ## build functions 
			  y0 <- (p.0 - p.0)
			  y1 <- (p.1 - aa) / (bb - aa)			
			  y2 <- (p.2 - cc) / (bb - cc)			
			  y4 <- (p.3 - p.3)
			  y <- c(y0, y1, y2, y4)
			  
			  return (y)
			}
			
			curve(f.y1, oo, mm, add = TRUE, col = "violet")
		  }
		  ## condition for trapezoid in left side
		  else if (var.mf[1, j] == 2){
			
			f.y2 <- function(x){
			  
			  p.1 <- x[x <= bb]
			  p.2 <- x[x > bb & x <= cc]
			  p.3 <- x[x > cc & x <= mm]
			  
			  y1 <- (p.1 - p.1 + 1)			
			  y2 <- (p.2 - cc) / (bb - cc)			
			  y3 <- (p.3 - p.3)
			  
			  y <- c(y1, y2, y3)
			  
			  return (y)
			}
			
			curve(f.y2, oo, mm, add=TRUE, col = "blue")
			
		  }
		  ## condition for trapezoid in right side
		  else if (var.mf[1, j] == 3){
			f.y3 <- function(x){
			  
			  p.0 <- x[x >= oo & x <= aa]
			  p.1 <- x[x > aa & x <= bb]
			  p.2 <- x[x > bb & x <= mm]
			  
			  y0 <- (p.0 - p.0)
			  y1 <- (p.1 - aa) / (bb - aa)
			  y2 <- (p.2 - p.2 + 1)			
			  y <- c(y0, y1, y2)
			  
			  return (y)
			}
			
			curve(f.y3, oo, mm, add = TRUE, col = "green")
			
		  }
		  ## condition for trapezoid in the middle
		  else if (var.mf[1, j] == 4){
			
			f.y4 <- function(x){
			  
			  p.0 <- x[x >= oo & x <= aa]
			  p.1 <- x[x > aa & x <= bb]
			  p.2 <- x[x > bb & x <= cc]
			  p.3 <- x[x > cc & x <= dd]
			  p.4 <- x[x > dd & x <= mm]
			  
			  y0 <- (p.0 - p.0)
			  y1 <- (p.1 - aa) / (bb - aa)
			  y2 <- (p.2 - p.2 + 1)			
			  y3 <- (p.3 - dd) / (cc - dd)
			  y4 <- (p.4 - p.4)
			  
			  y <- c(y0, y1, y2, y3, y4)
			  
			  return (y)
			}
			
			curve(f.y4, oo, mm, add = TRUE, col = "red")
			
		  }
		  ## condition for gaussian shape
		  else if (var.mf[1, j] == 5){
			
			f.y5 <- function(x){
			  y <- exp(- 0.5 * (x - aa) ^ 2 / bb ^ 2)
			  return (y)
			}
			## plot the functions
			curve(f.y5, oo, mm, add = TRUE, col = "gray")
		  }
		  ## condition for sigmoid/logistic
		  else if (var.mf[1, j] == 6){
			
			f.y6 <- function(x){
			  y <- 1 / (1 + exp(- aa * (x - bb)))
			  return (y)
			}
			#plot the functions
			curve(f.y6, oo, mm, add = TRUE, col = "black")
		  }
		  ## condition for generalized bell
		  else if (var.mf[1, j] == 7){
			
			f.y7 <- function(x){
			  y <- 1 / (1 + abs((x - cc)/aa) ^ (2 * bb))
			  return (y)
			}
			
			## plot the functions
			curve(f.y7, oo, mm, add = TRUE, col = "black")
		  }
		  
		  counter <- j + 1 
		  k <- k + 1
		  if (k > num.fvalinput[1, i]){
			k <- 1
			break
		  }
		}
	  }

	  par(op)
  } else {
	print("The plot is not supported by the used method, please have a look the documentation")
 }
}

# This function can be used to make representations of rules.
#
# @title The rule representing function
# 
# @param object an \code{\link{frbs-object}}.
# @export
rep.rule <- function(object){
	if(!inherits(object, "frbs")) stop("not a legitimate frbs model")

	## get names of variables
	colnames.var <- object$colnames.var
	
	## manipulate num.labels 
	if (object$type.model == "TSK")
		object$num.labels = cbind(object$num.labels, object$num.labels[1,1])
	
	## make description on rule
	if (any(object$method.type == c("WM", "HYFIS", "GFS.THRIFT", "ANFIS", "FIR.DM", "FS.HGD",
					"FRBCS.W", "FRBCS.CHI", "GFS.GCCL", "FH.GBML", "SLAVE", "GFS.LT.RS", "MANUAL"))){
		
		## get number of input variables and number of linguistic terms	
		num.varinput <- length(colnames.var) - 1
		num.labels <- object$num.labels
		
		## get names of linguistic terms
		if (!is.null(object$varinp.mf)){
			names.fTerms.inpvar <- colnames(object$varinp.mf)
		} else {names.fTerms.inpvar = NULL	}
		if (!is.null(object$varout.mf)){
			names.fTerms.outvar <- colnames(object$varout.mf)
		} else {names.fTerms.outvar = NULL }
		
		##	construct rule from rule.data.num (or rule in matrix format).
		if (!is.null(object$rule.data.num)){
			res <- generate.rule(object$rule.data.num, num.labels, 
			               names.fTerms.inpvar = names.fTerms.inpvar, names.fTerms.outvar = names.fTerms.outvar)
			rule <- res$rule
		} else {
			rule <- object$rule
		}
		
		## construct rule in IF ... THEN ... format
		new.rule <- matrix(nrow = nrow(rule), ncol = (ncol(rule) + 2 * (num.varinput + 1)))
		k <- 1
		
		for (j in 1 : num.varinput){				
			new.rule[, k] <- colnames.var[j] 
			new.rule[, k + 1] <- "is"
			new.rule[, k + 2] <- rule[, 2 * j - 1]
		
			if (j < num.varinput){
				## new.rule[, k + 3] <- "and"
				## A bug: when the boolean operator "or" (solved)
				new.rule[, k + 3] <- rule[, 2 * j]
			} else {
				new.rule[, k + 3] <- "THEN"
			}
			k <- k + 4
		}
		new.rule[, (ncol(new.rule) - 2)] <- colnames.var[num.varinput + 1] 
		new.rule[, (ncol(new.rule) - 1)] <- "is"
		new.rule[, ncol(new.rule)] <- rule[, ncol(rule)]
		
		## final checking for rules
		## TSK model
		if (object$type.model == "TSK"){
			rule <- new.rule[, 1 : (ncol(new.rule) - 3), drop = FALSE]
			if (object$method.type == "MANUAL"){
				rule <- cbind(rule, "THEN")
			}
		}
		
		## FRBCS model
		else if (object$type.model == "FRBCS" && any(object$method.type == c("FRBCS.W", "FRBCS.CHI", "MANUAL"))){
			rule <- new.rule[, 1 : (ncol(new.rule) - 3), drop = FALSE]
			rule <- cbind(rule, colnames.var[num.varinput + 1], "is", object$class)
		}
		
		## GFS.GCCL, FH.GBML, SLAVE methods
		else if (any(object$method.type == c("GFS.GCCL", "FH.GBML", "SLAVE"))){
			rule <- new.rule[, 1 : (ncol(new.rule) - 3), drop = FALSE]
			rule <- cbind(rule, colnames.var[num.varinput + 1], "is", object$rule.data.num[, ncol(object$rule.data.num)])
		}
		
		## Mamdani model
		else {
			rule <- new.rule
		}
		rule <- cbind("IF", rule)
	} else stop("It is not supported to create rule representation")
		
	return (rule)
} 

#' The purpose of this function is to generate a FRBS model from user-given 
#' input without a learning process.
#'
#' It can be used if rules have already been obtained manually, without employing the 
#' learning process. 
#' In the examples shown, we generate a fuzzy model using \code{frbs.gen} and generate the
#' fuzzy rule-based systems step by step manually. Additionally, the examples show several scenarios as follows.
#' \itemize{
#' \item Using \code{frbs.gen} for constructing the Mamdani model on a regression task. 
#' \item Using \code{frbs.gen} for constructing the Takagi Sugeno Kang model on a regression task.
#' \item Constructing the Mamdani model by executing internal functions such as \code{rulebase}, \code{fuzzifier},
#'       \code{inference}, and \code{defuzzifier} for the Mamdani model.
#' \item Using \code{frbs.gen} for constructing fuzzy rule-based classification systems (FRBCS) model.
#' }
#'
#' @title The frbs model generator
#' @param range.data a matrix (\eqn{2 \times n}) containing the range of the data, where \eqn{n} is the number of variables, and
#' first and second rows are the minimum and maximum values, respectively. 
#' @param num.fvalinput a matrix representing the number of linguistic terms of each input variables.
#' 
#' For example: \code{num.fvalinput <- matrix(c(3,2), nrow = 1)}
#' 
#' means that there are two variables where the first variable has three linguistic terms and the second one has two linguistic terms.
#' @param varinp.mf a matrix for constructing the shapes of the membership functions. See how to construct it in \code{\link{fuzzifier}}.
#' @param names.varinput a list containing names to the linguistic terms for input variables. See \code{\link{rulebase}}.
#' @param num.fvaloutput the number of linguistic terms of the output variable. This parameter is required for the Mamdani model only. 
#'
#' For example: \code{num.fvaloutput <- matrix(3, nrow = 1)}
#' 
#' means there are 3 linguistic terms for the output variable.
#' @param varout.mf a matrix for constructing the membership functions of the output variable. 
#' The form is the same as for the \code{varinp.mf} parameter. This parameter is required for the Mamdani model only. 
#' See \code{\link{fuzzifier}}.
#' @param names.varoutput a list giving names of the linguistic terms for the output variable. The form is the same as 
#' for the \code{names.varinput} parameter. This parameter is required for the Mamdani model only. See \code{\link{rulebase}}.
#' @param rule a list of fuzzy IF-THEN rules. There are some types of rule structures, for example: Mamdani, Takagi Sugeno Kang,
#' and fuzzy rule-based classification systems (FRBCS). If we use the Mamdani model then the consequent part is a linguistic term,
#' but if we use Takagi Sugeno Kang then we build a matrix representing linear equations in the consequent part.
#' e.g., "a1", "and", "b1, "->", "e1" means that 
#' "IF inputvar.1 is a1 and inputvar.2 is b1 THEN outputvar.1 is e1". 
#' Make sure that each rule has a "->" sign.
#' Furthermore, we are allowed to use linguistic hedges (e.g., "extremely", "slightly", etc), negation (i.e., "not"),
#' and the "dont_care" value representing degree of membership is always 1. 
#' For more detail, see \code{\link{rulebase}}. 
#' @param type.model the type of the model. There are three types available as follows. 
#' \itemize{
#' \item \code{MAMDANI} means we are using the Mamdani model. 
#' \item \code{TSK} means we are using the Takagi Sugeno Kang model.
#' \item \code{FRBCS} means we are using fuzzy rule-based classification systems (FRBCS).
#' }
#' @param type.defuz the type of the defuzzification method. It is used in the Mamdani model only. 
#'        See \code{\link{defuzzifier}}.
#' @param type.tnorm the type of the t-norm method. See \code{\link{inference}}.
#' @param type.snorm the type of the s-norm method. See \code{\link{inference}}.
#' @param func.tsk a matrix of parameters of the function on the consequent part using the Takagi Sugeno Kang model. 
#' This parameter must be defined when we are using Takagi Sugeno Kang. See \code{\link{rulebase}}.
#' @param colnames.var a list of names of input and output variables. 
#' @param type.implication.func a type of implication function. See \code{\link{WM}}.
#' @param name a name of the simulation.
#' @return The \code{\link{frbs-object}}. 
#' @examples 
#' 
#' #################################################
#' ## 1. The following codes show how to generate a fuzzy model 
#' ## using the frbs.gen function for regression tasks. 
#' ## The following are three scenarios:
#' ## 1a. Using the Mamdani model
#' ## 1b. Using the Takagi Sugeno Kang model
#' ## 1c. Using the Mamdani model and internal functions: fuzzifier, etc.
#' ## Note:
#' ## In the examples, let us consider four input variabels and one output variable.
#' ## Some variables could be shared together for other examples.  
#' #################################################
#'
#' ## Define shape and parameters of membership functions of input variables.
#' ## Please see the fuzzifier function to construct the matrix.
#' ## It can be seen that in this case we employ TRAPEZOID as the membership functions.
#' varinp.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
#'                       2, 0, 35, 75, NA, 3, 35, 75, 100, NA,
#'                       2, 0, 20, 40, NA, 1, 20, 50, 80, NA, 3, 60, 80, 100, NA,
#'                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
#'                       nrow = 5, byrow = FALSE)
#'
#' ## Define number of linguistic terms of the input variables.
#' ## Suppose, we have 3, 2, 3, and 3 numbers of linguistic terms 
#' ## for the first, second, third and fourth variables, respectively.
#' num.fvalinput <- matrix(c(3, 2, 3, 3), nrow=1)
#' 
#' ## Give the names of the linguistic terms of each input variables.
#' varinput.1 <- c("a1", "a2", "a3")
#' varinput.2 <- c("b1", "b2")
#' varinput.3 <- c("c1", "c2", "c3")
#' varinput.4 <- c("d1", "d2", "d3")
#' names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)
#'
#' ## Set interval of data.
#' range.data <- matrix(c(0,100, 0, 100, 0, 100, 0, 100, 0, 100), nrow=2)
#' 
#' ## Define inference parameters.
#' type.defuz <- "WAM"
#' type.tnorm <- "MIN"
#' type.snorm <- "MAX"
#' type.implication.func <- "ZADEH"
#'
#' ## Give the name of simulation.
#' name <- "Sim-0"
#'
#' ## Provide new data for testing. 
#' newdata <- matrix(c(15, 80, 85, 85, 45, 75, 78, 70), nrow = 2, byrow = TRUE)
#' ## the names of variables
#' colnames.var <- c("input1", "input2", "input3", "input4", "output1")
#'
#' ###################################################################
#' ## 1a. Using the Mamdani Model 
#' ####################################################################
#' ## Define number of linguistic terms of output variable.
#' ## In this case, we set the number of linguistic terms to 3.
#' num.fvaloutput <- matrix(c(3), nrow = 1)
#'
#' ## Give the names of the linguistic terms of the output variable.
#' varoutput.1 <- c("e1", "e2", "e3")
#' names.varoutput <- c(varoutput.1)
#'
#' ## Define the shapes and parameters of the membership functions of the output variables.
#' varout.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
#'                       nrow = 5, byrow = FALSE)
#'
#' ## Set type of model which is "MAMDANI" or "TSK" for Mamdani or 
#' ## Takagi Sugeno Kang models, respectively.
#' ## In this case, we choose the Mamdani model.
#' type.model <- "MAMDANI"
#' 
#' ## Define the fuzzy IF-THEN rules; In this case, we provide two scenarios using different operators:
#' rule.or <- matrix(c("a1", "or", "b1", "or", "c1", "or", "d1", "->", "e1",
#'                  "a2", "and", "b2", "and", "c2", "and", "d2", "->", "e2", 
#'                  "a3", "and", "b2", "and", "c2", "and", "d1", "->", "e3"), 
#'                  nrow = 3, byrow = TRUE) 
#'				  
#' ## Define the fuzzy IF-THEN rules; 
#' rule.and <- matrix(c("a1", "and", "b1", "and", "c1", "and", "d1", "->", "e1",
#'                  "a2", "and", "b2", "and", "c2", "and", "d2", "->", "e2", 
#'                  "a3", "and", "b2", "and", "c2", "and", "d1", "->", "e3"), 
#'                  nrow = 3, byrow = TRUE)  
#'
#' ## Generate a fuzzy model with frbs.gen.
#' object.or <- frbs.gen(range.data, num.fvalinput, names.varinput, 
#'                  num.fvaloutput, varout.mf, names.varoutput, rule.or, 
#'                  varinp.mf, type.model, type.defuz, type.tnorm, 
#'                  type.snorm, func.tsk = NULL, colnames.var, type.implication.func, name)
#'
#' object.and <- frbs.gen(range.data, num.fvalinput, names.varinput, 
#'                  num.fvaloutput, varout.mf, names.varoutput, rule.and, 
#'                  varinp.mf, type.model, type.defuz, type.tnorm, 
#'                  type.snorm, func.tsk = NULL, colnames.var, type.implication.func, name)
#' 
#' ## Plot the membership function.
#' plotMF(object.and)
#'
#' ## Predicting using new data.
#' res.or <- predict(object.or, newdata)$predicted.val
#' res.and <- predict(object.and, newdata)$predicted.val
#'
#' #####################################################################
#' ## 1b. Using the Takagi Sugeno Kang (TSK) Model 
#' #####################################################################
#' ## Define "TSK" for the Takagi Sugeno Kang model
#' type.model <- "TSK"
#' 
#' ## Define linear equations for consequent parts. 
#' ## The following command means that we have three equation related to the rules we have.
#' ## e.g., the first equation is 1*inputvar.1 + 1*inputvar.2 + 5*inputvar.3 + 2*inputvar.4 + 1, 
#' ## where inputvar.i is a value of the i-th input variable.
#' func.tsk <- matrix(c(1, 1, 5, 2, 1, 3, 1, 0.5, 0.1, 2, 1, 3, 2, 2, 2), 
#'             nrow = 3, byrow = TRUE)
#' 
#' ## Define the fuzzy IF-THEN rules; 
#' ## For TSK model, it isn't necessary to put linguistic term in consequent parts.
#' ## Make sure that each rule has a "->" sign. 
#' rule <- matrix(c("a1", "and", "b1", "and", "c1", "and", "d1", "->",
#'                  "a2", "and", "b2", "and", "c2", "and", "d2", "->",  
#'                  "a3", "and", "b2", "and", "c2", "and", "d1", "->"), 
#'                  nrow = 3, byrow = TRUE) 
#'				  
#' ## Generate a fuzzy model with frbs.gen.
#' ## It should be noted that for TSK model, we do not need to input: 
#' ## num.fvaloutput, varout.mf, names.varoutput, type.defuz.
#' object <- frbs.gen(range.data, num.fvalinput, names.varinput, 
#'              num.fvaloutput = NULL, varout.mf = NULL, names.varoutput = NULL, rule, 
#'				varinp.mf, type.model, type.defuz = NULL, type.tnorm, type.snorm, 
#'              func.tsk, colnames.var, type.implication.func, name)
#'				
#' ## Plot the membership function.
#' plotMF(object)
#'
#' ## Predicting using new data.
#' res <- predict(object, newdata)$predicted.val
#'
#' ######################
#' ## 1c. Using the same data as in the previous example, this example performs 
#' ## step by step of the generation of a fuzzy rule-based system
#' ######################
#' ## Using the Mamdani model.
#' type.model <- "MAMDANI"
#'
#' ## Construct rules.
#' rule <- matrix(c("a1", "and", "b1", "and", "c1", "and", "d1", "->", "e1",
#'                  "a2", "and", "b2", "and", "c2", "and", "d2", "->", "e2", 
#'                  "a3", "and", "b2", "and", "c2", "and", "d1", "->", "e3"), 
#'                  nrow = 3, byrow = TRUE) 
#' 
#' ## Check input data given by user.
#' rule <- rulebase(type.model, rule, func.tsk = NULL)
#' 
#' ## Fuzzification Module:
#' ## In this function, we convert crisp into linguistic values/terms
#' ## based on the data and the parameters of the membership function.
#' ## The output: a matrix representing the degree of the membership of the data
#' num.varinput <- ncol(num.fvalinput)
#' MF <- fuzzifier(newdata, num.varinput, num.fvalinput, varinp.mf)
#' 
#' ## Inference Module:
#' ## In this function, we will calculate the confidence factor on the antecedent for each rule
#' ## considering t-norm and s-norm.
#' miu.rule <- inference(MF, rule, names.varinput, type.tnorm, type.snorm)
#'
#' ## Defuzzification Module.
#' ## In this function, we calculate and convert the linguistic values back into crisp values. 
#' range.output <- range.data[, ncol(range.data), drop = FALSE]
#' result <- defuzzifier(newdata, rule, range.output, names.varoutput,
#'                   varout.mf, miu.rule, type.defuz, type.model, func.tsk = NULL)
#' 
#'
#' #################################################
#' ## 2. The following codes show how to generate a fuzzy model 
#' ## using the frbs.gen function for classification tasks using the Mamdani model. 
#' #################################################
#' ## define range of data.
#' ## Note. we only define range of input data. 
#' range.data.input <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1), nrow=2)
#' 
#' ## Define shape and parameters of membership functions of input variables.
#' ## Please see fuzzifier function to construct the matrix.
#' ## In this case, we are using TRIANGLE for membership functions.
#' varinp.mf <- matrix(c(1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA,
#'                       1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA,
#'                       1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA,
#'                       1, 0, 0, 0.5, NA, 1, 0, 0.5, 1, NA, 1, 0.5, 1, 1, NA),
#'                       nrow = 5, byrow = FALSE)
#'
#' ## Define number of linguistic terms of input variables.
#' ## Suppose, we have 3, 3, 3, and 3 numbers of linguistic terms 
#' ## for first up to fourth variables, respectively.
#' num.fvalinput <- matrix(c(3, 3, 3, 3), nrow=1)
#' 
#' ## Give the names of the linguistic terms of each input variable.
#' varinput.1 <- c("v.1_a.1", "v.1_a.2", "v.1_a.3")
#' varinput.2 <- c("v.2_a.1", "v.2_a.2", "v.2_a.3")
#' varinput.3 <- c("v.3_a.1", "v.3_a.2", "v.3_a.3")
#' varinput.4 <- c("v.4_a.1", "v.4_a.2", "v.4_a.3")
#' names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)
#'
#' ## Provide inference parameters.
#' type.tnorm <- "MIN"
#' type.snorm <- "MAX"
#' type.implication.func <- "ZADEH"
#' type.model <- "FRBCS"
#'
#' ## Give the name of simulation.
#' name <- "Sim-0"
#'
#' ## Provide new data for testing. 
#' newdata<- matrix(c(0.45, 0.5, 0.89, 0.44, 0.51, 0.99, 0.1, 0.98, 0.51,
#'                    0.56, 0.55, 0.5), nrow = 3, byrow = TRUE)
#'
#' ## the names of variables
#' colnames.var <- c("input1", "input2", "input3", "input4", "output1")
#' 
#' ## Construct rules.
#' ## It should be noted that on consequent parts we define categorical values instead of 
#' ## linguistic terms. 
#' rule <- matrix(
#'        c("v.1_a.2", "and", "v.2_a.2", "and", "v.3_a.3", "and", "v.4_a.2", "->", "3",
#'          "v.1_a.2", "and", "v.2_a.3", "and", "v.3_a.1", "and", "v.4_a.3", "->", "1",
#'          "v.1_a.2", "and", "v.2_a.2", "and", "v.3_a.2", "and", "v.4_a.2", "->", "2"), 
#'          nrow = 3, byrow = TRUE) 
#' 
#' ## Generate frbs object.
#' object <- frbs.gen(range.data = range.data.input, num.fvalinput, 
#'              names.varinput, num.fvaloutput = NULL, varout.mf = NULL, 
#'              names.varoutput = NULL, rule, varinp.mf, type.model, 
#'              type.defuz = NULL, type.tnorm, type.snorm, func.tsk = NULL, 
#'              colnames.var, type.implication.func, name)
#'				
#' ## Plot the shape of membership functions.
#' plotMF(object)
#'
#' ## Predicting using new data.
#' res <- predict(object, newdata)
#' 
#' ####################################################
#' ## 3. The following example shows how to convert 
#' ##    the frbs model into frbsPMML 
#' ####################################################
#' ## In this example, we are using the last object of FRBS.
#' ## Display frbsPMML in R
#' objPMML <- frbsPMML(object)
#' 
#' ## Write into a file with .frbsPMML extention
#' \dontrun{write.frbsPMML(objPMML, fileName="obj_frbsPMML")
#' 
#' ## Read the frbsPMML file into an R object of FRBS
#' obj <- read.frbsPMML("obj_frbsPMML.frbsPMML")}
#' @export
frbs.gen <- function (range.data, num.fvalinput, names.varinput, num.fvaloutput = NULL, varout.mf = NULL, names.varoutput = NULL, rule, varinp.mf,
                type.model = "MAMDANI", type.defuz = "WAM", type.tnorm = "MIN", type.snorm = "MAX", func.tsk = NULL, colnames.var = NULL, 
				type.implication.func = "ZADEH", name = "Sim-0"){

	## check linear eq. on consequent part and define num.labels 
	if (type.model == "TSK") {
		if (is.null(func.tsk))
			stop("Generating using this method, the consequent part should be given by linear equations as the Takagi Sugeno Model") 
		num.labels <- num.fvalinput
		num.varoutput <- 1
	}
	
	## check parameters required for Mamdani and define num.labels representing number of linguistic terms of input and output variabels
	else if (type.model == "MAMDANI"){
		if (is.null(num.fvaloutput) || is.null(varout.mf) || is.null(names.varoutput)) 
			stop("please complete the parameters needed for the Mamdani model")
		num.labels <- cbind(num.fvalinput, num.fvaloutput)
		var.mf <- cbind(varinp.mf, varout.mf)
		colnames(varout.mf) <- names.varoutput
		#get number of output variable
		num.varoutput <- ncol(num.fvaloutput)
	}
	
	## check parameters reguired for FRBCS
	else if (type.model == "FRBCS"){
		#class <- as.matrix(as.numeric(rule[, ncol(rule), drop = FALSE]))
		num.labels <- cbind(num.fvalinput, max(as.numeric(unique(rule[, ncol(rule)]))))
		num.varoutput <- 1
		class <- matrix(as.numeric(rule[, ncol(rule), drop = FALSE]), ncol = 1)
	}

	#get number of input variable
	num.varinput <- ncol(num.fvalinput)

	## keep colnames of training data into mod
	if (is.null(colnames.var)) {
		colnames.var <- paste("var", seq(1, ncol(range.data), sep = "."))
	}
	
	## define type of method as "MANUAL"
	method.type <- "MANUAL"
	colnames(varinp.mf) <- names.varinput		
	
	## define type of membership functions
	all.type.mf <- sort(unique(varinp.mf[1, ]))
	names.mf <- c("TRIANGLE", "TRAPEZOID", "TRAPEZOID", "TRAPEZOID", "GAUSSIAN", "SIGMOID", "BELL")
	type.mf.name <- matrix()
	for (i in 1 : length(all.type.mf))
		type.mf.name[i] <- names.mf[all.type.mf[i]]
	if (length(unique(type.mf.name)) > 1)
		type.mf <- "MIX"
	else
		type.mf <- c(unique(type.mf.name))
	
	mod <- list(num.labels = num.labels, varout.mf = varout.mf, rule = rule, varinp.mf = varinp.mf, range.data.ori = range.data,
			   type.model = type.model, type.tnorm = type.tnorm, type.implication.func = type.implication.func, type.mf = type.mf, 
			   type.defuz = type.defuz, type.snorm = type.snorm, func.tsk = func.tsk, 
			   method.type = method.type, name = name, colnames.var = colnames.var, class = class)

	## build into frbs class
	mod <- frbsObjectFactory(mod)
	
	## change rule format into IF ... THEN ...
	if (!is.null(mod$rule) && !is.null(mod$num.labels))
		mod$rule <- rep.rule(mod)
		
	return(mod)
}


#' This function is to transform from normalized data into real-valued data. 
#'
#' @title The data de-normalization
#' @param dt.norm a matrix (\eqn{n \times m}) of the normalized data.
#' @param range.data a matrix (\eqn{2 \times n}) containing the range of the data, where \eqn{n} is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
#' @param min.scale the minimum value within normalization.
#' @param max.scale the maximum value within normalization.
#' @seealso \code{\link{norm.data}}
#' @return the real-valued data
#' @export
denorm.data <- function(dt.norm, range.data, min.scale = 0, max.scale = 1){
	row.data <- nrow(dt.norm)
	col.data <- ncol(dt.norm)
	data.denorm <- matrix(nrow = row.data, ncol=col.data)
	
	# denormalize all data on each column 
	for (j in 1:col.data){
		min.data <- range.data[1, j]
		max.data <- range.data[2, j]

		for (i in 1:row.data){
			data.denorm[i, j] <- min.data + ((dt.norm[i, j] - min.scale)*(max.data - min.data))/ (max.scale - min.scale)
		}
	}
	
return(data.denorm)
}

#' This function is to transform from real-valued data into normalized data. 
#'
#' @title The data normalization
#' @param dt.ori a matrix (\eqn{n \times m}) of the original data.
#' @param range.data a matrix (\eqn{2 \times n}) containing the range of the data, where n is the number of variables, and
#' first and second rows are the minimum and maximum value, respectively. 
#' @param min.scale the minimum value within normalization.
#' @param max.scale the maximum value within normalization.
#' @seealso \code{\link{denorm.data}}
#' @return the normalized data
#' @export
norm.data <- function(dt.ori, range.data, min.scale = 0, max.scale = 1){
	row.data <- nrow(dt.ori)
	col.data <- ncol(dt.ori)
	data.norm <- matrix(nrow = row.data, ncol=col.data)
	
	# normalize all data on each column 
	for (j in 1:col.data){
		min.data <- range.data[1, j]
		max.data <- range.data[2, j]

		for (i in 1:row.data){
			data.norm[i, j] <- min.scale + (dt.ori[i, j] - min.data) * (max.scale - min.scale) / (max.data - min.data)
		}
	}
	
return(data.norm)
}


# This function is used to validate values of parameters. 
#
# @title The parameter validating function
# @param object a frbs object.
# @param newdata a new data.
# @return a valid object
# @export
validate.params <- function(object, newdata){
	if(!inherits(object, "frbs")) stop("not a legitimate frbs model")

	## in case, new data are not a matrix  
	if (class(newdata) != "matrix"){
		if (class(newdata) == "numeric")
			newdata <- matrix(newdata, nrow = 1)
		else
			newdata <- as.matrix(newdata)
	}  
	
	## Check range of new data
	for (i in 1 : ncol(newdata)){
		if (min(newdata[, i]) < object$range.data.ori[1, i])
			warning("There are your newdata which are out of the specified range")
		if (max(newdata[, i]) > object$range.data.ori[2, i])
			warning("There are your newdata which are out of the specified range")
	}
	
	## Check values of inference parameters
	object <- convert.params(object)
		
	## Check num.labels
	if (object$type.model == c("MAMDANI")){
		if (ncol(object$num.labels) != ncol(object$range.data.ori)) {
			stop("Please check your num.labels parameters")						
		}
	}
	else if (object$type.model == c("TSK")) {
		if (ncol(object$num.labels) != (ncol(object$range.data.ori) - 1)) {
			stop("Please check your num.labels parameters")						
		}
	}
	
	return(list(object = object, newdata = newdata))
}

# This function is used to convert numerical into string values of parameters. 
#
# @title The parameter converting function
# @param mod a frbs object.
# @return the normalized data
# @export
convert.params <- function(mod){
	## define a collection of values of parameters
	type.tnorm.str <- c("MIN", "HAMACHER", "YAGER", "PRODUCT", "BOUNDED")
	type.snorm.str <- c("MAX", "HAMACHER", "YAGER", "SUM", "BOUNDED")
	type.defuz.str <- c("WAM", "FIRST.MAX", "LAST.MAX", "MEAN.MAX", "COG")
	type.mf.str <- c("TRIANGLE", "TRAPEZOID", "GAUSSIAN", "SIGMOID", "BELL")
	type.model.str <- c("MAMDANI", "TSK", "FRBCS", "CLUSTERING", "APPROXIMATE")
	
	if (any(mod$type.model == c("MAMDANI", "TSK", "FRBCS", "APPROXIMATE"))){
		## check type.tnorm
		if (is.null(mod$type.tnorm)) {
			warning("type.tnorm is not defined, it will be assigned to 'MIN' ")
			mod$type.tnorm <- "MIN"
		}
		else if (class(mod$type.tnorm) == "numeric")
			mod$type.tnorm <- type.tnorm.str[mod$type.tnorm]
		else if (class(mod$type.tnorm) == "character")
			mod$type.tnorm <- toupper(mod$type.tnorm)
		
		## check type.snorm
		if (is.null(mod$type.snorm)){
			warning("type.snorm is not defined, it will be assigned to 'MAX' ")
			mod$type.snorm <- "MAX"
		}
		else if (class(mod$type.snorm) == "numeric")
			mod$type.snorm <- type.snorm.str[mod$type.snorm]
		else if (class(mod$type.snorm) == "character")
			mod$type.snorm <- toupper(mod$type.snorm)
		
		## check type.implication.func
		mod$type.implication.func <- toupper(mod$type.implication.func)
	}
	
		
	## check type.defuz
	if (!is.null(mod$type.defuz) && class(mod$type.defuz) == "numeric")
		mod$type.defuz <- type.defuz.str[mod$type.defuz]
	else if (!is.null(mod$type.defuz) && class(mod$type.defuz) == "character")
		mod$type.defuz <- toupper(mod$type.defuz)
		
	## check type.mf
	if (!is.null(mod$type.mf) && class(mod$type.mf) == "numeric")
		mod$type.mf <- type.mf.str[mod$type.mf]
	else if (!is.null(mod$type.mf) && class(mod$type.mf) == "character")
		mod$type.mf <- toupper(mod$type.mf)
	
	## check type.model
	if (is.null(mod$type.model))
		stop("please define type of model")
	else if (class(mod$type.model) == "numeric")
		mod$type.model <- type.model.str[mod$type.model]
	else if (class(mod$type.model) == "character")
		mod$type.model <- toupper(mod$type.model)
	
	## check type.defuz and MAMDANI
	if (mod$type.model == "MAMDANI" && is.null(mod$type.defuz)){
		warning("type.defuz is not defined, it will be assigned to 'WAM' ")
		mod$type.defuz <- "WAM"
	}
			
	return (mod)
}

# This function is to generate rule into string form
#
# @param rule.data.num a matrix of rules in integer form
# @param num.labels a matrix of the number of linguistic terms
# @param names.fTerms.inpvar a matrix of names of linguistic terms of input variables
# @param names.fTerms.outvar a matrix of names of linguistic terms of the output variable
generate.rule <- function(rule.data.num, num.labels, names.fTerms.inpvar = NULL, names.fTerms.outvar = NULL){
	
	## build the names of linguistic values	
	if (is.null(names.fTerms.inpvar) || is.null(names.fTerms.outvar)){
		fuzzyTerm <- create.fuzzyTerm(classification = FALSE, num.labels)
		names.inp.var <- fuzzyTerm$names.fvalinput
		names.out.var <- fuzzyTerm$names.fvaloutput
		names.variable <- c(names.inp.var, names.out.var)
	}
	else {
		names.inp.var <- names.fTerms.inpvar
		names.out.var <- names.fTerms.outvar
		names.variable <- c(names.inp.var, names.out.var)
	}
		
	## build the rule into list of string
	rule <- matrix(nrow = nrow(rule.data.num), ncol = 2 * ncol(rule.data.num) - 1)
	
	for (i in 1 : nrow(rule.data.num)){
		k <- 0
		for (j in 1 : ncol(rule.data.num)){
			k <- k + 1	
			if (j == ncol(rule.data.num) - 1){
				if (rule.data.num[i, j] == 0){
					rule[i, k] <- c("dont_care")
				} else {
					rule[i, k] <- c(names.variable[rule.data.num[i, j]])
				}
				rule[i, k + 1] <- c("->")
				k <- k + 1
			}
			else if (j == ncol(rule.data.num)){
				if (rule.data.num[i, j] == 0){
					rule[i, k] <- c("dont_care")
				} else {
					rule[i, k] <- c(names.variable[rule.data.num[i, j]])
				}
			}
			else{
				if (rule.data.num[i, j] == 0){
					rule[i, k] <- c("dont_care")
				} else {
					rule[i, k] <- c(names.variable[rule.data.num[i, j]])
				}
				rule[i, k + 1] <- c("and")
				k <- k + 1
			}
		}
	
	}
res <- list(rule = rule, names.varinput = names.inp.var, names.varoutput = names.out.var)
return (res)
}
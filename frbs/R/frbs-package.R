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
#' Fuzzy rule-based systems (FRBSs) are based on the fuzzy concept 
#' proposed by Zadeh in 1965, which represents the reasoning of human experts in production 
#' rules (a set of IF-THEN rules) to handle real-life problems from domains 
#' such as control, prediction and inference, data mining, bioinformatics data processing, 
#' robotics, and speech recognition. FRBSs are also known as fuzzy inference systems and 
#' fuzzy models. When applied to specific tasks, they may also be known under specific names 
#' such as fuzzy associative memories or fuzzy controllers.
#' In this package, we consider systems with multi-inputs and single-output (MISO), 
#' with real-valued data.
#' 
#' FRBSs are a competitive alternative to other classic models and algorithms in order to 
#' solve classification and regression problems. Generally, 
#' an FRBS consists of four functional parts: 
#' \itemize{
#' \item a fuzzification interface which transforms the crisp inputs into degrees 
#'       of membership functions of the linguistic term of each variable. 
#'       See \code{\link{fuzzifier}}.
#' \item a knowledge base consisting of a database (DB) and a rulebase (RB). While the database includes 
#'       the fuzzy set definitions, the rulebase contains fuzzy IF-THEN rules. 
#'       We will represent the knowledge as a set of rules. Each one has the following structure.
#'
#' \code{IF premise (antecedent) THEN conclusion (consequent)}
#' 
#'       See \code{\link{rulebase}}.
#' \item an inference engine which performs the inference operations on the fuzzy IF-THEN rules. 
#'       There are two kinds of inference for fuzzy systems based on linguistic rules: 
#'       The Mamdani and the Takagi Sugeno Kang model. See \code{\link{inference}}.
#' \item a defuzzification process to obtain the crisp values from linguistic values. There are several methods for 
#'       defuzzification such as the weighted average, centroid, etc. 
#'       See \code{\link{defuzzifier}}.
#' }
#' 
#' Since it may be difficult to obtain information from human experts in the form required,
#' an alternative and effective way to acquire the knowledge is to generate 
#' the fuzzy IF-THEN rules automatically from the numerical training data. 
#' In general, when modeling an FRBS, there are two important processes which should be conducted, 
#' namely structure identification and parameter estimation. 
#' Structure identification is a process to find appropriate fuzzy IF-THEN rules
#' and to determine the overall number of rules. 
#' Parameter estimation is applied to tune parameters of membership functions.
#' Many approaches have been proposed 
#' in order to perform this modeling such as a table-lookup scheme, heuristic procedures,
#' neuro-fuzzy techniques, clustering methods, genetic algorithms, least squares methods,
#' gradient descent, etc. In this package, the following approaches to generate 
#' fuzzy IF-THEN rules have been implemented: 
#' \enumerate{
#' \item FRBS based on space partition
#' \itemize{
#' \item Wang and Mendel's technique (\code{WM}): It is used to solve regression tasks. See \code{\link{WM}}.
#' \item Chi's technique (\code{FRBCS.CHI}): It is used to solve classification tasks. See \code{\link{FRBCS.CHI}}.
#' \item Ishibuchi's technique using weight factor (\code{FRBCS.W}): It is used to solve classification tasks. See \code{\link{FRBCS.W}}.
#' }
#' \item FRBS based on neural networks
#' \itemize{
#' \item The adaptive-network-based fuzzy inference system (\code{ANFIS}): 
#'              It is used to solve regression tasks. See \code{\link{ANFIS}}.
#' \item The hybrid neural fuzzy inference system (\code{HYFIS}): It is used to solve regression tasks. See \code{\link{HyFIS}}.
#' }
#' \item FRBS based on clustering approach
#' \itemize{
#' \item The subtractive clustering and fuzzy c-means (\code{SBC}): It is used to solve regression tasks. See \code{\link{SBC}}.
#' \item The dynamic evolving neural-fuzzy inference system (\code{DENFIS}): 
#'            It is used to solve regression tasks. See \code{\link{DENFIS}}.
#' }
#' \item FRBS based on genetic algorithms
#' \itemize{
#' \item The Thrift's method (\code{GFS.THRIFT}): It is used to solve regression tasks. See \code{\link{GFS.Thrift}}.
#' \item The Genetic fuzzy systems for fuzzy rule learning based on the MOGUL methodology (\code{GFS.FR.MOGUL}):
#'       It is used to solve regression tasks. See \code{\link{GFS.FR.MOGUL}}.
#' \item The Ishibuchi's method based on genetic cooperative-competitive learning (\code{GFS.GCCL}):
#'       It is used to solve classification tasks. See \code{\link{GFS.GCCL}}.
#' \item The Ishibuchi's method based on hybridization of genetic cooperative-competitive learning (GCCL) and Pittsburgh (\code{FH.GBML}):
#'       It is used to solve classification tasks. See \code{\link{FH.GBML}}.
#' \item The structural learning algorithm on vague environtment (\code{SLAVE}):
#'       It is used to solve classification tasks. See \code{\link{SLAVE}}.
#' \item The genetic for lateral tuning and rule selection of linguistic fuzzy system (\code{GFS.LT.RS}): 
#'       It is used to solve regression tasks. See \code{\link{GFS.LT.RS}}.
#' }
#' \item FRBS based on the gradient descent method
#' \itemize{
#' \item The FRBS using heuristics and gradient descent method (\code{FS.HGD}): 
#'        It is used to solve regression tasks. See \code{\link{FS.HGD}}
#' \item The fuzzy inference rules by descent method (\code{FIR.DM}): 
#'        It is used to solve regression tasks. See \code{\link{FIR.DM}}
#' }
#' }
#' The functions documented in the manual for the single methods are all called internally 
#' by \code{\link{frbs.learn}}, which is the central function of the package. 
#' However, in the documentation of each of the internal learning functions, 
#' we give some theoretical background and references to the original literature.
#' 
#'
#' \bold{Usage of the package:}
#' 
#' First of all, if you have problems using the package, find a bug, or have suggestions, 
#' please contact the package maintainer by email, instead of writing to the general R lists 
#' or to other internet forums and mailing lists.
#' 
#' The main functions of the package are the following:
#' \itemize{
#' \item The function \code{\link{frbs.learn}} allows to generate the model by 
#'       creating fuzzy IF-THEN rules or cluster centers from training data. 
#'       In other words, users just need to call this function to generate an FRBS model from 
#'       training data. The different algorithms mentioned above are all accessible through this function. 
#'       The outcome of the function is an \code{\link{frbs-object}}. 
#' \item Even though the main purpose of this package is to generate 
#'       the FRBS models from training data automatically, we provide the function \code{\link{frbs.gen}}, 
#'       which can be used to build a model manually without using a learning method.
#'       Moreover, we provide the following features: linguistic hedges, the "and" and "or" operators, 
#'       and the "dont_care" value for representing the degree of 1. The higher degree of interpretability 
#'       can also be achieved by using the "dont_care" value. If we want to define various length of rules
#'       in the rulebase, we can also define the "dont_care" value. See \code{\link{rulebase}}.    
#' \item The purpose of the function \code{\link{predict}} is to obtain predicted values 
#'       according to the testing data and the model (analogous to the \code{predict} function 
#'       that is implemented in many other R packages). 
#' \item There exist functions \code{\link{summary.frbs}} and \code{\link{plotMF}} to 
#'       show a summary about an \code{\link{frbs-object}}, and to plot the shapes of 
#'       the membership functions.
#' \item Exporting an FRBS model to the frbsPMML format can be done by executing \code{\link{frbsPMML}} and \code{\link{write.frbsPMML}}. 
#'       The frbsPMML format is a universal framework adopted from the Predictive Model Markup Language (PMML) format. Then, 
#'       in order to consume/import the frbsPMML format to an FRBS model, we call \code{\link{read.frbsPMML}}.  
#' }
#' 
#' To get started with the package, the user can have a look at some examples included in
#' the documentation of the functions \code{\link{frbs.learn}} and \code{\link{frbs.gen}} for generating models and
#' \code{\link{predict}} for the prediction phase.
#'  
#' Also, there are many demos that ship with the package. To get a list of them, type:
#' 
#' \code{demo()}
#' 
#' Then, to start a demo, type \code{demo(<demo_name_here>)}. All the demos are present as 
#' R scripts in the package sources in the \code{"demo"} subdirectory. Note that
#' some of them may take quite a long time which depends on specification hardwares.
#' 
#' Currently, there are the following demos available:
#' 
#' Regression using the Gas Furnance dataset: 
#' 
#' \code{demo(WM.GasFur)},
#' \code{demo(SBC.GasFur)},
#' \code{demo(ANFIS.GasFur)},
#' 
#' \code{demo(FS.HGD.GasFur)}, 
#' \code{demo(DENFIS.GasFur)},
#' \code{demo(HyFIS.GasFur)},
#' 
#' \code{demo(FIR.DM.GasFur)},
#' \code{demo(GFS.FR.MOGUL.GasFur)},
#'
#' \code{demo(GFS.THRIFT.GasFur)},
#' \code{demo(GFS.LT.RS.GasFur)}.
#'
#' Regression using the Mackey-Glass dataset:
#' 
#' \code{demo(WM.MG1000)},
#' \code{demo(SBC.MG1000)},
#' \code{demo(ANFIS.MG1000)},
#' 
#' \code{demo(FS.HGD.MG1000)}, 
#' \code{demo(DENFIS.MG1000)},
#' \code{demo(HyFIS.MG1000)},
#' 
#' \code{demo(GFS.THRIFT.MG1000)}, 
#' \code{demo(FIR.DM.MG1000)},
#'
#' \code{demo(GFS.FR.MOGUL.MG1000)},
#' \code{demo(GFS.LT.RS.MG1000)}.
#'
#' Classification using the Iris dataset:
#' 
#' \code{demo(FRBCS.W.Iris)},
#' \code{demo(FRBCS.CHI.Iris)},
#' \code{demo(GFS.GCCL.Iris)},
#'
#' \code{demo(FH.GBML.Iris)},
#' \code{demo(SLAVE.Iris)}.
#'
#' Generating FRBS model without learning process:
#' 
#' \code{demo(FRBS.Mamdani.Manual)},
#' \code{demo(FRBS.TSK.Manual)},
#' \code{demo(FRBS.Manual)}.
#'
#' Exporting/importing to/from frbsPMML:
#'
#' \code{demo(WM.GasFur.PMML)},
#' \code{demo(ANFIS.GasFur.PMML)},
#' \code{demo(GFS.GCCL.Iris.PMML)}.
#'
#' The Gas Furnance data and Mackey-Glass data are included in the package, 
#' please see \code{\link{frbsData}}. The Iris data is the standard Iris dataset that
#' ships with R.
#' 
#' Also have a look at the package webpage \url{http://sci2s.ugr.es/dicits/software/FRBS}, 
#' where we provide a more extensive introduction as well as additional explanations of 
#' the procedures.
#' 
#' @name frbs-package
#' @aliases frbs
#' @docType package
#' @title Getting started with the frbs package
#' @seealso \code{\link{frbs.learn}}, \code{\link{frbs.gen}}, \code{\link{frbsPMML}}, and \code{\link{predict}}. 
#' @references 
#' A. Guazzelli, M. Zeller, W.C. Lin, and G. Williams., 
#' "pmml: An open standard for sharing models", The R Journal, Vol. 1, No. 1, pp. 60-65 (2009).
#'
#' C.C. Lee, "Fuzzy Logic in control systems: Fuzzy logic controller part I", 
#' IEEE Trans. Syst., Man, Cybern., 
#' vol. 20, no. 2, pp. 404 - 418 (1990).
#' 
#' C.C. Lee, "Fuzzy Logic in control systems: Fuzzy logic controller part II",
#' IEEE Trans. Syst., Man, Cybern.,
#' vol. 20, no. 2, pp. 419 - 435 (1990).
#'
#' E.H. Mamdani and S. Assilian, "An experiment in linguistic synthesis with 
#' a fuzzy logic controller," International Journal of Man Machine Studies, vol. 7, no. 1, 
#' pp. 1 - 13 (1975).
#'
#' W. Pedrycz, "Fuzzy Control and Fuzzy Systems," New York: Wiley (1989).
#' 
#' L.S. Riza, C. Bergmeir, F. Herrera, and J.M. Benitez, 
#' "frbs: Fuzzy Rule-Based Systems for Classification and Regression in R,"
#' Journal of Statistical Software, vol. 65, no. 6, pp. 1 - 30 (2015).
#'
#' M. Sugeno and G.T. Kang, "Structure identification of fuzzy model," 
#' Fuzzy Sets Syst., vol. 28, pp. 15 - 33 (1988).
#'
#' T. Takagi and M. Sugeno, "Fuzzy identification of systems and its application to 
#' modelling and control", IEEE Transactions on Systems, Man and Cybernetics, vol. 15, no. 1, 
#' pp. 116 - 132 (1985).
#'
#' L.A. Zadeh, "Fuzzy sets", Information and Control, vol. 8, pp. 338 - 353 (1965).
#'
# @keywords package fuzzy rule based systems inference frbs regression classification
# @encoding UTF-8
# @encoding Latin-1
#' @author Lala Septem Riza \email{lala.s.riza@@decsai.ugr.es},
#' 
#' Christoph Bergmeir \email{c.bergmeir@@decsai.ugr.es}, 
#' 
#' Francisco Herrera \email{herrera@@decsai.ugr.es}, 
#' 
#' and Jose Manuel Benitez \email{j.m.benitez@@decsai.ugr.es}
#' 
#' DiCITS Lab, SCI2S group, DECSAI, University of Granada.
#' 
#' \url{http://sci2s.ugr.es/dicits/}, \url{http://sci2s.ugr.es}
#' 
#' @examples
#' ##################################
#' ## I. Regression Problem
#' ## In this example, we are using the gas furnace dataset that 
#' ## contains two input and one output variables.
#' ##################################
#'
#' ## Input data: Using the Gas Furnace dataset
#' ## then split the data to be training and testing datasets
#' data(frbsData)
#' data.train <- frbsData$GasFurnance.dt[1 : 204, ]
#' data.tst <- frbsData$GasFurnance.dt[205 : 292, 1 : 2]
#' real.val <- matrix(frbsData$GasFurnance.dt[205 : 292, 3], ncol = 1)
#'
#' ## Define interval of data
#' range.data <-apply(data.train, 2, range)
#'
#' ## Set the method and its parameters,
#' ## for example, we use Wang and Mendel's algorithm
#' method.type <- "WM"
#' control <- list(num.labels = 15, type.mf = "GAUSSIAN", type.defuz = "WAM", 
#'            type.tnorm = "MIN", type.snorm = "MAX", type.implication.func = "ZADEH",
#'                name = "sim-0") 
#'
#' ## Learning step: Generate an FRBS model
#' object.reg <- frbs.learn(data.train, range.data, method.type, control)
#'
#' ## Predicting step: Predict for newdata
#' res.test <- predict(object.reg, data.tst)
#' 
#' ## Display the FRBS model
#' summary(object.reg)
#' 
#' ## Plot the membership functions
#' plotMF(object.reg)
#'
#' ##################################
#' ## II. Classification Problem
#' ## In this example, we are using the iris dataset that 
#' ## contains four input and one output variables.
#' ##################################
#'
#' ## Input data: Using the Iris dataset
#' data(iris)
#' set.seed(2)
#'
#' ## Shuffle the data
#' ## then split the data to be training and testing datasets
#' irisShuffled <- iris[sample(nrow(iris)), ]
#' irisShuffled[, 5] <- unclass(irisShuffled[, 5])
#' tra.iris <- irisShuffled[1 : 105, ]
#' tst.iris <- irisShuffled[106 : nrow(irisShuffled), 1 : 4]
#' real.iris <- matrix(irisShuffled[106 : nrow(irisShuffled), 5], ncol = 1)
#'
#' ## Define range of input data. Note that it is only for the input variables.
#' range.data.input <- apply(iris[, -ncol(iris)], 2, range)
#' 
#' ## Set the method and its parameters. In this case we use FRBCS.W algorithm
#' method.type <- "FRBCS.W"
#' control <- list(num.labels = 7, type.mf = "GAUSSIAN", type.tnorm = "MIN", 
#'                type.snorm = "MAX", type.implication.func = "ZADEH")  
#'
#' ## Learning step: Generate fuzzy model
#' object.cls <- frbs.learn(tra.iris, range.data.input, method.type, control)
#'
#' ## Predicting step: Predict newdata
#' res.test <- predict(object.cls, tst.iris)
#' 
#' ## Display the FRBS model
#' summary(object.cls)
#' 
#' ## Plot the membership functions
#' plotMF(object.cls)
#'
#' #################################################
#' ## III. Constructing an FRBS model from human expert. 
#' ## In this example, we only consider the Mamdani model for regression. However, 
#' ## other models can be done in the same way.
#' ## Note:
#' ## In the examples, let us consider four input and one output variables.
#' #################################################
#'
#' ## Define a matrix representing shape and parameters of membership functions of input variables.
#' ## The matrix has 5 rows where the first row represent the type of the membership function whereas
#' ## others are values of its parameters.
#' ## Detailed explanation can be seen in the fuzzifier function to construct the matrix.
#' varinp.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA,
#'                       2, 0, 35, 75, NA, 3, 35, 75, 100, NA,
#'                       2, 0, 20, 40, NA, 1, 20, 50, 80, NA, 3, 60, 80, 100, NA,
#'                       2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
#'                       nrow = 5, byrow = FALSE)
#'
#' ## Define number of linguistic terms of input variables.
#' ## Suppose, we have 3, 2, 3, and 3 numbers of linguistic terms 
#' ## for the first, second, third and fourth variables, respectively.
#' num.fvalinput <- matrix(c(3, 2, 3, 3), nrow=1)
#' 
#' ## Give the names of the linguistic terms of each input variables.
#' varinput.1 <- c("low", "medium", "high")
#' varinput.2 <- c("yes", "no")
#' varinput.3 <- c("bad", "neutral", "good")
#' varinput.4 <- c("low", "medium", "high")
#' names.varinput <- c(varinput.1, varinput.2, varinput.3, varinput.4)
#'
#' ## Set interval of data.
#' range.data <- matrix(c(0, 100, 0, 100, 0, 100, 0, 100, 0, 100), nrow = 2)
#' 
#' ## Define inference parameters.
#' ## Detailed information about values can be seen in the inference function.
#' type.defuz <- "WAM"
#' type.tnorm <- "MIN"
#' type.snorm <- "MAX"
#' type.implication.func <- "ZADEH"
#'
#' ## Give the name of simulation.
#' name <- "Sim-0"
#'
#' ## Provide new data for testing. 
#' newdata<- matrix(c(25, 40, 35, 15, 45, 75, 78, 70), nrow = 2, byrow = TRUE)
#' ## the names of variables
#' colnames.var <- c("input1", "input2", "input3", "input4", "output1")
#'
#' ## Define number of linguistic terms of output variable.
#' ## In this case, we set the number of linguistic terms to 3.
#' num.fvaloutput <- matrix(c(3), nrow = 1)
#'
#' ## Give the names of the linguistic terms of the output variable.
#' varoutput.1 <- c("bad", "neutral", "good")
#' names.varoutput <- c(varoutput.1)
#'
#' ## Define the shapes and parameters of the membership functions of the output variables.
#' varout.mf <- matrix(c(2, 0, 20, 40, NA, 4, 20, 40, 60, 80, 3, 60, 80, 100, NA),
#'                       nrow = 5, byrow = FALSE)
#'
#' ## Set type of model which is "MAMDANI".
#' type.model <- "MAMDANI"
#' 
#' ## Define the fuzzy IF-THEN rules; 
#' ## In this example we are using the Mamdani model 
#' ## Note: e.g.,
#' ## "a1", "and", "b1, "->", "e1" means that 
#' ## "IF inputvar.1 is a1 and inputvar.2 is b1 THEN outputvar.1 is e1" 
#' ## Make sure that each rule has a "->" sign. 
#' rule <- matrix(
#'   c("low", "and", "yes", "and", "bad", "and", "low", "->", "bad",
#'     "medium", "and", "no", "and", "neutral", "and", "medium", "->", "neutral", 
#'     "high", "and", "no", "and", "neutral", "and", "low", "->", "good"), 
#'     nrow = 3, byrow = TRUE) 
#'
#' ## Generate a fuzzy model with frbs.gen.
#' object <- frbs.gen(range.data, num.fvalinput, names.varinput, 
#'                  num.fvaloutput, varout.mf, names.varoutput, rule, 
#'                  varinp.mf, type.model, type.defuz, type.tnorm, 
#'                  type.snorm, func.tsk = NULL, colnames.var, type.implication.func, name)
#' 
#' ## Plot the membership function.
#' plotMF(object)
#'
#' ## Predicting using new data.
#' res <- predict(object, newdata)$predicted.val
#' 
#' #################################################
#' ## IV. Specifying an FRBS model in the frbsPMML format.
#' ## other examples can be seen in the frbsPMML function.
#' #################################################
#' ## Input data
#' data(frbsData)
#' data.train <- frbsData$GasFurnance.dt[1 : 204, ]
#' data.fit <- data.train[, 1 : 2]
#' data.tst <- frbsData$GasFurnance.dt[205 : 292, 1 : 2]
#' real.val <- matrix(frbsData$GasFurnance.dt[205 : 292, 3], ncol = 1)
#' range.data<-matrix(c(-2.716, 2.834, 45.6, 60.5, 45.6, 60.5), nrow = 2)
#' 
#' ## Set the method and its parameters
#' method.type <- "WM"
#' control <- list(num.labels = 3, type.mf = "GAUSSIAN", type.defuz = "WAM", 
#'                 type.tnorm = "MIN", type.snorm = "MAX", 
#'                 type.implication.func = "ZADEH", name="sim-0") 
#' 
#' ## Generate fuzzy model
#' object <- frbs.learn(data.train, range.data, method.type, control)
#' 
#' ## 2. Constructing the frbsPMML format
#' frbsPMML(object)
NULL



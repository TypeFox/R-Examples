#' The spsann Package
#' 
#' Optimization of sample configurations using spatial simulated annealing
#' 
#' @section Introduction: 
#'
#' \pkg{spsann} is a package for the optimization of spatial sample configurations using spatial simulated 
#' annealing. It includes multiple objective functions to optimize spatial sample configurations for various 
#' purposes such as variogram estimation, spatial trend estimation, and spatial interpolation. Most of the 
#' objective functions were designed to optimize spatial sample configurations when a) multiple spatial 
#' variables must be modelled, b) we know very little about the model of spatial variation of those variables, 
#' and c) sampling is limited to a single phase.
#' 
#' Spatial simulated annealing is a well known method with widespread use to solve combinatorial optimization 
#' problems in the environmental sciences. This is mainly due to its robustness against local optima and 
#' easiness to implement. In short, the algorithm consists of randomly changing the spatial location of a 
#' candidate sampling point at a time and evaluating if the resulting spatial sample configuration is 
#' \emph{better} than the previous one with regard to the chosen quality criterion, i.e. an objective 
#' function. Sometimes a \emph{worse} spatial sample configuration is accepted so that the algorithm is able 
#' to scape from local optima solutions, i.e. those spatial sample configurations that are too good and appear 
#' to early in the optimization to be true. The chance of accepting a \emph{worse} spatial sample 
#' configuration reduces as the optimization proceeds so that we can get very close to the \emph{optimum} 
#' spatial sample configuration.
#' 
#' \pkg{spsann} also combines multiple objective functions so that spatial sample configurations can be
#' optimized regarding more than one modelling objective. Combining multiple objective functions gives rise to
#' a multi-objective combinatorial optimization problem (MOCOP). A MOCOP usually has multiple possible 
#' solutions. \pkg{spsann} finds a single solution by aggregating the objective functions using the 
#' weighted-sum method. With this method the relative importance of every objective function can be 
#' specified at the beginning of the optimization so that their relative influence on the resulting optimized
#' spatial sample configuration can be different. But this requires the objective functions first to be scaled 
#' to the same approximate range of values. The upper-lower bound approach is used for that end. In this 
#' approach, every objective function is scaled using as reference the respective minimum and maximum  
#' attainable objective function values, also known as the Pareto minimum and maximum.
#' 
#' @section Package Structure:
#'
#' \pkg{spsann} has a very simple structure composed of three families of functions. The first is the family of
#' \code{optim} functions. These are the functions that include the spatial simulated annealing algorithm, 
#' that is, the functions that perform the optimization regarding the chosen quality criterion (objective 
#' function). Every \code{optim} function is named after the objective function used as quality criterion. For
#' example, the quality criterion used by \code{\link[spsann]{optimMSSD}} is the \emph{mean squared shortest 
#' distance} (MSSD) between sample and prediction points. As the example shows, the name of the \code{optim}
#' functions is composed of the string \code{'optim'} followed by a suffix that indicates the respective 
#' objective  function. In the example this is \code{'MSSD'}.
#' 
#' There currently are nine function in the \code{optim} family: \code{\link[spsann]{optimACDC}}, 
#' \code{\link[spsann]{optimCLHS}}, \code{\link[spsann]{optimCORR}}, \code{\link[spsann]{optimDIST}}, 
#' \code{\link[spsann]{optimMSSD}}, \code{\link[spsann]{optimMKV}}, \code{\link[spsann]{optimPPL}}, 
#' \code{\link[spsann]{optimSPAN}}, and \code{\link[spsann]{optimUSER}}. The latter is a general purpose
#' function that enables to user to define his/her own objective function and plug it in the spatial simulated 
#' annealing algorithm.
#'
#' The second family of functions is the \code{obj} family. This family of functions is used to return the 
#' current objective function value of a spatial sample configuration. Like the family of \code{optim} 
#' functions, the name of the \code{obj} functions is composed of the string \code{'obj'} plus a suffix that
#' indicates the objective function being used. For example, \code{\link[spsann]{objMSSD}} computes the value
#' of the mean squared shortest distance between sample and prediction points of any spatial sample 
#' configuration. Accordingly, there is a \code{obj} function for every \code{optim} function, except for 
#' \code{\link[spsann]{optimUSER}}. A ninth \code{obj} function, \code{\link[spsann]{objSPSANN}}, returns the 
#' objective function value at any point of the optimization, irrespective of the objective function used.
#' 
#' The third family of functions implemented in \pkg{spsann} corresponds to a set of auxiliary functions. 
#' These auxiliary functions can be used for several purposes, such as organizing the information needed to
#' feed an \code{optim} function, retrieving information from an object of class 
#' \code{OptimizedSampleConfiguration}, i.e. an object containing an optimized sample configuration, 
#' generating plots of the spatial distribution an optimized sample configuration, and so on. These functions
#' are named after the purpose for which they have been designed. For example: \code{\link[spsann]{countPPL}}, 
#' \code{\link[spsann]{minmaxPareto}}, \code{\link[spsann]{scheduleSPSANN}}, \code{\link[spsann]{spJitter}},
#' and \code{\link[spsann]{plot}}.
#' 
#' Despite \pkg{spsann} functions are classified into three general family of functions defined according to 
#' the purpose for which they were designed, the documentation is constructed with regard to the respective 
#' objective functions. For example, every \pkg{spsann} function that uses as quality criterion the MSSD is 
#' documented in the same documentation page. The exception are the auxiliary functions, that generally are
#' documented separately.
#' 
#' @section Support:
#' \pkg{spsann} was initially developed as part of the PhD research project entitled \sQuote{Contribution to 
#' the Construction of Models for Predicting Soil Properties}, developed by Alessandro Samuel-Rosa under the 
#' supervision of Lúcia Helena Cunha dos Anjos \email{lanjos@@ufrrj.br} (Universidade Federal Rural do Rio de 
#' Janeiro, Brazil), Gustavo de Mattos Vasques \email{gustavo.vasques@@embrapa.br} (Embrapa Solos, Brazil), 
#' and Gerard B. M. Heuvelink \email{gerard.heuvelink@@wur.nl} (ISRIC -- World Soil Information, the 
#' Netherlands). The project was supported from March/2012 to February/2016 by the CAPES Foundation, Ministry 
#' of Education of Brazil, and the CNPq Foundation, Ministry of Science and Technology of Brazil.
#' 
#' @section Contributors:
#' Some of the solutions used to build \pkg{spsann} were found in the source code of other R-packages and 
#' scripts developed and published by other researchers. For example, the original skeleton of the 
#' optimization functions was adopted from the \pkg{intamapInteractive} package with the approval of the 
#' package authors, Edzer Pebesma \email{edzer.pebesma@@uni-muenster.de} and Jon Skoien 
#' \email{jon.skoien@@gmail.com}. The current skeleton is based on the later adoption of several solutions 
#' implemented in the script developed and published by Murray Lark \email{mlark@@bgs.ac.uk} as part of a
#' short course (\sQuote{Computational tools to optimize spatial sampling}) offered for the first time at the 
#' 2015 EGU General Assembly in Vienna, Austria.
#' 
#' A few small solutions were adopted from the packages \pkg{SpatialTools}, authored by Joshua French 
#' \email{joshua.french@@ucdenver.edu}, \pkg{clhs}, authored by Pierre Roudier 
#' \email{roudierp@@landcareresearch.co.nz}, and \pkg{spcosa}, authored by Dennis Walvoort 
#' \email{dennis.Walvoort@@wur.nl}, Dick Brus \email{dick.brus@@wur.nl}, and Jaap de Gruijter 
#' \email{Jaap.degruijter@@wur.nl}.
#' 
#' Major conceptual contributions were made by Gerard Heuvelink \email{gerard.heuvelink@@wur.nl}, Dick Brus 
#' \email{dick.brus@@wur.nl}, Murray Lark \email{mlark@@bgs.ac.uk}, and Edzer Pebesma 
#' \email{edzer.pebesma@@uni-muenster.de}.
#' 
# General information #########################################################################################
#' @author Author and Maintainer: Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}.
#' @name spsann-package
#' @aliases spsann-package spsann
#' @docType package
#' @useDynLib spsann
NULL

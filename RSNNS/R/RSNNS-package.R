#############################################################################
#
#   This file is part of the R package "RSNNS".
#
#   Author: Christoph Bergmeir
#   Supervisor: José M. Benítez
#   Copyright (c) DiCITS Lab, Sci2s group, DECSAI, University of Granada.
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Library General Public
#   License as published by the Free Software Foundation; either
#   version 2 of the License, or (at your option) any later version.
# 
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Library General Public License for more details.
# 
#   You should have received a copy of the GNU Library General Public License
#   along with this library; see the file COPYING.LIB.  If not, write to
#   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
#   Boston, MA 02110-1301, USA.
#
#############################################################################


# \tabular{ll}{
# Package: \tab RSNNS\cr
# Type: \tab Package\cr
# Version: \tab 0.3-2\cr
# Date: \tab 2010-11-30\cr
# License: \tab LGPL (>= 2)\cr
# LazyLoad: \tab yes\cr
# }
#
#' The Stuttgart Neural Network Simulator (SNNS) is a library containing many 
#' standard implementations of neural networks. This package wraps the SNNS 
#' functionality to make it available from within R.
#' 
#' If you have problems using RSNNS, find a bug, or have suggestions, please
#' contact the package maintainer by email, instead of writing to the general R
#' lists or contacting the authors of the original SNNS software.
#' 
#' If you use the package, please cite the following work in your publications:
#'
#' Bergmeir, C. and Benítez, J.M. (2012), Neural Networks in R Using the Stuttgart Neural Network Simulator: RSNNS. Journal of Statistical Software, 46(7), 1-26. \url{http://www.jstatsoft.org/v46/i07/}
#'
#' The package has a hierarchical architecture with three levels:
#' \itemize{
#' \item RSNNS high-level api (rsnns)
#' \item RSNNS low-level api (SnnsR)
#' \item The api of our C++ port of SNNS (SnnsCLib)
#' }
#' 
#' Many demos for using both low-level and high-level api of the package are
#' available. To get a list of them, type:
#' 
#' \code{library(RSNNS)}
#' 
#' \code{demo()}
#'  
#' It is a good idea to start with the demos of the high-level api (which is
#' much more convenient to use). E.g., to access the iris classification demo
#' type:
#' 
#' \code{demo(iris)}
#'  
#' or for the laser regression demo type:
#' 
#' \code{demo(laser)} 
#' 
#' As the high-level api is already quite powerful and flexible, you'll most
#' probably normally end up using one of the functions: \code{\link{mlp}},
#' \code{\link{dlvq}}, \code{\link{rbf}}, \code{\link{rbfDDA}},
#' \code{\link{elman}},  \code{\link{jordan}}, \code{\link{som}},
#' \code{\link{art1}}, \code{\link{art2}}, \code{\link{artmap}}, or
#' \code{\link{assoz}}, with some pre- and postprocessing. These S3 classes are
#' all subclasses of \code{\link{rsnns}}.
#' 
#' You might also want to have a look at the original SNNS program and the SNNS
#' User Manual 4.2, especially pp 67-87 for explications on all the parameters
#' of the learning functions, and pp 145-215 for detailed (theoretical)
#' explications of the methods and advice on their use. And, there is also the 
#' javaNNS, the sucessor of SNNS from the original authors. It makes the C core 
#' functionality available from a Java GUI.
#' 
#' Demos ending with "SnnsR" show the use of the low-level api. If you want to
#' do special things with neural networks that are currently not implemented in
#' the high-level api, you can see in this demos how to do it. Many demos are
#' present both as high-level and low-level versions.
#' 
#' The low-level api consists mainly of the class \code{\link{SnnsR-class}},
#' which internally holds a pointer to a C++ object of the class
#' \code{SnnsCLib}, i.e., an instance of the SNNS kernel. The class furthermore
#' implements a calling mechanism for methods of the \code{SnnsCLib} object, so
#' that they can be called conveniently using the "$"-operator. This calling
#' mechanism also allows for transparent masking of methods or extending the
#' kernel with new methods from within R. See
#' \code{\link{$,SnnsR-method}}. R-functions that are added by RSNNS to the
#' kernel are documented in this manual under topics beginning with
#' \code{SnnsRObject$}. Documentation of the original SNNS kernel user interface
#' functions can be found in the SNNS User Manual 4.2 pp 290-314.  A call to,
#' e.g., the SNNS kernel function \code{krui_getNoOfUnits(...)} can be done with
#' \code{SnnsRObject$getNoOfUnits(...)}. However, a few functions were excluded
#' from the wrapping for various reasons. Fur more details and other known
#' issues see the file /inst/doc/KnownIssues.
#' 
#' Another nice tool is the \code{NeuralNetTools} package, that can be used to
#' visualize and analyse the networks generated with RSNNS.
#' 
#' Most of the example data included in SNNS is also present in this package, see \code{\link{snnsData}}.
#' 
#' Additional information is also available at the project website: 
#' 
#' \url{http://sci2s.ugr.es/dicits/software/RSNNS}
#' 
#' @title Getting started with the RSNNS package
#' @name RSNNS-package
#' @aliases RSNNS
#' @docType package
#' @title Getting started with the RSNNS package
# @encoding UTF-8
# @encoding Latin-1
#' @author Christoph Bergmeir \email{c.bergmeir@@decsai.ugr.es} 
#' 
#' and José M. Benítez \email{j.m.benitez@@decsai.ugr.es}
#' 
#' DiCITS Lab, Sci2s group, DECSAI, University of Granada.
#' 
#' \url{http://dicits.ugr.es}, \url{http://sci2s.ugr.es}
#' 
#' @references 
#' 
#' Bergmeir, C. and Benítez, J.M. (2012), 'Neural Networks in R Using the Stuttgart Neural Network Simulator: RSNNS', Journal of Statistical Software, 46(7), 1-26.
#' \url{http://www.jstatsoft.org/v46/i07/}
#'
#' \emph{General neural network literature:}
#' 
#' Bishop, C. M. (2003), Neural networks for pattern recognition, University Press, Oxford.
#' 
#' Haykin, S. S. (1999), Neural networks :a comprehensive foundation, Prentice Hall, Upper Saddle River, NJ.
#' 
#' Kriesel, D. ( 2007 ), A Brief Introduction to Neural
#' Networks. http://www.dkriesel.com
#' 
#' Ripley, B. D. (2007), Pattern recognition and neural networks, Cambridge
#' University Press, Cambridge.
#' 
#' Rojas, R. (1996), Neural networks :a systematic introduction, Springer-Verlag, Berlin.
#' 
#' Rumelhart, D. E.; Clelland, J. L. M. & Group, P. R. (1986), Parallel
#' distributed processing :explorations in the microstructure of cognition, Mit,
#' Cambridge, MA etc..
#' 
#' \emph{Literature on the original SNNS software:}
#' 
#' Zell, A. et al. (1998), 'SNNS Stuttgart Neural Network Simulator User Manual, Version 4.2', IPVR, University of Stuttgart and WSI, University of Tübingen.
#' \url{http://www.ra.cs.uni-tuebingen.de/SNNS/}
#'
#' javaNNS, the sucessor of the original SNNS with a Java GUI: 
#' \url{http://www.ra.cs.uni-tuebingen.de/software/JavaNNS}
#'
#' Zell, A. (1994), Simulation Neuronaler Netze, Addison-Wesley.
#' 
#' \emph{Other resources:}
#' 
#' A function to plot networks from the \code{\link{mlp}} function:
#' \url{https://beckmw.wordpress.com/2013/11/14/visualizing-neural-networks-in-r-update/}
#' 

#' 
#' @keywords package neural networks SNNS
#' @seealso \code{\link{mlp}}, \code{\link{dlvq}}, \code{\link{rbf}}, \code{\link{rbfDDA}}, \code{\link{elman}}, 
#' \code{\link{jordan}}, \code{\link{som}}, \code{\link{art1}}, \code{\link{art2}}, \code{\link{artmap}}, \code{\link{assoz}}
#' @useDynLib RSNNS .registration=TRUE
#' @import methods
#' @import Rcpp
# @exportPattern "^[[:alpha:]]+"
# @examples
# \dontrun{mlp(...)}
NULL


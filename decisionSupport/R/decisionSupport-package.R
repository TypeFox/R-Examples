#
# file: decisionSupport-package.R
#
# This file is part of the R-package decisionSupport
# 
# Authors: 
#   Lutz Göhring <lutz.goehring@gmx.de>
#   Eike Luedeling (ICRAF) <eike@eikeluedeling.com>
#
# Copyright (C) 2015 World Agroforestry Centre (ICRAF) 
#	http://www.worldagroforestry.org
# 
# The R-package decisionSupport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# The R-package decisionSupport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with the R-package decisionSupport.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################################
#' @include decisionSupport.R
NULL
##############################################################################################
#' Quantitative Support of Decision Making under Uncertainty.
#' 
#' The \pkg{\code{decisionSupport}} package supports the quantitative analysis of 
#' welfare based decision making processes using Monte Carlo simulations. This
#' is an important part of the Applied Information Economics (AIE) approach
#' developed in Hubbard (2014). These decision making processes can be
#' categorized into two levels of decision making:
#' \enumerate{ 
#' 	\item The actual problem of interest of a policy maker
#' 	which we call the \emph{underlying welfare based decision} on how to influence an
#' 	ecological-economic system based on a particular information on the system
#' 	available to the decision maker and 
#' 	\item the \emph{meta decision} on how to allocate resources to reduce the
#' 	uncertainty in the underlying decision problem, i.e to increase the current
#' 	information to improve the underlying decision making process.
#' } 
#' The first problem, i.e. the underlying problem, is the problem of choosing
#' the decision which maximizes expected welfare. The welfare function can be
#' interpreted as a von Neumann-Morgenstern utility function. Whereas, the
#' second problem, i.e. the meta decision problem, is dealt with using the
#' \emph{Value of Information Analysis (VIA)}. Value of Information Analysis
#' seeks to assign a value to a certain reduction in uncertainty or,
#' equivalently, increase in information. Uncertainty is dealt with in a
#' probabilistic manner. Probabilities are transformed via Monte Carlo
#' simulations.
#' 
#' 
#' The functionality of this package is subdivided into three main parts: (i) the
#' welfare based analysis of the underlying decision, (ii) the meta decision of
#' reducing uncertainty and (iii) the Monte Carlo simulation for the
#' transformation of probabilities and calculation of expectation values. Furthermore, 
#' there is a wrapper function around these three parts which aims at providing an 
#' easy-to-use interface.
#' \subsection{Welfare based Analysis of the Underlying Decision Problem}{
#'				Implementation: \code{\link{welfareDecisionAnalysis}}
#'		}
#' \subsection{The Meta Decision of Reducing Uncertainty}{ 
#' 		The meta decision of how to allocate resources for uncertainty reduction can
#' 		be analyzed with this package in two different ways: via (i) Expected Value
#' 		of Information Analysis or (ii) via Partial Least Squares (PLS) analysis and
#' 		Variable Importance in Projection (VIP).
#' 		\subsection{Expected Value of Information (EVI)}{
#' 				Implementation: \code{\link{eviSimulation}}, \code{\link{individualEvpiSimulation}}
#' 		}
#' 		\subsection{Partial Least Squares (PLS) analysis and Variable Importance in Projection (VIP)}{
#' 				Implementation: \code{\link{plsr.mcSimulation}}, \code{\link[chillR:VIP]{VIP}}
#' 		}
#' }
#' \subsection{Solving the Practical Problem of Calculating Expectation Values by Monte Carlo Simulation}{
#' 		\subsection{Estimates}{
#' 			Implementation: \code{\link{estimate}}
#' 		}
#' 		\subsection{Multivariate Random Number Generation}{
#' 			Implementation: \code{\link{random.estimate}}
#' 		}
#' 		\subsection{Monte Carlo Simulation}{
#' 			Implementation: \code{\link{mcSimulation}}
#' 		}
#' }
#' \subsection{Integrated Welfare Decision and Value of Information Analysis: A wrapper function}{
#' 	The function \code{\link{decisionSupport}} integrates the most important features of this 
#' 	package into a single function. It is wrapped arround the functions 
#' 	\code{\link{welfareDecisionAnalysis}}, \code{\link{plsr.mcSimulation}}, 
#' 	\code{\link[chillR:VIP]{VIP}} and \code{\link{individualEvpiSimulation}}.
#' }
#'
#' 
#' @section Copyright \ifelse{latex}{\out{\copyright}}{\ifelse{html}{\out{&copy}}{(C)}}: 
#' 	  \href{http://www.worldagroforestry.org/}{World Agroforestry Centre (ICRAF)} 2015
#' @section License:
#'    The R-package \pkg{decisionSupport} is free software: you can redistribute it and/or modify
#'    it under the terms of the GNU General Public License as published by
#'    the Free Software Foundation, either version 3 of the License, or
#'    (at your option) any later version: 
#'    \href{http://www.gnu.org/licenses/gpl-3.0.html}{GNU GENERAL PUBLIC LICENSE, Version 3 (GPL-3)}
#' 
#'    The R-package \pkg{decisionSupport} is distributed in the hope that it will be useful,
#'    but WITHOUT ANY WARRANTY; without even the implied warranty of
#'    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#'    GNU General Public License for more details.
#' 
#'    You should have received a copy of the GNU General Public License
#'    along with the R-package decisionSupport.  If not, see \url{http://www.gnu.org/licenses/}.
#' @docType package
#' @name decisionSupport-package
#' @author Lutz \enc{Göhring}{Goehring} \email{lutz.goehring@@gmx.de},
#'	 Eike Luedeling (\href{http://www.worldagroforestry.org/}{ICRAF}) \email{eike@@eikeluedeling.com}	
#' @author Maintainer: Eike Luedeling \email{eike@@eikeluedeling.com}
#' @references Hubbard, Douglas W., \emph{How to Measure Anything? - Finding the Value of "Intangibles" in Business},
#'   John Wiley & Sons, Hoboken, New Jersey, 2014, 3rd Ed, \url{http://www.howtomeasureanything.com/}.
#'   
#'   Hugh Gravelle and Ray Rees, \emph{Microeconomics}, Pearson Education Limited, 3rd edition, 2004.
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline axis barplot hist par
#' @importFrom stats qnorm quantile rbeta rcauchy rchisq rexp rf rgamma rlnorm rlogis rnorm rt runif rweibull sd
#' @importFrom utils capture.output read.csv write.csv
#'   
#' @seealso \code{\link{welfareDecisionAnalysis}}, \code{\link{eviSimulation}}, \code{\link{mcSimulation}}	
NULL

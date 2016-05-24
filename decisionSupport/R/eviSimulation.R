#
# file: eviSimulation.R
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
#' @include welfareDecisionAnalysis.R
NULL
##############################################################################################
# eviSimulation(welfare, currentEstimate, prospectiveEstimate, numberOfModelRuns, functionSyntax)
##############################################################################################
#' Expected Value of Information (EVI) Simulation.
#' 
#' The Expected Value of Information (EVI) is calculated based on a Monte Carlo simulation of the
#' expected welfare (or values or benefits) of two different decision alternatives. The expected
#' welfare is calculated for the current estimate of variables determining welfare and a prospective
#' estimate of these variables. The prospective estimate resembles an improvement in information.
#' @param currentEstimate \code{\link{estimate}}: describing the distribution of the input variables
#'   as currently being estimated.
#' @param prospectiveEstimate \code{\link{estimate}} or \code{list} of \code{estimate} objects:
#'   describing the prospective distribution of the input variables which could hypothetically
#'   be achieved by collecting more information, viz. improving the measurement.
#' @inheritParams welfareDecisionAnalysis
#' @return An object of class \code{eviSimulation} with the following elements:
#'  \describe{
#' 			\item{\code{$current}}{
#' 			   \code{\link{welfareDecisionAnalysis}} object for \code{currentEstimate}
#' 			}
#' 			\item{\code{$prospective}}{
#' 			  \code{\link{welfareDecisionAnalysis}} object  for single \code{prospectiveEstimate} or a 
#' 			  list of \code{\link{welfareDecisionAnalysis}} objects for \code{prospectiveEstimate} being
#' 			  a list of \code{estimate}s.
#' 			}
#'  		\item{\code{$evi}}{
#'  		  Expected Value of Information(s) (EVI)(s)  gained by the prospective estimate(s) w.r.t. the 
#'  		  current estimate.
#'  		}
#'   }
#' @details
#'   \subsection{The Expected Value of Information (EVI)}{
#'     The Expected Value of Information is the decrease in the \eqn{\textrm{EOL}}{EOL} for an information
#'     improvement from the current (\eqn{\rho_X^{current}}{\rho_X_current}) to a better prospective (hypothetical)
#'     information (\eqn{\rho_X^{prospective}}{\rho_X_prospective}):
#'     \deqn{
#'       \textrm{EVI} := \textrm{EOL}(\rho_X^{current}) - \textrm{EOL}(\rho_X^{prospective}).
#'       }{
#'          EVI := EOL(\rho_X_current) - EOL(\rho_X_prospective).
#'       }
#'   }
#' @references Hubbard, Douglas W., \emph{How to Measure Anything? - Finding the Value of "Intangibles" in Business},
#'   John Wiley & Sons, Hoboken, New Jersey, 2014, 3rd Ed, \url{http://www.howtomeasureanything.com/}.
#'   
#'   Gravelle, Hugh and Ray Rees, \emph{Microeconomics}, Pearson Education Limited, 3rd edition, 2004.
#' @seealso \code{\link{welfareDecisionAnalysis}}, \code{\link{mcSimulation}}, \code{\link{estimate}},
#'   \code{\link{summary.eviSimulation}}
#' @examples 
#' #############################################################
#' # Example 1 Only one prospective estimate:
#' #############################################################
#' numberOfModelRuns=10000
#' # Create the estimate object:
#' variable=c("revenue","costs")
#' distribution=c("posnorm","posnorm")
#' lower=c(10000,  5000)
#' upper=c(100000, 50000)
#' currentEstimate<-as.estimate(variable, distribution, lower, upper)
#' prospectiveEstimate<-currentEstimate
#' revenueConst<-mean(c(currentEstimate$marginal["revenue","lower"],
#'                      currentEstimate$marginal["revenue","upper"]))
#' prospectiveEstimate$marginal["revenue","distribution"]<-"const"
#' prospectiveEstimate$marginal["revenue","lower"]<-revenueConst 
#' prospectiveEstimate$marginal["revenue","upper"]<-revenueConst 
#' # (a) Define the welfare function without name for the return value:
#' profit<-function(x){
#' 	x$revenue-x$costs
#' }
#' 
#' # Calculate the Expected Value of Information:
#' eviSimulationResult<-eviSimulation(welfare=profit,
#'                                    currentEstimate=currentEstimate,
#'                                    prospectiveEstimate=prospectiveEstimate,
#'                                    numberOfModelRuns=numberOfModelRuns,
#'                                    functionSyntax="data.frameNames")
#' # Show the simulation results:
#' print(summary(eviSimulationResult))
#' #############################################################
#' # (b) Define the welfare function with a name for the return value:
#' profit<-function(x){
#' 	list(Profit=x$revenue-x$costs)
#' }
#' # Calculate the Expected Value of Information:
#' eviSimulationResult<-eviSimulation(welfare=profit,
#'                                    currentEstimate=currentEstimate,
#'                                    prospectiveEstimate=prospectiveEstimate,
#'                                    numberOfModelRuns=numberOfModelRuns,
#'                                    functionSyntax="data.frameNames")
#' # Show the simulation results:
#' print(summary((eviSimulationResult)))
#' #############################################################
#' # (c) Two decision variables:
#' decisionModel<-function(x){
#'  list(Profit=x$revenue-x$costs,
#'       Costs=-x$costs)
#' }
#' # Calculate the Expected Value of Information:
#' eviSimulationResult<-eviSimulation(welfare=decisionModel,
#'                                    currentEstimate=currentEstimate,
#'                                    prospectiveEstimate=prospectiveEstimate,
#'                                    numberOfModelRuns=numberOfModelRuns,
#'                                    functionSyntax="data.frameNames")
#' # Show the simulation results:
#' print(summary((eviSimulationResult)))
#' #############################################################
#' # Example 2 A list of prospective estimates:
#' #############################################################
#' numberOfModelRuns=10000
#' #  Define the welfare function with a name for the return value:
#' profit<-function(x){
#'  list(Profit=x$revenue-x$costs)
#' }
#' # Create the estimate object:
#' variable=c("revenue","costs")
#' distribution=c("posnorm","posnorm")
#' lower=c(10000,  5000)
#' upper=c(100000, 50000)
#' currentEstimate<-as.estimate(variable, distribution, lower, upper)
#' perfectInformationRevenue<-currentEstimate
#' revenueConst<-mean(c(currentEstimate$marginal["revenue","lower"],
#'                      currentEstimate$marginal["revenue","upper"]))
#' perfectInformationRevenue$marginal["revenue","distribution"]<-"const"
#' perfectInformationRevenue$marginal["revenue","lower"]<-revenueConst 
#' perfectInformationRevenue$marginal["revenue","upper"]<-revenueConst
#' # (a) A list with one element
#' prospectiveEstimate<-list(perfectInformationRevenue=perfectInformationRevenue)
#' # Calculate the Expected Value of Information:
#' eviSimulationResult<-eviSimulation(welfare=profit,
#'                                    currentEstimate=currentEstimate,
#'                                    prospectiveEstimate=prospectiveEstimate,
#'                                    numberOfModelRuns=numberOfModelRuns,
#'                                    functionSyntax="data.frameNames")
#' # Show the simulation results:
#' print(summary(eviSimulationResult))
#' #############################################################
#' # (b) A list with two elements
#' perfectInformationCosts<-currentEstimate
#' costsConst<-mean(c(currentEstimate$marginal["costs","lower"], 
#'                    currentEstimate$marginal["costs","upper"]))
#' perfectInformationCosts$marginal["costs","distribution"]<-"const"
#' perfectInformationCosts$marginal["costs","lower"]<-costsConst 
#' perfectInformationCosts$marginal["costs","upper"]<-costsConst
#' prospectiveEstimate<-list(perfectInformationRevenue=perfectInformationRevenue,
#'                           perfectInformationCosts=perfectInformationCosts)
#' # Calculate the Expected Value of Information:
#' eviSimulationResult<-eviSimulation(welfare=profit,
#'                                    currentEstimate=currentEstimate,
#'                                    prospectiveEstimate=prospectiveEstimate,
#'                                    numberOfModelRuns=numberOfModelRuns,
#'                                    functionSyntax="data.frameNames")
#' # Show the simulation results:
#' print(summary(eviSimulationResult))
#' #############################################################
#' # Example 3 A list of prospective estimates and two decision variables:
#' #############################################################
#' numberOfModelRuns=10000
#' # Create the current estimate object:
#' variable=c("revenue","costs")
#' distribution=c("posnorm","posnorm")
#' lower=c(10000,  5000)
#' upper=c(100000, 50000)
#' currentEstimate<-as.estimate(variable, distribution, lower, upper)
#' # Create a list of two prospective estimates:
#' perfectInformationRevenue<-currentEstimate
#' revenueConst<-mean(c(currentEstimate$marginal["revenue","lower"],
#'                      currentEstimate$marginal["revenue","upper"]))
#' perfectInformationRevenue$marginal["revenue","distribution"]<-"const"
#' perfectInformationRevenue$marginal["revenue","lower"]<-revenueConst 
#' perfectInformationRevenue$marginal["revenue","upper"]<-revenueConst
#' perfectInformationCosts<-currentEstimate
#' costsConst<-mean(c(currentEstimate$marginal["costs","lower"], 
#'                    currentEstimate$marginal["costs","upper"]))
#' perfectInformationCosts$marginal["costs","distribution"]<-"const"
#' perfectInformationCosts$marginal["costs","lower"]<-costsConst 
#' perfectInformationCosts$marginal["costs","upper"]<-costsConst
#' prospectiveEstimate<-list(perfectInformationRevenue=perfectInformationRevenue,
#'                           perfectInformationCosts=perfectInformationCosts)
#' # Define the welfare function with two decision variables:
#' decisionModel<-function(x){
#'  list(Profit=x$revenue-x$costs,
#'       Costs=-x$costs)
#' }
#' # Calculate the Expected Value of Information:
#' eviSimulationResult<-eviSimulation(welfare=decisionModel,
#'                                    currentEstimate=currentEstimate,
#'                                    prospectiveEstimate=prospectiveEstimate,
#'                                    numberOfModelRuns=numberOfModelRuns,
#'                                    functionSyntax="data.frameNames")
#' # Show the simulation results:
#' print(sort(summary(eviSimulationResult)),decreasing=TRUE,along="Profit")
#' @export
eviSimulation<-function(welfare, currentEstimate, prospectiveEstimate, numberOfModelRuns, 
                        randomMethod="calculate", 
                        functionSyntax="data.frameNames",
                        relativeTolerance=0.05,
                        verbosity=0){
  # Return object:
  thisAnalysis<-NULL
  # Perform the current decision analysis:
  analysisCurrent<-welfareDecisionAnalysis( estimate=currentEstimate,
                                            welfare=welfare,
                                            numberOfModelRuns=numberOfModelRuns,
                                            randomMethod=randomMethod,
                                            functionSyntax=functionSyntax, 														
                                            relativeTolerance=relativeTolerance,
                                            verbosity=verbosity)
  
  # Perform the prospective decision analysis:
  if( class(prospectiveEstimate) == "estimate"){
    # Perform the decision analysis:
    analysisProspective<-welfareDecisionAnalysis( estimate=prospectiveEstimate,
                                                  welfare=welfare,
                                                  numberOfModelRuns=numberOfModelRuns,
                                                  randomMethod=randomMethod,
                                                  functionSyntax=functionSyntax, 														
                                                  relativeTolerance=relativeTolerance,
                                                  verbosity=verbosity)
    evi<-analysisCurrent$eol - analysisProspective$eol
  } else if ( is.list(prospectiveEstimate) ){
    analysisProspective<-lapply(X=prospectiveEstimate, 
                                FUN=function(estimate) welfareDecisionAnalysis(estimate=estimate,
                                                                               welfare=welfare,
                                                                               numberOfModelRuns=numberOfModelRuns,
                                                                               randomMethod=randomMethod,
                                                                               functionSyntax=functionSyntax, 														
                                                                               relativeTolerance=relativeTolerance,
                                                                               verbosity=verbosity)
    )
    evi<-lapply(X=analysisProspective, 
                FUN=function(x) analysisCurrent$eol - x$eol)
  } else {
    stop("prospectiveEstimate must be either an estimate or a list of estimates.")
  }	
  
  # Fill return object:
  thisAnalysis$call<-match.call()
  thisAnalysis$current<-analysisCurrent
  thisAnalysis$prospective<-analysisProspective
  thisAnalysis$evi<-as.data.frame(evi)
  class(thisAnalysis) <- "eviSimulation"
  return(thisAnalysis)
}
##############################################################################################
# summary.eviSimulation(object, ...)
##############################################################################################
#' Summarize EVI Simulation Results
#' 
#'  Produces result summaries of an Expected Value of Information (EVI) simulation obtained by 
#'  the function \code{\link{eviSimulation}}.
#' @param object An object of class \code{eviSimulation}.
#' @param ... Further arguments passed to \code{\link{summary.welfareDecisionAnalysis}}.
#' @inheritParams base::format
#' @return An object of class \code{summary.eviSimulation}.
#' @seealso \code{\link{eviSimulation}}, \code{\link{print.summary.eviSimulation}}, 
#' \code{\link{summary.welfareDecisionAnalysis}}, \ifelse{latex}{\cr}{ }
#'  \code{\link{sort.summary.eviSimulation}}
#' @export
summary.eviSimulation <- function(object,
                                  ...,
                                  digits = max(3, getOption("digits")-3)){	
  summaryList<-list(evi=format(x=object$evi, digits=digits),
                    current=summary(object$current, ..., digits=digits)$summary,
                    prospective=if( class(object$prospective)=="welfareDecisionAnalysis" ){
                      summary(object$prospective, ..., digits=digits)$summary
                    } else {
                      lapply(X=object$prospective, 
                             FUN=function(x) summary(x, ..., digits=digits)$summary
                      )
                    }
  )
  #	summaryList<-format(x=summaryList, digits=digits, ...)
  res<-list(summary=summaryList,
            call=object$call)
  
  class(res)<-"summary.eviSimulation"
  res
}
##############################################################################################
# sort.summary.eviSimulation(x, decreasing, ..., along)
##############################################################################################
#' Sort Summarized EVI Simulation Results..
#' 
#' Sort summarized EVI simulation results according to their EVI.
#' @param x An object of class \code{summary.eviSimulation}.
#' @param decreasing \code{logical}: if the EVI should be sorted in decreasing order.
#' @param ... currently not used
#' @param along \code{character}: the name of the valuation variable along which the EVI 
#'  should be sorted.
#' @return An object of class \code{summary.eviSimulation}.
#' @seealso \code{\link{eviSimulation}}, \code{\link{summary.eviSimulation}}, \code{\link[base]{sort}}
#' @export
sort.summary.eviSimulation <- function(x, decreasing=TRUE, ..., along=row.names(x$summary$evi)[[1]]){
  eviRanking<-order(x=as.numeric(x$summary$evi[along,]), decreasing=decreasing)
  eviRankingNames<-names(x$summary$evi)[eviRanking]
  x$summary$evi<-x$summary$evi[eviRankingNames]
  x$summary$prospective<-x$summary$prospective[eviRankingNames]
  x
}
##############################################################################################
# print.summary.eviSimulation(x, ...)
##############################################################################################
#' Print the Summarized EVI Simulation Results.
#' 
#' This function prints the summary of \code{eviSimulation} generated by 
#'\code{\link{summary.eviSimulation}}.
#' @param x An object of class \code{summary.eviSimulation}.
#' @param ... Further arguments to be passed to \code{\link{print.default}} and \ifelse{latex}{\cr}{ }
#'   \code{\link{print.summary.welfareDecisionAnalysis}}.
#' @seealso \code{\link{eviSimulation}}, \code{\link{print.summary.welfareDecisionAnalysis}}.
#' @export
print.summary.eviSimulation <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nSummary of EVI simulation:\n")
  cat("\nExpected Value of Information (EVI):\n")
  print(x$summary$evi,...)
  cat("\nUnderlying welfare decision analysis:\n")
  cat("Based on the current estimate:\n")
  print(x$summary$current, ...)
  cat("\nBased on the prospective estimate(s):\n")
  print(x$summary$prospective, ...)
}
##############################################################################################
# hist.eviSimulation(x, ...)
##############################################################################################
#' Plot Histograms of results of an EVI simulation
#' 
#' This function plots the histograms of the results of
#' \code{\link{eviSimulation}}.
#' @param x An object of class \code{eviSimulation}.
#' @param mainSuffix \code{character}: Suffix of the main titles of the histograms.
#' @inheritParams graphics::hist
#' @param ... Further arguments to be passed to \code{\link[graphics]{hist}}.
#' @param colorQuantile \code{character vector}: encoding the colors of the 
#'   quantiles defined in argument \code{colorProbability}.
#' @param colorProbability \code{numeric vector}: defines the quantiles that 
#'   shall be distinguished by the colors chosen in argument 
#'   \code{colorQuantile}. Must be of the same length as \code{colorQuantile}.
#' @param resultName \code{character}: indicating the name of the component of
#'   the simulation function (\code{model_function}) which results histogram
#'   shall be generated. If \code{model_function} is single valued, no name
#'   needs to be supplied. Otherwise, one valid name has to be specified.
#'   Defaults to \code{NULL}.
#' @return an object of class "\code{histogram}". For details see 
#'   \code{\link[graphics]{hist}}.
#' @seealso \code{\link{eviSimulation}}, \code{\link{hist}}. For a list of colors
#'   available in R see \code{\link[grDevices]{colors}}.
#' @export
hist.eviSimulation <- function(x, breaks=100, col=NULL, mainSuffix=" welfare simulation result", ...,
                               colorQuantile   =c("GREY", "YELLOW", "ORANGE", "DARK GREEN", "ORANGE", "YELLOW", "GREY"), 
                               colorProbability=c(1.00,    0.95,     0.75,     0.55,         0.45,     0.25,     0.05),
                               resultName=NULL){
  #numberOfPlots<- 1 + ifelse(class(x$prospective)=="welfareDecisionAnalysis", 
  #                           yes=1, 
  #                          no=length(x$prospective)
  #)
  #par(mfcol=c(numberOfPlots,1))
  #layout(1:numberOfPlots, respect=TRUE)
  #split.screen(figs=c(numberOfPlots,1))
  # Plot the distribution for the current information for the chosen component:
  hist(x$current, breaks=breaks, col=col,  
       main=paste("Current", mainSuffix), 
       ...,
       colorQuantile   =colorQuantile, 
       colorProbability=colorProbability,
       resultName=resultName)
  # Plot the distribution(s) for the prospective information for the chosen component:
  if( class(x$prospective)=="welfareDecisionAnalysis" ){
    hist(x$prospective, breaks=breaks, col=col, 
         main=paste("Prospective", mainSuffix), 
         ...,
         colorQuantile   =colorQuantile, 
         colorProbability=colorProbability,
         resultName=resultName)
  } else {
    for( i in names(x$prospective) )
      hist(x$prospective[[i]], breaks=breaks, col=col, 
           main=paste("Prospective (\"", i, "\")", mainSuffix), 
           ...,
           colorQuantile   =colorQuantile, 
           colorProbability=colorProbability,
           resultName=resultName)
  }
  
}

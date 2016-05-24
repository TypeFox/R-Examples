#
# file: individualEvpiSimulation.R
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
#' @include eviSimulation.R
NULL
##############################################################################################
# individualEvpiSimulation(welfare, currentEstimate, perfectProspectiveNames,perfectProspectiveValues,
#                           numberOfModelRuns, functionSyntax)
##############################################################################################
#' Individual Expected Value of Perfect Information Simulation
#' 
#' The Individual Expected Value of Perfect Information (Individual EVPI) is calculated based on a
#' Monte Carlo simulation of the values of two different decision alternatives.
#' @param currentEstimate \code{\link{estimate}}: describing the distribution of the input variables
#'   as currently being estimated.
#' @param perfectProspectiveNames \code{character vector}: input variable names that are assumed to be known perfectly with 
#' 				prospective information.
#' @param perfectProspectiveValues \code{numeric vector}: of the same length as \code{perfectProspectiveNames} with the corresponding
#' 				values assumed to be known perfectly.
#' @inheritParams welfareDecisionAnalysis
#' @return An object of class \code{eviSimulation} with the following elements:
#'  \describe{
#' 			\item{\code{$current}}{
#' 			   \code{\link{welfareDecisionAnalysis}} object for \code{currentEstimate}
#' 			}
#' 			\item{\code{$prospective}}{
#' 			  \code{\link{welfareDecisionAnalysis}} object  for single \code{perfectProspectiveNames} or a 
#' 			  list of \code{\link{welfareDecisionAnalysis}} objects for several \code{perfectProspectiveNames}.
#' 			}
#'  		\item{\code{$evi}}{
#'  		  Expected Value of Information(s) (EVI)(s) gained by the perfect knowledge of individual 
#'  		  variable(s) w.r.t. the current estimate.
#'  		}
#'   }
#' @details The Individual EVPI is defined as the EVI with respect to a prospective information 
#' that assumes perfect knowledge on one particular variable.
#'
#' @examples
#' # Number of running the underlying welfare model:
#' n=10000
#' # Create the current estimate from text:
#' estimateText<-"variable,  distribution, lower, upper
#'                revenue1,  posnorm,      100,   1000
#'                revenue2,  posnorm,      50,    2000
#'                costs1,    posnorm,      50,    2000
#'                costs2,    posnorm,      100,   1000"
#' currentEstimate<-as.estimate(read.csv(header=TRUE, text=estimateText, 
#'                           strip.white=TRUE, stringsAsFactors=FALSE))
#' # The welfare function:
#' profitModel <- function(x){
#'  list(Profit=x$revenue1 + x$revenue2 - x$costs1 - x$costs2)
#' }
#' # Calculate the Individual EVPI:
#' individualEvpiResult<-individualEvpiSimulation(welfare=profitModel,
#'                                                currentEstimate=currentEstimate,
#'                                                numberOfModelRuns=n,
#'                                                functionSyntax="data.frameNames")
#' # Show the simulation results:
#' print(sort(summary(individualEvpiResult)),decreasing=TRUE,along="Profit")
#' hist(individualEvpiResult, breaks=100)
#' @seealso \code{\link{eviSimulation}}, \code{\link{welfareDecisionAnalysis}}, \code{\link{mcSimulation}}, \code{\link{estimate}}
#' @export
individualEvpiSimulation <- function(welfare, currentEstimate, 
                                     perfectProspectiveNames=row.names(currentEstimate),
                                     #perfectProspectiveValues=colMeans(random(rho=currentEstimate, n=numberOfModelRuns, method=randomMethod, relativeTolerance=relativeTolerance)[,perfectProspectiveNames]),
                                     perfectProspectiveValues=colMeans(as.data.frame(random(rho=currentEstimate, n=numberOfModelRuns, method=randomMethod, relativeTolerance=relativeTolerance))[perfectProspectiveNames]),
                                     numberOfModelRuns,
                                     randomMethod="calculate",
                                     functionSyntax="data.frameNames",
                                     relativeTolerance=0.05,
                                     verbosity=0){
  prospectiveEstimate<-c()
  #print(perfectProspectiveValues)
  for( i in perfectProspectiveNames){
    # Set marginal information to certainty:
    prospectiveEstimate[[i]]<-currentEstimate
    prospectiveEstimate[[i]]$marginal[i,"distribution"]<-"const"
    prospectiveEstimate[[i]]$marginal[i,"lower"]<-perfectProspectiveValues[[i]]
    prospectiveEstimate[[i]]$marginal[i,"upper"]<-perfectProspectiveValues[[i]]
    # Remove correlation information, if exists, due to certainty:
    if( !is.null(prospectiveEstimate[[i]]$correlation_matrix) ){
      correlationMatrix<-corMat(prospectiveEstimate[[i]])[!(row.names(prospectiveEstimate[[i]]) %in% i),
                                                          !(row.names(prospectiveEstimate[[i]]) %in% i)]
      if( !is.matrix(correlationMatrix) )
        corMat(prospectiveEstimate[[i]])<-NULL
      else 
        corMat(prospectiveEstimate[[i]])<-correlationMatrix
    }
  }
  # Calculate the Expected Value of Perfect Information:
  evpiResult<-eviSimulation(welfare=welfare,
                            currentEstimate=currentEstimate,
                            prospectiveEstimate=prospectiveEstimate,
                            numberOfModelRuns=numberOfModelRuns,
                            randomMethod=randomMethod,
                            functionSyntax=functionSyntax, 														
                            relativeTolerance=relativeTolerance,
                            verbosity=verbosity)
  #	class(evpiResult)<-c("individualEvpiSimulation", class(evpiResult))
  evpiResult
}



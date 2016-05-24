#
# file: welfareDecisionAnalysis.R
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
#' @include mcSimulation.R
NULL
#############################################################
# welfareDecisionAnalysis(estimate, welfare, numberOfModelRuns, functionSyntax)
#############################################################
#' Analysis of the underlying welfare based decision problem.
#' 
#' The optimal choice between two different opportunities is calculated. Optimality is taken as 
#' maximizing expected welfare. Furthermore, the Expected Opportunity Loss (EOL) is calculated.
#' @param estimate \code{\link{estimate}} object describing the distribution of the input variables.
#' @param welfare either a \code{function} or a \code{list} with two \code{functions}, i.e.
#'   \code{list(p1,p2)}. In the first case the function is the net benefit (or welfare) of project approval (PA) vs.
#'   the status quo (SQ). In the second case the element \code{p1} is the function valuing the first
#'   project and the element \code{p2} valuing the second project, viz. the welfare function of \code{p1}
#'   and \code{p2} respectively.
#' @param numberOfModelRuns \code{integer}: The number of running the welfare model for the 
#' underlying Monte Carlo simulation.
#' @param randomMethod \code{character}: The method to be used to sample the distribution
#'   representing the input estimate. For details see option \code{method} in 
#'   \code{\link{random.estimate}}.
#' @param functionSyntax  \code{character}: function syntax used in the welfare function(s). For 
#'   details see \code{\link{mcSimulation}}.
#' @param relativeTolerance \code{numeric}: the relative tolerance level of deviation of the
#'   generated confidence interval from the specified interval. If this deviation is greater than
#'   \code{relativeTolerance} a warning is given.
#' @param verbosity \code{integer}: if \code{0} the function is silent; the larger the value the
#'   more verbose is output information.
#' @return An object of class \code{welfareDecisionAnalysis} with the following elements:
#'  \describe{
#'      \item{\code{$mcResult}}{The results of the Monte Carlo analysis of \code{estimate} 
#'      transformed by \code{welfare}} (cf. \code{\link{mcSimulation}}).
#'      \item{\code{$enbPa}}{Expected Net Benefit of project approval: ENB(PA)}
#' 			\item{\code{$elPa}}{Expected Loss in case of project approval: EL(PA)}
#' 			\item{\code{$elSq}}{Expected Loss in case of status quo: EL(SQ)}
#'  		\item{\code{$eol}}{Expected Opportunity Loss: EOL}
#'  		\item{\code{$optimalChoice}}{
#'  		    The optimal choice, i.e. either project approval (PA) or the status quo (SQ).
#'  		    }
#' }
#' @details 
#'   \subsection{The underlying decision problem and its notational framework}{ 
#'   We are considering a
#'   decision maker who can influence an ecological-economic system having two alternative decisions
#'   \eqn{d_1} and \eqn{d_2} at hand. We assume, that the system can be characterized by the 
#'   \eqn{n-}dimensional
#'   vector \eqn{X}. The characteristics \eqn{X}, are not necessarily known exactly to the decision maker.
#'   However, we assume furthermore that she is able to quantify this uncertainty which we call an
#'   \emph{\link{estimate}} of the characteristics. Mathematically, an estimate is a random variable with
#'   probability density \eqn{\rho_X}.
#'   
#'   Furthermore, the characteristics \eqn{X} determine the welfare \eqn{W(d)} according to the welfare
#'   function \eqn{w_d}: 
#'   \deqn{ 
#'       W_d = w_d (X) 
#'   }{
#'     W(d) = w_d (X) 
#'   } 
#'   Thus, the welfare of decision \eqn{d} is also a random
#'   variable whose probability distribution we call \eqn{\rho_{W_d}}{rho(W(d))}. The welfare function \eqn{w_d} values
#'   the decision \eqn{d} given a certain state \eqn{X} of the system. In other words, decision \eqn{d_2} is
#'   preferred over decision \eqn{d_1}, if and only if, the expected welfare of decision \eqn{d_2} is
#'   greater than the expected welfare (For a comprehensive
#'   discussion of the concept of social preference ordering and its representation by a welfare
#'   function cf. Gravelle and Rees (2004)) of decsion \eqn{d_1}, formally 
#'   \deqn{
#'     d_1 \prec d_2 \Leftrightarrow \textrm{E}[W_{d_1}] < \textrm{E}[W_{d_2}].
#'    }{
#'      d_1 < d_2 <==> E[W(d_1)] < E[W(d_2)].
#'    }
#'   This means the best decision \eqn{d^*}{d*} is the one which maximizes welfare: 
#'   \deqn{ 
#'     d^* := \arg \max_{d=d_1,d_2} \textrm{E}[W_d]
#'   }{
#'     d* := arg max (d=d_1,d_2) E[W(d)] 
#'   }
#'   
#'   This maximization principle has a dual minimization principle. We define the net benefit
#'   \eqn{\textrm{NB}_{d_1} := W_{d_1} - W_{d_2}}{NB(d_1) := W(d_1) - W(d_2)} as the difference 
#'   between the welfare of the two decision
#'   alternatives. A loss \eqn{L_d} is characterized if a decision \eqn{d} produces a negative net benefit.
#'   No loss occurs if the decision produces a positive net benefit. This is reflected in the formal
#'   definition 
#'   \deqn{
#'       L_d := - \textrm{NB}_d, \textrm{~if~} \textrm{NB}_d  < 0, \textrm{~and~} L_d := 0  
#'          \textrm{~otherwise}.
#'   }{
#'      L(d) :=  - NB(d) if NB(d)  < 0 and L(d) := 0 otherwise.
#'   }
#'   Using this notion it can be shown that the maximization of
#'   expected welfare is equivalent to the minimization of the expected loss 
#'   \eqn{\textrm{EL}_d := \textrm{E}[L_d]}{EL(d) := E[L(d)]}. 
#'   \subsection{The Expected Opportunity Loss (EOL)}{
#'     The use of this concept, here, is in line as described in Hubbard (2014). The Expected
#'     Opportunity Loss (\eqn{\textrm{EOL}}{EOL}) is defined as the expected loss for the best
#'     decision. The best decision minimizes the expected loss:
#'     \deqn{
#'       \textrm{EOL} := \min \left\{ \textrm{EL}_{d_1}, \textrm{EL}_{d_2}\right\}
#'      }{
#'       EOL := min \{ EL(d_1), EL(d_2) \}
#'      }
#'     
#'     The \eqn{\textrm{EOL}}{EOL} is always conditional on the available information which is
#'     characterized by the probability distribution of \eqn{X}
#'     \eqn{\rho_X}: \eqn{\textrm{EOL} = \textrm{EOL}(\rho_X)}{EOL = EOL(\rho_X)}.
#'   }
#'   \subsection{Special case: Status quo and project approval}{
#'     A very common actual binary decision problem is the question if a certain project shall be 
#'     approved versus continuing with business as usual, i.e. the status quo. It appears to be 
#'     natural to identify the status quo with zero welfare. This is a special case ( Actually, one
#'     can show, that this special case is equivalent to the discussion above.) of the binary
#'     decision problem that we are considering here. The two decision alternatives are
#'     \describe{
#'       \item{\eqn{d_1:}}{ project approval (PA) and }
#'       \item{\eqn{d_2:}}{ status quo (SQ),}
#'     }
#'     and the welfare of the approved project (or project outcome or yield of the project) is the
#'     random variable \eqn{W_\textrm{PA}}{W(PA)} with distribution 
#'     \eqn{P_{W_\textrm{PA}} = w_\textrm{PA}(P_X)}{P_W(PA) = w_PA(P_X)}:
#'       \deqn{
#'         W_\textrm{PA} \sim w_\textrm{PA}(P_X) = P_{W_\textrm{PA}}
#'        }{
#'        W(PA) ~ w_PA(P_X) = P_W(PA) 
#'        }
#'     and the welfare of the status quo serves as reference and is normalized to zero:
#'     \deqn{
#'       W_\textrm{SQ} \equiv 0,
#'     }{
#'       W(SQ) = 0
#'     }
#'     which implies zero expected welfare of the status quo:
#'     \deqn{
#'       \textrm{E}[W]_\textrm{SQ} 	= 0.
#'     }{
#'       E[W(SQ)] = 0.
#'     }
#'   }
#'   }
#'   
#' @references Hubbard, Douglas W., \emph{How to Measure Anything? - Finding the Value of "Intangibles" in Business},
#'   John Wiley & Sons, Hoboken, New Jersey, 2014, 3rd Ed, \url{http://www.howtomeasureanything.com/}.
#'   
#'   Gravelle, Hugh and Ray Rees, \emph{Microeconomics}, Pearson Education Limited, 3rd edition, 2004.
#' @seealso \code{\link{mcSimulation}}, \code{\link{estimate}}, \code{\link{summary.welfareDecisionAnalysis}}
#' @examples
#' #############################################################
#' # Example 1 (Creating the estimate from the command line):
#' #############################################################
#' # Create the estimate object:
#' variable=c("revenue","costs")
#' distribution=c("posnorm","posnorm")
#' lower=c(10000,  5000)
#' upper=c(100000, 50000)
#' costBenefitEstimate<-as.estimate(variable, distribution, lower, upper)
#' # (a) Define the welfare function without name for the return value:
#' profit<-function(x){
#'  x$revenue-x$costs
#' }
#' # Perform the decision analysis:
#' myAnalysis<-welfareDecisionAnalysis(estimate=costBenefitEstimate, 
#'                                     welfare=profit, 
#'                                     numberOfModelRuns=100000,
#'                                     functionSyntax="data.frameNames")
#' # Show the analysis results:
#' print(summary((myAnalysis)))
#' #############################################################
#' # (b) Define the welfare function with a name for the return value:
#' profit<-function(x){
#'  list(Profit=x$revenue-x$costs)
#' }
#' # Perform the decision analysis:
#' myAnalysis<-welfareDecisionAnalysis(estimate=costBenefitEstimate, 
#'                                     welfare=profit, 
#'                                     numberOfModelRuns=100000,
#'                                     functionSyntax="data.frameNames")
#' # Show the analysis results:
#' print(summary((myAnalysis)))
#' #############################################################
#' # (c) Two decsion variables:
#' welfareModel<-function(x){
#'  list(Profit=x$revenue-x$costs,
#'    Costs=-x$costs)
#' }
#' # Perform the decision analysis:
#' myAnalysis<-welfareDecisionAnalysis(estimate=costBenefitEstimate, 
#'                                     welfare=welfareModel, 
#'                                     numberOfModelRuns=100000,
#'                                     functionSyntax="data.frameNames")
#' # Show the analysis results:
#' print(summary((myAnalysis)))
#' @export
welfareDecisionAnalysis <- function(estimate, welfare, numberOfModelRuns, 
                                    randomMethod="calculate", 
                                    functionSyntax="data.frameNames",
                                    relativeTolerance=0.05,
                                    verbosity=0){
	# Auxiliary functions (ToDo: check!):
	# Expected loss of project approval
	elPa <- function(netBenefitSample){
		- mean( netBenefitSample*(netBenefitSample<0) )
	}

	# Expected loss of status quo
	elSq <- function(netBenefitSample){
		mean( netBenefitSample*(netBenefitSample>0) )
	}
	# Expected opportunity loss
	eol <- function(netBenefitSample){
		elPa_ <- elPa(netBenefitSample)
		elSq_ <- elSq(netBenefitSample)
		min(elPa_,elSq_)
	}
	# Return object:
	thisAnalysis<-NULL
	if ( is.function(welfare) ) {
		# Perform the Monte Carlo simulation:
		mcResult<-mcSimulation( estimate=estimate, 
														model_function=welfare, 
														numberOfModelRuns=numberOfModelRuns,
														randomMethod=randomMethod,
														functionSyntax=functionSyntax,
														relativeTolerance=relativeTolerance,
														verbosity=verbosity)
		# Expected net benefit of project approval:
		enbPa_<-colMeans(mcResult$y)
		# Expected loss for project aproval:
		elPa_<-apply(X=mcResult$y, MARGIN=2, FUN=elPa)
		# Expected loss for status quo:
		elSq_<-apply(X=mcResult$y, MARGIN=2, FUN=elSq)
		# Expected opportunity loss:
		eol_ <-pmin(elPa_,elSq_)
		# The optimal choice (either project aproval (PA) or the status quo (SQ)):
		optimalChoice_<-ifelse( eol_==elPa_, "PA", "SQ")
		# Fill return object:
		thisAnalysis$call<-match.call()
		thisAnalysis$mcResult<-mcResult
		thisAnalysis$enbPa<-enbPa_
		thisAnalysis$elPa<-elPa_
		thisAnalysis$elSq<-elSq_
		thisAnalysis$eol<-eol_
		thisAnalysis$optimalChoice<-optimalChoice_
	} else if ( is.list(welfare) ){
		stop("The general case of two welfare functions for project approval and status quo, 
				 respectively is not implemented, yet!")
	} else {
		stop("welfare must be either a function or a list of two functions.")
	}
	class(thisAnalysis) <- "welfareDecisionAnalysis"
	return(thisAnalysis)
}
##############################################################################################
# summary.welfareDecisionAnalysis(object, ...)
##############################################################################################
#' Summarize Welfare Decision Analysis results.
#' 
#'  Produce a summary of the results of a welfare decision analysis obtained by the function
#'  \ifelse{latex}{\cr}{ }\code{\link{welfareDecisionAnalysis}}.
#' @param object An object of class \code{welfareDecisionAnalysis}.
#' @param ... Further arguments passed to \code{\link{format}}.
#' @inheritParams base::format
#' @param probs \code{numeric vector}: quantiles that shall be displayed; if \code{=NULL} no 
#'   quantiles will be displayed.
#' @return An object of class \code{summary.welfareDecisionAnalysis}.
#' @seealso \code{\link{welfareDecisionAnalysis}}, 
#'  \code{\link{print.summary.welfareDecisionAnalysis}}, \code{\link{format}}
#' @export
summary.welfareDecisionAnalysis <- function(object,
																		 ...,
																		 digits = max(3, getOption("digits")-3),
																		 probs=c(0.05, 0.5, 0.95)){	
# 	summaryDf<-data.frame(enbPa=object$enbPa, 
# 												elPa=object$elPa, 
# 												elSq=object$elSq, 
# 												eol=object$eol, 
# 												optimalChoice=object$optimalChoice)	
summaryDf<-if(!is.null(probs)){
       data.frame(  t(apply(X=object$mcResult$y, MARGIN=2, FUN=quantile, probs=probs)),
                    enbPa=object$enbPa, 
                    elPa=object$elPa, 
                    elSq=object$elSq, 
                    eol=object$eol, 
                    optimalChoice=object$optimalChoice,
                    check.names=FALSE)
             }else{
	         data.frame(enbPa=object$enbPa, 
	                    elPa=object$elPa, 
	                    elSq=object$elSq, 
	                    eol=object$eol, 
	                    optimalChoice=object$optimalChoice)
  }
	
	summaryDf<-format(x=summaryDf, digits=digits, ...)
	# ToDo: combine this summary with summary(object$mcResult)
	res<-list(summary=summaryDf,
						call=object$call)
	
	class(res)<-"summary.welfareDecisionAnalysis"
	res
}
##############################################################################################
# print.summary.welfareDecisionAnalysis(x, ...)
##############################################################################################
#' Print the summarized Welfare Decision Analysis results.
#' 
#' This function prints the summary of a Welfare Decision Analysis generated by \ifelse{latex}{\cr}{ }
#' \code{\link{summary.welfareDecisionAnalysis}}.
#' @param x An object of class \code{summary.welfareDecisionAnalysis}.
#' @param ... Further arguments to \code{\link{print.data.frame}}.
#' @seealso \code{\link{welfareDecisionAnalysis}}, \code{\link{summary.welfareDecisionAnalysis}}, 
#'   \code{\link{print.data.frame}}.
#' @export
print.summary.welfareDecisionAnalysis <- function(x, ...){
	cat("Call:\n")
	print(x$call)
	cat("\nSummary of decision analysis:\n")
	print(x$summary,...)
}
##############################################################################################
# hist.welfareDecisionAnalysis(x, ...)
##############################################################################################
#' Plot Histogram of results of a Welfare Decision Analysis
#' 
#' This function plots the histograms of the results of
#' \code{\link{welfareDecisionAnalysis}}.
#' @param x An object of class \code{welfareDecisionAnalysis}.
#' @param xlab \code{character}: x label of the histogram. If it is not
#'   provided, i.e. equals \code{NULL} the name of the chosen variable by
#'   argument \code{resultName} is used.
#' @param main \code{character}: main title of the histogram.
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
#' @seealso \code{\link{welfareDecisionAnalysis}}, \code{\link{hist}}. For a list of colors
#'   available in R see \code{\link[grDevices]{colors}}.
#' @export
hist.welfareDecisionAnalysis <- function(x, breaks=100, col=NULL, xlab=NULL, main=paste("Histogram of " , xlab), ...,
                              colorQuantile   =c("GREY", "YELLOW", "ORANGE", "DARK GREEN", "ORANGE", "YELLOW", "GREY"), 
                              colorProbability=c(1.00,    0.95,     0.75,     0.55,         0.45,     0.25,     0.05),
                              resultName=NULL){

  hist(x$mcResult, breaks=breaks, col=col, xlab=xlab, main=main, ...,
       colorQuantile   =colorQuantile, 
       colorProbability=colorProbability,
       resultName=resultName)
}

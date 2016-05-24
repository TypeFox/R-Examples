################################################################
## Copyright 2014 Tracy Holsclaw.

## This file is part of NHMM.

## NHMM is free software: you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation, either version 3 of the License, or any later version.

## NHMM is distributed in the hope that it will be useful, but WITHOUT ANY
## WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
## A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## You should have received a copy of the GNU General Public License along with
## NHMM.  If not, see <http://www.gnu.org/licenses/>.
#############################################################



#' Calculates BIC, AIC, PLS, log-likelihood
#'
#' \code{OBIC} calculates BIC, AIC, approximate log-likelihood and plots the
#' log-likelihood for all iterations. The log-likelihood plot should 
#' be flat to show convergence to a stationary distribution. Minimize 
#' the AIC and BIC for the *best* model and maximize PLS.  The log likelihood is approximate in that
#' it is calculated by marginalizing over the current chain of hidden states instead
#' of using a recurrsive algorithm to compute it; every iterations produces an estimation
#' of the log-likelihood.  If yhold is provided the preditive log score (PLS) is also 
#' given. 
#' 
#' Predictive Log Score: mean(log( E(p(yhold|...)))  The expectation is over all of the
#' iterations of the algorithm.  And the mean is over the pT count of yhold. The scale of
#' the PLS is in the unit of t (usually days). 
#' 
#' @param nhmmobj an object created from the NHMM function
#' @param outfile a directory to put the .png plot
#' @return BIC  
#' @return output: AIC, BIC, PLS [if yhold data was provided], log-likelihood to the GUI and a plot of the log-likelihood
#' @export
#' @keywords AIC BIC log-likelihood PLS
#' @examples #OBIC(my.nhmm) 



OBIC=function(nhmmobj, outfile=NULL)
{
  T=nhmmobj$T
  J=nhmmobj$J
  K=nhmmobj$K
  B=nhmmobj$B
  A=nhmmobj$A
  iters=nhmmobj$iters
  burnin=nhmmobj$burnin 
  outboo=nhmmobj$outboo
  outdir=nhmmobj$outdir
  loglik=nhmmobj$loglik
  BICp=nhmmobj$BICp
  PLS=nhmmobj$PLS
  
  L=B+K
  
  
  BICf=-2*mean(loglik)+BICp*log(T)
  AICf=-2*mean(loglik)+2*BICp
  print(paste("Parameter count:  ",BICp,sep=""))
  print(paste("BIC:  ",round(BICf,2), " (favored method)",sep=""))
  print(paste("AIC:  ", round(AICf,2),sep=""))
  if(!is.null(PLS)){print(paste("PLS:  ",round(PLS,2)))}
  
  if(!is.null(outfile)){ png(paste(outfile,"BIC.png",sep=""), width=300, height=300)}
  plot(loglik, xlab="Iterations", ylab="Log-likelihood")
  if(!is.null(outfile)){ dev.off()}
  BICf
}  
  
  
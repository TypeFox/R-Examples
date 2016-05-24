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



#' Calculates the mean of the transition matrix 
#'
#' \code{OQQ} calculates the mean of the transition array which is a K by K by T.
#'  Each row is the state at time t-1 and the columns are time t. The Gui plot shows the 
#'  transition probability over time for each t-1 to t transition.
#' 
#' @param nhmmobj an object created from the NHMM function
#' @param outfile a directory to put the .png plot
#' @return QQmean the transition probabilities for each time step. [K by K by T] 
#' @return output: a plot where each row is the state at time t-1 and the columns are time t. 
#' The GUI plot shows the transition probability over time for each t-1 to 
#' t transition. If the columns are the same then the Markov property is probably weak.
#' @export
#' @keywords transition probabilities
#' @examples #OQQ(my.nhmm) 


OQQ=function(nhmmobj, outfile=NULL)
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
  L=B+K
 
   if(K==1){stop("K=1, there is no hidden state sequence")}
  
   meanQQ=nhmmobj$meanQQ
 
    if(!is.null(outfile)){ png(paste(outfile,"QQ.png",sep=""), width=K*150, height=K*150)}
     par(mfrow=c(K,K), mar=c(4,4,1,1))
     for(j  in 1:K)
     {   for(i in 1:K)
         {  plot(meanQQ[j,i,],ylim=c(0,1), ylab="Probability", xlab="T")
         }
     }
  if(!is.null(outfile)){ dev.off()}
  par(mfrow=c(1,1), mar=c(4,4,4,4))
  meanQQ
}

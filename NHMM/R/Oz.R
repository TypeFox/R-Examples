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



#' Most probable state (similar to Viterbi)
#'
#' \code{Oz} calculates the most probable state per time step (Viterbi like) with values from 1,...,K.
#'  The histogram of this sequence is displayed in the GUI output.   
#'  If there are ties for a given day then the lowest number state is chosen.
#' 
#' @param nhmmobj an object created from the NHMM function
#' @param outfile a directory to put the .png plot
#' @return zbest the most probable sequence from all iterations  
#' @return output: a plot of a histogram of the distribution of 
#' the most probable state sequence. If the number of states in 
#' the histogram is less than K, it probably means you should probably
#' re-run the model with smaller K as some of the states have disappeared.
#' @export
#' @keywords Viterbi
#' @examples #Oz(my.nhmm) 



Oz=function(nhmmobj, outfile=NULL)
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

    if(outboo==TRUE)  ##betasave was written to a file
    {      zsave=t(as.matrix(read.table(paste(outdir,"zsave.txt", sep=""))))  
    }else{   #zsave=matrix(0,T,iters)
           zsave=nhmmobj$zsave
    } 
    
    #### Most used z (like Viterbi)
    mode=function(vector_name)
    {  as.numeric(names(sort(table(vector_name),decreasing=TRUE))[1])  }

    zbest=numeric(T)
    for(t in 1:T){  zbest[t]=mode(zsave[t,])  }
    
    if(!is.null(outfile)){ png(paste(outfile,"z.png",sep=""), width=300, height=300)}
    hist(zbest, xlab="Most probable z")
    if(!is.null(outfile)){ dev.off()}
  
    zbest
}
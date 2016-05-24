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


#' Emission Parameters
#'
#' \code{Oemparams} calculates emission parameters 0.025, 0.05, mean, 0.50 (median),
#' 0.95, 0.975 quantiles from the iterations for each parameter. Each of the J sequences and K states
#' each has 2 parameters (Gamma and Normal) or 1 parameter (Poisson) per mixture component (nmix+delta)
#' 
#' 
#' @param nhmmobj an object created from the NHMM or HMM function
#' @param plots  TRUE/FALSE- default is FALSE because the plot window can grow quite large depending on the number of parameters.
#' [J by K*nmix] panes if outboo=TRUE and [K*nmix by J] panes if outboo=FALSE, 
#' where outboo is a parameter from NHMM which determines if output was written to a file (TRUE) or to a variable (FALSE).
#' If outfile is used then there will be K*J .png files containing dimension 2 parameters by nmix panes of trace plots.
#' Exception: NHMM_MVN, plots is always FALSE.
#' @param outfile a directory to put the .png plot
#' @return params [2 by nmix by K by J] by 6. There are six values returned: 0.025, 0.05, mean, 0.50 (median),
#' 0.95, 0.975 quantiles from the iterations. (0.025, 0.975) are used to construct 95% probability intervals (PIs),
#' likewise (0.05, 0.95) can be used to construct 90% PIs. The 2 parameter is for two-parameter distributions like Gamma
#' but the one-parameter distributions like the Exponential will only have meaningful data in the first row. 
#' Exception for NHMM_MVN then
#'  the parameters returned is the mean covariance matrix over all iterations [J by J by K]
#' @return output:  plot window can grow quite large depending on the number of parameters.
#' [J by K*nmix] panes if outboo=TRUE and [K*nmix by J] panes if outboo=FALSE, 
#' where outboo is a parameter from NHMM which determines if output was written to a file (TRUE) or to a variable (FALSE). 
#' Exception: NHMM_MVN object does not plot the covariance.
#'
#' @export
#' @keywords emission parameters
#' @examples 
#' #thetas=Oemparams(my.nhmm, FALSE); 
#' #thetas[,,,,3]  #mean values


Oemparams=function(nhmmobj, plots=FALSE, outfile=NULL)
{   T=nhmmobj$T
    J=nhmmobj$J
    K=nhmmobj$K
    B=nhmmobj$B
    A=nhmmobj$A
    iters=nhmmobj$iters
    burnin=nhmmobj$burnin 
    outboo=nhmmobj$outboo
    outdir=nhmmobj$outdir
    L=B+K
    
    if(iters<40){stop("Must have at least 40 iterations to get any results")}
    
    
    if(nhmmobj$fam==4)  #MVN
    {  if(outboo==TRUE)
       {   thetasave=array(0,dim=c(J,J,K))
           for(k in 1:K)
           { thetasave[,,k]=as.matrix(read.table(paste(outdir,"K",k,"-theta.txt", sep=""))) ## K*L by iters
           }
           thetasave
       }else{
         thetasave=nhmmobj$thetasave
         thetasave
       }
      
      
      
    }else{
      
     #not MVN
    if(outboo==TRUE)  ##psisave was written to a file
    {   first=as.matrix(read.table(paste(outdir,"K",1,"-J",1,"-theta1.txt", sep=""))) ## K*L by iters
        nmix=dim(first)[2]
        
        betar=array(0,dim=c(2,nmix,K,J,6))  #2,nmix,K,J,iters
        thetasaveKJ=matrix(0,nmix,iters)  #c(K+A,J,iters)
        
        
        for(j in 1:J)
        {  for(k in 1:K)
           {   if(!is.null(outfile)){ png(paste(outfile,"K",k,"J",j,"params.png",sep=""), height=2*220, width=nmix*220)}
               par(mfrow=c(2,nmix), mar=c(4,4,1,1))
               for(tt in 1:2)  #theta1 and theta2
               {  thetasaveKJ=as.matrix(read.table(paste(outdir,"K",k,"-J",j,"-theta",tt,".txt", sep=""))) ## K*L by iters
                  for(i in 1:nmix)
                  {  if(plots==TRUE)
                     {   plot(thetasaveKJ[,i], ylab=paste("K=",k," J=",j," Th=",tt,"nmix=",i, sep=""), xlab="Iterations")
                     } #######  .025, .05, mean, median, .95, .975
                      one=sort(thetasaveKJ[,i])
                      betar[tt,i,k,j,1]=one[.025*iters]
                      betar[tt,i,k,j,2]=one[.05*iters]
                      betar[tt,i,k,j,3]=mean(one)
                      betar[tt,i,k,j,4]=one[.50*iters]
                      betar[tt,i,k,j,5]=one[.95*iters]
                      betar[tt,i,k,j,6]=one[.975*iters]  
                   }
               }
               if(!is.null(outfile)){ dev.off()}
           }
        }
 
        
    }else{  #2,nmix,K,J,iters
      thetasave=nhmmobj$thetasave
      nmix=dim(thetasave)[2]
      
      if(plots==TRUE)
      { if(!is.null(outfile)){ png(paste(outfile,"emparams.png",sep=""),height=2*K*nmix*120,  width=J*120)}
        par(mfrow=c(2*K*nmix,J), mar=c(1,1,1,1))
        for(j in 1:J)
        {   for( i in 1:nmix)
            {   for(k in 1:K)
                {  plot(thetasave[1,i,k,j,], ylab=paste("1-K=",k," J=",j,"nmix=",i, sep=""),xlab="Iterations")
                   plot(thetasave[2,i,k,j,], ylab=paste("2-K=",k," J=",j,"nmix=",i, sep=""),xlab="Iterations")
                }
            }       
         }
       }
      if(!is.null(outfile)){ dev.off()}
      
      ##############################  95%  confidence for the Betas 
      #######  .025, .05, mean, median, .95, .975
      #####  95% and 90% confidence
      
      betar=array(0,dim=c(2,nmix,K,J,6))  #2,nmix,K,J,iters
      for(j in 1:J)
      {   for(k in 1:K)
          {   for(i in 1:nmix)
              {  for(tt in 1:2)
                  {  one=sort(thetasave[tt,i,k,j,])
                     betar[tt,i,k,j,1]=one[.025*iters]
                     betar[tt,i,k,j,2]=one[.05*iters]
                     betar[tt,i,k,j,3]=mean(one)
                     betar[tt,i,k,j,4]=one[.50*iters]
                     betar[tt,i,k,j,5]=one[.95*iters]
                     betar[tt,i,k,j,6]=one[.975*iters]  
                  }
              }
           }
       }
    } 
    par(mfrow=c(1,1), mar=c(4,4,4,4))
    betar  #2,nmix,K,J,6   (theta1, theta2)  
   }
}





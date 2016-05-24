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



#' Coefficients for X (transition inputs)
#'
#' \code{OXcoef} calculates transition coefficient 0.025, 0.05, mean, 0.50 (median),
#' 0.95, 0.975 quantiles from the iterations. Each of K-1 states (K-1 for identifiability)
#'  has K Markov coefficients and B input coeffients.
#'  
#' 
#' @param nhmmobj an object created from the NHMM function
#' @param plots  TRUE/FALSE- default is FALSE because the plot window can grow quite large depending on the number of parameters.
#' @param outfile a directory to put the .png plot
#' @return params [K-1 by K+B] by 6. There are six values returned: 0.025, 0.05, mean, 0.50 (median),
#' 0.95, 0.975 quantiles from the iterations. (0.025, 0.975) are used to construct 95% probability intervals (PIs),
#' likewise (0.05, 0.95) can be used to construct 90% PIs. 
#' @return output:  plot window can grow quite large depending on the number of states.
#' [K-1 by K+B] panes.
#' @return output: outputs statements of 90% and 95% significance of the Markov coefficients 
#' and each X input coefficients (if any of the K-1 coefficients for a variable are significant then that
#' variable is deemed signficant.)
#' @export
#' @keywords coefficients
#' @examples #thetas=OXcoef(my.nhmm, FALSE); 
#' #thetas[,,,,3]  #mean values


OXcoef=function(nhmmobj, plots=FALSE, outfile=NULL)
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
    
    
    if(outboo==TRUE)  ##betasave was written to a file
    {   betasaveK=matrix(0,L,iters)
        if(!is.null(outfile) && plots==TRUE){ png(paste(outfile,"beta.png",sep=""), height=(K-1)*120, width=L*120)}
        par(mfrow=c(K-1,L), mar=c(4,4,1,1))
        betar=array(0,dim=c(K-1,L,6))
        
        for(k in 1:(K-1))
        {  betasaveK=read.table(paste(outdir,"K-",k,"-beta.txt", sep="")) ## K*L by iters
            
           if(plots==TRUE)
           { for( i in 1:(L))
             {  if(i <= K)#markov X
               {  plot(betasaveK[,i], ylab=paste("K=",k," and Markov coeff ",i, sep=""),xlab="Iterations")
               }else{
                  plot(betasaveK[,i], ylab=paste("K=",k," and X ",i-K, sep=""),xlab="Iterations")
               }
               abline(h=0,col=2)
             }
           }
        
        ##############################  95%  confidence for the Betas 
        #######  .025, .05, mean, median, .95, .975
        #####  95% and 90% confidence

           for(j in 1:L)
           {     one=sort(betasaveK[,j])
                 betar[k,j,1]=one[.025*iters]
                 betar[k,j,2]=one[.05*iters]
                 betar[k,j,3]=mean(one)
                 betar[k,j,4]=one[.50*iters]
                 betar[k,j,5]=one[.95*iters]
                 betar[k,j,6]=one[.975*iters]    
            }
        }
        if(!is.null(outfile)){ dev.off()}
        
        
    }else{  #betasave=array(0,dim=c(K,L,iters))
      betasave=nhmmobj$betasave
      
      if(plots==TRUE)
      {  if(!is.null(outfile) && plots==TRUE){ png(paste(outfile,"beta.png",sep=""), height=(K-1)*120, width=L*120)}
         par(mfrow=c(K-1,L), mar=c(4,4,1,1))
         for(j in 1:(K-1))
         {   for( i in 1:(L))
             {  if(i <= K)#markov X
                { plot(betasave[j,i,], ylab=paste("K=",j," and Markov coeff ",i, sep=""),xlab="Iterations")
                }else{
                  plot(betasave[j,i,], ylab=paste("K=",j," and X ",i-K, sep=""),xlab="Iterations")
                }
                abline(h=0,col=2)
             }
          }
         if(!is.null(outfile)){ dev.off()}
      }
      ##############################  95%  confidence for the Betas 
      #######  .025, .05, mean, median, .95, .975
      #####  95% and 90% confidence
      betar=array(0,dim=c(K-1,L,6))
      for(j in 1:L)
      {   for(k in 1:(K-1))
          {   one=sort(betasave[k,j,])
              betar[k,j,1]=one[.025*iters]
              betar[k,j,2]=one[.05*iters]
              betar[k,j,3]=mean(one)
              betar[k,j,4]=one[.50*iters]
              betar[k,j,5]=one[.95*iters]
              betar[k,j,6]=one[.975*iters]    
          }
      }
    }
    
    
    
  
    
    #### each of the L sets of coefficients need to be examined for zeros
    signif1=matrix(TRUE,K-1,K)    #95% - first K are all about Markov and only need one indicator
    signif3=matrix(TRUE,K-1,K)    #90% - first K are all about Markov and only need one indicator
    
    #Is the Markov significant?
    for(i in 1:K) #first K coefficients
    {   for(j in 1:(K-1))
        {   if(betar[j,i,1] <0 & betar[j,i,6]>0)
            {  signif1[j,i]=FALSE
            }
            if(betar[j,i,2] <0 & betar[j,i,5]>0)
            {  signif3[j,i]=FALSE
            }
        }
    }
    ### Print Markov results 
    if(sum(signif1)>0)  #if at leasts 1 beta is significant then the Markov property is significant
    {   print("The Markov property of the model is significant with 95% probability.")   
    }else{
      print("The Markov property of the model is not significant with 95% probability.")
    }
    if(sum(signif3)>0)  #if at leasts 1 beta is significant then the Markov property is significant
    {   print("The Markov property of the model is significant with 90% probability.")   
    }else{
      print("The Markov property of the model is not significant at the 90% probability level.")
    }
    
    signif2=matrix(TRUE,K-1,L-K)  #95% - X variable significance
    signif4=matrix(TRUE,K-1,L-K)  #90% -  X variable significance
     
    for(i in (K+1):L)
    {   for(j in 1:(K-1))
        {   if(betar[j,i,1] <0 & betar[j,i,6]>0)
            {   signif2[j,i-K]=FALSE
            }
            if(betar[j,i,1] <0 & betar[j,i,6]>0)
            {   signif4[j,i-K]=FALSE
            }
        }
    }
    truthsL=apply(signif2,2,sum) 
    truthsL4=apply(signif4,2,sum) 
    for(k in 1:(L-K))
    {  if(truthsL[k]>0) #any of the K-1 coefficients for the variable are signif
       {  print(paste("Input",k,": X is significant with 95% probability",sep=""))
       }else{
       print(paste("Input",k,": X is not significant with 95% probability",sep=""))
       }
    }
    for(k in 1:(L-K))
    {  if(truthsL4[k]>0) #any of the K-1 coefficients for the variable are signif
       {  print(paste("Input",k,": X is significant with 90% probability",sep=""))
       }else{
          print(paste("Input",k,": X is not significant with 90% probability",sep=""))
       }
    }
    par(mfrow=c(1,1), mar=c(4,4,4,4))
    betar  #.025,.05, median (.50), mean, .95,.975   K by K+B by 6
    
}


#my.betar=piXcoef(my.nhmm)  


######################################################################

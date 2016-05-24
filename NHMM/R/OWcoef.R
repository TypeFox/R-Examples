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



#' Coefficients for W (emission inputs)
#'
#' \code{OWcoef} calculates emission coefficient 0.025, 0.05, mean, 0.50 (median),
#' 0.95, 0.975 quantiles from the iterations. These coeffiecients are 
#' used to determine the mixing weights of the mixture components of the emission distributions.
#'  There are K intercept coefficients and A input coeffients for each J.
#'  
#' 
#' @param nhmmobj an object created from the NHMM function
#' @param plots  TRUE/FALSE- default is FALSE because the plot window can grow quite large depending on the number of sequences.
#' @param outfile a directory to put the .png plot
#' @return params [K+A by J] by 6. There are six values returned: 0.025, 0.05, mean, 0.50 (median),
#' 0.95, 0.975 quantiles from the iterations. (0.025, 0.975) are used to construct 95% probability intervals (PIs),
#' likewise (0.05, 0.95) can be used to construct 90% PIs. 
#' @return output:  plot window can grow quite large depending on the number of sequences.
#' [K+A by J] panes.'
#' @return output: outputs statements of 90% and 95% significance of the intercept coefficients 
#' and each X input coefficients (if any of the J coefficients for a variable are significant then that
#' variable is deemed signficant.)
#' @export
#' @keywords coefficients
#' @examples #thetas=OWcoef(my.nhmm, FALSE); 
#' #thetas[,,,,3]  #mean values



OWcoef=function(nhmmobj, plots=FALSE, outfile=NULL)
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
    
    
    
    if(outboo==TRUE)  ##psisave was written to a file
    {   if(!is.null(outfile) && plots==TRUE){ png(paste(outfile,"psi.png",sep=""), height=J*120, width=(K+A)*120)}
        betar=array(0,dim=c(K+A,J,6))
        par(mfrow=c(J,K+A), mar=c(1,1,1,1))
     
        for(j in 1:J)
        {  psisaveJ=matrix(0,K+A,iters)  #c(K+A,J,iters)
           psisaveJ=as.matrix(read.table(paste(outdir,"J-",j,"-psi.txt", sep=""))) ## K*L by iters
           
           if(plots==TRUE)
           {  for(i in 1:(K+A))
              {   if(i <= K)   #intercept terms W
                  {   plot(psisaveJ[,i]  ,xlab="Iterations", ylab=paste("K=",i," and J=",j, sep=""))
                  }else{
                      plot(psisaveJ[,i], ylab=paste("K=",i-K," and J=",j, sep=""), xlab="Iterations")
                  }
                abline(h=0,col=2)
              }
           }
        
        ##############################  95%  confidence for the Betas 
        #######  .025, .05, mean, median, .95, .975
        #####  95% and 90% confidence
          for(k in 1:(K+A))
          {  one=sort(psisaveJ[,k])
             betar[k,j,1]=one[.025*iters]
             betar[k,j,2]=one[.05*iters]
             betar[k,j,3]=mean(one)
             betar[k,j,4]=one[.50*iters]
             betar[k,j,5]=one[.95*iters]
             betar[k,j,6]=one[.975*iters]    
          }
        }
        if(!is.null(outfile)){ dev.off()}
        
    }else{  #psisave=array(0,dim=c(K+A,J,iters))
        psisave=nhmmobj$psisave
      
        if(plots==TRUE)
        { if(!is.null(outfile)){ png(paste(outfile,"psi.png",sep=""), height=(K+A)*120, width=J*120)}
          par(mfrow=c(K+A,J), mar=c(1,1,1,1))
          for(j in 1:(K+A))
          {   for( i in 1:(J))
              {   if(i <= K)#markov X
                  {   plot(psisave[j,i,], ylab=paste("K=",j," and J= ",i, sep=""),xlab="Iterations")
                  }else{
                      plot(psisave[j,i,], ylab=paste("K=",j-K," and J= ",i, sep=""),xlab="Iterations")
                  }
                  abline(h=0,col=2)
              }
          }
        }
        ##############################  95%  confidence for the Betas 
        #######  .025, .05, mean, median, .95, .975
        #####  95% and 90% confidence
        
        betar=array(0,dim=c(K+A,J,6))
        for(j in 1:J)
        {   for(k in 1:(K+A))
        {    one=sort(psisave[k,j,])
             betar[k,j,1]=one[.025*iters]
             betar[k,j,2]=one[.05*iters]
             betar[k,j,3]=mean(one)
             betar[k,j,4]=one[.50*iters]
             betar[k,j,5]=one[.95*iters]
             betar[k,j,6]=one[.975*iters]    
        }
        }
        if(!is.null(outfile)){ dev.off()}
    } 
    

    

    
    #### each of the L sets of coefficients need to be examined for zeros
    signif1=matrix(TRUE,K,J)    #95% - first K are all about Markov and only need one indicator
    signif3=matrix(TRUE,K,J)    #90% - first K are all about Markov and only need one indicator
    
    #Is the Markov significant?
    for(i in 1:J) #first K coefficients
    {   for(j in 1:K)
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
    {   print("The intercepts of the model are significant with 95% probability.")   
    }else{
      print("The intercepts of the model are not significant with 95% probability.")
    }
    if(sum(signif3)>0)  #if at leasts 1 beta is significant then the Markov property is significant
    {   print("The intercepts of the model are significant with 90% probability.")   
    }else{
      print("The intercepts of the model are not significant at the 90% probability level.")
    }
    
    
    if(A>0)
    {  signif2=matrix(TRUE,A,J)  #95% - X variable significance
       signif4=matrix(TRUE,A,J)  #90% -  X variable significance
    
       for(i in 1:J)
       {   for(j in (K+1):(K+A))
           {   if(betar[j,i,1] <0 & betar[j,i,6]>0)
               {   signif2[j-(K+1),i]=FALSE
               }
               if(betar[j,i,1] <0 & betar[j,i,6]>0)
               {   signif4[j-(K+1),i]=FALSE
               }
           }
       }
       truthsL=apply(signif2,2,sum) 
       truthsL4=apply(signif4,2,sum) 
       for(k in 1:A)
       {  if(truthsL[k]>0)
          {  print(paste(k,". X is significant with 95% probability",sep=""))
          }else{
            print(paste(k,". X is not significant with 95% probability",sep=""))
          }
       }
       for(k in 1:(L-K))
       {  if(truthsL4[k]>0) 
          {  print(paste(k,". X is significant with 90% probability",sep=""))
          }else{
             print(paste(k,". X is not significant with 90% probability",sep=""))
          }
       }
    }
    
    par(mfrow=c(1,1), mar=c(4,4,4,4))
    betar  #.025,.05, median (.50), mean, .95,.975   K+A by J by 6
    
}

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



#' Bayesian Non-homogeneous Markov Model with Multivariate Normal emission distribution (NHMMMVN)
#'
#' \code{NHMM_MVN} calculates an NHMM for multiple sequences of data. 
#' The sequences can actually be short sets of equal length sequences (subseq).
#' The traditional input variables (X) influence the non-homogenous transition probabiities
#' of the model. An additional set of input variables (W) can be included to influence the mean
#' of the emission distribution. All parameters are sampled
#' via Gibbs steps (latent variables such that no tuning is needed.) The X variable coefficients
#' are sampled through an unordered Multinomial logit model Polya-Gamma formulation (Polson, Scott, Windle). 
#' The hidden states are sampled through a blocked Gibbs sampler.  
#' 
#'
#' @param y   T by J matrix of data (J=1 is sufficient)
#'          -  missing data is denoted with NA
#' @param subseq  [optional] if y is actually a set of subsequences then give the length of those 
#'         sequences (122 for JJAS) (365 is not it!). Default is subseq=T.
#'
#' @param X       B by T matrix for the transition input data (B different inputs)  
#'              Missing values are not allowed. If there are no Xs then use HMM function.
#' @param betapriorm [optional]  default=NULL which is reference prior. Or a [K+B by K matrix] for
#'            the mean of the Normal prior for the beta coefficients.
#' @param betapriorp [optional]  default=NULL which is reference prior. Or a [K+B by K+B by K array] for
#'            the precision prior(1/sig^2) of the Normal prior for the beta coefficients.
#'
#' @param K       number of states (default=2)
#' @param iters   number of iterations to keep after burn in  (default=1000)
#' @param burnin the number of burn in  (default=200)
#' 
#' @param W    [optional] is an A by T by J array of emission input data (A different inputs), 
#'          missing values are not allowed, do not include an intercept term here
#'          The mean function depends on W. 
#' @param psipriorm [optional]  default=NULL which is reference prior. Or a [K+A by J matrix] for
#'           the mean of the Normal prior for the beta coefficients.
#' @param psipriorp [optional]  default=NULL which is reference prior. Or a [K+A by J matrix] for
#'           the precision prior(1/sig^2) of the Normal prior for the beta coefficients.


#' @param priors1  [optional] scale parameter for each state (vector of length K) of the Wishart prior for Sigma 
#'                 priors1=NULL is set to 1 
#' @param priors2  [optional] Covariance matrix [J by J by K] for the MVN for each state, second parameter
#'                 of the Wishart for the Sigma matrix.
#'                  priors2=NULL  will be a diagonal matrix with diagonal set to 1
#'            

#' @param outdir [optional] can output each set of parameters to output files in a directory 
#'         use this with larger dimension data sets or with large number of iterations 
#'         the output will be written line by line and not overburden the memory limit
#'         outdir needs to end with a slash or double slash depending on OS.
#'
#'
#' @param ymiss [optional-TRUE/FALSE] if outdir is specified then draws for any missing
#'                data points will be saved to ymiss-J*.txt. There will be one data point for
#'                each sequence of length Tmiss (the amount missing from that sequence).
#'                Each row of the output file will be an iteration of the algorithm after burnin is
#'                removed.    
#'                
#' @param  yrep  [optional] number (ie. 100,200,or 500) of output replicate data sets 
#'         to print to outdir. The replicates will be the same dimension as y and start
#'         after burn in period. default is zero. These replicate data sets are 
#'         generated from the same input variables.
#'         - must be shorter than iters
#'        
#' @param ypred  number of predictive sets (ie. 0,100, 200, 500)
#'         - must be shorter than iters
#'        Predicted chains: will produce ypred set of predictives [pT by J] to print to outdir.
#'        yrep uses the same inputs to make replications but ypred uses new input 
#'        values of X and W to make predictions over a different time span (pT) 
#'        Also outputs a set of predictive z values or length (pT by ypred)
#' @param  Xp    predictive set of Xs of length pT  [B by pT]
#'        -missing values is not allowed
#'        -ensure that the Xp inputs are in the same order as X
#' @param  Wp    predictive set of Ws of length pT [A by pT by J]
#'        -missing valus are not allowed
#'        -ensure that the Wp inputs are in the same order as W
#'@param yhold [optioinal] a sequence of y observed values [pT by J], held out data that is 
#'           of length ypred that is used to compute  the predictive log score (PLS)
#'           which is a metric like BIC (ie. hold out last 10% of data or do 5-fold CV). 
#'           missing data values are filled in with mean PLS value. 
#'           PLS is the average PLS across sequences.      
#'
#' @return my.nhmm object  
#' @export
#' @keywords Bayesian NHMM
#' @examples ## Multivariate Normal data
#' \dontrun{ 
#' data(NHMMdata)
#' attach(NHMMdata)
#' 
#' my.nhmm1=NHMM_MVN(y=ymvn,  X=tX, W=tW2, K=3, iters=50, burnin=10, 
#'                  priors1=rep(2,3))
#' OBIC(my.nhmm1)
#' zz=Oz(my.nhmm1)  #compare with the truth zgamma
#' qq=OQQ(my.nhmm1)
#' bb=OXcoef(my.nhmm1)
#' pp=OWcoef(my.nhmm1,FALSE)
#' tt=Oemparams(my.nhmm1,FALSE)  #just Sigma matrix for MVN, returns mean of Sigma
#'  
#'  
#'  #filelocation="C:\\Users\\iamrandom\\Desktop\\here\\"
#'  #my.nhmm6=NHMM_MVN(y=ymvn[1:1800,],   X=matrix(tX[,1:1800],1,1800), 
#'  #          W=array(tW2[,1:1800,],dim=c(2,1800,15)), K=3, iters=50, 
#'  #            burnin=10,outdir=filelocation, ymiss=TRUE, yrep=10, 
#'  #            Xp=matrix(tX[,1801:2000],1,200), 
#'  #           Wp=array(tW2[,1801:2000,],dim=c(2,200,15)), ypred=10, 
#'  #           yhold=ymvn[1801:2000,])
#'  #OBIC(my.nhmm6)
#'  
#'  #Could try it with K=1, to compare K=1 to K=3
#'  }




NHMM_MVN=function(y, subseq=NULL, X=NULL, betapriorm=NULL, betapriorp=NULL, K=2, iters=1000, burnin=200,  W=NULL, psipriorm=NULL,psipriorp=NULL, priors1=NULL, priors2=NULL, outdir=NULL, ymiss=FALSE, yrep=0 ,ypred=0, Xp=NULL, Wp=NULL, yhold=NULL)
{  
  

    Xnull=is.null(X)
    
    if(is.null(outdir) && yrep>0)
    { stop("Please specify outdir if you want yrep=TRUE.")}
    
    

    
    if(iters>20000){stop("Use outdir parameter, not enough memory to store all of the output")}
    if(!is.null(outdir)){outboo=TRUE}else{outboo=FALSE}  #there is an output file

    if(!is.matrix(y)==TRUE){stop("y needs to be a matrix")}
    T=dim(y)[1]
    if(T<=0){stop("No y data added")}
    J=dim(y)[2]
   
    if(T<=4){stop("y sequences are too short")}
    if(T<=15){warning("y sequences are quiet short")}  #15 is an arbitrary number
    
  

    ### missing data
    yboo=is.na(y)  #finds missing values
    y[yboo]=mean(y)  #simple imputation of missings to mean
    if(sum(yboo)/(T*J) > .50){ warning("Over 50% of your data is missing. Results may be questionable.")}
    y[yboo]=mean(na.omit(y))  #fill in missingness


    
    if(K >= 0)
    {  if(K==0){K=1}                   #K=0 and K=1 are treated the same
    }else{stop("Invalid choice for K")}
    
    if(T < K){stop("Data sequence is shorter than the number of states")}
    
    
    #######  X
    if(is.null(X))
    {  stop("No X inputs provided, use HMM or MVHMM instead; transitions are homogeneous")
    }
    if(sum(is.na(X))>0){stop("X has NA values; this is not allowed")}
    
    B=dim(X)[1]
    if(B < 1){stop("Need to provide X data")}
    
    
    ######## W 
    if(is.null(W)){  A=0;               
    }else{  
            if(dim(W)[1]<=0){stop("Provide W data or leave it set to NULL")}
            if(dim(W)[2] != T || dim(W)[3] != J){stop("W is not the right dimension")}
            if(sum(is.na(W))>0){stop("W has NA values; this is not allowed")}
            A=dim(W)[1]
    }
    
    ###########################  subseq #########################################
    if(is.null(subseq)){subseq=T}
    if(subseq > 0)
    {  if(T%%subseq !=0)
    {  stop("subseq is incorrect and the length of y cannot be split evenly into subseq pieces")
    }else{  subseqy=rep(1:(T/subseq),each=subseq)
    }
    }
    subboo=rep(0,T)
    for(t in 1:(T-1))
    {  if(subseqy[t]!=subseqy[t+1]){ subboo[t]=1}
    }
    if(subseq <=4){stop("subsequences are too short")}
    if(subseq <=15){warning("subsequences are quite short")}
    if(subseq<K){stop("subsequences are shorter than K")}
    ############################################################################
    

    if(burnin < 0){stop("Use different amount of burnin")}
    if(iters <0){stop("Use different iterations")}
    
    ### Set the priors for betas (B > 0)
    if(is.null(betapriorm) )         #reference prior is all zeros 
    {  betapriorm=matrix(0,K+B,K)  #mean for the B input variables and  K states
    }
    if(is.null(betapriorp))         #reference prior is all zeros 
    {  betapriorp=array(0, dim=c(K+B,K+B,K))  #precision for the B input variables each have
    }
    if(is.null(psipriorm) )         #reference prior is all zeros 
    {  psipriorm=matrix(0,A+K,J)  #mean for the A input variables
    }  
    if( is.null(psipriorp))         #reference prior is all zeros 
    {   psipriorp=matrix(0,A+K,J)  #precision for the A input variables 
    }  
    ###################################################  
    ### set up hidden states (latent variable z) ###########
    z=numeric(T)    ### initialize z by roughly sorting y into K bins
    yy=apply(y,1,sum)  
    z=cutree(hclust(dist(yy)), k=K)
    sumy=numeric(K)
    for(k in 1:K){ sumy[k]=sum(z==k)}
    if(is.element(0,sumy)==0){z[order(yy)]=c(rep(1:K,each=floor(T/K)),rep(K,T%%K))}  
    
    

    
    

    
    theta=array(0,dim=c(J,J,K))
    for(k in 1:K)
    { theta[,,k]=var(y) #initialize
    }
                    
                    
        if(is.null(priors1))   #non-informative prior
        {   priors1=rep(1,K)
        }
        if(is.null(priors2))   #non-informative prior
        {      priors2=array(0,dim=c(J,J,K))
               for(k in 1:K)
               {  priors2[,,k]=diag(J)*1  #not precision
               }
        }
   
      if( dim(priors2)[1]!=J){stop("Wrong first dimension of *priors2* input")}
      if( dim(priors2)[2]!=J){stop("Wrong second dimension of *priors2* input")}
      if( dim(priors2)[3]!=K){stop("Wrong first dimension of *priors2* input")}
    
        if(sum(priors1<= 0) >0)
        {   stop("priors must be greater than zero")}
        
    
   


   
    if(ypred> iters ||  yrep > iters){stop("iters must be larger than both ypred or yrep.")}
    if(ypred>0 || !is.null(yhold))  #Wp: A pT J,  Xp:  B pT
    { pT=0

      if(!is.null(Xp)){ pT=dim(Xp)[2]}else{stop("For ypred greater than zero, there must be an Xp.")}
      if(pT<1){stop("Xp is not long enough (must be greater than zero)")}
      if(!is.null(Wp))
      { if(dim(Wp)[2] != dim(Xp)[2]){stop("Wp and Xp are not the same length")} 
        if(A != dim(Wp)[1]){stop("A for W and A for Wp are different")}
        if(J != dim(Wp)[3]){stop("J for Wp is wrong")} 
        if(sum(is.na(Wp))>0){stop("W has NA values; this is not allowed")}
      }
      if(B != dim(Xp)[1]){stop("B for X and B for Xp are different")}
      if(sum(is.na(Xp))>0){stop("W has NA values; this is not allowed")}
      if(!is.null(Xp)){pT=dim(Xp)[2]}else{stop("Xp is needs data")}
    }
    
    if(!is.null(yhold))
    {  if(dim(yhold)[1]!=pT | dim(yhold)[2]!=J){stop("yhold is not the correct dimensions.")}
    } 
    
    
    if(ymiss==TRUE && is.null(outdir)){stop("must specify outdir to have ymiss=TRUE")}
    
    

    NHMM_MVNmain( z, theta, y, yboo, subseqy, subboo, X, betapriorm, betapriorp, K, iters,  burnin,   W, psipriorm, psipriorp,  priors1, priors2, outdir,outboo,  yrep, Xp, Wp, ypred, yhold, ymiss)
    
}





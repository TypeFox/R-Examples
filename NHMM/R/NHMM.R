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


#' Bayesian Non-homogeneous Markov Model (NHMM)
#'
#' \code{NHMM} calculates an NHMM for multiple sequences of data. 
#' The sequences can actually be short sets of equal length sequences (subseq).
#' The traditional input variables (X) influence the non-homogenous transition probabiities
#' of the model. An additional set of input variables (W) can be included to influence the mixture 
#' proportions of the emission distributions. The NHMM follows the general weather
#' state formulation of Hughes and Guttorp but in a Bayesian fashion. All parameters are sampled
#' via Gibbs steps (latent variables such that no tuning is needed.) The W variable coefficients
#' are sampled through an ordered Mulinomial probit (Albert and Chib). The X variable coefficients
#' are sampled through an unordered Multinomial logit model Polya-Gamma formulation (Polson, Scott, Windle). 
#' The hidden states are sampled through a blocked Gibbs sampler.  
#' 
#'
#' @param y        T by J matrix of data (J=1 is sufficient)
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

#' @param emdist   emission distribution:  "normal", "poisson", "gamma"  
#'           actual choices are Normal, Poisson,  Gamma,  Exponential,  
#'          or finite mixtures of mixtures or zero inflated version of any of these
#' @param nmix    [optional] number of mixture components for emdist, default is one (do not include delta)
#' @param delta   [TRUE/FALSE] TRUE-if we are using a zero inflated distribution
#'         (adds a delta function at zero as the first mixture component)
#' @param W       [optional] is an A by T by J array of emission input data (A different inputs), 
#'          missing values are not allowed, do not include an intercept term here
#'           The mixture components of emission depend on W.
#' @param psipriorm [optional]  default=NULL which is reference prior. Or a [K+A by J matrix] for
#'           the mean of the Normal prior for the beta coefficients.
#' @param psipriorp [optional]  default=NULL which is reference prior. Or a [K+A by J matrix] for
#'           the precision prior(1/sig^2) of the Normal prior for the beta coefficients.


#' @param priors  priors for emission components, each state can have a different prior
#'         dimension 5 by nmix by K by J  (some of the 5 dimensions are filled with zeros for some distributions)
#'         Normal(mu,sig2) or reference prior: mu~Normal(a.mu, 1/b.mu^2), sig2~IG(a.sig,b.sig)
#'            -the first 4 rows are used [a.mu, b.mu, a.sig, b.sig] but the 5th row is not used but must be present. 
#'            -reference prior: priors=NULL  
#'            
#'   	    Poisson(lambda):  lambda~Gamma(a.lam,b.lam)  
#'            -the first 2 rows are used [a.lam, b.lam] but the 3rd-5th rows are not used but must be present. 
#'
#'         Gamma(alpha,beta): alpha~Gamma(a.alpha,b.alpha) beta~Gamma(a.beta,b.beta) 
#'               -if the second row is supplied as NA then the first row supplies fixed first parameters for the Gamma
#'                 this is the suggested parametrization of the Gamma. Estimating the first parameter of the Gamma distribution
#'                 switches to an MCMC that is not as stable. 
#'               -to create an Exponential distribution the first row is set to 1 and the second is set to NA
#'                 then the third and fourth rows contain the (a.lam, b.lam) priors for the parameter of the Exponential
#'                 if prior=NULL then these are set to ones. setting these parameters to zeros is not a good idea.
#'               -a non-informative prior may not work, use a matrix of ones for a lowly informative prior instead of zeros.
#'               -the 5th row is for the tuning parameter for the MCMC (try 0.5) if a two parameter Gamma is desired.
#'                 A 2 parameter Gamma needs an informative prior and results may be more unstable.
#'               -default prior=NULL is the Exponential distribution with lowly informative prior of ones
#'         

#' @param outdir [optional] can output each set of parameters to output files in a directory 
#'         use this with larger dimension data sets or with large number of iterations 
#'         the output will be written line by line and not overburden the memory limit
#'         outdir needs to end with a slash or double slash depending on OS. 
#'         
#'  @param ymiss [optional-TRUE/FALSE] if outdir is specified then draws for any missing
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
#'         - must be shorter than iters.
#'          Predicted chains: will produce ypred sets of predictives [pT by J] to print to outdir.
#'        yrep uses the same inputs to make replications, this uses new inputs
#'        to make predictions over a different time span (pT)
#'        Also outputs a set of predictive z values or length (pT by ypred)
#' @param  Xp    predictive set of Xs of length pT  [B by pT]
#'        -missing values is not allowed
#'        -ensure that the Xp inputs are in the same order as X
#' @param  Wp    predictive set of Ws of length pT [A by pT by J]
#'        -missing valus are not allowed
#'        -ensure that the Wp inputs are in the same order as W
#'        
#' @param yhold [optioinal] a sequence of y observed values [pT by J], held out data that is 
#'           of length pT that is used to compute  the predictive log score (PLS)
#'           which is a metric like BIC (ie. hold out last 10% of data or do 5-fold CV). 
#'           missing data values are filled in with mean PLS value. 
#'           PLS is the average PLS across sequences.
#'             
#' @return my.nhmm object  
#' @export
#' @keywords Bayesian NHMM
#' @examples ## Gamma or Exponential
#' ### if "priors" is not specified, this is an Exponetial distribution
#' data(NHMMdata)
#' attach(NHMMdata)
#' 
#' ## Set to iters=40 for example only this should be in the thousands
#' my.nhmm=NHMM(y=ygamma[1:200,1:3],  X=matrix(tX[,1:200],1,200),
#'     K=3, iters=40, burnin=2, emdist="gamma", nmix=3, delta=TRUE)
#'           
#' OBIC(my.nhmm)
#' Oz(my.nhmm)  #compare with the truth: tz1
#' OQQ(my.nhmm) #transition probabilities 
#' \dontrun{
#' bb=OXcoef(my.nhmm)
#' pp=OWcoef(my.nhmm,FALSE)
#' tt=Oemparams(my.nhmm,FALSE)
#'  
#'  ## Normal - X is not used to create this data, so it should not be significant
#'  my.nhmm2=NHMM(y=ynormal, subseq=1000, X=tX, K=3, iters=100, 
#'            burnin=10, emdist="normal", nmix=2, delta=FALSE)
#'  OBIC(my.nhmm2)
#'  
#'  ## Poisson
#'  my.nhmm3=NHMM(y=ypoisson, X=tX, K=3, iters=100, burnin=10,
#'               emdist="poisson", nmix=2, delta=FALSE)
#'  OBIC(my.nhmm3)
#'  
#'  ## Predictive estimation - make 15 predictive data sets (new X) and 20 replicate data sets (same X)
#'  #filelocation="C:\\Users\\iamrandom\\Desktop\\here\\"
#'  #my.nhmm4=NHMM(y=ygamma,  X=tX, K=3, iters=100, burnin=10, 
#'  #              emdist="gamma", nmix=3, delta=TRUE, 
#'  #              outdir=filelocation, yrep=20, Xp=Xp1, ypred=15)
#'  #OBIC(my.nhmm4)  #needed larger burnin
#'  #tt=Oemparams(my.nhmm4,TRUE,filelocation)
#'  
#'  ## Exponential with W variable
#'  #filelocation="C:\\Users\\iamrandom\\Desktop\\here\\"
#'  #my.nhmm5=NHMM(y=ygamma,  X=tX, W=tW1, K=3, iters=50, burnin=10,
#'  #               emdist="gamma", nmix=3, delta=TRUE, 
#'  #              outdir=filelocation, yrep=20, Xp=Xp1, Wp=Wp1,ypred=35)
#'  #OBIC(my.nhmm5)
#'  #pp=OWcoef(my.nhmm5,filelocation)
#'  
#'  ## Gamma with fixed first variables nmix=2
#'  nmix=2; K=3; J=dim(ygamma)[2]
#'  prior1=array(1,dim=c(5,nmix,K,J));  prior1[1,1,,]=1;  
#'                 prior1[1,2,,]=2;  prior1[2,,,]=NA
#'  my.nhmm6=NHMM(y=ygamma,  X=tX, priors=prior1, K=3, iters=100, 
#'                 burnin=10, emdist="gamma", nmix=2, delta=TRUE)
#'  
#'  
#'  ### Compare my.nhmm6 (K=3) and my.nhmm7 (K=1) using both BIC 
#'  ###     and PLS (yhold is the last 10% of the data)
#'  #ygamma2=ygamma
#'  #ygamma2[1600,10]=NA  #add some missingness
#'  #ygamma2[1840,10]=NA  #add some missingness to yhold
#'  #filelocation="C:\\Users\\iamrandom\\Desktop\\here\\"
#'  #my.nhmm7=NHMM(y=ygamma2[1:1800,],  X=matrix(tX[,1:1800],1,1800), 
#'  #              W=array(tW[,1:1800,],dim=c(1,1800,15)), 
#'  #              K=3, iters=50, burnin=10, emdist="gamma", nmix=3, 
#'  #               delta=TRUE, outdir=filelocation, ymiss=TRUE, yrep=10, 
#'  #              Xp=matrix(tX[,1801:2000],1,200), Wp=array(tW[,1801:2000,],dim=c(1,200,15)), 
#'  #              ypred=10, yhold=ygamma2[1801:2000,])
#'  #OBIC(my.nhmm7)
#'  
#'  #compare K=1 and K=3
#'  }









###  # mixprior  REMOVED --- W=NULL will include only intercept terms
###          [optional] Dirichlet prior parameters used for the mixing weights if nmix > 1
###           and if W not specified. [nmix by K by J]  or [nmix+1 by K by J] if delta=TRUE
###   		
### # dirprior  [optional] prior for Dirichlet prior on the rows of the transition matrix
###            only for the HMM and MHMM, must be size KxK. If not supplied a flat prior is used.  
######################################################################################

# subseq=NULL; betapriorm=NULL; betapriorp=NULL;  iters=100; burnin=10; nmix=2; delta=FALSE; W=NULL; psipriorm=NULL; psipriorp=NULL; priors=NULL; outdir=NULL; ymiss=FALSE; yrep=0 ;ypred=0; Xp=NULL; Wp=NULL; yhold=NULL
#y=yt; X=X1t; W=W1t; K=4; iters=100; burnin=50; priors=prior1; emdist="normal"; nmix=2; delta=TRUE


NHMM=function(y, subseq=NULL, X=NULL, betapriorm=NULL, betapriorp=NULL, K=2, iters=1000, burnin=200, emdist="normal", nmix=1, delta=FALSE, W=NULL, psipriorm=NULL,psipriorp=NULL, priors=NULL, outdir=NULL, ymiss=FALSE, yrep=0 ,ypred=0, Xp=NULL, Wp=NULL, yhold=NULL)
{  
  
   # library(BayesLogit) 
   # library(msm)

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
    if(is.null(nmix)){nmix=1}	 
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
    
    ### prior checking
    #if(is.null(dirprior)){dirprior=matrix(1,K,K)} #flat prior for each row
    
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
    # emdist, nmix, delta, W(y/n), prior (inf/ninf),  fixed parameters
 
    ### set up hidden states (latent variable z) ###########
    z=numeric(T)    ### initialize z by roughly sorting y into K bins
    yy=apply(y,1,sum)  
    z=cutree(hclust(dist(yy)), k=K)
    sumy=numeric(K)
    for(k in 1:K){ sumy[k]=sum(z==k)}
    if(is.element(0,sumy)==0){z[order(yy)]=c(rep(1:K,each=floor(T/K)),rep(K,T%%K))}  
    
    
    if(!is.null(priors))#5,nmix,K,J
    { if( dim(priors)[1]!=5){stop("Wrong first dimension of *priors* input")}
      if( dim(priors)[2]!=nmix){stop("Wrong second dimension of *priors* input")}
      if( dim(priors)[3]!=K){stop("Wrong third dimension of *priors* input")}
      if( dim(priors)[4]!=J){stop("Wrong fourth dimension of *priors* input")}
    }
    
  #if(is.null(mixprior)) {  mixprior=array(1,dim=c(nmix,K,J)) }
    emcode=0
    
    theta=array(1,dim=c(2,nmix,K,J)) 
   ################################## NORMAL ########################################
    if(emdist=="normal" )
    {   emcode=1
        Rgettheta=get("RgetNormaltheta")
        if(delta==TRUE){   stop("A point mass with a Normal distribution??? try again")}
                            
                     for(k in 1:K)  #mu and sig2
                    {  for(j in 1:J)
                       {  if(sum(z==k)!=0){ if(mean(y[z==k,j])>0){theta[1,,k,j]=mean(y[z==k,j])}}
                       }   
                    }
        if(is.null(priors))   #non-informative prior
        {   priors=array(0,dim=c(5,nmix,K,J))
        }
        if(sum(priors<0, na.rm=TRUE) >0)
        {   stop("priors must be greater than or equal to zero")}
        
    }	
    
   ################################### GAMMA ##################################### 
#Gamma(alpha,beta): alpha~Gamma(a.alpha,b.alpha) beta~Gamma(a.beta,b.beta) [4 by nmix by K by J]
    if(emdist=="gamma")
    {   
        emcode=2 
        Rgettheta=get("RgetGammatheta")
        if(is.null(priors))  #set to exponential distributions- low weight prior
        {   priors=array(1,dim=c(5,nmix,K,J))
            priors[1,,,]=1   
            priors[2,,,]=NA
        }
        if(sum(priors<0, na.rm=TRUE) >0)
        {   stop("priors must be greater than or equal to zero")}
      

        for(k in 1:K) 
        {  for(j in 1:J)
           {  if(sum(z==k)!=0)  #MLE estimates to start theta
              {  theta[1,,k,j]=1  ## start at an exponential
                 if((mean(y[z==k,j])/theta[1,1,k,j]) >0){theta[2,,k,j]=mean(y[z==k,j])/theta[1,1,k,j]}
              }
           }
        }
                        	
    }
    if(emdist=="poisson")
    {   emcode=3
        Rgettheta=get("RgetPoissontheta")
       for(k in 1:K) 
         {  for(j in 1:J)
            {  if(sum(z==k)!=0)  
                {  theta[1,,k,j]=1  ## start at an exponential
                   if((mean(y[z==k,j])/theta[1,1,k,j])>0){theta[2,,k,j]=mean(y[z==k,j])/theta[1,1,k,j]}
                }
            }
         }
         if(is.null(priors))  #low weight prior
         {   priors=array(1,dim=c(5,nmix,K,J))  
         }
         
    }
   

if(emcode==0){stop("emdist is not correctly specified.")}
    
if(dim(priors)[1]!=5 || dim(priors)[2]!=nmix || dim(priors)[3]!=K || dim(priors)[4]!=J)
{   stop("priors should be a 5 by nmix by K by J matrix")  }


    if(delta==FALSE & sum(y==0)>0 & emcode==2)  #there are zeros in the Gamma data set
    {  stop("Too many zeros in the data set, use delta=TRUE for a zero inflated distribution.")
    }


    if(ypred> iters ||  yrep > iters){stop("iters must be larger than both ypred or yrep.")}


   if(ypred>0 || !is.null(yhold))  #Wp: A pT J,  Xp:  B pT
   { pT=0

     if(!is.null(Xp)){ pT=dim(Xp)[2]}else{stop("For ypred greater than zero or a yhold, there must be an Xp.")}
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

    #if(K==1){stop("K cannot equal one for now. Special case...")}
    NHMMmain(Rgettheta, z, theta, y, yboo, subseqy, subboo, X, betapriorm, betapriorp, K, iters, emdist, burnin, nmix,  W, psipriorm, psipriorp,  priors, outdir,outboo,  delta, yrep, Xp, Wp, ypred, yhold, ymiss)
    
}





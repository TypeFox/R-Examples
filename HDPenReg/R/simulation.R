#' Simulate copy number data for a case-control study.
#'
#' @param n Number of individuals.
#' @param nbSNP Size of the DNA sequence.
#' @param probCas Probability to be a case individual.
#' @param nbSeg Number of causal segments.
#' @param meanSegmentSize The mean size of anormal segment.
#' @param prob A 2*2 matrix containing probabilities:
#' 
#' prob[1,1]=probability to have an anomaly to a SNP given the person does not have the disease and the SNP is causal.
#' 
#' prob[1,2]=probability to have an anomaly to a SNP given the person does not have the disease and the SNP is not causal.
#' 
#' prob[2,1]=probability to have an anomaly to a SNP given the person has the disease and the SNP is causal.
#' 
#' prob[2,2]=probability to have an anomaly to a SNP given the person has the disease and the SNP is not causal.
#' 
#' @param alpha Parameter of the beta(alpha,alpha).
#' @return a list containing: 
#' \item{data}{A matrix of size n*nbSeg, containing values of the copy-number signal.}
#' \item{response}{A vector of size n containing the cas/control status.}
#' \item{causalSNP}{A vector of size nbSeg containing the center of causal segments.}
#' @author Quentin Grimonprez, Serge Iovleff
#' @examples 
#' data=simul(50,10000,0.4,10,150,matrix(c(0.1,0.8,0.001,0.001),nrow=2))
#' @export 

simul=function(n,nbSNP,probCas,nbSeg,meanSegmentSize,prob,alpha=15)
{
  if(missing(n))
    stop("n is missing.")
  if(missing(nbSNP))
    stop("nbSNP is missing.")
  if(missing(probCas))
    stop("probCas is missing.")
  if(missing(nbSeg))
    stop("nbSeg is missing.")
  if(missing(meanSegmentSize))
    stop("meanSegmentSize is missing.")
  if(missing(prob))
    stop("prob is missing.")
  
  .checkSim(n,nbSNP,probCas,nbSeg,meanSegmentSize,prob,alpha)
  
  #choose the causal SNPs
  SNPcaus=sort(sample(1:nbSNP,nbSeg))
  
  #generation of the response
  y=rbinom(n,1,probCas)
  
  #initialization of the copy-number signal with a normal signal
  X=matrix(rbeta(n*nbSNP,alpha,alpha)+1.5,n,nbSNP)
  
  #decide if the anomaly at a causal SNP is a gain or a loss
  anoSNPcaus=2*rbinom(nbSeg,1,0.5)+0.5
  
  anomalie=rep(0,nbSNP)
  for(i in 1:n)
  {
    #where the individuals has an anomaly?
    anomalie[SNPcaus]=rbinom(nbSeg,1,prob[y[i]+1,1])
    anomalie[-SNPcaus]=rbinom(nbSNP-nbSeg,1,prob[y[i]+1,2])
    
    #size of anomaly
    segSize=rpois(which(anomalie==1),meanSegmentSize)
    segSize[segSize==0]=1
    
    compteur=1
    for(j in which(anomalie==1))
    { 
      sequence=max(1,j-floor(segSize[compteur]/2)):min(nbSNP,j+floor(segSize[compteur]/2))
      if(j %in% SNPcaus)
        X[i,sequence]=anoSNPcaus[which(SNPcaus==j)]+rbeta(length(sequence),alpha,alpha)
      else
        X[i,sequence]=2*rbinom(1,1,0.5)+0.5+rbeta(length(sequence),alpha,alpha)
      compteur=compteur+1
    }
    
  }
  
  return(list(data=X,response=y,causalSNP=SNPcaus)) 
}

# check arguments from simulation function
.checkSim=function(n,nbSNP,probCas,nbSeg,meanSegmentSize,prob,alpha)
{
  
  ## n
  if(!.is.wholenumber(n))
    stop("n must be a positive integer")
  if(n<=0)
    stop("n must be a positive integer")
  
  ## nbSNP
  if(!.is.wholenumber(nbSNP))
    stop("nbSNP must be a positive integer")
  if(n<=0)
    stop("nbSNP must be a positive integer")
  
  ## probCas
  if(!is.double(probCas))
    stop("probCas must be a real between 0 and 1")  
  if( (probCas<0) || (probCas>1) )
    stop("probCas must be a real between 0 and 1")		
  
  ## nbSeg
  if(!.is.wholenumber(nbSeg))
    stop("nbSeg must be a positive integer smaller than nbSNP")
  if( (nbSeg<0) || (nbSeg> nbSNP) )
    stop("nbSeg must be a positive integer smaller than nbSNP")
  
  ## meanSegmentSize
  if(!.is.wholenumber( meanSegmentSize))
    stop(" meanSegmentSize must be a positive integer smaller than nbSNP")
  if( (meanSegmentSize<=0) || (meanSegmentSize>nbSNP) )
    stop(" meanSegmentSize must be a positive integer smaller than nbSNP")
  
  ## prob
  if(!is.numeric(prob) || !is.matrix(prob))
    stop("prob must be a 2*2 matrix of probability")
  if( (ncol(prob)!=2) || (nrow(prob)!=2) || any(prob<0) || any(prob>1) )
    stop("prob must be a 2*2 matrix of probability")
  
  ## alpha
  if(!is.double(alpha))
    stop("alpha must be a real greater than 1")	
  if( alpha<1 )
    stop("alpha must be a real greater than 1")		
}
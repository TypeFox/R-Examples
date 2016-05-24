#' This function computes the Kullback-Leibler divergence between two mixtures of multidimensional ISR distributions.
#' @title Kullback-Leibler divergence 
#' @author Quentin Grimonprez
#' @param proportion1,proportion2 vectors (which sums to 1) containing the K mixture proportions.
#' @param pi1,pi2 matrices of size K*p, where K is the number of clusters and p the number of dimension, containing the probabilities of a good comparaison of the model (dispersion parameters).
#' @param mu1,mu2 matrices of size K*sum(m), containing the modal ranks. Each row contains the modal rank for a cluster. In the case of multivariate ranks, the reference rank for each dimension are set successively on the same row.
#' @param m a vector containing the size of ranks for each dimension.
#' @return a real, the Kullback-Leibler divergence. 
#' @references 
#' http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
#' @examples
#' proportion1=c(0.4,0.6)
#' pi1=matrix(c(0.8,0.75),nrow=2)
#' mu1=matrix(c(1,2,3,4,4,2,1,3),nrow=2,byrow=TRUE)
#' proportion2=c(0.43,0.57)
#' pi2=matrix(c(0.82,0.7),nrow=2)
#' mu2=matrix(c(1,2,3,4,4,2,1,3),nrow=2,byrow=TRUE)
#' dK=kullback(proportion1,pi1,mu1,proportion2,pi2,mu2,4)
#' @export
kullback <-function(proportion1,pi1,mu1,proportion2,pi2,mu2,m)
{
  if(missing(proportion1))
    stop("proportion1 is missing")
  if(missing(mu1))
    stop("mu1 is missing")
  if(missing(pi1))
    stop("pi1 is missing")
  if(missing(proportion2))
    stop("proportion2 is missing")
  if(missing(mu2))
    stop("mu2 is missing")
  if(missing(pi2))
    stop("pi2 is missing")
  
  eps=1e-10
  #proportion1
  if(!is.vector(proportion1,mode="numeric"))
    stop("proportion1 must be a vector of positive real whose sum equal 1")
  if(min(proportion1)<0)
    stop("proportion1 must be a vector of positive real whose sum equal 1")
  if(abs(1-sum(proportion1))>eps)
    stop("proportion1 must be a vector of positive real whose sum equal 1")
  
  #proportion2
  if(!is.vector(proportion2,mode="numeric"))
    stop("proportion2 must be a vector of positive real whose sum equal 1")
  if(min(proportion2)<0)
    stop("proportion2 must be a vector of positive real whose sum equal 1")
  if(abs(1-sum(proportion2))>eps)
    stop("proportion2 must be a vector of positive real whose sum equal 1")
  if(length(proportion1)!=length(proportion2))
    stop("proportion1 and proportion2 must have the same length.")
  
  #m
  if(!is.vector(m,mode="numeric"))
    stop("m must be a (vector of) integer strictly greater than 1")
  if(length(m)!=length(m[m>1]))
    stop("m must be a (vector of) integer strictly greater than 1")
  if(!min(m==round(m)))
    stop("m must be a (vector of) integer strictly greater than 1")
  if( (length(m)!=ncol(pi1)) || (length(m)!=ncol(pi2)) )
    stop("The number of column of pi1 or pi2 and m don't match.")
  if( (sum(m)!=ncol(mu1)) || (sum(m)!=ncol(mu2)) )
    stop("The number of column of mu1 or mu2 and sum(m) don't match.")
  
  #pi1
  if(!is.numeric(pi1) || !is.matrix(pi1))
    stop("pi1 must be a matrix of probabilities")
  if( (min(pi1)<0) && (max(pi1)>1) )
    stop("pi1 must be a matrix of probabilities")
  
  #pi2
  if(!is.numeric(pi2) || !is.matrix(pi2))
    stop("pi2 must be a matrix of probabilities")
  if( (min(pi2)<0) && (max(pi2)>1) )
    stop("pi2 must be a matrix of probabilities")
  if(length(pi1)!=length(pi2))
    stop("pi1 and pi2 must have the same size.")
  if( (nrow(pi1)!=length(proportion1)) || (nrow(pi1)!=nrow(mu1)) )
    stop("The number of rows of pi1 doesn't match with the others parameters.")
  if( (nrow(pi2)!=length(proportion2)) || (nrow(pi2)!=nrow(mu2)) )
    stop("The number of rows of pi2 doesn't match with the others parameters.") 
  
  #mu1
  if(!is.numeric(mu1) || !is.matrix(mu1))
    stop("mu1 must be a matrix of positive integer")
  if(min(mu1)<1)
    stop("mu1 must be a matrix of positive integer")
  if(nrow(mu1)!=length(proportion1))
    stop("The number of rows of mu1 and the length of proportion1 don't match.")
  if(nrow(mu1)!=nrow(pi1))
    stop("The number of rows of mu1 and pi1 doesn't match.")
  
  #mu2
  if(!is.numeric(mu2) || !is.matrix(mu2))
    stop("mu2 must be a matrix of positive integer")
  if(min(mu2)<1)
    stop("mu2 must be a matrix of positive integer")
  if(nrow(mu2)!=length(proportion2))
    stop("The number of rows of mu2 and the length of proportion2 don't match.")
  if(nrow(mu2)!=nrow(pi2))
    stop("The number of rows of mu2 ans pi2 doesn't match.")
  
  #check if mu contains ranks
  for(i in 1:length(m))
  {
    if(sum(apply(mu1[,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])],1,checkRank,m[i]))!=nrow(mu1))
      stop("mu1 is not correct")
    if(sum(apply(mu2[,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])],1,checkRank,m[i]))!=nrow(mu2))
      stop("mu2 is not correct")
  }
  
  a=t(pi1)
  b=t(pi2)
  dKL=.Call("kullback",m,mu1,mu2,a,b,proportion1,proportion2,PACKAGE="Rankcluster")
  
  return(dKL)
}



#' This function computes the p-value of the khi2 adequation test (only for univariate data).
#' @title khi2 test
#' @author Quentin Grimonprez
#' @param data a matrix in which each row is a rank of size m.
#' @param proportion a vector (which sums to 1) containing the K mixture proportion.
#' @param pi a vector of size K, where K is the number of clusters, containing the probabilities of a good paired comparaison of the model (dispersion parameters).
#' @param mu a matrix of size K*m, where m is the size of a rank, containing the modal rankings of the model (position parameters).
#' @param nBoot number of bootstrap iterations used to estimate the khi2 adequation test p-value.
#' @return a real, the p-value of the khi2 adequation test. 
#' @examples
#' proportion=c(0.4,0.6)
#' pi=c(0.8,0.75)
#' mu=matrix(c(1,2,3,4,4,2,1,3),nrow=2,byrow=TRUE)
#' #simulate a data set with declared parameters.
#' data=rbind(simulISR(proportion[1]*100,pi[1],mu[1,]),
#' simulISR(proportion[2]*100,pi[2],mu[2,]))
#' pval=khi2(data,proportion,mu,pi)
#' @export
khi2 <-function(data,proportion,mu,pi,nBoot=1000)
{
  if(missing(proportion))
    stop("proportion is missing")
  if(missing(mu))
    stop("mu is missing")
  if(missing(pi))
    stop("pi is missing")
  if(missing(data))
    stop("data is missing")
  
  eps=1e-10
  #proportion
  if(!is.vector(proportion,mode="numeric"))
    stop("proportion must be a vector of positive real whose sum equal 1")
  if(min(proportion)<0)
    stop("proportion must be a vector of positive real whose sum equal 1")
  if(abs(1-sum(proportion))>eps)
    stop("proportion must be a vector of positive real whose sum equal 1")
  
  #pi
  if(!is.vector(pi,mode="numeric"))
    stop("pi must be a vector of probabilities")
  if( (min(pi)<0) && (max(pi)>1) )
    stop("pi must be a vector of probabilities")
  
  #mu
  if(!is.numeric(mu) || !is.matrix(mu))
    stop("mu must be a matrix of positive integer")
  if(min(mu)<1)
    stop("mu must be a matrix of positive integer")
  if(nrow(mu)!=length(proportion))
    stop("The number of rows of mu and the length of proportion don't match.")
  if(nrow(mu)!=length(pi))
    stop("The number of rows of mu and the length of pi don't match.")
  
  #data
  if(missing(data))
    stop("data is missing")
  if(!is.numeric(data) || !is.matrix(data))
    stop("X must be a matrix of positive integer")
  if(length(data[data>=0])!=length(data))
    stop("data must be a matrix of positive integer")
  
  if(ncol(data)!=ncol(mu))
    stop("mu and data must have the same number of columns.")
  
  
  #nBoot
  if(!is.numeric(nBoot) )
    stop("nBoot must be a positive integer.")
  if(length(nBoot)!=1 )
    stop("nBoot must be a positive integer.")
  if( (nBoot<0) || (nBoot!=round(nBoot)) )
    stop("nBoot must be a positive integer.")
  
  #check if mu and data are rank
  if(sum(apply(data,1,checkPartialRank))!=nrow(data))
    stop("Data are not correct")
  if(sum(apply(mu,1,checkPartialRank))!=nrow(mu))
    stop("mu is not correct")
  
  pval=.Call("adkhi2partial",data,pi,proportion,mu,nBoot,PACKAGE="Rankcluster")
  
  return(pval)
}

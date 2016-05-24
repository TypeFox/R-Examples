#' This function estimates the loglikelihood of a mixture of multidimensional ISR model, as well as the BIC and ICL model selection criteria.
#' @title criteria estimation
#' @author Quentin Grimonprez
#' @param data a matrix in which each row is a rank (partial or not; for partial rank, missing elements of a rank are put to 0 ).
#' @param proportion a vector (which sums to 1) containing the K mixture proportions.
#' @param pi a matrix of size K*p, where K is the number of clusters and p the number of dimension, containing the probabilities of a good comparaison of the model (dispersion parameters).
#' @param mu a matrix of size K*sum(m), containing the modal ranks. Each row contains the modal rank for a cluster. In the case of multivariate ranks, the reference rank for each dimension are set successively on the same row.
#' @param m a vector containing the size of ranks for each dimension.
#' @param Ql number of iterations of the Gibbs sampler used for the estimation of the log-likelihood.
#' @param Bl burn-in period of the Gibbs sampler.
#' @param IC number of run of the computation of the loglikelihood.
#' @param nb_cpus number of cpus for parallel computation
#' @return a list containing:
#'   \item{ll}{the estimated log-likelihood.}
#'   \item{bic}{the estimated BIC criterion.}
#'   \item{icl}{the estimated ICL criterion.}
#' @examples
#' data(big4)
#' res=rankclust(big4$data,m=big4$m,K=2,Ql=100,Bl=50,maxTry=2)
#' if(res@@convergence)
#' 	crit=criteria(big4$data,res[2]@@proportion,res[2]@@pi,res[2]@@mu,big4$m,Ql=200,Bl=100)
#' @export
criteria <-function(data,proportion,pi,mu,m,Ql=500,Bl=100,IC=1, nb_cpus=1)
{
  if(missing(proportion))
    stop("proportion is missing")
  if(missing(mu))
    stop("mu is missing")
  if(missing(pi))
    stop("pi is missing")
  if(missing(m))
    stop("m is missing")
  
  #data
  if(missing(data))
    stop("data is missing")
  if(!is.numeric(data) || !is.matrix(data))
    stop("data must be a matrix of positive integer")
  if(length(data[data>=0])!=length(data))
    stop("data must be a matrix of positive integer")
  
  #proportion
  if(!is.vector(proportion,mode="numeric"))
    stop("proportion must be a vector of positive real whose sum equal 1")
  if(min(proportion)<0)
    stop("proportion must be a vector of positive real whose sum equal 1")
  if(abs(1-sum(proportion))>1e-10)
    stop("proportion must be a vector of positive real whose sum equal 1")
  
  #m
  if(!is.vector(m,mode="numeric"))
    stop("m must be a (vector of) integer strictly greater than 1")
  if(length(m)!=length(m[m>1]))
    stop("m must be a (vector of) integer strictly greater than 1")
  if(!min(m==round(m)))
    stop("m must be a (vector of) integer strictly greater than 1")
  if( (length(m)!=ncol(pi)) )
    stop("The number of column of p and m don't match.")
  if(sum(m)!=ncol(mu)) 
    stop("The number of column of mu and sum(m) don't match.")
  
  #p
  if(!is.numeric(pi) || !is.matrix(pi))
    stop("pi must be a matrix of probabilities")
  if( (min(pi)<0) && (max(pi)>1) )
    stop("pi must be a matrix of probabilities")
  if( (nrow(pi)!=length(proportion)) || (nrow(pi)!=nrow(mu)) )
    stop("The number of rows of pi doesn't match with the others parameters.")
  
  #Ql
  if(!is.numeric(Ql) || (length(Ql)>1))
    stop("Ql must be a strictly positive integer")
  if( (Ql!=round(Ql)) || (Ql<=0))
    stop("Ql must be a strictly positive integer")
  
  #IC
  if(!is.numeric(IC) || (length(IC)>1))
    stop("IC must be a strictly positive integer")
  if( (IC!=round(IC)) || (IC<=0))
    stop("IC must be a strictly positive integer")
  
  #nb_cpus
  if(!is.numeric(nb_cpus) || (length(nb_cpus)>1))
    stop("nb_cpus must be a strictly positive integer")
  if( (nb_cpus!=round(nb_cpus)) || (nb_cpus<=0))
    stop("nb_cpus must be a strictly positive integer")
  
  #Bl
  if(!is.numeric(Bl) || (length(Bl)>1))
    stop("Bl must be a strictly positive integer lower than Ql")
  if( (Bl!=round(Bl)) || (Bl<=0) || (Bl>=Ql))
    stop("Bl must be a strictly positive integer lower than Ql")
  
  #mu
  if(!is.numeric(mu) || !is.matrix(mu))
    stop("mu must be a matrix of positive integer")
  if(min(mu)<1)
    stop("mu must be a matrix of positive integer")
  if(nrow(mu)!=length(proportion))
    stop("The number of rows of mu and the length of proportion don't match.")
  if(nrow(mu)!=nrow(pi))
    stop("The number of rows of mu and pi doesn't match.")
  
  
  #check if mu contains ranks
  for(i in 1:length(m))
  {
    if(sum(apply(mu[,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1]),drop=FALSE],1,checkRank,m[i]))!=nrow(mu))
      stop("mu is not correct")
  }
  
  #check data
  for(i in 1:length(m))
  {
    if(sum(apply(data[,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])],1,checkTiePartialRank,m[i]))!=nrow(data))
      stop("Data are not correct")
  }
  
  
  a=t(pi)
  
  LL=.Call("loglikelihood",data,mu,a,proportion,m,Ql,Bl,IC,nb_cpus,PACKAGE="Rankcluster")
  
  if(LL$ll[1]=="pb")
    stop("Data are not correct.")
  return(LL)
}



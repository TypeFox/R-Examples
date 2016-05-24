#' This functions estimates a clustering of ranking data, potentially multivariate, partial and containing tied, based on a mixture of multivariate ISR model [2].
#' By specifying only one cluster, the function performs a modelling of the ranking data using the multivariate ISR model.
#' The estimation is performed thanks to a SEM-Gibbs algorithm in the general case.
#'
#' @title model-based clustering for multivariate partial ranking
#' @author Quentin Grimonprez
#' @param data a matrix in which each row is a ranking (partial or not; for partial ranking, 
#' missing elements must be 0. Tied are replaced by the lowest position they share). For multivariate rankings, the rankings of each dimension are 
#' placed end to end in each row. The data must be in ranking notation (see Details or 
#' \link{convertRank} functions).
#' @param m a vector composed of the sizes of the rankings of each dimension (default value is the number of column of the matrix data).
#' @param K an integer or a vector of integer with the number of clusters.
#' @param criterion criterion "bic" or "icl", criterion to minimize for selecting the number of clusters.
#' @param Qsem the total number of iterations for the SEM algorithm (defaut value=40).
#' @param Bsem burn-in period for SEM algorithm (default value=10).
#' @param RjSE a vector containing, for each dimension, the number of iterations of the Gibbs sampler 
#' used both in the SE step for partial rankings and for the presentation orders generation (default value=mj(mj-1)/2).
#' @param RjM a vector containing, for each dimension, the number of iterations of the Gibbs sampler used in the  M step (default value=mj(mj-1)/2)
#' @param Ql number of iterations of the Gibbs sampler 
#' for estimation of log-likelihood (default value=100).
#' @param Bl burn-in period for estimation of log-likelihood (default value=50).
#' @param maxTry maximum number of restarts of the SEM-Gibbs algorithm in the case of non convergence (default value=3).
#' @param run number of runs of the algorithm for each value of K.
#' @param detail boolean, if TRUE, time and others informations will be print during the process (default value FALSE).
#' @return An object of class rankclust (See \code{\link{Output-class}} and \code{\link{Rankclust-class}}).
#'
#' For example :
#' res=rankclust(data,K=1:2,m=m)
#'
#' You can access the result by res[number of groups]@@slotName where slotName is an element of the class Output.
#' @references 
#' [1] C.Biernacki and J.Jacques (2013), A generative model for rank data based on sorting algorithm, Computational Statistics and Data Analysis, 58, 162-176.
#'
#'[2] J.Jacques and C.Biernacki (2012), Model-based clustering for multivariate partial ranking data, Inria Research Report n 8113.
#'
#' @examples
#' data(big4)
#' result=rankclust(big4$data,K=2,m=big4$m,Ql=200,Bl=100,maxTry=2)
#' 
#' @details
#' 
#' The ranks have to be given to the package in the ranking notation (see \link{convertRank} function), with the following convention :
#' 
#' - missing positions are replaced by 0
#' 
#' - tied are replaced by the lowest position they share"
#' 
#' 
#'   The ranking representation r=(r_1,...,r_m) contains the
#' ranks assigned to the objects, and means that the ith
#' object is in r_ith position.
#' 
#' The ordering representation o=(o_1,...,o_m) means that object
#' o_i is in the ith position.
#' 
#' Let us consider the following example to illustrate both
#' notations: a judge, which has to rank three holidays
#' destinations according to its preferences, O1 =
#'   Countryside, O2 =Mountain and O3 = Sea, ranks first Sea,
#' second Countryside, and last Mountain. The ordering
#' result of the judge is o = (3, 1, 2) whereas the ranking
#' result is r = (2, 3, 1).
#' 
#' 
#' @seealso See \code{\link{Output-class}} and \code{\link{Rankclust-class}} for available output.
#' 
#' @export
#' 
rankclust<-function(data,m=ncol(data),K=1,criterion="bic",Qsem=100,Bsem=20,RjSE=m*(m-1)/2,RjM=m*(m-1)/2,Ql=500,Bl=100,maxTry=3,run=1,detail=FALSE)
{
  
  .checkArgRankclust(data,m,K,criterion,Qsem,Bsem,RjSE,RjM,Ql,Bl,detail,maxTry,run)
  
  
  result=c()
  
  G=c()
  for(k in K)
  {
    ## first run
    res=mixtureSEM(data,k,m,Qsem,Bsem,Ql,Bl,RjSE,RjM,maxTry,run,detail)
    if(res@convergence)
    {	
      G=c(G,k)
      result=c(result,list(res))
    }
    else
    {
      cat("\n for K=",k,"clusters, the algorithm has not converged (a proportion was equal to 0 during the process), please retry\n")
    }	
  }
  
  
  
  if(length(G)==0)
  {
    resultat=new("Rankclust",convergence=FALSE)
    cat("No convergence for all values of K (a proportion was equal to 0 during the process). Please retry")
  }
  else
  {
    colnom=c()
    for(i in 1:length(m))
      colnom=c(colnom,paste0("dim",i),rep("",m[i]-1))
    
    colnames(data)=colnom
    
    resultat=new("Rankclust",K=G,criterion=criterion,results=result,data=data,convergence=TRUE)
  }
  
  
  return(resultat)
}


.checkArgRankclust=function(data,m,K,criterion,Qsem,Bsem,RjSE,RjM,Ql,Bl,detail,maxTry,run)
{
  ##################check the arguments
  #data
  if(missing(data))
    stop("data is missing")
  if(!is.numeric(data) || !is.matrix(data))
    stop("data must be a matrix of positive integer")
  if(length(data[data>=0])!=length(data))
    stop("data must be a matrix of positive integer")
  
  
  #m
  if(!is.vector(m,mode="numeric"))
    stop("m must be a (vector of) integer strictly greater than 1")
  if(length(m)!=length(m[m>1]))
    stop("m must be a (vector of) integer strictly greater than 1")
  if(!min(m==round(m)))
    stop("m must be a (vector of) integer strictly greater than 1")
  if(sum(m)!=ncol(data))
    stop("The number of column of data and m don't match.")
  
  
  #K
  if(!is.vector(K,mode="numeric"))
    stop("K must be a (vector of) integer strictly greater than 0")
  if(length(K)!=length(K[K>0]))
    stop("K must be a (vector of) integer strictly greater than 0")
  if(!min(K==round(K)))
    stop("K must be a (vector of) integer strictly greater than 0")
  
  
  #criterion
  if(criterion!="bic" && criterion!="icl")
    stop("criterion must be \"bic\" or \"icl\" ")	
  
  #Qsem
  if(!is.numeric(Qsem) || (length(Qsem)>1))
    stop("Qsem must be a strictly positive integer")
  if( (Qsem!=round(Qsem)) || (Qsem<=0))
    stop("Qsem must be a strictly positive integer")
  
  #Bsem
  if(!is.numeric(Bsem) || (length(Bsem)>1))
    stop("Bsem must be a strictly positive integer lower than Qsem")
  if( (Bsem!=round(Bsem)) || (Bsem<=0) || (Bsem>= Qsem))
    stop("Bsem must be a strictly positive integer lower than Qsem")
  
  #RjM
  if(!is.numeric(RjM) || (length(RjM)!=length(m)))
    stop("RjM must be a vector of strictly positive integer")
  if( (RjM!=round(RjM)) || (RjM<=0))
    stop("RjM must be a vector of strictly positive integer")
  
  #RjSE
  if(!is.numeric(RjSE) || (length(RjSE)!=length(m)))
    stop("RjSE must be a vector of strictly positive integer")
  if( (RjSE!=round(RjSE)) || (RjSE<=0))
    stop("RjSE must be a vector of strictly positive integer")
  
  #Ql
  if(!is.numeric(Ql) || (length(Ql)>1))
    stop("Ql must be a strictly positive integer")
  if( (Ql!=round(Ql)) || (Ql<=0))
    stop("Ql must be a strictly positive integer")
  
  #Bl
  if(!is.numeric(Bl) || (length(Bl)>1))
    stop("Bl must be a strictly positive integer lower than Ql")
  if( (Bl!=round(Bl)) || (Bl<=0) || (Bl>=Ql))
    stop("Bl must be a strictly positive integer lower than Ql")
  
  #maxTry
  if(!is.numeric(maxTry) || (length(maxTry)!=1))
    stop("maxTry must be a positive integer")
  if( (maxTry!=round(maxTry)) || (maxTry<=0))
    stop("maxTry must be a positive integer")
  
  #run
  if(!is.numeric(run) || (length(run)!=1))
    stop("run must be a positive integer")
  if( (run!=round(run)) || (run<=0))
    stop("run must be a positive integer")
  
  #detail
  if(!is.logical(detail))
    stop("detail must be a logical.")
}

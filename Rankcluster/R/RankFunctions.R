convertInter=function(x)
{
  xb=rep(0,length(x))
  ind=x[x>0]
  for(i in ind)
    xb[i]=which(x==i)
  
  return(xb)
}


#'convertRank converts a rank from its ranking representation to its ordering representation, and vice-versa. The function does not work with partial ranking.
#'The transformation to convert a rank from ordering to ranking representation is the same that from ranking to ordering representation, there is no need to precise the representation of rank x.
#'
#'The ranking representation r=(r_1,...,r_m) contains the ranks assigned to the objects,
#'and means that the ith object is in r_ith position.
#'
#'The ordering representation o=(o_1,...,o_m) means that object o_i is in the ith position. 
#'
#'
#'Let us consider the following example to illustrate both notations: a judge, which has to rank three holidays destinations according to its preferences, O1 = Countryside, O2 =Mountain and O3 = Sea, ranks first Sea, second Countryside, and last Mountain. The ordering result of the judge is o = (3, 1, 2) whereas the ranking result is r = (2, 3, 1).
#' @useDynLib Rankcluster
#' @title change the representation of a rank
#' @author Julien Jacques
#' @param x a rank (vector) datum either in its ranking or ordering representation.
#' @return a rank (vector) in its ordering representation if its ranking representation has been given in input of convertRank, and vice-versa.
#' @examples
#' x=c(2,3,1,4,5)
#' convertRank(x)
#' @export
convertRank <- function(x)
{
  if(is.matrix(x))
    return(t(apply(x,1,convertInter ) ) )
  else
    return(convertInter(x))
}

# '
# ' This function checks if data are correct.
# ' 
# '
# ' @title Check the data
# ' @param X a matrix containing ranks
# ' @param m a vector composed of the sizes of the rankings of each dimension (default value is the number of column of the matrix data).
# ' 
# ' @return a list containing for each dimension, numbers of rows with problem.
# ' 
# ' @examples
# ' data(big4)
# ' #add a bad rank
# ' big4$data[1,1:4] = c(1,5,2,3)
# ' 
# ' res=checkData(big4$data,big4$m)
# ' print(res)
# ' 
# ' @export
# checkData = function(X,m=length(X))
# {
#   if(!is.matrix(X))
#   {
#     X=as.matrix(X)
#   }
#   if(sum(m)!=ncol(X))
#     stop("the number of columns of X does not match with m.")
#   
#   d=length(m)
#   n=nrow(X)
#   cm=cumsum(c(0,m))
#   pb=list()
#   for(i in 1:d)
#   {
#     check=apply(X[,(1+cm[i]):cm[i+1],drop=FALSE],1,checkTiePartialRank,m[i])
#     if(sum(check)!=n)
#     {
#       indfalse=which(check==0)
#       pb[[i]]=indfalse
#     }
#   }
#   
#   return(pb)
# }

# checkRank  check if a vector is a rank
checkRank <- function(x,m=length(x))
{
  if(sum(sort(x)==(1:m))==m)
    return(TRUE)
  else
    return(FALSE)	
}

# checkPartialRank check if a vector is a partial rank
checkPartialRank <- function(x,m=length(x))
{
  if((length(x[x<=m])==m)&& (length(x[x>=0])==m) && (length(unique(x[x!=0]))==length(x[x!=0])))
    return(TRUE)
  else
    return(FALSE)
}


# checkPartialRank check if a vector is a partial rank
checkTiePartialRank <- function(x,m=length(x))
{
  if((length(x[x<=m])==m) && (length(x[x>=0])==m) )
    return(TRUE)
  else
    return(FALSE)
}

# completeRank complete partial that have only one missing element
completeRank <-function(x)
{
  if(length(x[x==0])==1)
  {	
    m=length(x)
    a=1:m
    a[x[x!=0]]=0
    x[x==0]=a[a!=0]
  }
  return(x)
}

# check if a number is an integer
# @param x number
# @param tol tolerance
#
# @ return TRUE if the number is an integer, FALSE else
#
is.wholenumber=function(x, tol = .Machine$double.eps^0.5)  
{
  #if(!is.double(x))
  #  return(FALSE)
  
  abs(x - round(x)) < tol
}

#' This function takes in input a matrix containing all the observed ranks (a rank can be repeated) and returns a matrix containing all the different observed ranks with their observation frequencies (in the last column).
#' @title Convert data storage
#' @author Quentin Grimonprez
#' @param X a matrix containing ranks.
#' @param m a vector with the size of ranks of each dimension.
#' @return A matrix containing each different observed ranks with its observation frequencies in the last column.
#' @examples
#' X=matrix(1:4,ncol=4,nrow=5,byrow=TRUE)
#' Y=frequence(X)
#' Y
#' @export
frequence <-function(X,m=ncol(X))
{
  if(missing(X))
    stop("X is missing")
  if(!is.numeric(X) || !is.matrix(X))
    stop("X must be a matrix of positive integer")
  if(length(X[X>=0])!=length(X))
    stop("X must be a matrix of positive integer")
  if(!is.vector(m,mode="numeric"))
    stop("m must be a (vector of) integer strictly greater than 1")
  if(length(m)!=length(m[m>1]))
    stop("m must be a (vector of) integer strictly greater than 1")
  
  if(length(m)==1)
  {
    if(m!=ncol(X))
    {
      print(paste0("You put m=",m,", but X has ",ncol(X)," columns(rank of size ",ncol(X)-1," and 1 for the frequence)."))
      print(paste0("The algorithm will continue with m=",ncol(X)-1))
    }
  }
  
  res=.Call("freqMultiR",X,m,PACKAGE="Rankcluster")
  
  data=matrix(0,ncol=length(res$data[[1]])+1,nrow=length(res$data))
  for(i in 1:nrow(data))
    data[i,]=c(res$data[[i]],res$freq[[i]])
  
  
  return(data)
  
}

#' This function simulates univariate rankings data (ordering representation) according to the ISR(pi,mu).
#' @title simulate a sample of ISR(pi,mu)
#' @author Julien Jacques
#' @param n size of the sample.
#' @param pi dispersion parameter: probability of correct paired comparaison according to mu.
#' @param mu position parameter: modal ranking in ordering representation.
#' @return a matrix with simulated ranks.
#' 
#' @details
#' The ranking representation r=(r_1,...,r_m) contains the
#' ranks assigned to the objects, and means that the ith
#' object is in r_ith position.
#' 
#' The ordering representation o=(o_1,...,o_m) means that
#' object o_i is in the ith position.
#' 
#' Let us consider the following example to illustrate both
#' notations: a judge, which has to rank three holidays
#' destinations according to its preferences, O1 =
#'   Countryside, O2 =Mountain and O3 = Sea, ranks first Sea,
#' second Countryside, and last Mountain. The ordering
#' result of the judge is o = (3, 1, 2) whereas the ranking
#' result is r = (2, 3, 1).
#' 
#' You can see the \link{convertRank} function to convert the simualted ranking drom ordering to ranking representation.
#' 
#' @references 
#' [1] C.Biernacki and J.Jacques (2013), A generative model for rank data based on sorting algorithm, Computational Statistics and Data Analysis, 58, 162-176.
#' @examples
#' x=simulISR(30,0.8,1:4)
#' @export
simulISR <-function(n,pi,mu)
{
  if(missing(n))
    stop("n is missing")
  if(missing(mu))
    stop("mu is missing")
  if(missing(pi))
    stop("pi is missing")
  
  if(!is.numeric(n) || (length(n)>1))
    stop("n must be a strictly positive integer")
  if( (n!=round(n)) || (n<=0))
    stop("n must be a strictly positive integer")
  
  if(!is.numeric(pi) || (length(pi)>1))
    stop("pi must be a real between 0 and 1")
  if( (pi>1) || (pi<0))
    stop("pi must be a real between 0 and 1")
  
  if(!is.vector(mu,mode="numeric"))
    stop("mu must be a complete rank")
  if(!checkRank(mu))
    stop("mu must be a complete rank")
  
  
  
  res=.Call("simulISRR",n,length(mu),mu,pi,PACKAGE="Rankcluster")
  
  return(res)
}

#' This function takes in input a matrix in which the m first columns are the different observed ranks and the last column contains the observation frequency, and returns a matrix containing all the ranks (ranks with frequency>1 are repeated).
#' @title Convert data
#' @param data a matrix containing rankings and observation frequency.
#' @return a matrix containing all the rankings.
#' @examples
#' data(quiz)
#' Y=unfrequence(quiz$frequency)
#' Y
#' @export
unfrequence=function(data)
{
  X=matrix(ncol=ncol(data)-1,nrow=sum(data[,ncol(data)]))
  colnames(X)=colnames(data)[-ncol(data)]
  compteur=1
  for(i in 1:nrow(data))
    for(j in 1:data[i,ncol(data)])
    {
      X[compteur,]=data[i,-ncol(data)]
      compteur=compteur+1
    }
  return(X)    
}

#' This function computes the probability of x according to a multivariate ISR o parameter mu and pi.
#' @title Probability computation
#' @author Quentin Grimonprez
#' @param x a vector or a matrix of 1 row containing the rankings in ranking notation (see Details or \link{convertRank} function). The rankings of each dimension are 
#' placed end to end. x must contain only full ranking (no partial or tied).
#' @param pi a vector of size p, where p is the number of dimension, containing the probabilities of a good comparaison of the model (dispersion parameters).
#' @param mu a vector of length sum(m) or a matrix of size 1*sum(m), containing the modal ranks in ranking notation (see Details or \link{convertRank} function). The rankings of each dimension are 
#' placed end to end. mu must contain only full ranking (no partial or tied).
#' @param m a vector containing the size of ranks for each dimension.
#' @return the probability of x according to a multivariate ISR o parameter mu and pi.
#' 
#' @details
#' 
#'  The ranks have to be given to the package in the ranking notation (see \link{convertRank} function), with the following convention :
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
#' @examples
#' m=c(4,5)
#' x=mu=matrix(nrow=1,ncol=9)
#' x[1:4] = c(1,4,2,3)
#' x[5:9] = c(3,5,2,4,1)
#' mu[1:4] = 1:4
#' mu[5:9] = c(3,5,4,2,1)
#' pi=c(0.75,0.82)
#' 
#' prob=probability(x,mu,pi,m)
#' prob
#' @export
#' 
probability = function(x,mu,pi,m=length(x))
{
  ### check parameters
  if(missing(x))
    stop("x is missing.")
  if(missing(mu))
    stop("mu is missing.")
  if(missing(pi))
    stop("pi is missing.")
  
  #x
  if(!(is.vector(x) || is.matrix(x)))
    stop("x must be either a matrix or a vector.")
  if(is.vector(x))
    x=t(as.matrix(x))
  if(!is.numeric(x))
    stop("x must be either a matrix or a vector of integer.")
  
  #mu
  if(!(is.vector(mu) || is.matrix(mu)))
    stop("mu must be either a matrix or a vector.")
  if(is.vector(mu))
    mu=t(as.matrix(mu))
  if(!is.numeric(mu))
    stop("mu must be either a matrix or a vector of integer.")
  
  #pi
  if(!is.numeric(pi))
    stop("pi must be a vector of probabilities.")
  if(!is.vector(pi))
    stop("pi must be a vector of probabilities.") 
  if((min(pi)<0) || max(pi)>1)
    stop("pi must be a vector of probabilities.")
  
  #m
  if(!is.numeric(m))
    stop("m must be a vector of integer.")
  if(!is.vector(m))
    stop("m must be a vector of integer.") 
  if(sum(unlist(lapply(m,is.wholenumber)))!=length(m))
    stop("m contains non integer.")
  if(sum(m)!=length(x))
    stop("sum(m) and the length of x do not match.")
  if(sum(m)!=length(mu))
    stop("sum(m) and the length of mu do not match.")
  if(length(m)!=length(pi))
    stop("the length of pi and m do not match.")
  
  #check if mu contains ranks
  for(i in 1:length(m))
  {
    if(!checkRank(mu[,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])],m[i]))
      stop("mu is not correct.")
    if(!checkRank(x[,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])],m[i]))
      stop("x is not correct.")
  }
  
  #convert to ordering
  for(i in 1:length(m))
  {
    mu[1,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])] = convertRank(mu[1,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])])
    x[1,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])] = convertRank(x[1,(1+cumsum(c(0,m))[i]):(cumsum(c(0,m))[i+1])])
  }
  
  
  prob=.Call("computeProba",x,mu,pi,m,PACKAGE = "Rankcluster")
  
  return(prob)
}


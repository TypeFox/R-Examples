#' Trim probability distribution to unique events with positive probability
#'
#' @param dist list with named numeric vectors \code{x} and \code{fx}, denoting respectively the events and probabilities of the discrete distribution.
#' @details The function reduces \code{x} to the unique values and sums the corresponding elements from \code{fx}.
#' @return list with named numeric vectors \code{x} and \code{fx}, denoting respectively the events and probabilities of the discrete distribution.
#' @examples
#' dist.unique.events(list(x=c(0,1,1,2),fx=c(0.2,0.25,0.15,0.4)))
#' @export
dist.unique.events <- function(dist){
  x0 <- sort(unique(dist$x)) #retain unique vals
  fx0 <- as.vector(tapply(dist$fx,match(dist$x,x0),FUN=sum,simplify=TRUE)) # sum prs by unique vals
  list(x=x0[fx0>0],fx=fx0[fx0>0])
}
NULL
#' Distribution of product of several discrete random variables as product X*Y
#'
#' @param dists a list of distributions
#' @param n.max maximum number of mass points of discrete distribution used in the process
#' @param appr if TRUE, then the distributions are shrunken (approximated), if necessary, to not exceed n.max
#' @param appr.method integer: 1 (merge mass points to lower bound); 2 (merge to upper bound)
#' @param n.max.appr maximum number of mass points of shrunken distributions
#' @param r0 numeric, relative tolerance used in first step of shrinking the distributions
#' @param R numeric, \code{r0} is multiplied with \code{R} until the number of mass points is at most \code{n.max.appr}
#' @return list with named sublists: 
#' \itemize{
#'  \item cumdist1: a list with vectors \code{x}, \code{Fx}
#'  \item dist2: a list with vectors \code{x}, \code{fx}
#' }
#' @examples
#' 
#' data(freqsNLngm)
#' 
#' set.seed(123)
#' x <- sample.profiles(1,freqsNLngm)
#' 
#' # per locus distribution of kinship index
#' dists <- ki.dist(x,hyp.1="FS",hyp.2="UN",hyp.true="UN")
#' 
#' n <- sapply(dists,function(x) length(x$fx))
#' prod(n) # too many outcomes to store!

#' # but, for two subsets of the loci, the distribution can be obtained
#' pair <- dists.product.pair(dists)
#' str(pair) # with these, we can compute exceedance probabilities quickly
#' 
#' # obtain the cdf as a function
#' cdf <- dist.pair.cdf(pair)
#' cdf(1)
#' 
#' # plot the cdf
#' x0 <- seq(from=-10,to=5,length=50)
#' plot(x0,cdf(10^x0),type="l",xlab="x",ylab="Fn(x)")
#' @export
dists.product.pair <- function(dists,n.max=1e6,appr=FALSE,appr.method=1L,n.max.appr=1e3,r0=1e-2,R=1.05){
  # if possible, computes the dist of two partial products of the rv's,
  # s.t. both have max. n (defaults to 1e7) events
  # in this way, a product with at most n.max^2 events can be studied
  # if this is not possible, then obtain shrunken distributions
  if (length(dists)<2) stop("Supply at least two dists to obtain the distribution of the product")
  
  nn <- sapply(dists,function(x) length(x$x))
  dists.subsets <- Zfind.subsets.with.max.product(nn,max.product=n.max)
  if (length(dists.subsets)>2){
    if (!appr){
      stop("Problem is too big. Could not find subsets of the distributions such that all products have at most n.max events without approximating the distribution. Increase n.max or set appr=TRUE.")
    }else{
      # approximate the distributions until the product X*Y can be obtained
      repeat{
        n <- sapply(dists,function(x) length(x$x)) # events of remaining variables
        dists <- dists[order(n)] # sort from few to many events
        n <- sapply(dists,function(x) length(x$x)) # events of remaining variables
        dists.subsets <- Zfind.subsets.with.max.product(n,max.product=n.max)
        
        if (length(dists.subsets)<=2) break
        
        if (length(dists.subsets)<length(dists)){
          # reduce the number of variables by computing the product distribution
          dists <- lapply(dists.subsets,function(s) dists.product(dists[s]))
        }else{
          # we can not reduce the number of variables by computing the product
          # we have to reduce the number of events by approximation
          dists <- lapply(dists,Zdistapprox,maxn=n.max.appr,r0=r0,R=R,method=appr.method) 
          if (!(length(Zfind.subsets.with.max.product(sapply(dists,function(x) length(x$x)),max.product=n.max))<length(dists.subsets)))
            stop("Approximation does not decrease number of events such that n^2<n.max. Ensure n.max is at least n.max.app^2")      
        } 
      }
      # create a pair (cumdist of a variable, dist of other variable)
      return(dists.product.pair(dists,n.max=n.max))
    }
  }
      

  if (length(dists.subsets)==1){
    # number of events of product is smaller than n.max -> fit all but the last marginal in the cdf
    dists.subsets <- list(1:(length(dists)-1),length(dists))
  }
  list(cumdist1=dists.product(dists[dists.subsets[[1]]],n.max=n.max,return.cumdist=TRUE),
       dist2= dists.product(dists[dists.subsets[[2]]],n.max=n.max,return.cumdist=FALSE))
}
NULL
#' Distribution of product of several discrete random variables
#'
#' @param dists a list of distributions
#' @param n.max memory limit, distribution will not be obtained when number of mass points exceeds n.max
#' @param return.cumdist when TRUE, returns the cumulative distribution
#' @return list with named numeric vectors \code{x} and \code{fx} (or \code{Fx} when \code{return.cumdist==TRUE}), denoting respectively the events and probabilities of the discrete distribution.
#' @examples
#' 
#' # x is a fair die, i.e. has probability distribution:
#' die <- list(x=1:6,fx=rep(1/6,6))
#' 
#' # we have 5 of them
#' dists <- replicate(5,die,simplify=FALSE)
#' 
#' # what is the distribution of x1*x2*x3*x4*x5?
#' prod.dist <- dists.product(dists)
#' 
#' plot(prod.dist$x,prod.dist$fx,xlab="x1*x2*x3*x4*x5",ylab="fx",type="h")
#' @export
dists.product <- function(dists,n.max=1e8,return.cumdist=FALSE){
  # computes the dist of a product of nonnegative rv's with given dists
  # actual work is done in a not-exported c++ function, which requires some preprocessing
  # since that function expects strictly positive and finite values
  
  tmp <- Zpdists.properties(dists)
#   dists.pr.0 <- sapply(dists,function(f) ifelse(f$x[1]==0,f$fx[1],0) )
#   dists.pr.inf <- sapply(dists,function(f) ifelse(f$x[length(f$x)]==Inf,f$fx[length(f$x)],0) )
#   i0 <- as.integer(dists.pr.0>0)
#   n0 <- sapply(dists,function(f) length(f$x)) - as.integer(dists.pr.inf>0)
#   # only process the values from i0 to n0 and add the events 0, +Inf (when pr>0)
#   prod.pr.0 <- 1-prod(1-dists.pr.0)
#   prod.pr.inf <- 1-prod(1-dists.pr.inf)
#   prod.N <- prod(n0-i0)+(prod.pr.0>0)+(prod.pr.inf>0) # number of events of the product
  
  if (tmp$prod.N>n.max) stop("Distribution of product has possibly more than n.max events. Increase n.max to proceed.")
  
  # filter the markers for which no events are left after removing the x=0 event:
  ind <- tmp$i0!=tmp$n0
  Zproductdist(x=Zdiststomatrix.X(dists=dists[ind]),prob=Zdiststomatrix.P(dists=dists[ind]),i=tmp$i0[ind],n=tmp$n0[ind],N=tmp$prod.N,pr0=tmp$prod.pr.0,prinf=tmp$prod.pr.inf,returncumdist=return.cumdist)  
}

Zpdists.properties <- function(dists){
  # finds out properties of the product of the rv's dists[[1]]*dists[[2]]*...*dists[[n]]
  # such as pr(product=0), pr(product=Inf) and number of events of the product
  
  dists.pr.0 <- sapply(dists,function(f) ifelse(f$x[1]==0,f$fx[1],0) )
  dists.pr.inf <- sapply(dists,function(f) ifelse(f$x[length(f$x)]==Inf,f$fx[length(f$x)],0) )
  i0 <- as.integer(dists.pr.0>0)
  n0 <- sapply(dists,function(f) length(f$x)) - as.integer(dists.pr.inf>0)
  # only process the values from i0 to n0 and add the events 0, +Inf (when pr>0)
  prod.pr.0 <- 1-prod(1-dists.pr.0)
  prod.pr.inf <- 1-prod(1-dists.pr.inf)
  prod.N <- prod(n0-i0)+(prod.pr.0>0)+(prod.pr.inf>0) # number of events of the product
  
  list(pr.0=dists.pr.0,pr.inf=dists.pr.inf,prod.pr.0=prod.pr.0,prod.pr.inf=prod.pr.inf,i0=i0,n0=n0,prod.N=prod.N)
}
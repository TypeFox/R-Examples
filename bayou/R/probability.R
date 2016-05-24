#' Conditional Poisson distribution
#' 
#' \code{cdpois} calculates the probability density of a value \code{k} from a Poisson distribution with a maximum \code{kmax}. \code{rdpois} draws random numbers from a conditional Poisson distribution.
#' 
#' @rdname cdpois
#' @param k random variable value
#' @param n number of samples to draw
#' @param kmax maximum value of the conditional Poisson distribution
#' @param log log transformed density
#' @param lambda rate parameter of the Poisson distribution
#' @param ... additional parameters passed to \code{dpois} or \code{rpois}
#' @export
#' @examples
#' cdpois(10,1,10)
#' cdpois(11,1,10)
#' #rdpois(5,10,10)
cdpois <- function(k,lambda,kmax,log=TRUE){
  if(kmax < lambda) stop("lambda is too high relative to kmax")
  kmax <- ceiling(kmax)
  i <- 0:kmax
  R <- sum(dpois(i,lambda))
  num <- ifelse(k<=kmax, dpois(k,lambda), 0)
  if(log){
    log(num/R)
  } else {num/R }
}
#' @rdname cdpois
#' @export
rdpois <- function(n,lambda,kmax, ...){
  kmax <- ceiling(kmax)
  i=rep(kmax+1,n)
  j=0
  while(any(i>kmax)){
    i[i>kmax] <- rpois(sum(i>kmax),lambda, ...)
    j <- j+1
    if(j>100){stop ("Lambda too high relative to kmax")}
  }
  return(i)
}
#' Probability density functions for bayou
#' 
#' \code{dsb} calculates the probability of a particular arrangement of shifts for a given set of assumptions. 
#' 
#' @rdname dsb
#' @param sb A vector giving the branch numbers (for a post-ordered tree)
#' @param ntips The number of tips in the phylogeny
#' @param bmax A single integer or a vector of integers equal to the number of branches in the phylogeny indicating the
#' maximum number of shifts allowable in the phylogeny. Can take values 0, 1 and Inf.
#' @param prob A single value or a vector of values equal to the number of branches in the phylogeny indicating the probability that
#' a randomly selected shift will lie on this branch. Can take any positive value, values need not sum to 1 (they will be scaled to sum to 1)
#' @param log A logical indicating whether the log probability should be returned. Default is 'TRUE'
#' @param k The number of shifts to randomly draw from the distribution
#' 
#' @description This function provides a means to specify the prior for the location of shifts across the phylogeny. Certain combinations are not
#' allowed. For example, a maximum shift number of Inf on one branch cannot be combined with a maximum shift number of 1 on another. Thus, bmax must be
#' either a vector of 0's and Inf's or a vector of 0's and 1's. Also, if bmax == 1, then all probabilities must be equal, as bayou cannot sample unequal 
#' probabilities without replacement. 
#' 
#' @return The log density of the particular number and arrangement of shifts.
#' 
#' @examples
#' n=10
#' tree <- sim.bdtree(n=n)
#' tree <- reorder(tree, "postorder")
#' nbranch <- 2*n-2
#' sb <- c(1,2, 2, 3)
#' 
#' # Allow any number of shifts on each branch, with 
#' # probability proportional to branch length
#' dsb(sb, ntips=n, bmax=Inf, prob=tree$edge.length)
#' 
#' # Disallow shifts on the first branch, returns -Inf 
#' # because sb[1] = 1
#' dsb(sb, ntips=n, bmax=c(0, rep(1, nbranch-1)), 
#'                              prob=tree$edge.length)
#' 
#' # Set maximum number of shifts to 1, returns -Inf 
#' # because two shifts are on branch 2
#' dsb(sb, ntips=n, bmax=1, prob=1)
#' 
#' #Generate a random set of k branches
#' rsb(5, ntips=n, bmax=Inf, prob=tree$edge.length)
#' @export
dsb <- function(sb, ntips=ntips, bmax=1, prob=1, log=TRUE){
  if(any(!(bmax %in% c(0,1,Inf)))) stop("Number of shifts allowed per branch must be 0, 1, or Inf") 
  if(length(bmax)==1) bmax <- rep(bmax, 2*ntips-2)
  if(length(bmax)!=(2*ntips-2)) stop ("bmax not a multiple of the number of branches")
  sbt <- table(sb)
  if(any(sbt > bmax[as.numeric(names(sbt))])){
    dens <- 0
    if(log) return(log(0)) else 0
  } else {
    if(max(bmax)==1){
      if(length(prob)>1) warning("cannot sample unequal probabilities without replacement, assuming equal probabilities for each branch")
      dens <- 1/choose(sum(bmax),sum(sbt))
      if(log) return(log(dens)) else return(dens)
    } else {
      if(any(!(bmax %in% c(0,Inf)))) stop("Cannot sample unequal probabilities without replacement")
      if(length(prob)==1) prob <- rep(1,2*ntips-2)
      if(length(prob)!=2*ntips-2) stop("Number of probabilities provided must equal number of branches")  
      prob[bmax==0] <- 0
      sbp.all <- prob/sum(prob)
      sbp <- c(sbp.all[as.numeric(names(sbt))],1-sum(sbp.all[as.numeric(names(sbt))]))
      if(log) return(dmultinom(c(sbt,0),prob=sbp,log=TRUE)) else return(dmultinom(c(sbt,0),prob=sbp))
    }
  }    
}
#' @rdname dsb
#' @export
rsb <- function(k, ntips=ntips, bmax=1, prob=1, log=TRUE){
  if(any(!(bmax %in% c(0,1,Inf)))) stop("Number of shifts allowed per branch must be 0, 1, or Inf") 
  if(length(bmax)==1) bmax <- rep(bmax, 2*ntips-2)
  if(length(bmax)!=(2*ntips-2)) stop ("bmax not a multiple of the number of branches")
  if(max(bmax)==1){
    if(length(prob)>1) warning("cannot sample unequal probabilities without replacement, assuming equal probabilities for each branch")
    sb <- .sample((1:(2*ntips-2))[bmax==1], k, replace=FALSE)
    return(sb)
  } else {
    if(any(!(bmax %in% c(0,Inf)))) stop("Cannot sample unequal probabilities without replacement")
    if(length(prob)==1) prob <- rep(1,2*ntips-2)
    if(length(prob)!=2*ntips-2) stop("Number of probabilities provided must equal number of branches")  
    prob[bmax==0] <- 0
    sbp.all <- prob/sum(prob)
    sb <- suppressWarnings(.sample((1:(2*ntips-2)), k, prob=sbp.all, replace=TRUE))
    return(sb)
  }
}    
  
#' Probability density function for the location of the shift along the branch
#' 
#' \code{dloc} calculates the probability of a shift occuring at a given location along the branch assuming a uniform distribution of unit length
#' \code{rloc} randomly generates the location of a shift along the branch
#' 
#' @param loc The location of the shift along the branch
#' @param min The minimum position on the branch the shift can take
#' @param max The maximum position on the branch the shift can take
#' @param log A logical indicating whether the log density should be returned
#' @param k The number of shifts to return along a branch
#' 
#' @description Since unequal probabilities are incorporated in calculating the density via \code{dsb}, all branches are assumed to be of unit length. 
#' Thus, the \code{dloc} function simply returns 0 if \code{log=TRUE} and 1 if \code{log=FALSE}. 
#' @rdname dloc
#' @export
dloc <- function(loc,min=0,max=1,log=TRUE) if(log) return (rep(0,length(loc))) else return(rep(1,length(loc)))
#' @rdname dloc
#' @export
rloc <- function(k,min=0,max=1){
  return(runif(k))
}

#' Half cauchy distribution taken from the R package LaplacesDemon (Hall, 2012).
#' 
#' \code{dhalfcauchy} returns the probability density for a half-Cauchy distribution
#' 
#' @param x A parameter value for which the density should be calculated
#' @param scale The scale parameter of the half-Cauchy distributoin
#' @param log A logical indicating whether the log density should be returned
#' @param q A vector of quantiles
#' @param p A vector of probabilities
#' @param n The number of observations
#' 
#' @rdname dhalfcauchy
#' 
#'@export
dhalfcauchy <- function(x, scale=25, log=FALSE)
{
  x <- as.vector(x); scale <- as.vector(scale)
  if(any(scale <= 0)) stop("The scale parameter must be positive.")
  NN <- max(length(x), length(scale))
  x <- rep(x, len=NN); scale <- rep(scale, len=NN)
  dens <- log(2*scale) - log(pi*{x*x + scale*scale})
  if(log == FALSE) dens <- exp(dens)
  return(dens)
}
#' @rdname dhalfcauchy
#' @export
phalfcauchy <- function(q, scale=25)
{
  q <- as.vector(q); scale <- as.vector(scale)
  if(any(scale <= 0)) stop("The scale parameter must be positive.")
  NN <- max(length(q), length(scale))
  q <- rep(q, len=NN); scale <- rep(scale, len=NN)
  z <- {2/pi}*atan(q/scale)
  return(z)
}
#' @rdname dhalfcauchy
#' @export
qhalfcauchy <- function(p, scale=25)
{
  p <- as.vector(p); scale <- as.vector(scale)
  if(any(p < 0) || any(p > 1)) stop("p must be in [0,1].")
  if(any(scale <= 0)) stop("The scale parameter must be positive.")
  NN <- max(length(p), length(scale))
  p <- rep(p, len=NN); scale <- rep(scale, len=NN)
  q <- scale*tan({pi*p}/2)
  return(q)
}
#' @rdname dhalfcauchy
#' @export
rhalfcauchy <- function(n, scale=25)
{
  scale <- rep(scale, len=n)
  if(any(scale <= 0)) stop("The scale parameter must be positive.")
  p <- runif(n, 0, 1)
  x <- scale*tan({pi*p}/2)
  return(x)
}

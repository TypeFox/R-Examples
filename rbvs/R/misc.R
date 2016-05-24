#' @title Standardise data
#' @description Standardises the columns of a numeric matrix \code{x} (similar to R-function \code{scale}).
#' If \code{x} is a vector, it is treated as a 1-column matrix.
#' @details This function is much faster than \code{scale}.  
#' @param x A numeric matrix (or vector).
#' @param scale A logical; if \code{TRUE} each column of \code{x} is divided by the square root of the sum of its centred squares.
#' @return Matrix with centred (and optionally scaled) columns.
#' @useDynLib rbvs standardise_matrix_r
#' @examples
#' x <- matrix(rnorm(100*10), nrow = 100, ncol = 10)
#' x <- standardise(x)
#' standard.deviations <- apply(x,2,sd)
#' print(standard.deviations)
#' @export 

standardise <- function(x,scale=TRUE){
  
  x <- as.matrix(x)
  
  storage.mode(x) <- "double"
  scale <- as.integer(scale)
  
  .Call("standardise_matrix_r", x, scale)
  
} 

#' @title Generates subsamples. 
#' @details Generates \code{m}-element subsamples drawn \eqn{\lfloor\frac{{n}}{{m}}\rfloor}{floor(n/m)} times  
#' from \code{1,...,n} independently without replacement; such subsampling is repeated \code{B} times.
#' @param n The sample size.
#' @param m Subsample size (an integer lower or equal than \code{n}).
#' @param B Number of sample splits.
#' @return Matrix with the indices of the subsamples drawn in each column. 
#' @references R. Baranowski, P. Fryzlewicz (2015), Ranking Based Variable Selection, in submission (\url{http://personal.lse.ac.uk/baranows/rbvs/rbvs.pdf)}).  
#' @examples
#' subsample(10,5,2)
#' subsample(10,3,10)
#' @export 

subsample <- function(n,m,B){
	
	r <- floor(n/m);
	r.t.m <- r*m
	
	res <- matrix(0,nrow=r.t.m,ncol=B)
	
	for(i in 1:B) res[,i] <- sample.int(n,r.t.m)
	
	matrix(as.integer(res),nrow=m,byrow=FALSE)
}

		
#' @title Evaluate rankings
#' @description Returns the non-increasing order of the values in the columns of \code{x}. Ties are solved at random.
#' @param x Numeric matrix. 
#' @param k.max Integer. Indices of k.max largest elements are returned.
#' @return Matrix with the indices corresponding to the \code{k.max} largest values in \code{x}. 
#' @examples
#' omega <- abs(matrix(rnorm(100*5), nrow = 10, ncol = 5))
#' rankings(omega, k.max = 10)
#' @useDynLib rbvs rankings_r
#' @export 

rankings <- function(x,k.max){
  
  
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  if(any(is.na(x))) stop("'x' cannot contain NA's.")
  
  k.max <- as.integer(k.max)
  p <- nrow(x)
  
  if(length(k.max) > 1){
      warning("'k.max' cannot contain more than 1 elements. Taking the first one.")
      k.max <- k.max[1]
  }
  
  if(length(k.max) < 1) stop("'k.max' cannot be NULL")
  if(is.na(k.max)) stop("'k.max' cannot be NA")
  
  
  
  return(.Call("rankings_r", 
               as.matrix(x),
               as.integer(k.max)))

	
}

#' @title Find k-top-ranked sets
#' @description Finds k-top-ranked sets defined in Baranowski and Fryzlewicz (2015). This routine is used inside \code{\link{rbvs}}; it typically will be not called directly by the user. 
#' @details Uses Portable qsort_r / qsort_s library ( Turner (2013)).
#' @param rankings Matrix with rankings in each column. 
#' @param k.max Positive integer. 
#' @param min.max.freq Maximum frequency.
#' @param active A vector with previously found active variables.
#' @return List containing the following fields.
#' \item{frequencies}{Frequencies corresponding to the most frequent subsests at the top of the rankings.}
#' \item{subsets}{The moost frequent subsets.}
#' @useDynLib rbvs k_top_ranked_sets_r
#' @references R. Baranowski, P. Fryzlewicz (2015), Ranking-Based Variable Selection, in submission (\url{http://personal.lse.ac.uk/baranows/rbvs/rbvs.pdf)}).\cr
#' I. Turner (2013), Portable qsort_r / qsort_s, GitHub repository (\url{https://github.com/noporpoise/sort_r}).
#' @export 

top.ranked.sets <- function(rankings, k.max, min.max.freq=1,active=NULL){

  B <- ncol(rankings)
	p <- nrow(rankings)
	
	if(missing(k.max)) k.max <- p
	else k.max <- min(k.max,p)
	
	k.max <- as.integer(k.max)
	min.max.freq <- as.integer(min.max.freq)
	active <- unique(as.integer(active))
	
  rankings <- as.matrix(rankings)
  if(!is.integer(x[1])) storage.mode(x) <- "integer"
  
  #Check for NA's
  
  if(any(is.na(rankings))) stop("'rankings' cannot contain NA's")
  if(any(is.na(k.max))) stop("'k.max' cannot contain NA's")
  if(any(is.na(min.max.freq))) stop("'min.max.freq' cannot contain NA's")
  if(any(is.na(active))) stop("'active' cannot contain NA's")
  
  #check if parameters have proper lenghts
  
  if(length(k.max) < 0) stop("k.max cannot be NULL")
  if(length(min.max.freq) < 0) stop("min.max.freq cannot be NULL")
  
  #check if active is correct
  if(length(active) >0 ){
    min.active <- min(active)
    if(min.active <= 0 ) stop("Indices in 'active' must be > 0.")
  }

  
  
  
	
	res <- .Call("k_top_ranked_sets_r",
	    rankings,
      k.max,
	    min.max.freq,
	    active)
	
	return(res)
}

#' @title Estimate the size of the top-ranked set
#' @description Estimates the number of elements in the top-ranked set.
#' @details See Baranowski and Fryzlewicz (2015).
#' @param prob Vector with probabilities.
#' @return A list with the following fields:
#' \item{scores}{Vector with the values of the criterion.}
#' \item{s.hat}{The estimate of the number of important covariates.}
#' @references R. Baranowski, P. Fryzlewicz (2015), Ranking Based Variable Selection, in submission (\url{http://personal.lse.ac.uk/baranows/rbvs/rbvs.pdf)}).
#' @export

s.est.quotient <- function(prob){
  
  prob <- as.numeric(prob)
	if(any(is.na(prob))) stop("'prob' cannot contain NA's.")
  if(any(prob<0 || prob>1)) stop("Values in 'prob' must be between 0 and 1.")
	  
	scores <- prob/c(1,prob[-length(prob)])
	s.hat <- which.min(scores)-1
	return(list(scores=scores,s.hat=s.hat))
	
}

#' @title  Generate factor model design matrix.
#' @description This function enables a quick generation of random design matrices (see details).
#' @param n Number of independent realisations of the factor model.
#' @param p Number of covariates.
#' @param n.factors Number of factors.
#' @param sigma Standard deviation for the normal distribution (see details).
#' @return \code{n} by \code{p} matrix with independent rows following factor model (see details).
#' @useDynLib rbvs factor_model_r
#' @export
#' @details The elements of the matrix returned by this routine satisfy \eqn{X_{ij} = \sum_{l=1}^{n.factors} f_{ijl} \varphi_{il} + \theta_{ij}}{X_{ij} = \sum_{l=1}^{K} f_{ijl} \varphi_{il} + \theta_{ij}}
#' with \eqn{f_{ijl}}{f_{ijl}}, \eqn{\varphi_{il}}{\varphi_{il}}, \eqn{\theta_{ij}}{\theta_{ij}}, \eqn{\varepsilon_{i}}{\varepsilon_{i}}  i.i.d. \eqn{\mathcal{N}(0,(sigma)^2)}{\mathcal{N}(0,(sigma)^2)}.
 

factor.model.design <- function(n,p, n.factors, sigma=1){
  
  n <- as.integer(n)
  p <- as.integer(p)
  n.factors <- as.integer(n.factors)
  sigma <- as.double(sigma)
  
  if(any(is.na(n))) stop("'n' cannot be NA.")
  if(any(is.na(p))) stop("'p' cannot be NA.")
  if(any(is.na(n.factors))) stop("'n.factors' cannot be NA.")
  if(any(is.na(sigma))) stop("'n.factors' cannot be NA.")
  
  if(length(n) < 1) stop("'n' cannot be empty.")
  if(length(p) < 1) stop("'p' cannot be empty.")
  if(length(n.factors) < 1) stop("'n.factors' cannot be empty.")
  if(length(sigma) < 1) stop("'sigma' cannot be empty.")
  
  if(length(n) > 1){
    warning("'n' should be a scalar. Only first element will be used.")
    n <- n[1]
  } 

  if(length(p) > 1) {
    warning("'p' should be a scalar. Only first element will be used.")
    p <- p[1]
  }
  
  if(length(n.factors) > 1){
    warning("'n.factors' should be a scalar. Only first element will be used.")
    n.factors <- n.factors[1]
  } 
  if(length(sigma) > 1) {
    warning("'sigma' should be a scalar. Only first element will be used.")
    sigma <- sigma[1]
  }
  
  if(n <= 0)  stop("'n' must be > 0.")
  if(p <= 0)  stop("'p' must be > 0.")
  if(n.factors < 0)  stop("'n.factors' must be >= 0.")
  if(sigma <= 0)  stop("'sigma' must be > 0.")
  

  .Call("factor_model_r", n, p, n.factors, sigma)
}


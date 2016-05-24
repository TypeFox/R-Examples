#' @title Ranking-Based Variable Selection
#' @description The package implements the Ranking-Based Variable Selection 
#' algorithm proposed in Baranowski and Fryzlewicz (2015) for variable selection in high-dimensional data.
#' @details The main routine of the package is \code{\link{rbvs}}.
#' @docType package
#' @references R. Baranowski, P. Fryzlewicz (2015), Ranking-Based Variable Selection, in submission (\url{http://personal.lse.ac.uk/baranows/rbvs/rbvs.pdf)}).
#' @name rbvs-package

NULL

#' @title Ranking-Based Variable Selection
#' @description Performs Rankings-Based Variable Selection using various measures of the dependence between the predictors and the response.
#' @details Currently supported measures are: Pearson correlation coefficient (\code{measure="pc"}), Distance Correlation (\code{measure="dc"}), the regression coefficients estimated via Lasso (\code{measure="lasso"}), the regression coefficients estimated via MC+ (\code{measure="mcplus"}).
#' @param x Matrix with \code{n} observations of \code{p} covariates in each row.
#' @param y Response vector with \code{n} observations.
#' @param m Subsample size used in the RBVS algorithm.
#' @param B Number of sample splits.
#' @param measure Character with the name of the method used to measure the association between 
#' the response and the covariates. See Details below.
#' @param fun Function used to evaluate the measure given in \code{measure}. It is required when
#' method=="user". Must have at least three arguments: \code{x} (covariates matrix),  \code{.y} (response vector), \code{subsamples} (a matrix, each row contains indices of the observations to be used); return a vector of the same length as 
#' the number of covariates in \code{.x}. See for example \code{\link{pearson.cor}} or \code{\link{lasso.coef}}.   
#' @param s.est Function used to estimate the number of important covariates based on the RBVS path. Must accept \code{probs} (a vector with probabilities) as an argument. See \code{\link{s.est.quotient}} and Details below.
#' @param iterative Logical variable indicating the type of the procedure. If \code{TRUE}, an iterative extension of the RBVS algorithm is launched.  
#' @param use.residuals Logical. If true, the impact of the previously detected variables is removed from the response in the IRBVS procedure.
#' @param min.max.freq Positive integer. Optional parameter - the algorithm stops searching for the most frequent set when the frequencies reach this value.
#' @param max.iter Maximum number of iterations fot the IRBVS algorithm.
#' @param k.max Maximum size of the subset of important variables..
#' @param verbose Logical indicating wheter the progress of the algorithm should be reported.
#' @param ... Other parameters that may be passed to \code{fun} and\code{s.est}. 
#' @return Object of class rbvs with the following fields
#' \item{measure}{Character indicating type of measure used.}
#' \item{score}{List with scores at each iteration.}
#' \item{subsets}{A list with subset candidates at each iteration.}
#' \item{frequencies}{A list with observed frequencies at each iteration.}
#' \item{ranks}{Rankings evaluated (for the last iteration \code{iterative=TRUE})}
#' \item{s.hat}{Vector with the number of the covariates selected at each iteration.}
#' \item{active}{Vector with the selected covariates.}
#' \item{timings}{Vector reporting the amount of time the (I)RBVS algorithm took at each iteration.}
#' @export rbvs
#' @rdname rbvs
#' @examples
#' set.seed(1) 
#' 
#' x <- matrix(rnorm(200*1000),200,1000)
#' active <- 1:4
#' beta <- c(3,2.5,-1.7,-1)
#' y <- 1*rnorm(200) +x[,active]%*%beta
#' #RBVS algorithm
#' rbvs.object <- rbvs(x,y, iterative=FALSE)
#' rbvs.object$active
#' rbvs.object$subsets[[1]][[4]]
#' #IRBVS algorithm
#' rbvs.object <- rbvs(x,y)
#' rbvs.object$active
#' @references R. Baranowski, P. Fryzlewicz (2015), Ranking-Based Variable Selection, in submission (\url{http://personal.lse.ac.uk/baranows/rbvs/rbvs.pdf)}).
 
rbvs <- function(x,y,...) UseMethod("rbvs")

#' @rdname rbvs
#' @method rbvs default 
#' @useDynLib rbvs matrix_projection_in_place_r
#' @export

rbvs.default <- function(x, y, m, B=500, measure=c("pc","dc","lasso","mcplus", "user"), fun=NULL, s.est=s.est.quotient, iterative=TRUE, use.residuals=TRUE, k.max,min.max.freq=0,max.iter=10, verbose = TRUE, ...){
	
	## VERIFY THE INPUT ARGUMENT
	## data 
  
  x <- as.matrix(x)
  y <- as.vector(y)
  
	n <- nrow(x)
	p <- ncol(x)
	
	if(n!=length(y)) stop("The number of rows in 'x' should be equal to the length of 'y'.")
	if(p<10) stop("'x' should be a matrix with at least 10 columns.")
	
	if(any(is.na(x)))  stop("x must not contain NA values")
  if(any(is.na(y))) stop("y must not contain NA values")
	
	#veryfing B
	
	B <- as.integer(B)
	
	if(length(B) == 0){
	  warning("'B' cannot be empty. Setting B to 500.")
	  B <- 500
	}
	
	if(length(B) > 1){
	  warning("'B' must be a single value. Taking B <- B[1].")
	  B <- B[1]
	}
  

  #veryfing m
  if(B <= 3){
    warning("'B' must be a single integer. Setting B <- 500")
    B <- 500
  }

	if(missing(m)) m <- floor(n/2)
	else {
	  
	  m <- as.integer(m)
	  if(length(m) > 1){
	    warning("'m' must be a single integer. Taking m <- m[1].")
	    m <- m[1]
	  }

  if(length(m) == 0){
    warning("'m' must be a single integer. Taking m <- n/2.")
    m <- floor(n/2)
   }
	  
	 if( m<=3 || m >= n) stop("'m' must be between 4 and length(n)-1.")
	}
	
	#veryfing k.max
	if(missing(k.max)) k.max <- min(p,m)
	else {
	  
	  k.max <- as.integer(k.max)
	  if(length(k.max) > 1){
	    warning("'k.max' must be a single integer. Taking k.max <- k.max[1].")
	    k.max <- k.max[1]
	  }
	    
    if(length(k.max) == 0){
      warning("'k.max' must be a single integer. Taking k.max <- min(p,m).")
      k.max <- min(p,m)
    }
	    
	  if( k.max<=5 || k.max > p) stop("'k.max' must be between 6 and ncol(x).")
	}
	
	#veryfing min.max.freq

	min.max.freq <- as.integer(min.max.freq)
  if(length(min.max.freq) > 1){
    warning("'min.max.freq' must be a single integer. Taking min.max.freq <- min.max.freq[1].")
    min.max.freq <- min.max.freq[1]
  }
  
  if(length(min.max.freq) == 0){
    warning("'min.max.freq' must be a single integer. Taking min.max.freq <- min(p,m).")
    min.max.freq <- min(p,m)
  }
	
	#veryfing max.iter
	
	max.iter <- as.integer(max.iter)
	if(length(max.iter) > 1){
	  warning("'max.iter' must be a single integer. Taking max.iter <- max.iter[1].")
	  max.iter <- max.iter[1]
	}
	
	if(length(max.iter) == 0){
	  warning("'max.iter' must be a single integer. Taking max.iter <- 10.")
	  min.max.freq <- 10
	}  
	
	#veryfing verbose
	verbose <- as.logical(verbose)
	
	if(length(verbose) > 1){
	  warning("'verbose' must be a single integer. Taking verbose <- verbose[1].")
	  verbose <- verbose[1]
	}
	
	if(length(verbose) == 0){
	  warning("'verbose' must be a single integer. Taking verbose <- TRUE")
	  min.max.freq <- TRUE
	}
  
	## measure function
	measure <- match.arg(measure)

	
	fun <- switch(measure,
	                      pc=pearson.cor,
	                      dc=distance.cor,
	                      lasso=lasso.coef,
	                      mcplus=mcplus.coef,
	                      user=fun)

	fun.args <- names(as.list(args(fun)))

	if(length(intersect(fun.args,c("x","y","subsamples")))<3) stop("fun should be of form function(x,y,subsamples,...)")
	
	## create rbvs object
	rbvs.object <- list()
  class(rbvs.object) <- "rbvs"
  rbvs.object$measure <- as.character(measure)
	
	## initial iterators
	no.iterations <- 1;	
	i <- 1;
	
	rbvs.object$active <- rbvs.object$s.hat <- rbvs.object$timings <- c()
	rbvs.object$frequencies <- rbvs.object$subsets <- rbvs.object$scores <- rbvs.object$ranks <- list()
  
  total.samples <-  floor(n/m) * B
  
  x.current <- x
  y.current <- y 
	
	while(i <= min(no.iterations,max.iter)){
    
	  tic <- proc.time()
	  
	  k.max <- min(k.max, p - length(rbvs.object$active))    
	  rbvs.object$ranks[[i]] <- matrix(integer(1), k.max, total.samples)

    
    if(verbose) cat(paste0("Iteration ",i,": ", "evaluating omega\n"))      
    omega <- fun(x.current, y.current, subsample(n,m,B),  verbose=verbose, ...)

    if(verbose) cat(paste0("Iteration ",i,": ", "evaluating rankings\n"))      
    rbvs.object$ranks[[i]]  <- rankings(omega,k.max=k.max)

		if(verbose) cat(paste0("Iteration ",i,": ", "selecting variables\n"))
		tmp <- top.ranked.sets(rbvs.object$ranks[[i]],k.max,active=rbvs.object$active,min.max.freq=min.max.freq,...)

    
		rbvs.object$frequencies[[i]] <- tmp$frequencies
		rbvs.object$subsets[[i]] <- tmp$subsets
		
		tmp <- s.est(rbvs.object$frequencies[[i]]/ncol(rbvs.object$ranks[[i]]),...)
	
		
		rbvs.object$scores[[i]] <- tmp$scores
		
		rbvs.object$s.hat <- cbind(rbvs.object$s.hat,tmp$s.hat)
		
		if(rbvs.object$s.hat[i]>0) rbvs.object$active <- c(rbvs.object$active,as.vector(rbvs.object$subsets[[i]][[rbvs.object$s.hat[i]]]))
		
    if(verbose) cat(paste0("Iteration ",i,": ",rbvs.object$s.hat[i], " variables selected\n"))

		if(iterative & tmp$s.hat>0) {
		  if(verbose) cat(paste0("Iteration ",i,": ",rbvs.object$s.hat[i], " removing impact of the selected variables\n"))
    
			no.iterations <- no.iterations +1;

			Q <- qr.Q(qr(x[,rbvs.object$active],LAPACK=TRUE))
      
			PI <- diag(n)-tcrossprod(Q,Q)	
  
      if(use.residuals) y.current <- PI%*%y
			else y.current <- y
			
      # y.current <- y
      rm(x.current)
      x.current <- .Call("matrix_projection_in_place_r", x, PI, rbvs.object$active)
			
		}
    toc <- proc.time()
    
    rbvs.object$timings <- c(rbvs.object$timings, (toc-tic)[3])
    
		i <- i+1
	}
	
	return(rbvs.object)
	
}
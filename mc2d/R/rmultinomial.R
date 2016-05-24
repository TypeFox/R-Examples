#<<BEGIN>>
dmultinomial <- function (x, size = NULL, prob, log = FALSE)
#TITLE The Vectorized Multinomial Distribution
#KEYWORDS distribution
#DESCRIPTION
#Generate multinomially distributed random number vectors and compute multinomial probabilities.
#INPUTS
#{x}<<vector or matrix of length (or ncol) K of integers in \samp{0:size}.>>
#{n}<<number of random vectors to draw.>>
#{size}<<a vector of integers, say N, specifying the total number of objects that are put
#into K boxes in the typical multinomial experiment. For \samp{dmultinom}, it defaults to \samp{sum(x)}.
#The first element correspond to the vector \samp{prob} or the first row of \samp{prob}, ...>>
#{prob}<<Numeric non-negative vector of length K, or matrix of size \samp{(x x K)}
#specifying the probability for the K classes; is internally normalized to sum 1.>>
#{log}<<Logical; if TRUE, log probabilities are computed.>>
#EXAMPLE
#x <- c(100,200,700)
#x1 <- matrix(c(100,200,700,200,100,700,700,200,100), byrow=TRUE, ncol=3)
#p <- c(1,2,7)
#p1 <- matrix(c(1,2,7,2,1,7,7,2,1),byrow=TRUE,ncol=3)
#dmultinomial(x1,prob=p) 
### is equivalent to 
#c(	dmultinom(x1[1,],prob=p),
#	dmultinom(x1[2,],prob=p),
#	dmultinom(x1[3,],prob=p))
#
#dmultinomial(x1,prob=p1, log=TRUE) 
### is equivalent to 
#c(	dmultinom(x1[1,],prob=p1[1,], log=TRUE),
#	dmultinom(x1[2,],prob=p1[2,], log=TRUE),
#	dmultinom(x1[3,],prob=p1[3,], log=TRUE))
#
#dmultinomial(x,prob=p1,log=TRUE)
### is equivalent to 
#c(	dmultinom(x,prob=p1[1,], log=TRUE),
#	dmultinom(x,prob=p1[2,], log=TRUE),
#	dmultinom(x,prob=p1[3,], log=TRUE))
#
#prob <- c(1,2,7)
#rmultinomial(4,1000,prob)
#rmultinomial(4,c(10,100,1000,10000),prob)
#
### rmultinomial used with mcstoc
### (uncertain size and prob)
#s <- mcstoc(rpois,"U",lambda=50)
#p <- mcstoc(rdirichlet,"U",nvariates=3,alpha=c(4,10,20))
#mcstoc(rmultinomial,"VU",nvariates=3,size=s, prob=p)
#DETAILS
#These functions are the vectorized versions of \code{\link{rmultinom}} and \code{\link{dmultinom}}.
#Recycling is permitted.
#--------------------------------------------
{ 
    if(length(x) == 0) return(numeric(0))
	if(is.vector(prob)) prob <- t(prob)
    if(is.vector(x)) x <- t(x)
    K <- ncol(prob)
    maxn <- max(nx <- nrow(x), np <- nrow(prob))
	if (ncol(x) != K) stop("x[] and prob[] must be equal length vectors or equal col matrix.")
    if(nx != maxn) x <- matrix(t(x), ncol=K, nrow=maxn, byrow=TRUE)
	if(np != maxn) prob <- matrix(t(prob), ncol=K, nrow=maxn, byrow=TRUE)
	
	N <- apply(x,1,sum)
	if (is.null(size)) size <- N
    size <- rep(size, length = maxn)
	if (any(size != N)) stop("at least one size != sum(x), i.e. one is wrong")
    
	if (any(prob < 0) || any((s <- apply(prob,1,sum)) == 0))
        stop("probabilities cannot be negative nor all 0.")
    prob <- prob/s
    if (any(as.integer(x + 0.5) < 0)) stop("'x' must be non-negative")
    
    r <- sapply(1:maxn, function(y) {
										xx <- as.integer(x[y,]+.5)
										pp <- prob[y,]
										i0 <- pp == 0
										if (any(i0)) {
											if (any(xx[i0] != 0)) return(0)
											xx <- xx[!i0]
											pp <- pp[!i0]
													}
										return(lgamma(size[y] + 1) + sum(xx * log(pp) - lgamma(xx + 1)))
										}
				)
		
    if (log) return(r) else return(exp(r))
}

#<<BEGIN>>
rmultinomial <- function(n, size, prob)
#ISALIAS dmultinomial
#--------------------------------------------
{
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) can not be negative in rmultinomial")
  
  if(is.vector(prob) || (dim(prob)[1]) == 1) {
    if(length(size)==1) return(t(rmultinom(n,size,prob)))      # classical
    prob <- matrix(prob,nrow=1)
    }

  nrp <- nrow(prob)
  mnr <- min( max(nrp ,length(size)), n)
  ss  <- rep(size,length.out=mnr)              # recycling size
  
  if(nrp != mnr) prob <- matrix(t(prob),ncol=ncol(prob),nrow=mnr,byrow=TRUE)    # recycling prob

  n1 <- n%/%mnr
  n2 <- n%%mnr

  res <- sapply(1:mnr,function(x) rmultinom(n1,ss[x],prob[x,]))
  res <- matrix(res,ncol=ncol(prob),byrow=TRUE)
  index <- as.vector(matrix(1:(mnr*n1),ncol=mnr,byrow=TRUE))
  res <- res[index,]

  if (n2 != 0){
    res <- rbind(res,t(sapply(1:n2,function(x) rmultinom(1,ss[x],prob[x,]))))
    }
  return(res)
}


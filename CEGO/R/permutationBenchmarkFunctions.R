#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Create Quadratic Assignment Problem (QAP) Benchmark
#'
#' Creates a benchmark function for the Quadratic Assignment Problem.
#'
#' @param a distance matrix
#' @param b flow matrix
#'
#' @return the function of type cost=f(permutation)
#'
#' @examples
#' set.seed(1)
#' n=5
#' #ceate a flow matrix
#' A <- matrix(0,n,n) 
#' for(i in 1:n){
#' 	for(j in i:n){
#' 		if(i!=j){
#' 			A[i,j] <- sample(100,1)
#' 			A[j,i] <- A[i,j]
#'	 	}
#' 	}
#' }
#' #create a distance matrix
#' locations <- matrix(runif(n*2)*10,,2)
#' B <- as.matrix(dist(locations))
#' #create QAP objective function 
#' fun <- benchmarkGeneratorQAP(A,B)
#' #evaluate
#' fun(1:n)
#' fun(n:1)
#'
#' @seealso \code{\link{benchmarkGeneratorFSP}}, \code{\link{benchmarkGeneratorTSP}}, \code{\link{benchmarkGeneratorWT}}
#' @export
###################################################################################
benchmarkGeneratorQAP <- function(a, b) { # Generator function. 
	a 
	b #lazy evaluation fix, faster than force()
	function(x){
		bx<-b[x,x] 
		sum(a*bx) # divide by 2 if exact cost required
	}
}
#other example: http://www.neos-guide.org/content/qap4
#A <- matrix(c(0,22,53,53,22,0,40,62,53,40,0,55,53,62,55,0),4,4)
#B <- matrix(c(0,3,0,2,3,0,0,1,0,0,0,4,2,1,4,0),4,4)
#test <- benchmarkGeneratorQAP(A,B)
#test(c(2,4,1,3))/2
#test(c(3,4,1,2))/2


###################################################################################
#' Create Flow shop Scheduling Problem (FSP) Benchmark
#'
#' Creates a benchmark function for the Flow shop Scheduling Problem.
#'
#' @param a matrix of processing times for each step and each machine
#' @param n number of jobs
#' @param m number of machines
#'
#' @return the function of type cost=f(permutation)
#'
#' @examples
#' n=10
#' m=4
#' #ceate a matrix of processing times
#' A <- matrix(sample(100,replace=TRUE),n,m) 
#' #create FSP objective function 
#' fun <- benchmarkGeneratorFSP(A,n,m)
#' #evaluate
#' fun(1:n)
#' fun(n:1)
#' 
#' @seealso \code{\link{benchmarkGeneratorQAP}}, \code{\link{benchmarkGeneratorTSP}}, \code{\link{benchmarkGeneratorWT}}
#' @export
###################################################################################
benchmarkGeneratorFSP <- function(a, n, m) { # Generator function. see Reeves1995
	a 
	n #lazy evaluation fix, faster than force()
	m
	function(x){
		C=matrix(NA,n,m)
		ax <- a[x,]
		C[,1]<-as.numeric(cumsum(ax[,1]))
		C[1,]<-as.numeric(cumsum(ax[1,]))
		for(i in 2:n){
			for(j in 2:m){
				C[i,j]=max(C[i-1,j],C[i,j-1])+ax[i,j]
			}
		}		
		C[n,m]
	}
}

###################################################################################
#' Create (Asymmetric) Travelling Salesperson Problem (TSP) Benchmark
#'
#' Creates a benchmark function for the (Asymmetric) Travelling Salesperson Problem.
#' Path (Do not return to start of tour. Start and end of tour not fixed.) 
#' or Cycle (Return to start of tour). Symmetry depends on supplied distance matrix.
#'
#' @param distanceMatrix Matrix that collects the distances between travelled locations.
#' @param type Can be "Cycle" (return to start, default) or "Path" (no return to start).
#'
#' @return the function of type cost=f(permutation)
#'
#' @examples
#' set.seed(1)
#' #create 5 random locations to be part of a tour
#' n=5
#' cities <- matrix(runif(2*n),,2)
#' #calculate distances between cities
#' cdist <- as.matrix(dist(cities))
#' #create objective functions (for path or cycle)
#' fun1 <- benchmarkGeneratorTSP(cdist, "Path") 
#' fun2 <- benchmarkGeneratorTSP(cdist, "Cycle") 
#' #evaluate
#' fun1(1:n)
#' fun1(n:1)
#' fun2(n:1)
#' fun2(1:n)
#'
#' @seealso \code{\link{benchmarkGeneratorQAP}}, \code{\link{benchmarkGeneratorFSP}}, \code{\link{benchmarkGeneratorWT}}
#' @export
###################################################################################
benchmarkGeneratorTSP <- function(distanceMatrix, type="Cycle") { # Generator function
	distanceMatrix #lazy evaluation fix, faster than force()
	
	if(type=="Path"){
		f <- function (x){
			x1 <- x[-1] #without return to start point: path.
			x <- x[-length(x)]
			sum(distanceMatrix[cbind(x,x1)])
		}
	}else{
		f <- function (x){
			x1 <- c(x[-1],x[1]) #with return to start point: cycle.
			sum(distanceMatrix[cbind(x,x1)])
		}
	}
	return(f)
}

###################################################################################
#' Create single-machine total Weighted Tardiness (WT) Problem Benchmark
#'
#' Creates a benchmark function for the single-machine total Weighted Tardiness Problem.
#'
#' @param p processing times 
#' @param w weights
#' @param d due dates
#'
#' @return the function of type cost=f(permutation)
#'
#' @examples
#' n=6
#' #processing times
#' p <- sample(100,n,replace=TRUE)
#' #weights
#' w <- sample(10,n,replace=TRUE)
#' #due dates
#' RDD <- c(0.2, 0.4, 0.6,0.8,1.0)
#' TF <- c(0.2, 0.4, 0.6,0.8,1.0)
#' i <- 1
#' j <- 1
#' P <- sum(p)
#' d <- runif(n,min=P*(1-TF[i]-RDD[j]/2),max=P*(1-TF[i]+RDD[j]/2))
#' #create WT objective function
#' fun <- benchmarkGeneratorWT(p,w,d)
#' fun(1:n)
#' fun(n:1)	
#'
#' @seealso \code{\link{benchmarkGeneratorQAP}}, \code{\link{benchmarkGeneratorTSP}}, \code{\link{benchmarkGeneratorFSP}}
#' @export
###################################################################################
benchmarkGeneratorWT <- function(p, w, d) { # Generator function
	p 
	w #lazy evaluation fix, faster than force()
	d
	n= length(p)
	function(x){
		px <- p[x]
		dx <- d[x]
		wx <- w[x]
		s=c(0,cumsum(px[-n]))
		Ti=pmax(s+px-dx,0)*wx
		return(sum(Ti))
	}
}
#   Copyright (c) 2014-2015 by Martin Zaefferer, Cologne University of Applied Sciences

###################################################################################
#' Interchange Distance for Permutations
#' 
#' The interchange distance is an edit-distance, counting how many edit operation (here: interchanges, i.e., transposition of two arbitrary elements) have to be
#' performed to transform permutation x into permutation y.
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @references Schiavinotto, Tommaso, and Thomas Stuetzle. "A review of metrics on permutations for search landscape analysis." Computers & operations research 34.10 (2007): 3143-3153.
#'
#' @examples
#' x <- 1:5
#' y <- c(1,4,3,2,5)
#' distancePermutationInterchange(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationInterchange)
#'
#' @export
#' @useDynLib CEGO
###################################################################################
distancePermutationInterchange <- function(x, y){
	N<-length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	#x<-y[order(x)]
	result <- .Call("permutationDistanceInterchange", as.integer(x),as.integer(y), PACKAGE="CEGO")
	(N-result) / (N-1) 
}

###################################################################################
#' Longest Common Subsequence Distance for Permutations
#' 
#' DEPRECATED, see \code{\link{distancePermutationInsert}}. 
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @export
#' @keywords internal
###################################################################################
distancePermutationLCSeq<- function(x, y){
	.Deprecated("distancePermutationInsert")
	#N <- length(x)
	#if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	#result <- .Call("permutationDistanceLongestCommonSubsequence", as.integer(x),as.integer(y), PACKAGE="CEGO")
	#(N-result)/(N-1) # N-1 is for PERMUTATIONS only, because two permutations have always at least a common string/sequence of length one. 
                  # this will be different for strings where two strings may have no single character in common. permutations always contain all 
                  # of N characters.
	### 
	## insert is identical but faster
	###
	distancePermutationInsert(x,y)									
}

###################################################################################
#' Longest Common Substring Distance for Permutations
#' 
#' Distance of permutations. Based on the longest string of adjacent elements that two permutations have in common.
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' #' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @references Hirschberg, Daniel S. "A linear space algorithm for computing maximal common subsequences." Communications of the ACM 18.6 (1975): 341-343.
#'
#' @examples
#' x <- 1:5
#' y <- c(5,1,2,3,4)
#' distancePermutationLCStr(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationLCStr)
#'
#' @export
#' @useDynLib CEGO
###################################################################################
distancePermutationLCStr <- function(x, y){
	N <- length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	result <- .Call("permutationDistanceLongestCommonSubstring",  as.integer(x),as.integer(y), PACKAGE="CEGO") 
	(N-result)/(N-1)
}

###################################################################################
#' Levenshtein Distance for Permutations
#' 
#' Levenshtein Distance, often just called "Edit Distance". The number of insertions, substitutions or deletions to turn one permutation (or string of equal length) into another.
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @references Levenshtein, Vladimir I. "Binary codes capable of correcting deletions, insertions and reversals." Soviet physics doklady. Vol. 10. 1966.
#'
#' @examples
#' x <- 1:5
#' y <- c(1,2,5,4,3)
#' distancePermutationLevenshtein(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationLevenshtein)
#'
#' @export
#' @useDynLib CEGO
###################################################################################
distancePermutationLevenshtein <- function(x, y){
	N <- length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	.Call("permutationDistanceLevenshtein",  as.integer(x),as.integer(y), PACKAGE="CEGO") / N
}

###################################################################################
#' Swap-Distance for Permutations
#' 
#' The swap distance is an edit-distance, counting how many edit operation (here: swaps, i.e., transposition of two adjacent elements) have to be
#' performed to transform permutation x into permutation y.
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @references Schiavinotto, Tommaso, and Thomas Stuetzle. "A review of metrics on permutations for search landscape analysis." Computers & operations research 34.10 (2007): 3143-3153.
#'
#' @examples
#' x <- 1:5
#' y <- c(1,2,3,5,4)
#' distancePermutationSwap(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationSwap)
#'
#' @export
#' @useDynLib CEGO
###################################################################################
distancePermutationSwap <- function(x, y){
	N <- length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	result <- .Call("permutationDistanceSwap", as.integer(x),as.integer(y),  PACKAGE="CEGO")
	2*result / (N^2 - N)
}

###################################################################################
#' R-Distance for Permutations
#' 
#' R distance or unidirectional adjacency distance. Based on count of number of times that a two element sequence in x also occurs in y, in the same order.
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @references Sevaux, Marc, and Kenneth Soerensen. "Permutation distance measures for memetic algorithms with population management." Proceedings of 6th Metaheuristics International Conference (MIC'05). 2005.
#' @references Reeves, Colin R. "Landscapes, operators and heuristic search." Annals of Operations Research 86 (1999): 473-490.
#'
#' @examples
#' x <- 1:5
#' y <- c(1,2,3,5,4)
#' distancePermutationR(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationR)
#'
#' @export
#' @useDynLib CEGO
###################################################################################
distancePermutationR <- function(x, y){
	N<-length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	result <- .Call("permutationDistanceR", as.integer(x), as.integer(y), PACKAGE="CEGO")
	result / (N - 1)
}

###################################################################################
#' Adjacency Distance for Permutations
#' 
#' Bi-directional adjacency distance for permutations, depending on how often two elements are neighbours in both permutations x and y.
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @references Sevaux, Marc, and Kenneth Soerensen. "Permutation distance measures for memetic algorithms with population management." Proceedings of 6th Metaheuristics International Conference (MIC'05). 2005.
#' @references Reeves, Colin R. "Landscapes, operators and heuristic search." Annals of Operations Research 86 (1999): 473-490.
#'
#' @examples
#' x <- 1:5
#' y <- 5:1
#' distancePermutationAdjacency(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationAdjacency)
#'
#' @export
#' @useDynLib CEGO
###################################################################################
distancePermutationAdjacency <- function(x, y){
	N<-length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	result <- .Call("permutationDistanceAdjacency", as.integer(x), as.integer(y), PACKAGE="CEGO")
	(N-result-1) / (N - 1)
}

###################################################################################
#' Position Distance for Permutations
#' 
#' Position distance (or Spearmans Correlation Coefficient), scaled to values between 0 and 1. 
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @references Schiavinotto, Tommaso, and Thomas Stuetzle. "A review of metrics on permutations for search landscape analysis." Computers & operations research 34.10 (2007): 3143-3153.
#' @references Reeves, Colin R. "Landscapes, operators and heuristic search." Annals of Operations Research 86 (1999): 473-490.
#'
#' @examples
#' x <- 1:5
#' y <- c(1,3,5,4,2)
#' distancePermutationPosition(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationPosition)
#'
#' @export
#' @useDynLib CEGO
###################################################################################
distancePermutationPosition <- function(x, y){
	N<-length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	dis <- .Call("permutationDistancePosition", as.integer(x), as.integer(y), PACKAGE="CEGO")
	if(N%%2) #scale to [0;1] in case of odd N
		dis <- dis / ((N^2-1)/ 2)
	else #scale to [0;1] in case of even N
		dis <- dis / (N^2 / 2)
	dis	
}

###################################################################################
#' Squared Position Distance for Permutations
#' 
#' Squared position distance (or Spearmans Footrule), scaled to values between 0 and 1. 
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @references Schiavinotto, Tommaso, and Thomas Stuetzle. "A review of metrics on permutations for search landscape analysis." Computers & operations research 34.10 (2007): 3143-3153.
#' @references Reeves, Colin R. "Landscapes, operators and heuristic search." Annals of Operations Research 86 (1999): 473-490.
#'
#' @examples
#' x <- 1:5
#' y <- c(1,3,5,4,2)
#' distancePermutationPosition2(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationPosition2)
#'
#' @export
#' @useDynLib CEGO
###################################################################################
distancePermutationPosition2 <- function(x, y){
	N<-length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	result <- .Call("permutationDistancePosition2", as.integer(x), as.integer(y), PACKAGE="CEGO")
	result / ((N^3-N)/3) 
}

###################################################################################
#' Hamming Distance for Permutations
#' 
#' Hamming distance for permutations, scaled to values between 0 and 1.
#' That is, the number of unequal elements of two permutations, divided by the permutations length.
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @examples
#' x <- 1:5
#' y <- c(5,1,2,3,4)
#' distancePermutationHamming(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationHamming)
#'
#' @export
###################################################################################
distancePermutationHamming <- function(x, y){
	N <- length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	sum(x != y)/N
}
#distancePermutationHamming <- function(x, y){
#	N <- length(x)
#	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
#	.Call("permutationDistanceHamming", as.integer(x), as.integer(y),PACKAGE="CEGO") / N #this works,but is actually slower than a pure R implementation
#}

###################################################################################
#' Euclidean Distance for Permutations
#' 
#' Euclidean distance for permutations, scaled to values between 0 and 1:
#' \deqn{d(x,y) = frac(1){r} sqrt(\sum_{i=1}^n (x_i - y_i)^2) }{ d(x,y) = 1/r * sqrt(\sum_{i=1}^n  (x_i - y_i)^2)}
#' where n is the length of the permutations x and y, and scaling factor \eqn{r=sqrt(2*4*n*(n+1)*(2*n+1)/6)} (if n is odd)
#' or \eqn{r=sqrt(2*n*(2*n-1)*(2*n+1)/3)} (if n is even).
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @examples
#' x <- 1:5
#' y <- c(5,1,2,3,4)
#' distancePermutationEuclidean(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationEuclidean)
#'
#' @export
###################################################################################
distancePermutationEuclidean <- function(x, y){
	N<-length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	#result <- .Call("permutationDistanceEuclidean", as.integer(x), as.integer(y), PACKAGE="CEGO") #this works,but is actually slower than a pure R implementation, as below. at least for permutation length 30
	#mdis <- sqrt(result)
	mdis <- sqrt(sum((x-y)^2))
	if(N%%2){ #scale to [0;1] in case of odd N
		n=(N-1)/2
		dis <- mdis / sqrt(8*n*(n+1)*(2*n+1)/6)
	}else{ #scale to [0;1] in case of even N
		n=N/2
		dis <- mdis / sqrt(2*n*(4*n^2-1)/3)
	}	
	dis
}


###################################################################################
#' Manhattan Distance for Permutations
#' 
#' Manhattan distance for permutations, scaled to values between 0 and 1:
#' \deqn{d(x,y) = frac(1){r} \sum_{i=1}^n |x_i - y_i| }{ d(x,y) = 1/r * \sum_{i=1}^n  |x_i - y_i|}
#' where n is the length of the permutations x and y, and scaling factor \eqn{r=(n^2-1)/2} (if n is odd)
#' or \eqn{r=((n^2)/2} (if n is even).
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @examples
#' x <- 1:5
#' y <- c(5,1,2,3,4)
#' distancePermutationManhattan(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationManhattan)
#'
#' @export
###################################################################################
distancePermutationManhattan <- function(x, y){
	N <- length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	mdis<-sum(abs(x-y))
	if(N%%2) #scale to [0;1] in case of odd N
		dis <- mdis / ((N^2-1)/ 2)
	else #scale to [0;1] in case of even N
		dis <- mdis / (N^2 / 2)
	dis # TODO is the scaling equivalent to "dis/round(N^2/2)" ? if yes, what is faster?
}

###################################################################################
#' Chebyshev Distance for Permutations
#' 
#' Chebyshev distance for permutations. Specific to permutations is only the scaling to values of 0 to 1:
#' \deqn{d(x,y) = \frac{max(|x - y|) }{ (n-1) } }{d(x,y) = max(|x - y|) / (n-1)}
#' where n is the length of the permutations x and y.
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @examples
#' x <- 1:5
#' y <- c(5,1,2,3,4)
#' distancePermutationChebyshev(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationChebyshev)
#'
#' @export
###################################################################################
distancePermutationChebyshev <- function(x, y){
	N <- length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	max(abs(x-y)) / (N-1)
}

###################################################################################
#' Lee Distance for Permutations
#' 
#' Usually a string distance, with slightly different definition.
#' Adapted to permutations as: 
#' \deqn{d(x,y) = \sum_{i=1}^n min(|x_i - y_i|), n- |x_i - y_i|)}{ d(x,y) = \sum_{i=1}^n  min(|x_i - y_i|), n- |x_i - y_i|)}
#' where n is the length of the permutations x and y.
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @references Lee, C., "Some properties of nonbinary error-correcting codes," Information Theory, IRE Transactions on, vol.4, no.2, pp.77,82, June 1958
#'
#' @examples
#' x <- 1:5
#' y <- c(5,1,2,3,4)
#' distancePermutationLee(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationLee)
#'
#' @export
#' @useDynLib CEGO
###################################################################################
distancePermutationLee <- function(x,y){
	N <- length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	mdis <- .Call("permutationDistanceLee", as.integer(x), as.integer(y), PACKAGE="CEGO") #this works,but is actually slower than a pure R implementation, as below. at least for permutation length 30
	#mdis <- sum(pmin(abs(x-y),N-abs(x-y)))
	if(N%%2) #scale to [0;1] in case of odd N
		dis <- mdis / ((N^2-1)/ 2)
	else #scale to [0;1] in case of even N
		dis <- mdis / (N^2 / 2)
	dis
}

###################################################################################
#' Insert Distance for Permutations
#' 
#' The Insert Distance is an edit distance. It counts the minimum number of delete/insert operations
#' required to transform one permutation into another. A delete/insert operation shifts one element to a new position.
#' All other elements move accordingly to make place for the element. E.g., the following shows a single delete/insert move that
#' sorts the corresponding permutation: 1 4 2 3 5 -> 1 2 3 4 5. This distance is also called Ulam's metric and can as well
#' be interpreted to be based on the longest common subsequence of two permutations. 
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @references Schiavinotto, Tommaso, and Thomas Stuetzle. "A review of metrics on permutations for search landscape analysis." Computers & operations research 34.10 (2007): 3143-3153.
#' @references Wikipedia contributors, "Longest increasing subsequence", Wikipedia, The Free Encyclopedia, 12 November 2014, 19:38 UTC, <http://en.wikipedia.org/w/index.php?title=Longest_increasing_subsequence&oldid=633565014> [accessed 13 November 2014] 
#'
#' @examples
#' x <- 1:5
#' y <- c(5,1,2,3,4)
#' distancePermutationInsert(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationInsert)
#'
#' @export
#' @useDynLib CEGO
###################################################################################
distancePermutationInsert <- function(x, y){
	N<-length(x) 
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	#x<-order(y)[x] #Important note: this line means "p1 * p2^{-1}". Schiavinotto say it should be  "p2^{-1} * p1" . That seems to be a typo/error. (of course it does not matter which permutation is is inverted. but the inverted permutation should be at the right.). Is now calculated in c code
	result <- .Call("permutationDistanceInsert", as.integer(x), as.integer(y), PACKAGE="CEGO")
	(N-result)/(N-1)
}

###################################################################################
#' Cosine Distance for Permutations
#' 
#' The Cosine distance for permutations is derived from the Cosine similarity measure
#' which has been applied in fields like text mining.
#' It is based on the scalar product of two vectors (here: permutations).
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @references Singhal, Amit (2001)."Modern Information Retrieval: A Brief Overview". Bulletin of the IEEE Computer Society Technical Committee on Data Engineering 24 (4): 35-43
#'
#' @examples
#' x <- 1:5
#' y <- c(5,1,2,3,4)
#' distancePermutationCos(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationCos)
#'
#' @export
###################################################################################
distancePermutationCos <- function(x,y){
	N=length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	as.numeric(1-((x%*%y)-(N*(N+1)*(N+2)/6))/((N^3-N)/6)) # |a| = |b| because permutation -> |a| * |a| = sum_{i=1}^n i^2 = (N*(N+1)*(2*N+1)/6)
}#the minus factor is tetrahedral pyramid




###################################################################################
#' Lexicographic permutation distance
#' 
#' This function calculates the lexicographic permutation distance. That is the difference of positions
#' that both positions would receive in a lexicographic ordering. Note, that this distance
#' measure can quickly become inaccurate if the length of the permutations grows too large, due
#' to being based on the factorial of the length. In general, permutations longer than 100 elements should
#' be avoided.
#'
#' @param x first permutation (integer vector)
#' @param y second permutation (integer vector)
#'
#' @return numeric distance value \deqn{d(x,y)}, scaled to values between 0 and 1 (based on the maximum possible distance between two permutations)
#'
#' @seealso \code{\link{lexicographicPermutationOrderNumber}}
#'
#' @examples
#' x <- 1:5
#' y <- c(1,2,3,5,4)
#' distancePermutationLex(x,y)
#' p <- replicate(10,sample(1:5),simplify=FALSE)
#' distanceMatrix(p,distancePermutationLex)
#'
#' @export
###################################################################################
distancePermutationLex <- function(x,y){
	N <- length(x)
	if(N!=length(y)|!is.numeric(x)|!is.numeric(y)) stop("Incorrect input to distance function, only permutations of same length are allowed.")
	abs(lexicographicPermutationOrderNumber(x) - lexicographicPermutationOrderNumber(y)) / (factorial(N)-1)
}

###################################################################################
#' Lexicographic order number
#' 
#' This function returns the position-number that a permutation would receive in a lexicographic ordering.
#' It is used in the lexicographic distance measure.
#'
#' @param x permutation (integer vector)
#'
#' @return numeric value giving position in lexicographic order.
#'
#' @seealso \code{\link{distancePermutationLex}}
#'
#' @examples
#' lexicographicPermutationOrderNumber(1:5)
#' lexicographicPermutationOrderNumber(c(1,2,3,5,4))
#' lexicographicPermutationOrderNumber(c(1,2,4,3,5))
#' lexicographicPermutationOrderNumber(c(1,2,4,5,3))
#' lexicographicPermutationOrderNumber(c(1,2,5,3,4))
#' lexicographicPermutationOrderNumber(c(1,2,5,4,3))
#' lexicographicPermutationOrderNumber(c(1,3,2,4,5))
#' lexicographicPermutationOrderNumber(5:1)
#' lexicographicPermutationOrderNumber(1:7)
#' lexicographicPermutationOrderNumber(7:1)
#'
#' @export
#' @useDynLib CEGO
###################################################################################
lexicographicPermutationOrderNumber <- function(x){
	N <- length(x)
	if(!is.numeric(x)) stop("Incorrect input to lexicographicPermutationOrderNumber, should be permutation.")
	result <- .Call("lexPermOrder", as.integer(x), PACKAGE="CEGO")
	sum(result*factorial((N-1):0))+1
}

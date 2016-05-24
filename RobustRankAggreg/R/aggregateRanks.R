

#' Create rank matrix
#' 
#' Convert a set of ranked lists into a rank matrix
#' 
#' The lists are converted to a format that is used by aggregateRanks. If partial
#' rankings are given to the function, all the missing values are subtituted by the
#' maximum rank N, which can be specified manually. This parameter has a very strong
#' influence on the performance of RRA algorithm, therfore it should be reasonably
#' accurate. If the N is different for the gene lists, it can be also given as a vector. 
#' 
#' Parameter full is used, when full rankings are given, but the sets of ranked elements
#' do not match perfectly. Then the structurally missing values are substituted with
#' NA-s.
#'
#' @param glist list of preference lists
#' @param N number of all rankable elements
#' @param full logical showing if the given rankings are complete
#' @return A matrix, with as many columns as input rankings and rows as unique elements
#' in all the rankings combined.
#' @author  Raivo Kolde \email{rkolde@@gmail.com}
#' @examples
#' # Make sample input data
#' glist <- list(sample(letters, 4), sample(letters, 10), sample(letters, 12))
#' 
#' r = rankMatrix(glist)
#' r = rankMatrix(glist, full = TRUE)
#' 
#' # Use real data
#' data(cellCycleKO)
#' r = rankMatrix(cellCycleKO$gl, N = cellCycleKO$N)
#' 
#' @export
rankMatrix <-  function(glist, N = NA, full = FALSE){
	u = unique(c(glist, recursive = TRUE))
	
	if(all(is.na(N))){
		N = length(u)
	}
	if(!full){
		rmat = matrix(1, nrow = length(u), ncol = length(glist), dimnames = list(u, names(glist)))
		if(length(N) == 1){
			N = rep(N, ncol(rmat))
		}
	}
	else{
		rmat = matrix(NA, nrow = length(u), ncol = length(glist), dimnames = list(u, names(glist)))
		N = unlist(lapply(glist, length))
	}
	for(i in 1:length(glist)){
		rmat[glist[[i]], i] = (1:length(glist[[i]])) / N[i]
	}
	return(rmat)
}


# Output function
formatOutput <- function(scores, score.names, ordering = "ascending"){
	res = data.frame(Name = score.names, Score = scores)
	if(ordering == "ascending"){
		res = res[order(res$Score), ]
	}
	else{
		res = res[order(res$Score, decreasing = TRUE), ]
	}
	return(res)
} 


# Stuart-Aerts method helper functions
sumStuart <- function(v, r){
	k = length(v)
	l_k = 1:k
	ones = (-1)**(l_k + 1)
	f = factorial(l_k)
	p = r ** l_k
	return(ones %*% (rev(v) * p / f))
}

qStuart <- function(r){
	N = sum(!is.na(r))
	v = rep(1, N + 1)
	for(k in 1:N){
		v[k + 1] = sumStuart(v[1:k], r[N - k + 1])
	}
	return(factorial(N) * v[N + 1])
}

stuart <- function(rmat){
	rmat <- t(apply(rmat, 1, sort, na.last = TRUE))
	return(apply(rmat, 1, qStuart))
}

# RRA helper functions
 
#' Calculate beta scores
#' 
#' Calculate the beta scores for normalized rank vector.
#' 
#' Takes in a vector with values in [0, 1]. It sorts the values to get the order
#' statistics and calculates p-values for each of the order statistics. These are based 
#' on their expected distribution under the null hypothesis of uniform distribution. 
#'
#' In RRA algorithm context the inputs are supposed to be normalized ranks. However, 
#' p-values in general follow the uniform distribution, therefore it can be used with any 
#' kind of p-value vectors, to see if there are more small values than expected. 
#' 
#' The NA values are removed before calculation and all results take into account only 
#' existing values.  
#' 
#' @param r vector of values in [0, 1]
#' @return  The functions returns a vector of p-values, that correspond to the sorted
#' input vector. The NA-s are pushed to the end.
#' 
#' @references  Kolde et al "Robust Rank Aggregation for gene list integration and 
#' meta-analysis" (in preparation)
#' 
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#'  betaScores(c(runif(15)))
#'  betaScores(c(runif(10), rbeta(5, 1, 50)))
#' 
#' @export
betaScores <- function(r){
	n <- sum(!is.na(r))
	p <- rep(1, n)
	r <- sort(r, na.last = TRUE)
	p <- pbeta(r, 1:n, n - 1:n + 1)
	return(p)
} 



thresholdBetaScore <- function(r, k = seq_along(r), n = length(r), sigma = rep(1,n)){
	if(length(sigma) != n) stop("The length of sigma does not match n")
	if(length(r) != n) stop("The length of pvalues does not match n")
	if(min(sigma)< 0 || max(sigma) > 1) stop("Elements of sigma are not in the range [0,1]")
   if(any(!is.na(r) & r > sigma)) stop("Elements of r must be smaller than elements of sigma")

	x <- sort(r[!is.na(r)])
	sigma <- sort(sigma, decreasing = TRUE)
	beta <- rep(NA, length(k))
	for(i in seq_along(k))
	{
		if(k[i] > n)
		{
			beta[i] <- 0
			next;
		}
		if(k[i] > length(x))
		{
			beta[i] <- 1
			next;
		}
		if(sigma[n] >= x[k[i]])
		{
			beta[i] <- pbeta(x[k[i]], k[i], n + 1 - k[i])
			next;
		}
		
		# Non-trivial cases
		# Find the last element such that sigma[n0] <= x[k[i]]
		n0 <- which(sigma < x[k[i]])[1] - 1
		
		# Compute beta score vector beta(n,k) for n = n0 and k = 1..k[i] 
		if(n0 == 0) 
		{ 
			B <- c(1, rep(0, k[i]))
		} else if(k[i] > n0)
		{
			B <- c(1, pbeta(x[k[i]], 1 : n0, n0 : 1), rep(0, k[i] - n0))
		} else 
		{ 
			B <- c(1, pbeta(x[k[i]], 1 : k[i], n0 + 1 - c(1 : k[i])))
		}

		# In the following update steps sigma < x[k[i]] 
		z <- sigma[(n0 + 1) : n]
		for(j in seq_len(n - n0))
		{
  			B[2 : (k[i] + 1)] <- (1 -z[j]) * B[2 : (k[i] + 1)] + z[j] * B[1 : k[i]] 
		}
		
		beta[i] <- B[k[i]+1]
	}
	
	names(beta) <- k
	return(beta)
}


correctBetaPvalues <- function(p, k){
	p <- min(p * k, 1) 
	
	return(p)
}

correctBetaPvaluesExact <- function(p, k){
	rm = 1 - t(sapply(p, qbeta, 1:k, k - 1:k + 1))
	res = 1 - stuart(rm)
	return(res)
}

#' Calculate rho scores
#' 
#' Calculate Rho score for normalized rank vector
#' 
#' Takes in a vector with values in [0, 1]. Applies \code{\link{betaScores}} to the vector, takes the minimum of the beta scores and converts it to a valid p-value. 
#'
#' @param r vector of values in [0, 1]
#' @param topCutoff a vector of cutoff values used to limit the number of elements in the 
#' input lists
#' @param exact indicator if exact p-values should be calculated (Warning: it is computationally unstable and does ot give considerable gain)
#' @references  Kolde et al "Robust Rank Aggregation for gene list integration and 
#' meta-analysis" (in preparation)
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#'  rhoScores(c(runif(15)))
#'  rhoScores(c(runif(10), rbeta(5, 1, 50)))
#' 
#' @export
rhoScores <- function(r, topCutoff = NA, exact = F){
	if(is.na(topCutoff[1])){
		x <- betaScores(r)
	}
	else{
		r <- r[!is.na(r)]
		r[r == 1] <- NA
		x <- thresholdBetaScore(r, sigma = topCutoff)
	}
	if(exact){
		rho <- correctBetaPvaluesExact(min(x, na.rm = T), k = sum(!is.na(x)))
	}
	else{
		rho <- correctBetaPvalues(min(x, na.rm = T), k = sum(!is.na(x)))
	}
	
	return(rho)
}

# The dynamic algorithm for more accurate BetaScore calculation


#' Aggregate ranked lists
#' 
#' Method implementing various gene list aggregation methods, most notably Robust Rank 
#' Aggregation.
#' 
#' All the methods implemented in this function make an assumtion that the number of
#' ranked items is known. This assumption is satisfied for example in the case of 
#' gene lists (number of all genes known to certain extent), but not when aggregating 
#' results from google searches (there are too many web pages). This parameter N can be 
#' set manually and has strong influence on the end result. The p-values from RRA 
#' algorithm can be trusted only if N is close to the real value.
#' 
#' The rankings can be either full or partial. Tests with the RRA algorithm show that one 
#' does not lose too much information if only top-k rankings are used. The missing values 
#' are assumed to be equal to maximal value and that way taken into account 
#' appropriately. 
#' 
#' The function can handle also the case when elements of the different rankings do not 
#' overlap perfectly. For example if we combine resutls from different microarray 
#' platforms with varying coverage. In this case these structurally missing values are 
#' substituted with NA-s and handled differently than omitted parts of the rankings. 
#' The function accepts as an input either list of rankings or rank matrix based on them. 
#' It converts the list to rank matrix automatically using the function 
#' \code{\link{rankMatrix}}. For most cases the ranking list is more convenient. Only 
#' in complicated cases, for example with top-k lists and structural missing values one 
#' would like to construct the rank matrix manually.  
#' 
#' When the number of top elements included into input is specified in advance, for 
#' example some lists are limited to 100 elements, and the lengths of these lists differ 
#' significantly, we can use more sensitive and accurate algorithm for the score 
#' calculation. Then one has to specify in the input also the parameter topCutoff, which 
#' is a vector defining an cutoff value for each input list. For example if we have three 
#' lists of 1000 elements but first is limited to 100, second 200 and third to 900 
#' elements, then the topCutoff parameter should be c(0.1, 0.2, 0.9).
#' 
#' @param glist list of element vectors, the order of the vectors is used as the ranking.
#' @param rmat the rankings in matrix format. The glist is by default converted to this 
#' format.
#' @param N the number of ranked elements, important when using only top-k ranks, by 
#' default it is calculated as the number of unique elements in the input.
#' @param method rank aggregation method, by defaylt \code{'RRA'}, other options are 
#' \code{'min'}, \code{'geom.mean'}, \code{'mean'}, \code{'median'} and \code{'stuart'} 
#' @param full indicates if the full rankings are given, used if the the sets of ranked 
#' elements do not match perfectly
#' @param exact indicator showing if exact p-value will be calculated based on rho score (Default: if number of lists smaller than 10, exact is used)
#' @param topCutoff a vector of cutoff values used to limit the number of elements in the 
#' input lists
#' elements do not match perfectly
#' @return  Returns a two column dataframe with the element names and associated scores 
#' or p-values.
#' @references  Kolde et al "Robust Rank Aggregation for gene list integration and 
#' meta-analysis" (in preparation)
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#' # Make sample input data
#' glist <- list(sample(letters, 4), sample(letters, 10), sample(letters, 12))
#' 
#' # Aggregate the inputs
#' aggregateRanks(glist = glist, N = length(letters))
#' aggregateRanks(glist = glist, N = length(letters), method = "stuart")
#' 
#' # Since we know the cutoffs for the lists in advance (4, 10, 12) we can use
#' # the more accurate algorithm with parameter topCutoff
#' 
#' # Use the rank matrix instead of the gene lists as the input
#' r = rankMatrix(glist)
#' 
#' aggregateRanks(rmat = r)
#' 
#' # Example, when the input lists represent full rankings but the domains do not match 
#' glist <- list(sample(letters[4:24]), sample(letters[2:22]), sample(letters[1:20]))
#' r = rankMatrix(glist, full = TRUE) 
#' head(r)
#' 
#' aggregateRanks(rmat = r, method = "RRA")
#' 
#' # Dataset representing significantly changed genes after knockouts 
#' # of cell cycle specific trancription factors
#' data(cellCycleKO)
#' r = rankMatrix(cellCycleKO$gl, N = cellCycleKO$N)
#' ar = aggregateRanks(rmat = r)
#' head(ar)
#' 
#' @aliases RobustRankAggreg
#' 
#' @export
aggregateRanks <-  function(glist, rmat = rankMatrix(glist, N, full = full), N = NA, method = "RRA", full = FALSE, exact = F, topCutoff = NA){

	if(!(method %in% c("mean", "min", "median", "geom.mean", "RRA", "stuart"))){
		stop("method should be one of:  'min', 'geom.mean', 'mean', 'median', 'stuart' or 'RRA' ")
	}
	
	if(is.na(N)){
		N <- nrow(rmat)
	}
	
	if(is.null(rownames(rmat))){
		rownames(rmat) <- 1:nrow(rmat)
	}
	
	if(method == "min"){
		a <-  apply(rmat, 1, min, na.rm = TRUE)
		return(formatOutput(scores = a, score.names = names(a), ordering = "ascending"))
	}
	if(method == "median"){
		a <- apply(rmat, 1, median, na.rm = TRUE)
		return(formatOutput(scores = a, score.names = names(a), ordering = "ascending"))
	}
	if(method == "geom.mean"){
		a <-  apply(rmat, 1, function(x) exp(mean(log(x), na.rm = TRUE)))
		return(formatOutput(scores = a, score.names = names(a), ordering = "ascending"))
	}
	if(method == "RRA"){
		a = apply(rmat, 1, rhoScores, topCutoff = topCutoff, exact = exact)
		names(a) <- rownames(rmat)
		return(formatOutput(scores = a, score.names = names(a), ordering = "ascending"))
	}
	if(method == "mean"){
		a <- apply(rmat, 1, mean, na.rm = TRUE)
		n <- apply(rmat, 1, function(x) sum(!is.na(x)))
		b <- pnorm(a, 0.5, sqrt(1/12/n))
		return(formatOutput(scores = b, score.names = names(a), ordering = "ascending"))
	}
	if(method == "stuart"){
		a <- stuart(rmat)
		return(formatOutput(scores = a, score.names = names(a), ordering = "ascending"))
	}
}
# 
# require(foreach)
# space = paste("X", 1:10000, sep = "")
# gl = foreach(i = 1:10) %do% {sample(space)[1:1000]} 
# ar = aggregateRanks(gl, exact = T, N = 10000)
# ar2 = aggregateRanks(gl, exact = F, N = 10000)
# quartz(); hist(ar[, 2])

#' A dataset based on Reimand \emph{et al} and Hu \emph{et al}. It contains lists  
#' yeast genes that were most influenced by 12 cell cycle related transcription factor 
#' knockouts. 
#' The dataset is a list with 3 slots
#' \enumerate{ 
#' 	\item \code{gl} - set of gene lists in a format suitable for 
#' \code{\link{aggregateRanks}};
#' 	\item \code{N} - number of yeast genes; 
#' 	\item \code{ref} - reference list of cell cycle related genes taken from  de 
#' Lichtenberg \emph{et al}.
#' }
#'
#' @name cellCycleKO
#' @docType data
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @references Reimand, J., Vaquerizas, J. M., Todd, A. E., Vilo, J., and Luscombe, N. M.
#' (2010). "Comprehensive reanalysis of transcription factor knockout expression data
#' in saccharomyces cerevisiae reveals many new targets. Nucleic Acids Res."
#' 
#' Hu, Z., Killion, P. J., and Iyer, V. R. (2007). "Genetic reconstruction of 
#' a functional transcriptional regulatory network." Nat. Genet., 39(5), 683-7
#' 
#' de Lichtenberg, U., Jensen, L. J., Fausboll, A., Jensen, T. S., Bork, P., 
#' and Brunak, S. (2005). "Comparison of computational methods for the identification of 
#' cell cycle- regulated genes. Bioinformatics, 21(7), 1164-71."
#' @keywords data
NULL






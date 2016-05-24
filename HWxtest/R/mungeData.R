# Functions for arranging data used in HW test
# (c) William R. Engels, 2014
#' 

#' @title 
#' mungeData
#' @description
#' Utility functions for handling genotype counts and arranging data
#' 

#' @details
#' Interconvert between different formats for genotype counts. 
#' 
#' Let \code{k} be the number of alleles:
#' \itemize{

	#' \item \code{clearUpper} fills the upper-right half of the \eqn{k x k} matrix with \code{NA}
	#' \item \code{fillUpper} makes the \eqn{k x k} matrix symmetrical by filling the upper-right half with numbers from the lower half.
	#' \item \code{vec.to.matrix} converts genotype counts in vector form and returns a matrix. The vector must have \eqn{k(k+1)/2} non-negative integers.
	#' \item  \code{matrix.to.vec} converts a \eqn{k x k} matrix of genotype counts to a vector of length \eqn{k(k+1)/2}
	#' \item \code{alleleCounts} returns a vector of length \eqn{k} containing the numbers of each allele. The sum of this vector will be twice the number of diploids in the sample.
	#' \item \code{remove.missing.alleles} returns a matrix with no \code{0}'s for allele counts
	#' \item \code{df.to.matrices} converts a data frame to a list of genotype count matrices. The data frame should be of the kind produced in the package \code{adegenet} with \code{genind2df}
#' }
#' 
#' 
#' @param gmat a matrix of non-negative integers representing genotype counts. In a matrix of genotype counts, \code{a[i,j]} and \code{a[j,i]} both represent the same heterozygote. Only the lower-left half of \code{gmat} is used. Numbers along the diagonal represent counts of the homozygotes.
#' 
#' @param gvec vector containing \code{k(k+1)/2} genotype counts. All non-negative integers.  Genotype counts should be in the order: \code{a11, a21, a22, a31, a32, ..., akk}
#' 
#' @param alleleNames an optional list of names for the alleles. The length should be \eqn{k}
#' @param df a dataframe containing individual genotypes. Each row represents an individual. The first column, named \dQuote{pop} names the population. Each other column is named for a particular locus. The genotypes are as \dQuote{123/124}
#' @param sep For a dataframe, this is the separator character. typically \dQuote{/}
#' 
#' 
#' @examples
#' gvec <- c(0,3,1,5,18,1,3,7,5,2)
#' gmat <- vec.to.matrix(gvec, alleleNames=letters[1:4])
#' alleleCounts(gmat)
#' 
#' 
#' @rdname mungeData
#' @export
fillUpper <- 
function(gmat){
	if(!(is.matrix(gmat) && (nrow(gmat)==ncol(gmat)))) stop("Must be square matrix at least 2x2")
	k <- nrow(gmat);
	if(k<2) return(gmat)
	for(j in 2:k) {gmat[1:(j-1), j] <- gmat[j,1:(j-1)]};
	gmat
}


#' @rdname mungeData 
#' @export
alleleCounts <- 
function(gmat) {
	if(class(gmat)!="matrix") gmat <- vec.to.matrix(gmat)
	t <- fillUpper(gmat);
	k <- dim(t)[1];
	m <- integer(k);
	for(i in 1:k) {m[i] <- sum(t[i,]) + t[i,i]};
	if(!is.null(rownames(gmat))) names(m) <- rownames(gmat)
	m
}


#' @title vec.to.matrix
#' @rdname mungeData
#' @export
vec.to.matrix <- 
function(gvec, alleleNames=""){
	if(!(is.vector(gvec) && is.numeric(gvec))) stop("\nMust be a vector")
	nGenotypes <- length(gvec)
	nAlleles <- as.integer((sqrt(8*nGenotypes + 1) - 1)/2)
	if(nGenotypes != nAlleles*(nAlleles + 1)/2) stop("\nWrong number of genotype counts")
	t <- matrix(NA, nAlleles, nAlleles)
	for(i in 1:nAlleles){t[i, 1:i] <- gvec[(i*(i-1)/2 + 1):(i*(i+1)/2)]}
	if(length(alleleNames)>=nAlleles) {
		rownames(t) <- alleleNames;
		colnames(t) <- alleleNames;
	}
	t	
}
#' 
#' @title remove.missing.alleles
#' @description remove missing alleles
#' @details none
#' @rdname mungeData
#' @export
remove.missing.alleles <- 
function(gmat) {
	if(class(gmat)!="matrix") gmat <- vec.to.matrix(gmat)
	m <- alleleCounts(gmat)
	zm <- which(m==0)
	if(length(zm)) gmat <- gmat[-zm,-zm]
	gmat
}

#' 
#' @title matrix.to.vec
#' @description converts matrix to vector
#' @details none
#' @rdname mungeData
#' @export
matrix.to.vec <- 
function(gmat){
	if(!(is.matrix(gmat) && (nrow(gmat)==ncol(gmat)))) stop("Must be square matrix")
	v <- c();
	k <- nrow(gmat)
	for(i in 1:k){v <- append(v,gmat[i,1:i])}
	names(v) <- NULL;
	v	
}


#' @title clearUpper
#' @description Clears upper-right of matrix
#' @rdname mungeData
#' @export
clearUpper <- 
function(gmat){
	if(!(is.matrix(gmat) && (nrow(gmat)==ncol(gmat)))) stop("Must be square matrix")
	k <- nrow(gmat);
	for(j in 2:k) {gmat[1:(j-1), j] <- NA};
	gmat
}


#' @title df.to.matrices
#' @rdname mungeData
#' @export
df.to.matrices <- function(df, sep="/"){
	if(!is.data.frame(df)) stop("Must be a data frame")
	pnames <- levels(factor(df$pop)) # population names
	gnames <- colnames(df)[-1] #locus names (with "pop" removed)
	res <- list()
	for(pn in pnames){
		popn <- list()
		popsel <- df$pop==pn
		for(gn in gnames){
			anames <- levels(factor(unlist(strsplit(df[[gn]][popsel], sep))))
			k <- length(anames)
			ta <- table(factor(df[[gn]][popsel]))
			gtypes <- names(ta)
			t <- matrix(0, k,k)
			rownames(t) <- anames
			colnames(t) <- anames
			if(k > 1) {
			for(g in gtypes) {
				diploid <- unlist(strsplit(g,sep))
				t[diploid[[1]], diploid[[2]]] <- ta[[g]]
			}}
			t <- t + t(t)
			diag(t)  <- diag(t)/2
			popn[[gn]] <- NA
			if(k > 1) popn[[gn]] <- clearUpper(t)
		} # for gn
		res[[pn]] <- popn
	} #for pn
	res
}
NULL


# #' Convert a list of genotypes into a genotype count matrix
# #' 
# #' Genotype lists as are used by packages `genetics` and `adegenet` are converted to an array of genotype counts.
# #' This function requires package `genetics`
# #' 
# #' @param g List of text objects indicating genotypes. Alleles are separated by \dQuote{/}
# #' 
# #' @return matrix of \eqn{k x k} genotype counts
# #' 
# #' @examples
# #' g <- c(rep("a/b",3),
# #'		"b/b",
# #'		rep("a/c", 5),
# #'		rep("b/c", 18),
# #'		"c/c",
# #'		rep("a/d",3),
# #'		rep("b/d", 7),
# #'		rep("c/d", 5),
# #'		rep("d/d", 2))
# #' genotypeList.to.matrix(g)

# #' @export
# genotypeList.to.matrix <- 
# function(g){
	# if(!require(genetics)) stop("\ngenetics package is required")
	# x <- genotype(g);
	# tab <- table(factor(allele(x, 1), levels = allele.names(x)), factor(allele(x, 2), levels = allele.names(x)));
	# t(tab)
# }


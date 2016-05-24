#' Determines coherence
#'
#' This function determines the number of embedded absences in an interaction
#' matrix, and compares this value against null simulated matrices. Species
#' ranges should be coherent along the ordination axis, as this axis represents
#' a latent environmental gradient. A negative value of coherence (empirical
#' matrix has more embedded absences than null matrices) indicates a
#' 'checkerboard' pattern (Leibold & Mikkelson 2002). Nonsignificance has been
#' historically interpreted as being indicative of a 'random' pattern, though
#' this may be seen as accepting the null hypothesis, as nonsignificance cannot
#' be used to infer a process.
#'
#' 'method' is an argument handed to functions in the 'vegan' package. Leibold
#' & Mikkelson advocated the use of equiprobable rows and columns (provided
#' that rows and columns had at least one entry). This method is called 'r00'.
#' Methods maintaining row (site) frequencies include 'r0','r1' & 'r2'. The
#' default method argument is 'r1', which maintains the species richness of a
#' site (row totals) and fills species ranges (columns) based on their marginal
#' probabilities. Arguably the most conservative null algorithm is the fixed
#' row - fixed column total null, which is implemented as 'fixedfixed'. See the
#' help file for 'commsimulator' or Wright et al. 1998 for more information.
#'
#' If 'order' is FALSE, the interaction matrix is not ordinated, allowing the
#' user to order the matrix based on site characteristics or other biologically
#' relevant characteristics.
#'
#' This function can either be used as a standalone, or can be used through the
#' 'metacommunity()' function, which determines all 3 elements of metacommunity
#' structure (coherence, boundary clumping, & turnover) (Leibold & Mikkelson
#' 2002)
#'
#' @param comm community data in the form of a presence absence matrix
#' @param method null model randomization method used by 'nullmaker'. See
#' details below (and the help file of fucntion 'nullmaker') for more
#' information.
#' @param sims number of simulated null matrices to use in analysis
#' @param scores axis scores to ordinate matrix. 1: primary axis scores
#' (default) 2: secondary axis scores
#' @param order logical argument indicating whether to ordinate the interaction
#' matrix or not. See details.
#' @param allowEmpty logical argument indicating whether to allow null
#' matrices to have empty rows or columns
#' @param binary logical argument indicating whether to ordinate the community
#' matrix based on abundance or binary (default) data.
#' @param verbose Logical. Prints a graphical progress bar that tracks the
#' creation of null matrices. Useful for conservative null models on large
#' and/or sparse data.
#' @return A vector containing the number of embedded absences (embAbs), z-score (z), p-value (pval), mean (simulatedMean) and variance (simulatedVariance) of simulations, and null model randomization method (method).
#'
#' @author Tad Dallas
#' @export
#' @references Leibold, M.A. and G.M. Mikkelson. 2002. Coherence, species
#' turnover, and boundary clumping: elements of meta-community structure. Oikos
#' 97: 237 - 250.
#'
#' Wright, D.H., Patterson, B.D., Mikkelson, G.M., Cutler, A. & Atmar, W.
#' (1998). A comparative analysis of nested subset patterns of species
#' composition. Oecologia 113, 1-20.
#' @examples
#'
#' #define an interaction matrix
#' data(TestMatrices)
#' intmat=TestMatrices[[7]]
#'
#' #determine coherence of interaction matrix
#' coh.intmat <- Coherence(intmat, method='r1', sims=100, scores=1, order=TRUE, binary=TRUE)
#'
#' #return results
#' coh.intmat
#'
Coherence <-function(comm, method='r1', sims=1000, scores=1, order=TRUE, allowEmpty=FALSE, binary=TRUE, verbose=FALSE){

coherence <- function(web){
	zeros <- which(web==0, arr.ind=TRUE)
  ret <- matrix(0, ncol=2)
	uncols <- which(colSums(web)>1)
	for(i in 1:length(uncols)){
		temp <- zeros[which(zeros[,2] == uncols[i]),]
		tempmin <- min(which(web[,uncols[i]]==1))
		tempmax <- max(which(web[,uncols[i]]==1))
		if(length(temp) < 3){
			if(temp[1] %in% tempmin:tempmax){ret <- rbind(ret,as.vector(temp))}
		  }else{
		  temp=temp[which(temp[,1] %in% tempmin:tempmax),]
		  ret=rbind(ret,temp)
	    }
	}

	unrows <- which(rowSums(web)>1)
	for(j in 1:length(unrows)){
		temp <- zeros[which(zeros[,1]==unrows[j]),]
		tempmin <- min(which(web[unrows[j], ] == 1))
		tempmax <- max(which(web[unrows[j], ] == 1))
		if(length(temp) < 3) {
			if(temp[1] %in% tempmin:tempmax){ret=rbind(ret,as.vector(temp))}
		}else{
		temp <- temp[which(temp[,2] %in% tempmin:tempmax),]
  	ret <- rbind(ret,temp)
	  }
	}
ret <- ret[-1,]
ret <- unique(ret)

nrow(ret)
}

  if(order==TRUE){comm <- OrderMatrix(comm, scores=scores, binary=binary)
    }else{comm <- comm}
	statistic <- coherence(comm)
	nulls <- NullMaker(comm=comm, sims=sims, method=method, ordinate=order, allowEmpty=allowEmpty, verbose=verbose)

	simstat <- as.numeric(lapply(nulls,coherence))
	varstat <- sd(simstat)
	z <- (mean(simstat)-statistic)/(varstat)
	pval <- 2*pnorm(-abs(z))
	return(list(embAbs = statistic, z = z, pval = pval, simulatedMean = mean(simstat), simulatedVariance = varstat, method = method))
}

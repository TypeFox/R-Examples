#' computes logLik of two blurGenotypes 
#' 
#' Compute logLik of genotypes of person 1 
#' given genotypes of person 2 and assumed relation 
#' between the two persons (expressed with transition 
#' probability matrix; as returned with 
#' 'makeTransitionMatrix').
#' 
#' @param bGenotype1 blurred genotype of person 1
#' @param bGenotype2 blurred genotype of person 2
#' @param TransitionMatrix transition probability matrix
#' @param q vector of effect allele frequencies
#' 
#' @author Yurii Aulchenko
#' 
#' @examples 
#' data(srdta)
#' # select 10 first SNPs
#' df <- srdta[,1:10]
#' # compute effect allele freq
#' EAF <- summary(gtdata(df))$"Q.2"
#' # get genotypes of first 2 people
#' g1 <- as.numeric(df[1:2,])
#' # blur all genotypes of person 1; use HWE to infer missing
#' bg1 <- blurGenotype(g1[1,],q=EAF)
#' # blur all genotypes of person 2; use HWE to infer missing
#' bg2 <- blurGenotype(g1[2,],q=EAF)
#' # generate sib-sib transision matrices
#' trss <- makeTransitionMatrix(EAF,nmei=c(2,2))
#' getLogLikelihoodGivenRelation(bg1,bg2,trss,EAF)
#' 
getLogLikelihoodGivenRelation <- function(bGenotype1,bGenotype2,TransitionMatrix,q) {
# sanity checks
	if (class(bGenotype1) != "matrix")
		stop("bGenotype1 must be a matrix with 3 rows")
	if (class(bGenotype2) != "matrix")
		stop("bGenotype2 must be a matrix with 3 rows")
	if (any(dim(bGenotype1) != dim(bGenotype2)))
		stop("dim(bGenotype1) != dim(bGenotype2)")
	if (nrow(bGenotype1)!=3)
		stop("bGenotype1 must be a matrix with 3 rows")
	if (nrow(bGenotype2)!=3)
		stop("bGenotype2 must be a matrix with 3 rows")
	if (any(is.na(bGenotype1)))
		stop("missing values not allowed in bGenotype1")
	if (any(is.na(bGenotype2)))
		stop("missing values not allowed in bGenotype2")
	nloci <- dim(bGenotype1)[2]
	if (nloci == 1) {
		dmTM <- dim(TransitionMatrix)
		if (    (length(dmTM)!=2 && length(dmTM)!=3) || 
				(length(dmTM)==2 && any(dmTM!=c(3,3))) ||
				(length(dmTM)==3  && any(dmTM!=c(3,3,1))) )
			stop("for single locus, TransitionMatrix must be 3x3 matrix or 3x3x1 array")
	} 
	
# compute likelihood 
# P(call1,call2|nmeioses) = SUM_{g1,g2} P(call1|g1) P(call2|g2) 
#                                       P(g1|g2,nmeioses) P(g2)	
	mm <- matrix(ncol=nloci,nrow=3)
	mm[1,] <- bGenotype1[1,] * TransitionMatrix[1,1,] + 
			bGenotype1[2,] * TransitionMatrix[2,1,] + 
			bGenotype1[3,] * TransitionMatrix[3,1,]
	mm[2,] <- bGenotype1[1,] * TransitionMatrix[1,2,] + 
			bGenotype1[2,] * TransitionMatrix[2,2,] + 
			bGenotype1[3,] * TransitionMatrix[3,2,]
	mm[3,] <- bGenotype1[1,] * TransitionMatrix[1,3,] + 
			bGenotype1[2,] * TransitionMatrix[2,3,] + 
			bGenotype1[3,] * TransitionMatrix[3,3,]
	x <- mm * bGenotype2
	x <- apply(x,FUN=sum,MARGIN=2)
	out <- list(locusLik = x, logLik = sum(log(x)))
	out
}

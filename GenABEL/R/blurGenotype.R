#' blur genotype calls into probabilites
#' 
#' 'blurs' genotype calls into probabilities: translates 
#' single genotype g2, into probability distribution P(g1|g2), 
#' that is probability that true genotype is g1 
#' given g2 is the observed 'called' genotype and error rate is 
#' epsilon. Probability that 'true' genotype is called genotype 
#' is set to (1-epsilon)^2, the probability that true genotype 
#' differs at 1 allele is set to epsilon*(1-epsilon), and 
#' both allelel differ = epsilon^2.
#' 
#' @param g vector of genotypes for a particular person 
#' (at locus 1, locus 2, etc., coded as 0, 1, 2 (corresponding 
#' to genotypes AA, AB, and BB, respectively) and NA. 
#' @param q (optional) vector of coded allele freqeuncies for
#' locus 1, locus 2, etc.
#' @param epsilon error rate
#' 
#' @return matrix with columns corresponding to SNPs 
#' and rows corresponding to 'g0', 'g1', 'g2'. For 
#' a particular SNP, a vale in cell 'gK' is the 
#' probability that true genotype is 'K', given 
#' thw original call and error-rate. 
#' 
#' @author Yurii Aulchenko
#' 
#' @examples 
#' data(srdta)
#' # select 10 first SNPs
#' df <- srdta[,1:10]
#' # compute effect allele freq
#' EAF <- summary(gtdata(df))$"Q.2"
#' EAF
#' # get genotypes of first 5 people
#' g1 <- as.numeric(df[1:5,])
#' g1
#' # blur the genotype of person 1, snp 1
#' blurGenotype(g1[1,1])
#' # blur all genotypes of person 2; assume no info for missing
#' blurGenotype(g1[2,])
#' # blur all genotypes of person 2; use HWE to infer missing
#' blurGenotype(g1[2,],q=EAF)
#' 
#' 
blurGenotype <- function(g,q=NULL,epsilon=0.01) {
# sanity checks
	if (is.vector(g)) {
		g1 <- matrix(g,nrow=1);
		if (!is.null(names(g))) colnames(g1) <- names(g)
		g <- g1
		rm(g1)
	}
	if (is.matrix(g) && nrow(g)!=1) 
		stop("g should be a vector or an one-row matrix");
	if (!is.null(q)) {
		if (length(q) != dim(g)[2]) 
			stop("length of q should be the same as length of g");
		if (!is(q,"numeric")) 
			warning("q should be real (scalar or vector) between 0 and 1")
		if (any(q<0) | any(q>1)) 
			stop("(some) q<0 or q>1")
	}
# set P(g1|g2): probability that true genotype is g1 
# if g2 is observed and error rate is epsilon
# as a function of IBS (0,1,2)
	delta <- c((1-epsilon)^2,epsilon*(1-epsilon),epsilon^2)
# function to convert calls to probabilities assuming 
# error epsilon
	getNondegenerate <- function(glocal,delta) {
		if (!is.na(glocal)) {
			trma <- matrix(c(0,1,2,1,0,1,1,2,1,2,1,0),ncol=4)+1
			out <- delta[trma[glocal+1,]]
			out[2] <- out[2] + out[3]
			out <- out[-3]
			return(out)
		} else {
			return(c(1/4,1/2,1/4))
		}
	}
# setting missing to HWE
	setToHWE <- function(qlocal) {
		plocal <- 1 - qlocal
		out <- matrix(c(plocal*plocal,2*plocal*qlocal,qlocal*qlocal),ncol=3)
		out <- t(out)
		out
	}
# blur genotypes
	gx3 <- apply(g,MARGIN=2,FUN=getNondegenerate,delta=delta)
# apply HWE to missing
	if (!is.null(q)) 
		gx3[,is.na(g)] <- setToHWE(q[is.na(g)])
# done
	if (!is.null(colnames(g))) colnames(gx3) <- colnames(g)
	rownames(gx3) <- c("g0","g1","g2")
	return(gx3)
}

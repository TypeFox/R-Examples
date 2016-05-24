#' csSamWrapper function - performs entire functionality
#' 
#' csSamWrapper function - performs entire functionality
#' 
#' @param G Matrix of gene expression, columns ordered in the same order at the
#' cell-frequency matrix (n by g, n samples, g genes)
#' @param cc Matrix of cell-frequency. (n by k, n samples, k cell-types)
#' @param y A numeric vector of group association of each sample. Either 1 or 2.
#' @param nperms The number of permutations to perform.
#' @param alternative two.sided less greater
#' @param standardize Standardize sample or not. Default is TRUE.
#' @param medianCenter Median center rhat distributions. Default is TRUE.
#' @param logRm Exponentiate data for deconvolution stage. Default is FALSE
#' @param logBase Base of logaritm used to determine exponentiation factor.
#' Default is 2
#' @param nonNeg For single channel arrays. Set any cell-specific expression
#' estimated as negative, to a ceiling of 0. It is conservative in its study of
#' differential expression. Default is FALSE.
#' @param fileName PDF file containing plots of FDR vs. number of genes called
#' for whole tissue comparison (via SAM) as well as each cell-type (by csSAM)
#' @return Returns a list containing:
#' \item{deconv}{A list object containing a fit (cell-type specfic
#' expression) for each group. Each element in the list is an object returned by
#' csFit.}
#' \item{fdr.csSAM}{A list output of the fdrCsSAM function.}
#' \item{fdr.SAM}{A list output of the fdrSAM function.}
#' \item{sigGene.csSAM}{A list of significant genes.}
#' \item{fileName}{The filename into whcih the FDR plots are dumped.}
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @seealso
#' \code{\link{csfit}},\code{\link{csSAM}},\code{\link{fdrCsSAM}},\code{\link{plotCsSAM}}
#' @cite Shen-Orr2010
csSamWrapper <- function(G,cc,y,nperms = 200,alternative = 'two.sided',standardize=TRUE,medianCenter=TRUE, logRm =FALSE,logBase = 2,nonNeg=TRUE,fileName='csSAMout.pdf') {

	deconv <- list()

	numset = length(unique(y))
	numgene = ncol(G)
	numcell = ncol(cc)
	cellNames = colnames(cc)
	
	n <- vector(mode = "logical", length  = numset)
	for (i in 1:numset) n[i] =sum(y==i)

	for (curset in 1:numset) {
		deconv[[curset]]= csfit(cc[y==curset,], G[y==curset,],logRm,logBase)
	}
	rhat <- array(dim = c(numcell,numgene))
	rhat[, ] <- csSAM(deconv[[1]]$ghat, deconv[[1]]$se,
				  n[1], deconv[[2]]$ghat, deconv[[2]]$se, n[2],
				  standardize, medianCenter, nonNeg)  
	tt.sam <- runSAM(G, y)
	fdr.csSAM <- fdrCsSAM(G,cc,y,n,numcell,numgene, rhat,
					nperms = nperms,alternative,standardize,
					medianCenter, logRm,logBase,nonNeg)
	fdr.sam <- fdrSAM(G, y, nperms=nperms, tt.sam, alternative)
	sigGene <- findSigGene(G, cc, y, rhat, fdr.csSAM)

	plotCsSAM(fdr.csSAM, fdr.sam,alternative,cellNames,numcell, fileName = fileName)
	return(list(deconv = deconv, fdr.csSAM = fdr.csSAM, fdr.SAM = fdr.sam,sigGene.csSAM = sigGene,fileName = fileName))
}

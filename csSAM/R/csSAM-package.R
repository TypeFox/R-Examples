

#' Internal csSAM functions
#' 
#' These functions are not to be called by the user.
#' 
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' @keywords internal
#' @name csSAM-internal
NULL





#' Cell-specific Differential Expression (csSAM)
#' 
#' SAM for Cell-specific Differential Expression SAM.
#' 
#' \tabular{ll}{ Package: \tab csSAM\cr Type: \tab Package\cr Version: \tab
#' 1.2\cr Date: \tab 2011-10-08\cr License: \tab LGPL\cr LazyLoad: \tab yes\cr
#' 
#' Tissues are often made up of multiple cell-types. Each with its own
#' functional attributes and molecular signature.  Yet, the proportions of any
#' given cell-type in a sample can vary markedly. This results in a significant
#' loss of sensitivity in gene expression studies and great difficulty in
#' identifying the cellular source of any perturbations. Here we present a
#' statistical methodology (cell-type specific Significance Analysis of
#' Microarrays or csSAM) which, given microarray data from two groups of
#' biological samples and the relative cell-type frequencies of each sample,
#' estimates in a virtual manner the gene expression data for each cell-type at
#' a group level, and uses these to identify differentially expressed genes at a
#' cell-type specific level between groups.
#' 
#' The lower limit for the number of samples needed for deconvolving the
#' cell-specific expression of N cell-types is N+1. For a singe color array -
#' the result could be interperted as the avg. expression level of a given gene
#' in a cell-type of that group. Multiplied by the frequecy of a given cell-type
#' in an individual in the group, it is the amount contributed by that cell type
#' to the overall measured expression on the array.  \cr Key functions for this
#' package:\cr csSamWrapper - Single wrapper function performs all
#' functionality. csfit: For deconvolving the average cell-type specific
#' expression for each cell-type in a given group.\cr csSAM: For calculating the
#' constrast between every pair of cells being compared between the two
#' groups.\cr fdrCsSAM: Estimate the false discovery rate for each cell-type
#' specific comparison.\cr findSigGenes:Identifies the list of differentially
#' expressed genes in a given cell-type at a given FDR cutoff.\cr
#' plotCsSAM:Plots a fdr plot of ther results.\cr } Additional functions exists
#' (runSAM and fdrSAM to contrast csSAM with the tissue heterogeneity ignorant
#' SAM).
#' 
#' @name csSAM-package
#' @docType package
#' @author Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' 
#' Maintainer: Shai Shen-Orr <shenorr@@stanford.edu>
#' @bibliography ~/Documents/articles/library.bib
#' @cite Shen-Orr2010
#' @exportPattern ^[^\\.]
#' @examples
#' 
#' library("csSAM")
#' ##
#' ## Generate random dataset
#' ##
#' set.seed(143)
#' k <- 5 # number of cell types 
#' ng <- 500 # number of genes
#' p <- 20 # number of samples
#' ndiff <- 100 # number of genes differentially expressed
#' 
#' # true cell-specific signatures
#' H1 <- matrix(rnorm(5*ng), ncol=ng)
#' H2 <- H1
#' # create differential expression for 3rd cell type
#' H2[3,1:ndiff] <- H2[3,1:ndiff] + 5 
#' 
#' # cell frequency matrix per sample
#' cc <- matrix(runif(p*k), ncol=k) 
#' cc <- t(scale(t(cc), center=FALSE, scale=rowSums(cc))) 
#' colnames(cc) <- paste('cellType', 1:ncol(cc), sep="")
#' 
#' # global expression matrix
#' G <- rbind(cc[1:10, ] %*% H1, cc[11:p, ] %*%H2 ) + matrix(rnorm(p*ng), ncol=ng)
#' # sample classes (2 groups) 
#' y <- gl(2, p/2)  
#' 
#' fileName = "Example File.pdf";
#' \dontshow{ on.exit(unlink(filename)) }
#' 
#' # Now run, either using the wrapper 
#' # NB: more permutations would be needed for real data
#' deconvResults = csSamWrapper(G, cc, y, nperms = 50, alternative = "two.sided"
#' 								, standardize = TRUE
#' 								, medianCenter = TRUE
#' 								, fileName = fileName)
#' 
#' # Or by calling each function independently: 
#' # this is useful if you want to perform only cell-specific expression 
#' # without differential expression.
#' \dontrun{
#' numset = nlevels(y)
#' n <- summary(y, maxsum=Inf) # number of samples in each class
#' numgene = ncol(G)
#' numcell = ncol(cc)
#' geneID = colnames(G)
#' cellID = colnames(cc)
#' 	
#' deconv <- list()
#' # run analysis
#' for (curset in levels(y))
#' 	deconv[[curset]]= csfit(cc[y==curset,], G[y==curset,])
#' 
#' rhat <- array(dim = c(numcell,numgene))
#' rhat[, ] <- csSAM(deconv[[1]]$ghat, deconv[[1]]$se,
#'                   n[1], deconv[[2]]$ghat, deconv[[2]]$se, n[2],
#'                   standardize=TRUE, medianCenter=TRUE, nonNeg=TRUE)  
#' 
#' tt.sam <- runSAM(G, y)
#' falseDiscovR <- fdrCsSAM(G,cc,y,n,numcell,numgene, rhat,
#'                     nperms = 200,standardize=TRUE,alternative='two.sided',
#'                     medianCenter=TRUE, nonNeg=TRUE)
#' falseDiscovRSAM <- fdrSAM(G, y, nperms=200, alternative = 'two.sided',tt.sam)
#' sigGene <- findSigGene(G, cc, y, rhat, falseDiscovR)
#' 
#' plotCsSAM(falseDiscovR, falseDiscovRSAM,alternative='two.sided',cellID,numcell, fileName)
#' print (falseDiscovR$fdr.g[ , ] )
#' }
#' 
NULL




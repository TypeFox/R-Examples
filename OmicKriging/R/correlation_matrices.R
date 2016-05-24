#' Compute gene expression correlation matrix.
#'
#' This function computes a gene expression correlation matrix given a file of
#' transcript expression levels for each sample in the study. It returns a
#' correlation matrix with rownames and colnames as sample IDs.
#'
#' @param expFile Path to gene expression file.
#' @param gxmFilePrefix File path prefixes for outputting GCTA style binary
#'   correlation matrices.
#' @param idfile Path to file containing family IDs and sample IDs with header FID and IID. 
#'
#' @return Returns a correlation matrix of (N-samples x N-samples), with
#'   rownames and colnames as sample IDs.
#'
#' @importFrom utils read.delim
#' @importFrom utils read.table
#' @importFrom stats cor 
#' @export
#' @examples
#'  ## load gene expression values from vignette
#'  expressionFile <- system.file(package = "OmicKriging",
#'                      "doc/vignette_data/ig_gene_subset.txt.gz")
#'  ## compute correlation matrix
#'  geneCorrelationMatrix <- make_GXM(expressionFile)
make_GXM <- function(expFile = NULL, gxmFilePrefix = NULL, idfile = NULL) {

  ## data input
  genedata <- read.delim(expFile, sep=" ", as.is=T, header=T)
  if(!is.null(idfile)) {
    iddata <- read.table(idfile,header=T,as.is=T)
    genedata <- merge(iddata,genedata,by.x=c("FID","IID"),by.y=c("FID","IID"))
  }
  cordata.id <- genedata[c("FID","IID")]
  genemat <- as.matrix(genedata[,!(names(genedata) %in% c("FID","IID"))])
  
  ## center and scale genemat
  genemat <- scale(genemat, center = TRUE, scale = TRUE)
  
  ## impute mean for missing values
  genemat[is.na(genemat)] <- 0.0
   
  ## compute cor mat
  cormat <- cor(t(genemat))
  
  ## give it row names
  colnames(cormat) <- cordata.id[,2]
  rownames(cormat) <- cordata.id[,2]

  ## write to disk
  if(!is.null(gxmFilePrefix)) {
    write_GRMBin(X = cormat, prefix = gxmFilePrefix, n.snps = length(genemat[1,])) 
  }
  return(cormat)
}

#' Run Principal Component Analysis (PCA) using base R svd() function.
#'
#' A simple wrapper around the base R svd() function which returns the top N
#' eigenvectors of a matrix. Use this function to generate covariates for use
#' with the \code{\link{okriging}} or \code{\link{krigr_cross_validation}}
#' functions. This wrapper preserves the rownames of the original matrix.
#'
#' @param X A correlation matrix.
#' @param n.top Number of top principal compenents to return 
#'
#' @return A matrix of Principal Components of dimension (# of samples) x
#'   (n.top). As expected, eigenvectors are ordered by eigenvalue. Rownames
#'   are given as sample IDs.
#'
#' @keywords covariate, PCA, GRM
#' @export
#' @examples
#'  ## compute PC's using the  gene expression correlation matrix from vignette
#'  ## load gene expression values from vignette
#'  expressionFile <- system.file(package = "OmicKriging",
#'                      "doc/vignette_data/ig_gene_subset.txt.gz")
#'  ## compute correlation matrix
#'  geneCorrelationMatrix <- make_GXM(expressionFile)
#'  ## find top ten PC's of this matrix using SVD
#'  topPcs <- make_PCs_svd(geneCorrelationMatrix, n.top=10) 
make_PCs_svd <- function(X, n.top = 2) {
  res <- La.svd(X, nu = n.top)
  rownames(res["u"][[1]]) <- rownames(X)
  return(res["u"][[1]])
}

#' Run Principal Component Analysis (PCA) using the irlba package.
#'
#' A simple wrapper around the irlba() function which computes a partial SVD
#' efficiently. This function's run time depends on the number of eigenvectors
#' requested but scales well. Use this function to generate covariates for use
#' with the \code{\link{okriging}} or \code{\link{krigr_cross_validation}}
#' functions.
#'
#' @param X A correlation matrix.
#' @param n.top Number of top principal compenents to return 
#'
#' @return A matrix of Principal Components of dimension (# of samples) x
#'   (n.top). As expected, eigenvectors are ordered by eigenvalue. Rownames
#'   are given as sample IDs.
#'
#' @keywords covariate, PCA, GRM
#'
#' @references library(irlba)
#'
#' @import irlba
#' @export
#' @examples
#'  ## compute PC's using the  gene expression correlation matrix from vignette
#'  ## load gene expression values from vignette
#'  expressionFile <- system.file(package = "OmicKriging",
#'                      "doc/vignette_data/ig_gene_subset.txt.gz")
#'  ## compute correlation matrix
#'  geneCorrelationMatrix <- make_GXM(expressionFile)
#'  ## find top ten PC's of this matrix using SVD
#'  topPcs <- make_PCs_irlba(geneCorrelationMatrix, n.top=10) 
make_PCs_irlba <- function(X, n.top = 2) {
  res <- irlba(X, nu = n.top)
  rownames(res$u) <- rownames(X)
  return(res$u)
}


#' Read the GRM binary file.
#'
#' Function provided by GCTA maintainers (modified slightly) for accessing their
#' recently introduced binary GRM format. The GRM is stored as a vector of numerics
#' which correspond to the lower triangular elements including the diagonal. We simply read these, pull
#' the diagonal elements, and inflate them into a full symmetric matrix. We add
#' sample IDs to colnames and rownames for compatibility with other Kriging 
#' functions.
#'
#' Note that the GRM is described by three files, and this function assumes that all
#' have a common prefix that is passed in.
#'
#' @param prefix The file path prefix to GRM binary files (e.g., test.grm.bin, test.grm.N.bin, test.grm.id.)
#' @param size The length (in bytes) of each value in the raw GRM vector. Default is 4, and matches GRM writen by GCTA 1.11.
#'
#' @return GRM of dim (N.samples x N.samples) with rownames and colnames as sample ID.
#'
#' @references http://www.complextraitgenomics.com/software/gcta/estimate_grm.html
#'
#' @importFrom utils read.table
#' @export
#' @examples
#'   ## read binary Genetic Relatedness Matrix (GRM) generated by GCTA
#'   grmFile <- system.file(package = "OmicKriging",
#'                          "doc/vignette_data/ig_genotypes.grm.bin")
#'   grmFileBase <- substr(grmFile,1, nchar(grmFile) - 4)
#'   GRM <- read_GRMBin(grmFileBase)
read_GRMBin <- function(prefix, size = 4){
  sum_i <- function(i){
    return(sum(1:i))
  }

  ## open file connections and read in data
  BinFileName <- paste(prefix,".bin",sep="")
  NFileName <- paste(prefix,".N.bin",sep="")
  IDFileName <- paste(prefix,".id",sep="")
  id <- read.table(IDFileName)
  n <- dim(id)[1]
  BinFile <- file(BinFileName, "rb")
  grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile <- file(NFileName, "rb")

  ## read in the number of SNPs used to calculate the GRM (does not appear to work)
  N <- readBin(NFile, n=1, what=numeric(0), size=size)
  i <- sapply(1:n, sum_i)
  
  ## clean up file connections
  close(BinFile)
  close(NFile)

  ## pull apart diagonal and lower triagular elements
  diag.elem <- grm[i]
  off.diag.elem <- grm[-i]

  ## create the full symmetric correlation matrix
  X <- diag(diag.elem)
  X[ upper.tri(X, diag = FALSE) ] <- off.diag.elem
  X <- X + t(X) - diag(diag(X)) 

  ## add sample IDs to rownames and colnames
  rownames(X) <- id$V2
  colnames(X) <- id$V2

  return(X)
}

#' Write GRM binary files.
#'
#' Function to write a binary GRM format recently introduced by GCTA. It takes
#' a correlation matrix as used by other Kriging functions, and writes three
#' files: binary file for storing the diagonal + lower triangular elements, a
#' text file for sample IDs, and a binary file storing the number of SNPs used
#' in the correlation matrix calculation.
#'
#' @param X Correlation matrix with rownames and colnames as sample IDs.
#' @param prefix Base file path and names for the three output files.
#' @param n.snps Number of SNPs used in correlation matrix calculation. Default is 0.0.
#' @param size Number of bytes to write for each value. Default is 4
#'
#' @return None. Though side effects are writing three files as described above.
#'
#' @references http://www.complextraitgenomics.com/software/gcta/estimate_grm.html
#'
#' @importFrom utils write.table
#' @export
#' @examples
#'   
#'   ## create a random genotype matrix
#'   nSamples <- 10
#'   mMarkers <- 100
#'   X <- matrix(rbinom(n=100, size=2, prob=0.5), nrow=nSamples)
#'   ## compute the Genetric Relatedness Matrix
#'   grm <- cor(X)
#'   ## write a Genetic Relatedness Matrix (GRM)
#'   ## NOTE: to following is not run here -- not writing any files in examples
#'   #write_GRMBin(grm, n.snps=mMarkers, prefix="grm.out")
write_GRMBin <- function(X, n.snps = 0.0, prefix, size = 4) {

  sum_i <- function(i){
    return(sum(1:i))
  }

  ## file connections
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  BinFileName <- paste(prefix,".grm.bin",sep="")

  ## pull sample ids and dimension of GRM
  id <- rownames(X)
  n <- length(id)

  ## pull diagonal elements
  diag.elem <- diag(X)

  ## pull lower triangular elements
  off.diag.elem <- X[lower.tri(X, diag=FALSE)]

  ## collapse GRM into vector
  i <- sapply(1:n, sum_i)
  collapsed.grm <- vector(mode="numeric", length = n*(n+1)/2)
  collapsed.grm[i] <- diag.elem
  collapsed.grm[-i] <- off.diag.elem

  ## write binary files
  BinFile <- file(BinFileName, "wb")
  NFile <- file(NFileName, "wb")
  writeBin(con = BinFile, collapsed.grm, size = size)
  writeBin(con = NFile, n.snps, size = size )
  close(BinFile)
  close(NFile)

  ## write sample ID file -- we are dropping sample family IDs here
  write.table(cbind(id, id), file = IDFileName)

}


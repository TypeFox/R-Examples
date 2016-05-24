#' Fragments per kilobase per million
#' 
#' Takes a count matrix and a vector of gene lengths and returns an optionally \code{log2}-transformed FPKM matrix. Modified from edgeR.
#' 
#' @author Davis McCarthy
#' @author Gordon Smyth
#'   
#' @param x a matrix of counts
#' @param length a vector of length \code{nrow(x)} giving length in bases
#' @param log logical, if \code{TRUE}, then \code{log2} values are returned.
#' @param prior.count average count to be added to each observation to avoid 
#'   taking log of zero. Used only if \code{log=TRUE}.
#'   
#' @return A matrix of FPKM values.
#' 
#' @examples
#' \dontrun{
#' library(readr)
#' library(dplyr)
#' countdata <- read_csv("http://files.figshare.com/2439061/GSE37704_featurecounts.csv") 
#' counts <- countdata %>% select(countdata, starts_with("SRR")) %>% as.matrix
#' counts2fpkm(counts, countdata$length)
#' }
#'   
#' @export
counts2fpkm <- function(x, length, log=FALSE, prior.count=.25) {
    # sanity checks
    if (class(x)!="matrix") stop("x must be a matrix")
    if (nrow(x)!=length(length)) stop("dimensions of count matrix and gene lengths don't match")
    # library size is sum of reads in each sample
    lib.size <- colSums(x)
    if (log) {
        # If you're log scaling, you'll have to add something to the zeros, so
        # adjust the library sizes accordingly. prior.count is the average count to
        # be added to each observation to avoid taking log of zero.
        prior.count.scaled <- lib.size/mean(lib.size) * prior.count
        lib.size <- lib.size + 2 * prior.count.scaled
    }
    # Per million
    lib.size <- 1e-06 * lib.size
    if (log) {
        cpm <- log2(t((t(x) + prior.count.scaled)/lib.size))
    } else {
        cpm <- t(t(x)/lib.size)
    }
    # per kilobase
    length.kb <- length/1000
    if (log) {
        fpkm <- cpm-log2(length.kb)
    } else {
        fpkm <- cpm/length.kb
    }
    return(fpkm)
}

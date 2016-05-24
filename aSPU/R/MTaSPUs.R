#' The SPU and aSPU tests for multiple traits - single SNP association with GWAS summary statistics.
#'
#' SNP based adaptive association test for multiple phenotypes with GWAS summary statistics.
#'
#' @param Z matrix of summary Z-scores, SNPs in rows and traits in columns. Or a vector of summary Z-scores for a single snp
#'
#' @param v estimated correlation matrix based on the summary Z-scores (output of estcov)
#'
#' @param pow power used in SPU test. A vector of the powers.
#'
#' @param transform if TRUE, the inference is made on transformed Z
#' 
#' @param B number of Monte Carlo samples simulated to compute p-values, the maximum number of MC simulations is 1e8
#'
#' @return compute p-values for SPU(gamma) i.e. pow=1:8, and infinity
#'             aSPU, based on the minimum p-values over SPU(power)
#'             each row for single SNP
#'

#' @author Junghi Kim, Yun Bai and Wei Pan 
#'
#' @references
#' Junghi Kim, Yun Bai and Wei Pan (2015) An Adaptive Association Test for Multiple Phenotypes with GWAS Summary Statistics, Genetic Epidemiology, 8:651-663
#'
#' @examples
#'
#' # -- n.snp: number of SNPs
#' # -- n.trait: number of traits
#' # -- n.subject: number of subjects
#'
#' n.snp <- 100
#' n.traits <- 10
#' n.subjects <- 1000
#' traits <- matrix(rnorm(n.subjects*n.traits), n.subjects, n.traits)
#' v <- cov(traits)
#' allZ <- rmvnorm(n.snp, sigma=v)
#' colnames(allZ) <- paste("trait", 1:n.traits, sep="")
#' rownames(allZ) <- paste("snp", 1:n.snp, sep="")
#'
#' 
#' r <- estcov(allZ)
#' MTaSPUs(Z = allZ, v = r, B = 100, pow = c(1:4, Inf), transform = FALSE)
#' MTaSPUs(Z = allZ[1,], v = r, B = 100, pow = c(1:4, Inf), transform = FALSE)
#' minP(Zi= allZ[1,], r = r)
#'
#' @seealso \code{\link{minP}} \code{\link{estcov}}

MTaSPUs <- function(Z, v, B, pow, transform = FALSE){
   # -- Z: matrix of summary Z-scores, SNPs in rows and traits in columns  
   # -- Or a vector of summary Z-scores for a single snp
   # -- v: output of estcov
   # --    estimated estimated correlation matrix  based on the summary Z-scores
   # -- transform: if TRUE, the inference is made on transformed Z
   # -- B: number of Monte Carlo samples simulated to compute p-values max(B)=1e8
   # -- results: compute p-values for SPU(gamma) i.e. pow=1:8, and infinity
   # --          aSPU, based on the minimum p-values over SPU(power)
   # --          each row for single SNP    
	if (B < 1e8) p <- MTaSPUsmallB(Z, v, B, pow, transform) else  p <- MTaSPUsB1e8(Z, v, pow, transform)
	return(p)
}



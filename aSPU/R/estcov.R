#' estcov
#'
#' Estimate the covariance matrix of multiple traits based on their (null) summary Z-scores.
#'
#' @param allZ matrix of summary Z-scores for all SNP. each row for SNP; each column for single trait.
#'
#' @return estimated correlation matrix.
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
#' @seealso \code{\link{MTaSPUs}} \code{\link{minP}}

estcov <- function(allZ){
    n.snp <- dim(allZ)[1]
    n.traits <- dim(allZ)[2]
    minZ = rowMins(2*pnorm(-abs(allZ)))
    nullSNPs = which(minZ>(0.05/(n.snp/n.traits)))
    return( cor(allZ[nullSNPs,]))
}

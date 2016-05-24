#' minP test.
#'
#' Return exact minP test p-value for multiple traits - single SNP association. 
#'
#' @param Zi a vector of summary Z-scores for single SNP
#'
#' @param r estimated correlation matrix based on the summary Z-scores (output of estcov)
#'
#' @return return exact minP test
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
#' @seealso \code{\link{estcov}} \code{\link{MTaSPUs}}

minP <- function(Zi, r){
    n <- dim(r)[1]
    x <- as.numeric(max(abs(Zi)))
    return(as.numeric(1 - pmvnorm(lower=c(rep(-x,n)), upper=c(rep(x,n)), mean=c(rep(0,n)), sigma=r)))
}


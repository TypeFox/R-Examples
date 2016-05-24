#' @title Fu's Fs
#' @description Calculate Fu's Fs for a set of sequences to test 
#'   for selection.
#' 
#' @param x set of DNA sequences or a haploid \linkS4class{gtypes} 
#'   object with sequences.
#' 
#' @references Fu, Y-X. 1997. Statistical tests of neutrality of mutations 
#'   against population growth, hitchiking and background selection. 
#'   Genetics 147:915-925.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @importFrom copula Stirling1
#' @export
#' 
fusFs <- function(x) {
  dna <- if(inherits(x, "gtypes")) {
    sequences(x, as.haplotypes = FALSE)
  } else {
    as.multidna(x)
  }
  
  haps <- lapply(getSequences(dna, simplify = FALSE), labelHaplotypes)
  
  sapply(haps, function(h) {
    if(is.null(h)) return(NA)
    h$hap.seqs <- as.matrix(h$hap.seqs)
    h$haps <- na.omit(h$haps)
    n <- length(h$haps)
    k0 <- length(unique(h$haps))
    
    pws.diff <- dist.dna(h$hap.seqs, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE)
    pws.diff <- pws.diff[h$haps, h$haps]
    theta.pi <- mean(pws.diff[lower.tri(pws.diff)])
    Sn.theta.pi <- prod(theta.pi + 0:(n - 1)) # potential typo on page 916 that implies prod(theta.pi - 0:(n + 1))
    Sk.theta.k <- sapply(k0:n, function(k) abs(Stirling1(n, k)) * (theta.pi ^ k))
    s.prime <- sum(Sk.theta.k / Sn.theta.pi)
    
    log(s.prime / (1 - s.prime))
  })
}

#' @rdname popStructStat
#' @export
#' 
Hstats <- function(g) {
  if(ploidy(g) < 2) {
    result <- matrix(NA, nrow = 3, ncol = nLoc(g))
    rownames(result) <- c("Ho", "Hs", "Ht")
    colnames(result) <- locNames(g)
    return(result)
  }
  
  result <- Hstats_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    as.numeric(strata(g)) - 1,
    ploidy(g)
  )
  rownames(result) <- c("Ho", "Hs", "Ht")
  colnames(result) <- locNames(g)
  
#   nloc <- ncol(g@loci)
#   result <- matrix(0, nrow = 3, ncol = nloc)
#   rownames(result) <- c("Ho", "Hs", "Ht")
#   colnames(result) <- colnames(g@loci)
#   for(i in 1:nloc) {
#     ## Estimate Ho (frequency of all heterozygotes): Equation 5, page 254
#     locus <- g@loci[, i]
#     loc.mat <- matrix(locus, ncol = g@ploidy)
#     genotype <- sapply(1:nrow(loc.mat), function(i) {
#       alleles <- loc.mat[i, ]
#       if(any(is.na(alleles))) return(NA)
#       if(length(unique(alleles)) == 1) alleles[1] else "het"
#     })
#     hom.freq <- prop.table(table(g@strata, genotype), 1)
#     hom.freq <- hom.freq[, colnames(hom.freq) != "het", drop = FALSE]
#     Ho <- if(ncol(hom.freq) == 0) 0 else {
#       1 - sum(colMeans(hom.freq))
#     }
# 
#     ## Estimate Hs (expected heterozygosity within strata): Equation 9, page 255
#     allele.freq <- table(strata, locus)
#     strata.freq <- rowSums(allele.freq) / g@ploidy
#     hom.freq <- prop.table(allele.freq, 1)
#     mean.het <- mean(1 - rowSums(hom.freq ^ 2))
#     harm.n <- harmonic.mean(strata.freq)
#     Hs <- (harm.n / (harm.n - 1)) * (mean.het - (Ho / 2 / harm.n))
# 
#     ## Estimate Ht (expected heterozygosity overall): Equation 11, page 256
#     mean.het <- 1 - sum(colMeans(hom.freq) ^ 2)
#     harm.n.s <- harm.n * sum(allele.freq)
#     Ht <- mean.het + (Hs / harm.n.s) - (Ho / 2 / harm.n.s)
#     result[, i] <- c(Ho, Hs, Ht)
#   }
  
  result
}
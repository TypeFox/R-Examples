#' @rdname popStructStat
#' @importFrom swfscMisc harmonic.mean 
#' @export
#' 
statJostD <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {
  if(is.null(strata.mat)) g <- g[, , strataNames(g)]
  if(ploidy(g) < 2 | nStrata(g) == 1) {
    return(list(
      stat.name = "D", 
      result = c(estimate = NA, p.val = NA),
      null.dist = NULL
    ))
  }
  
  strata.mat <- .checkStrataMat(strata.mat, g, nrep)
  
  result <- statJostD_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    strata.mat, ploidy(g)
  )
  
  if(length(result) == 1) keep.null <- FALSE
  p.val <- if(length(result) == 1) NA else mean(result >= result[1], na.rm = TRUE) 
  return(list(
    stat.name = "D", 
    result = c(estimate = result[1], p.val = p.val),
    null.dist = if(keep.null) result[-1] else NULL
  ))
  
#   allele.freqs <- alleleFreqs(g, by.strata = TRUE)
#   terms <- sapply(locNames(g), function(x) {
#     num.strata <- dim(allele.freqs[[x]])[3]
#     i.terms <- sapply(dimnames(allele.freqs[[x]])[[1]], function(a) {
#       j.terms <- sapply(dimnames(allele.freqs[[x]])[[3]], function(s) {
#         Nj <- sum(allele.freqs[[x]][, "freq", s])
#         Nij <- allele.freqs[[x]][a, "freq", s]
#         a.term1 <- Nij / Nj
#         a.term2 <- a.term1 ^ 2
#         b.term <- Nij * (Nij - 1) / (Nj * (Nj - 1))
#         c(a.term1 = a.term1, a.term2 = a.term2, b.term = b.term)
#       })
#       a.term1 <- sum(j.terms["a.term1", ]) ^ 2
#       a.term2 <- sum(j.terms["a.term2", ])
#       c(a = (a.term1 - a.term2) / (num.strata - 1), 
#         b = sum(j.terms["b.term", ])
#       )
#     })
#     
#     c(a = sum(i.terms["a", ]), b = sum(i.terms["b", ]))
#   })
#   d.by.locus <- 1 - terms["a", ] / terms["b", ]
#   d.by.locus <- ifelse(d.by.locus < 0, 0, d.by.locus)
#   est <- harmonic.mean(d.by.locus)
# 
#   c(D = est)
}
#' @rdname popStructStat
#' @export
#' 
statChi2 <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {
  if(is.null(strata.mat)) g <- g[, , strataNames(g)]
  strata.mat <- .checkStrataMat(strata.mat, g, nrep)
  
  result <- statChi2_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    strata.mat, ploidy(g)
  )
  
  if(length(result) == 1) keep.null <- FALSE
  p.val <- if(length(result) == 1) NA else mean(result >= result[1], na.rm = TRUE) 
  return(list(
    stat.name = "Chi2", 
    result = c(estimate = result[1], p.val = p.val),
    null.dist = if(keep.null) result[-1] else NULL
  ))
  
#   if(nStrata(g) == 1) return(Chi2 = NA)
#   strata <- if(is.null(strata)) {
#     strata(g)
#   } else {
#     rep(strata, length.out = nInd(g))
#   }
#   if(!is.factor(strata)) strata <- factor(strata)
#   
#   if(any(is.na(strata))) {
#     toUse <- !is.na(strata)
#     strata <- strata[toUse]
#     g <- g[toUse, , ]
#   }
#   
#   est <- statChi2_C(
#     sapply(loci(g), function(x) as.numeric(x) - 1), 
#     as.numeric(strata) - 1,
#     ploidy(g)
#   )
  
#   est <- vector("numeric", ncol(g@loci))
#   for(i in 1:length(est)) {
#     obs.freq <- table(g@loci[, i], strata, useNA = "no")
#     if(nrow(obs.freq) > 0 & ncol(obs.freq) > 1) {
#       exp.freq <- outer(rowSums(obs.freq), colSums(obs.freq)) / sum(obs.freq)
#       est[i] <- sum((obs.freq - exp.freq) ^ 2 / exp.freq, na.rm = TRUE)
#     }
#   }
  
#  c(Chi2 = est)
}
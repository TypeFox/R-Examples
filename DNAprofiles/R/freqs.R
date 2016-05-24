#' Recode allelic frequencies with different levels
#'
#' @param freqs list with named numeric vectors \code{x} and \code{fx}, denoting respectively the events and probabilities of the discrete distribution.
#' @param along.with second list of allelic frequencies from which the levels are taken
#' @details Profiles are stored with integers corresponding to the corresponding index of the names attribute of the allelic frequencies. This funciton recodes a set of frequencies to include all names of a larger set of allelic frequencies.
#' @return list with named numeric vectors \code{x} and \code{fx}, denoting respectively the events and probabilities of the discrete distribution.
#' @export
#' @examples
#'# The following example demonstrates the effect on LRs of uncertainty from a small reference db
#'# We sample a small reference db and compute LRs with estimated frequencies
#'
#'data(freqsNLsgmplus)
#'
#'fr.true <- freqsNLsgmplus
#'
#'# sample a small db
#'n <- 1e2
#'x <- sample.profiles(n,fr.true)
#'x.char <- profiles.to.chars(x,two.cols.per.locus = TRUE)
#'
#'# estimate frequencies
#'fr.hat <- list()
#'for (L in names(fr.true)){
#'  f0.tab <- table(x.char[,paste(L,c(".1",".2"),sep = "")])
#'  fr.hat[[L]] <- setNames(as.vector(f0.tab)/(2*n),names(f0.tab))  
#'}
#'
#'# not all alleles are seen in the sample, so the ladders don't line up
#'# rmp(x,fr.hat) # reallly wrong!
#'fr.hat <- recode.freqs(freqs = fr.hat,along.with = fr.true) # fix
#'
#'plot(log10(rmp(x)),log10(rmp(x,freqs = fr.hat)))
#'abline(a=0,b=1)
#' @export
recode.freqs <- function(freqs,along.with){
  f1 <- freqs[names(along.with)]
  f2 <- along.with
  
  f1.lev <- lapply(names(f1), function(L) names(f1[[L]]))
  names(f1.lev) <- names(f1)
  f2.lev <- lapply(names(f1), function(L) names(f2[[L]]))
  names(f2.lev) <- names(f1)
  
  if (!all(sapply(names(f1),function(L) all(names(f1[[L]]) %in% names(f2[[L]]))))){
    stop("Not all allele names of freqs are found in along.with")
  }
  
  ret <- list()
  for (L in names(f1)){
    f <- numeric(length(f2.lev[[L]]))
    f[match(f1.lev[[L]],table=f2.lev[[L]])] <- as.vector(f1[[L]])
    names(f) <- f2.lev[[L]]    
    ret[[L]] <- f
  }
  ret  
}
NULL
Zfreqs.to.mat <- function(freqs,markers=names(freqs)){
  freqs.L <- sapply(freqs[markers],length)
  ret <- matrix(0.,nrow = max(freqs.L),ncol = length(markers))
  for(m in seq_along(markers)) ret[seq_len(freqs.L[m]),m] <- as.vector(freqs[[markers[m]]])
  ret
}
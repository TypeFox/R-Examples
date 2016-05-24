#' A new function for calculating allelic richness by resampling
#' 
#' Kevin Keenan, 2015
#' 
# define a function for calculating allelic richness 
arSample <- function(dat, nrep, ci = FALSE, alpha = 0.05){
  ps <- sapply(dat$genos, function(x) dim(x)[1])
  n <- min(ps)
  # create resample indexes
  idx <- lapply(1:nrep, function(i){
    lapply(ps, function(x) sample(x, n, replace = TRUE))
  })
  
  Rcpp::sourceCpp("allCount.cpp")
  
  # wrapper function for calculating number of alleles
  allCountWrap <- function(genos, idx) {
    apply(genos[idx,,], 2, allCount)
  }
  
  # calculate the number of alleles per sub-sample
  nall <- lapply(idx, function(x){
    mapply(allCountWrap, genos = dat$genos, idx = x, SIMPLIFY = FALSE)
  })

  # organise allele counts per population sample
  nall <- lapply(1:length(ps), function(i) sapply(nall, "[[", i))
  
  # calculate allelic richness
  out <- list(ar = data.frame(sapply(nall, rowMeans, na.rm = TRUE)))
  out$ar <- cbind(dat$locs, out$ar)
  colnames(out$ar) <- c("Locus", gsub(",", "", sapply(dat$indnms, "[", 1)))
  # calculate CI if called for
  if(ci){
    cis <- lapply(nall, function(x){
      apply(x, 1, quantile, prob = c(alpha/2, 1-(alpha/2)))
    })
    out$lo <- do.call(cbind, lapply(cis, function(x) return(x[1,])))
    out$lo <- data.frame(dat$locs, out$lo)
    colnames(out$lo) <- c("Locus", gsub(",", "", sapply(dat$indnms, "[", 1)))
    out$hi <- do.call(cbind, lapply(cis, function(x) return(x[2,])))
    out$hi <- data.frame(dat$locs, out$hi)
    colnames(out$hi) <- c("Locus", gsub(",", "", sapply(dat$indnms, "[", 1)))
    out$mean_ci <- lapply(nall, function(x){
      quantile(colMeans(x), prob = c(alpha/2, 1-(alpha/2)))
    })
  }
  return(out)
}
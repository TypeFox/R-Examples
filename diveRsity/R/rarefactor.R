#' allelic richness calculations
#' 
#' Kevin Keenan, 2015

rarefactor <- function(dat){
  nl <- length(dat$af)
  count_convert <- function(af, ps){
    t(apply(af, 1, function(x){ x * (ps*2)}))
  }
  all_count <- mapply(count_convert, af = dat$af, ps = dat$ps, SIMPLIFY = FALSE)
  rf <- function(all_count, ps){
    n <- min(ps)*2
    apply(all_count, 2, function(x){
      ns <- sum(x)
      samp <- choose(ns - x, n)/choose(ns, n)
      samp[is.na(samp)] <- 0
      sum(1-samp)
    })
  }
  
  out <- data.frame(locus  = dat$locs)
  ar <- t(mapply(rf, all_count = all_count, ps = dat$ps, SIMPLIFY = TRUE))
  colnames(ar) <- gsub(",", "", sapply(dat$indnms, "[", 1))
  out <- cbind(out, ar)
  return(out)
}
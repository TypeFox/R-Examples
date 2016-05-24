specCons <- function(x, log) {
  
  if ( nrow(x) > 1 ){
    if ( !missing(log) ){
      acc <- do.call(rbind, lapply(rownames(x), splitGiTaxon))
      slog("\n ", acc[1, 2], "-", paste(acc[, "gi"], collapse = ", "), 
           file = log)
    }
    obsvalue <- as.raw(c(136, 40, 72, 24))
    tableRaw <- function(x, obsvalue){
      sapply(obsvalue, function(o, x) length(which(x %in% o)), x = x)
    }
    ## corresponds to method 'profile' in seqinr::consensus
    x <- apply(x, 2, tableRaw, obsvalue = obsvalue)
    return(apply(x, 2, function(x, o) o[which.max(x)], o = obsvalue))
  } else {
    return(as.vector(x))
  }
}
# prepSNPInfo
prepSNPInfo <- function(snpinfo, snpNames, aggregateBy, wt1=NULL, wt2=NULL) {
  
  if(!is.character(snpinfo[ , snpNames])) {
    snpinfo[ , snpNames] <- as.character(snpinfo[ , snpNames])
    warning("Converting snpNames column to character")
  }
  if(!is.character(snpinfo[ , aggregateBy])) {
    snpinfo[ , aggregateBy] <- as.character(snpinfo[ , aggregateBy])
    warning("Converting aggregateBy column to character")
  }
  
  dups <- duplicated(snpinfo[ , c(aggregateBy, snpNames)])
  snp_na <- is.na(snpinfo[ , snpNames])
  gene_na <- is.na(snpinfo[ , aggregateBy])
  
  idx <- which(!(dups | snp_na | gene_na))
  
  cols <- c(aggregateBy, snpNames)
  if (is.character(wt1) && length(wt1) > 0L) {
    cols <- c(cols, wt1)
  }
  if (is.character(wt2) && length(wt2) > 0L) {
    cols <- c(cols, wt2)
  }
  
  if (length(cols) > 2L) {
    si_unique <- unique(snpinfo[!(snp_na | gene_na) , cols])
    for (i in 3:length(cols)) {
      if (any(is.na(si_unique[ , i]))) {
        stop("Cannot have missing weights.  Please impute and re-run")
      }
    }           
    if (nrow(si_unique) != length(idx)) {
      stop("Non-unique weight for snp-gene pair")
    }    
  } 
  
  nmiss <- sum(snp_na | gene_na)
  if (nmiss > 0L) {
    warning("Removed ", nmiss, " missing snps and/or genes from SNPinfo file") 
  } 
  
  ndups <- sum(dups)
  if (ndups > 0L) {
    warning("Removed ", ndups," duplicate snp-gene pairs from SNPinfo file") 
  }   
  
  snpinfo[idx, cols]
}


# is_monomorphic
is_monomorphic <- function(x, tolerance = .Machine$double.eps ^ 0.5) {
 
  caf <- colMeans(x, na.rm=TRUE)
  caf[is.na(caf)] <- -1
  
  caf1 <- rep(FALSE, ncol(x))
  names(caf1) <- colnames(x)
  caf1_idx <- which(caf > 1-tolerance & caf < 1+tolerance)
  if (length(caf1_idx) > 0) {
    tmp <- apply(x[ , caf1_idx, drop=FALSE], 2, function(x) all(stats::na.omit(x) == 1))
    caf1[names(tmp)] <- tmp    
  }
  
  (caf == 0 | caf1 | caf == 2)
}

monomorphic_snps <- function(x, tolerance = .Machine$double.eps ^ 0.5) {
  monos <- is_monomorphic(x, tolerance)
  names(monos[monos])
}


is_mono <- function(Z) {
  apply(Z, 2, function(x) {
    y <- stats::na.omit(x)
    if (length(y) == 0) {
      FALSE
    } else {
      (all(y == 0) | all(y == 1) | all(y == 2))}
  })
}

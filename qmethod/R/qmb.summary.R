qmb.summary <- function(qmboots) {
  # Basic info from the analysis
  nfactors <- qmboots$orig.res$brief$nfactors
  nstat    <- qmboots$orig.res$brief$nstat
  nqsorts  <- qmboots$orig.res$brief$nqsorts
  
  #-------------------------------------------------------
  # Gather results of Q-sorts
  obj.loa <- as.array(paste0("qmboots$loa.stats$factor", 1:nfactors))
  
  loa.std <- qmboots$orig.res$loa
  loa.bts <- apply(obj.loa, 1,
                   function(x) eval(parse(text=paste0(x, "[,c('mean','sd')]"))))
  loa.frq <- apply(obj.loa, 1, 
                   function(x) eval(parse(text=paste0(x, "[,c('flag_freq')]"))))
  # Give appropriate column names
  dimnames(loa.frq) <- list(rownames(loa.bts[[1]]), paste0("flag.freq", 1:nfactors))
  colnames(loa.std) <- paste0("f", 1:nfactors, ".std")
  for (i in 1:nfactors) {
    colnames(loa.bts[[i]]) <- paste0("f", i, c(".loa", ".SE"))
  }
  
  # Reorder rows in standard results (bootstrap reorders Q-sorts alphab.)
  loa.std <- loa.std[rownames(loa.frq),]
  
  # Calculate estimate of bias
  loa.bts.est <- apply(obj.loa, 1, 
                       function(x) eval(parse(text=paste0(x, "[,'mean']"))))
  loa.bias <- loa.std - loa.bts.est 
  names(loa.bias) <- paste0("f", 1:nfactors, ".bias")
  
  # Bind together
  qs <- data.frame(loa.std, do.call("cbind", loa.bts), loa.frq, loa.bias)
  
  #-------------------------------------------------------
  # Gather results of statements
  obj.zsc <- as.array(paste0("qmboots$'zscore-stats'$factor", 1:nfactors))
  
  zsc.std <- qmboots$orig.res$zsc
  zsc.bts <- apply(obj.zsc, 1,
                   function(x) eval(parse(text=paste0(x, "[,c('mean','sd')]"))))
  # Appropriate column names
  colnames(zsc.std) <- paste0("f", 1:nfactors, ".zsc.std")
  for (i in 1:nfactors) {
    colnames(zsc.bts[[i]]) <- paste0("f", i, c(".zsc.bts", ".SE"))
  }
  
  # And factor scores
  zscn.std <- qmboots$orig.res$zsc_n
  zscn.bts <- qmboots$'zscore-stats'$'Bootstraped factor scores'
  colnames(zscn.bts) <- paste0("fsc.bts.", 1:nfactors)
  
  # Calculate estimate of bias for z-scores
  zsc.bts.est <- apply(obj.zsc, 1, 
                       function(x) eval(parse(text=paste0(x, "[,'mean']"))))
  zsc.bias <- zsc.std - zsc.bts.est 
  names(zsc.bias) <- paste0("f", 1:nfactors, ".bias")
  
  # Calculate estimate of bias for factor scores
  zscn.bias <- zscn.std - zscn.bts 
  names(zscn.bias) <- paste0("f", 1:nfactors, ".fsc.bias")
  
  # Bind together
  st <- data.frame(zsc.std, do.call("cbind", zsc.bts), zsc.bias, 
                   zscn.std, zscn.bts, zscn.bias)
  qmb <- list(qs, st)
  names(qmb) <- c("qsorts", "statements")
  return(qmb)
}
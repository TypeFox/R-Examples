
update.sum.stat <- function(sum.stat, exc.snps){
  
  if(length(exc.snps) == 0){
    return(sum.stat)
  }
  
  stat <- list()
  n <- length(sum.stat$stat)
  fid <- 0
  snps.in.study <- NULL
  lambda <- NULL
  nsamples <- list()
  ncases <- list()
  ncontrols <- list()
  for(i in 1:n){
    id <- which(!(sum.stat$stat[[i]][, "SNP"] %in% exc.snps))
    
    if(length(id) == 0){
      next
    }
    
    fid <- fid + 1
    stat[[fid]] <- sum.stat$stat[[i]][id, ]
    snps.in.study <- unique(c(snps.in.study, stat[[fid]][, "SNP"]))
    lambda[fid] <- sum.stat$lambda[i]
    nsamples[[fid]] <- sum.stat$nsamples[[i]]
    ncases[[fid]] <- sum.stat$ncases[[i]]
    ncontrols[[fid]] <- sum.stat$ncontrols[[i]]
    
  }
  
  if(length(stat) == 0){
    msg <- "All SNPs excluded in update.sum.stat"
    stop(msg)
  }
  
  snps.in.study <- unique(snps.in.study)
  SNP.sample.size <- sum.stat$SNP.sample.size[snps.in.study, ]
  
  rm(sum.stat)
  gc()
  
  list(stat = stat, snps.in.study = snps.in.study, lambda = lambda, 
       nsamples = nsamples, ncases = ncases, ncontrols = ncontrols, 
       SNP.sample.size = SNP.sample.size)
  
}


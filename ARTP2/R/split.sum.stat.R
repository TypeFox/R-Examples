
split.sum.stat <- function(sum.stat, sub.pathway){
  
  ngrp <- length(sub.pathway)
  
  ss <- list()
  nstudy <- length(sum.stat$stat)
  for(i in 1:ngrp){
    snps <- sub.pathway[[i]]$SNP
    ss[[i]] <- sum.stat
    for(j in 1:nstudy){
      ss[[i]]$stat[[j]] <- ss[[i]]$stat[[j]][ss[[i]]$stat[[j]]$SNP %in% snps, ]
    }
    ss[[i]]$snps.in.study <- snps
    ss[[i]]$SNP.sample.size <- ss[[i]]$SNP.sample.size[snps, ]
  }
  
  ss
  
}

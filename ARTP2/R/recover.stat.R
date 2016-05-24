
recover.stat <- function(sum.stat, pathway, ref.geno, allele.info, options){
  
  msg <- paste("Recovering test statistics:", date())
  if(options$print) message(msg)
  
  nsamples <- sum.stat$nsamples
  sum.info <- sum.stat$stat
  all.snp <- sort(sum.stat$snps.in.study)
  rm(sum.stat)
  gc()
  
  nstudy <- length(sum.info)
  
  chr <- unique(pathway$Chr)
  ngrp <- length(unique(chr))
  
  ref.cor <- list()
  for(i in 1:ngrp){
    grp.snps <- pathway$SNP[pathway$Chr == chr[i]]
    grp.snps <- unique(grp.snps)
    suppressWarnings(ref.cor[[i]] <- cor(ref.geno[, grp.snps, drop = FALSE], use = "pairwise.complete.obs", method = "pearson"))
    id <- which(colnames(ref.geno) %in% grp.snps)
    ref.geno <- ref.geno[, -id, drop = FALSE]
    ref.cor[[i]][is.na(ref.cor[[i]])] <- 0
  }
  names(ref.cor) <- as.character(chr)
  
  score0 <- list()
  ra <- list()
  rownames(allele.info) <- allele.info$SNP
  for(i in 1:ngrp){
    grp.snps <- colnames(ref.cor[[i]])
    nsnp <- length(grp.snps)
    score0[[i]] <- rep(0, nsnp)
    ra[[i]] <- rep(NA, nsnp)
    names(score0[[i]]) <- grp.snps
    names(ra[[i]]) <- grp.snps
    
    for(rs in grp.snps){
      for(k in 1:nstudy){
        if(rs %in% sum.info[[k]][, "SNP"]){
          nr <- which(sum.info[[k]][, "SNP"] == rs)
          
          if(is.na(ra[[i]][rs])){
            ra[[i]][rs] <- sum.info[[k]][nr, "RefAllele"]
          }
          
          sgn <- ifelse(ra[[i]][rs] == sum.info[[k]][nr, "RefAllele"], 1, -1)
          score0[[i]][rs] <- score0[[i]][rs] + sgn * sum.info[[k]][nr, "BETA"] / sum.info[[k]][nr, "SE"]^2
        }
      }
    }
  }
  names(score0) <- as.character(chr)
  names(ra) <- as.character(chr)
  
  for(i in 1:ngrp){
    grp.snps <- colnames(ref.cor[[i]])
    for(rs in grp.snps){
      if(allele.info[rs, "RefAllele"] != ra[[i]][rs]){
        ref.cor[[i]][rs, ] <- -ref.cor[[i]][rs, ]
        ref.cor[[i]][, rs] <- -ref.cor[[i]][, rs]
      }
    }
  }
  
  wt <- list()
  V <- list()
  for(i in 1:ngrp){
    grp.snps <- colnames(ref.cor[[i]])
    nsnp <- length(grp.snps)
    wt[[i]] <- matrix(0, nsnp, nsnp)
    rownames(wt[[i]]) <- grp.snps
    colnames(wt[[i]]) <- grp.snps
    
    for(k in 1:nstudy){
      id <- which(grp.snps %in% sum.info[[k]][, "SNP"])
      if(length(id) == 0){
        next
      }
      ks <- grp.snps[id]
      
      rownames(sum.info[[k]]) <- sum.info[[k]][, "SNP"]
      tmp <- strsplit(sum.info[[k]][ks, "Direction", drop = TRUE], '')
      tmp <- sapply(tmp, as.integer)
      
      if(is.vector(tmp)){
        tmp <- matrix(tmp, nrow = 1)
      }
      colnames(tmp) <- sum.info[[k]][ks, "SNP"]
      
      es <- t(tmp) %*% (tmp * nsamples[[k]])
      se <- sum.info[[k]][ks, "SE"]
      
      wt[[i]][ks, ks] <- wt[[i]][ks, ks] + es * outer(diag(es), diag(es), "*")^(-.5) / outer(se, se, "*")
    }
    
    V[[i]] <- wt[[i]] * ref.cor[[i]]
    score0[[i]] <- score0[[i]][colnames(V[[i]])]
  }
  
  max.total.N <- sum(unlist(nsamples))
  
  for(i in 1:length(V)){
    rs <- sort(names(score0[[i]]))
    score0[[i]] <- score0[[i]][rs] / sqrt(max.total.N)
    V[[i]] <- V[[i]][rs, rs, drop = FALSE] / max.total.N
  }
  
  names(V) <- as.character(chr)
  names(score0) <- as.character(chr)
  
  list(V = V, score0 = score0)
  
}


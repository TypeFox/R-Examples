
filter.raw.geno <- function(raw.geno, pathway, options, control.id = NULL, print = FALSE){
  
  # initialize them in case all filters are turned off
  deleted.snps <- data.frame(SNP = NULL, reason = NULL, comment = NULL, stringsAsFactors = FALSE)
  deleted.genes <- data.frame(Gene = NULL, reason = NULL, stringsAsFactors = FALSE)
  rs <- colnames(raw.geno)
  
  pathway <- pathway[pathway$SNP %in% colnames(raw.geno), ]
  
  if(options$snp.miss.rate < 1){
    msg <- paste("Removing SNPs with high missing rate:", date())
    if(options$print && print) message(msg)
    
    # use for loop to avoid high memory consumption
    snp.miss.rate <- rep(NA, ncol(raw.geno))
    for(i in 1:ncol(raw.geno)){
      snp.miss.rate[i] <- mean(is.na(raw.geno[, i]))
      if((i %% 50) == 0){
        gc()
      }
    }
    gc()
    names(snp.miss.rate) <- rs
    #snp.miss.rate <- apply(raw.geno, 2, function(x){mean(is.na(x))})
    id <- which(snp.miss.rate > options$snp.miss.rate)
    if(length(id) > 0){
      exc.snps <- colnames(raw.geno)[id]
      del.snps <- data.frame(SNP = exc.snps, reason = "SNP_MISS_RATE", comment = snp.miss.rate[id], stringsAsFactors = FALSE)
      raw.geno <- raw.geno[, -id, drop = FALSE]
      rs <- setdiff(rs, exc.snps)
      if(length(rs) == 0){
        msg <- "No SNPs were left due to SNP_MISS_RATE"
        #stop(msg)
        return(list(deleted.snps = deleted.snps, deleted.genes = deleted.genes))
      }
      deleted.snps <- rbind(deleted.snps, del.snps)
      pathway <- pathway[!(pathway$SNP %in% exc.snps), ]
    }
  }
  
  ######
  
  if(options$maf > 0){
    msg <- paste("Removing SNPs with low MAFs:", date())
    if(options$print && print) message(msg)
    
    maf <- rep(NA, ncol(raw.geno))
    for(i in 1:ncol(raw.geno)){
      m <- mean(raw.geno[, i], na.rm = TRUE)/2
      if((i %% 50) == 0){
        gc()
      }
      maf[i] <- min(m, 1-m)
    }
    gc()
    names(maf) <- colnames(raw.geno)
    #maf <- apply(raw.geno, 2, function(x){m <- mean(x, na.rm = TRUE)/2; pmin(m, 1-m)})
    id <- which(maf < options$maf)
    if(length(id) > 0){
      exc.snps <- colnames(raw.geno)[id]
      del.snps <- data.frame(SNP = exc.snps, reason = "SNP_LOW_MAF", comment = maf[id], stringsAsFactors = FALSE)
      raw.geno <- raw.geno[, -id, drop = FALSE]
      rs <- setdiff(rs, exc.snps)
      if(length(rs) == 0){
        msg <- "No SNPs were left due to SNP_LOW_MAF"
        #stop(msg)
        return(list(deleted.snps = deleted.snps, deleted.genes = deleted.genes))
      }
      deleted.snps <- rbind(deleted.snps, del.snps)
      pathway <- pathway[!(pathway$SNP %in% exc.snps), ]
    }
  }
  
  ######
  
  if(TRUE){
    msg <- paste("Removing constant SNPs:", date())
    if(options$print && print) message(msg)
    
    SD <- rep(NA, ncol(raw.geno))
    for(i in 1:ncol(raw.geno)){
      SD[i] <- sd(raw.geno[, i], na.rm = TRUE)
      if((i %% 50) == 0){
        gc()
      }
    }
    gc()
    names(SD) <- colnames(raw.geno)
    #SD <- apply(raw.geno, 2, function(x){sd(x, na.rm = TRUE)})
    id <- which(SD == 0)
    if(length(id) > 0){
      exc.snps <- colnames(raw.geno)[id]
      del.snps <- data.frame(SNP = exc.snps, reason = "SNP_CONST", comment = "", stringsAsFactors = FALSE)
      raw.geno <- raw.geno[, -id, drop = FALSE]
      rs <- setdiff(rs, exc.snps)
      if(length(rs) == 0){
        msg <- "No SNPs were left due to SNP_CONST"
        #stop(msg)
        return(list(deleted.snps = deleted.snps, deleted.genes = deleted.genes))
      }
      deleted.snps <- rbind(deleted.snps, del.snps)
      pathway <- pathway[!(pathway$SNP %in% exc.snps), ]
    }
  }
  
  #######
  
  if(options$HWE.p > 0){
    msg <- paste("Removing SNPs fail to pass HWE test:", date())
    if(options$print && print) message(msg)
    
    if(is.null(control.id) || length(control.id) == 0){
      hwe.pval <- rep(NA, ncol(raw.geno))
      for(i in 1:ncol(raw.geno)){
        hwe.pval[i] <- HWE.exact(raw.geno[, i])
        if((i %% 50) == 0){
          gc()
        }
      }
      gc()
      names(hwe.pval) <- colnames(raw.geno)
      #hwe.pval <- apply(raw.geno, 2, HWE.exact)
    }else{
      if(max(control.id) > nrow(raw.geno)){
        msg <- "Fail to perform HWE test"
        stop(msg)
      }
      hwe.pval <- rep(NA, ncol(raw.geno))
      for(i in 1:ncol(raw.geno)){
        hwe.pval[i] <- HWE.exact(raw.geno[control.id, i])
        if((i %% 50) == 0){
          gc()
        }
      }
      gc()
      names(hwe.pval) <- colnames(raw.geno)
      #hwe.pval <- apply(raw.geno[control.id, ], 2, HWE.exact)
    }
    id <- which(hwe.pval < options$HWE.p)
    if(length(id) > 0){
      exc.snps <- colnames(raw.geno)[id]
      del.snps <- data.frame(SNP = exc.snps, reason = "SNP_HWE", comment = "", stringsAsFactors = FALSE)
      raw.geno <- raw.geno[, -id, drop = FALSE]
      rs <- setdiff(rs, exc.snps)
      if(length(rs) == 0){
        msg <- "No SNPs were left due to SNP_HWE"
        #stop(msg)
        return(list(deleted.snps = deleted.snps, deleted.genes = deleted.genes))
      }
      deleted.snps <- rbind(deleted.snps, del.snps)
      pathway <- pathway[!(pathway$SNP %in% exc.snps), ]
    }
  }
  
  #######
  
  if(options$gene.R2 < 1){
    msg <- paste("Removing high LD SNPs within genes:", date())
    if(options$print && print) message(msg)
    
    gene <- unique(pathway$Gene)
    exc.snps <- NULL
    exc.snps2 <- NULL
    comment <- NULL
    for(g in gene){
      snps.in.gene <- pathway$SNP[pathway$Gene == g]
      snps.in.gene <- intersect(snps.in.gene, rs)
      snps.in.gene <- sort(snps.in.gene)
      if(length(snps.in.gene) <= 1){
        next
      }
      
      rg <- raw.geno[, snps.in.gene, drop = FALSE]
      suppressWarnings(cor2 <- cor(rg, use = "pairwise.complete.obs", method = "pearson")^2)
      cor2[is.na(cor2)] <- 0 # Sometimes two SNPs are approximately independent due to the location of missing, see rs4673651 and rs16847776 in 1000 Genomes EUR as an example
      diag(cor2) <- -1
      
      tmp <- apply(cor2, 2, function(x){max(x, na.rm = TRUE) > options$gene.R2})
      if(!any(tmp)){
        next
      }
      rg <- rg[, tmp, drop = FALSE]
      cor2 <- cor2[tmp, tmp, drop = FALSE]
      snps.in.gene <- snps.in.gene[tmp]
      
      maf <- apply(rg, 2, function(x){m <- mean(x, na.rm = TRUE)/2; pmin(m, 1-m)})
      rm(rg)
      gc()
      names(maf) <- snps.in.gene
      
      if(length(snps.in.gene) > options$huge.gene.size && options$trim.huge.chr){
        tmp <- order(maf)
        cor2 <- cor2[tmp, tmp]
        cor2[lower.tri(cor2)] <- -1
        maf <- maf[tmp]
        snps.in.gene <- snps.in.gene[tmp]
        count <- apply(cor2, 1, function(x){any(x > options$huge.gene.R2, na.rm = TRUE)})
        if(!any(count)){
          next
        }
        rm.snps <- snps.in.gene[count]
        exc.snps2 <- c(exc.snps2, rm.snps)
      }else{
        while(1){
          if(nrow(cor2) == 1){
            break
          }
          max.r2 <- max(cor2, na.rm = TRUE)
          if(max.r2 <= options$gene.R2){
            break
          }
          id <- which(cor2 == max.r2, arr.ind = TRUE)
          s <- unique(rownames(id))
          snp.lower.maf <- names(which.min(maf[s]))
          k <- which(colnames(cor2) == snp.lower.maf)
          exc.snps <- c(exc.snps, snp.lower.maf)
          tmp <- c(cor2[snp.lower.maf, ], cor2[, snp.lower.maf])
          tmp <- tmp[tmp > 0]
          cc <- paste(names(tmp)[which(tmp == max.r2)[1]], round(max.r2, 3), sep = "_")
          comment <- c(comment, cc)
          cor2 <- cor2[-k, -k, drop = FALSE]
        }
      }
      rm(cor2)
      gc()
    }
    
    if(!is.null(exc.snps)){
      del.snps <- data.frame(SNP = exc.snps, reason = "GENE_R2", comment = comment, stringsAsFactors = FALSE)
      id <- which(colnames(raw.geno) %in% exc.snps)
      raw.geno <- raw.geno[, -id, drop = FALSE]
      rs <- setdiff(rs, exc.snps)
      if(length(rs) == 0){
        msg <- "No SNPs were left due to GENE_R2"
        #stop(msg)
        return(list(deleted.snps = deleted.snps, deleted.genes = deleted.genes))
      }
      deleted.snps <- rbind(deleted.snps, del.snps)
      pathway <- pathway[!(pathway$SNP %in% exc.snps), ]
    }
    
    if(!is.null(exc.snps2)){
      exc.snps2 <- unique(exc.snps2)
      if(!is.null(exc.snps)){
        exc.snps2 <- setdiff(exc.snps2, exc.snps)
      }
      if(!is.null(exc.snps2)){
        del.snps2 <- data.frame(SNP = exc.snps2, reason = "HUGE_GENE_R2", comment = "", stringsAsFactors = FALSE)
        id <- which(colnames(raw.geno) %in% exc.snps2)
        raw.geno <- raw.geno[, -id, drop = FALSE]
        rs <- setdiff(rs, exc.snps2)
        if(length(rs) == 0){
          msg <- "No SNPs were left due to HUGE_GENE_R2"
          #stop(msg)
          return(list(deleted.snps = deleted.snps, deleted.genes = deleted.genes))
        }
        deleted.snps <- rbind(deleted.snps, del.snps2)
        pathway <- pathway[!(pathway$SNP %in% exc.snps2), ]
      }
    }
  }
  
  ######
  
  if(options$chr.R2 < 1){
    msg <- paste("Removing SNPs in high LD within chromosomes:", date())
    if(options$print && print) message(msg)
    
    chr <- unique(pathway$Chr)
    exc.snps <- NULL
    comment <- NULL
    for(c in chr){
      snps.in.chr <- pathway$SNP[pathway$Chr == c]
      snps.in.chr <- intersect(snps.in.chr, rs)
      snps.in.chr <- sort(snps.in.chr)
      if(length(snps.in.chr) <= 1){
        next
      }
      rg <- raw.geno[, snps.in.chr, drop = FALSE]
      suppressWarnings(cor2 <- cor(rg, use = "pairwise.complete.obs", method = "pearson")^2)
      cor2[is.na(cor2)] <- 0
      diag(cor2) <- -1
      
      tmp <- apply(cor2, 2, function(x){max(x, na.rm = TRUE) > options$chr.R2})
      if(!any(tmp)){
        next
      }
      rg <- rg[, tmp, drop = FALSE]
      cor2 <- cor2[tmp, tmp, drop = FALSE]
      snps.in.chr <- snps.in.chr[tmp]
      
      maf <- apply(rg, 2, function(x){m <- mean(x, na.rm = TRUE)/2; pmin(m, 1-m)})
      rm(rg)
      gc()
      names(maf) <- snps.in.chr
      while(1){
        if(nrow(cor2) == 1){
          break
        }
        max.r2 <- max(cor2, na.rm = TRUE)
        if(max.r2 <= options$chr.R2){
          break
        }
        id <- which(!is.na(cor2) & (cor2 == max.r2), arr.ind = TRUE)
        s <- unique(rownames(id))
        snp.lower.maf <- names(which.min(maf[s]))
        k <- which(colnames(cor2) == snp.lower.maf)
        exc.snps <- c(exc.snps, snp.lower.maf)
        tmp <- c(cor2[snp.lower.maf, ], cor2[, snp.lower.maf])
        tmp <- tmp[tmp > 0]
        cc <- paste(names(tmp)[which.max(tmp)], round(max.r2, 3), sep = "_")
        comment <- c(comment, cc)
        cor2 <- cor2[-k, -k, drop = FALSE]
      }
      rm(cor2)
      gc()
    }
    
    if(!is.null(exc.snps)){
      del.snps <- data.frame(SNP = exc.snps, reason = "CHR_R2", comment = comment, stringsAsFactors = FALSE)
      id <- which(colnames(raw.geno) %in% exc.snps)
      raw.geno <- raw.geno[, -id, drop = FALSE]
      rs <- setdiff(rs, exc.snps)
      if(length(rs) == 0){
        msg <- "No SNPs were left due to CHR_R2"
        #stop(msg)
        return(list(deleted.snps = deleted.snps, deleted.genes = deleted.genes))
      }
      deleted.snps <- rbind(deleted.snps, del.snps)
      pathway <- pathway[!(pathway$SNP %in% exc.snps), ]
    }
  }
  
  #########
  
  chr.name <- unique(pathway$Chr)
  chr.size <- NULL
  for(chr in chr.name){
    chr.size <- c(chr.size, length(unique(pathway$SNP[pathway$Chr == chr])))
  }
  
  large.chr <- chr.name[chr.size > options$huge.chr.size]
  
  if(length(large.chr) > 0 && options$trim.huge.chr){
    
    if(options$huge.gene.R2 < 1){
      msg <- paste("Removing high LD SNPs within genes of huge chromosome:", date())
      if(options$print && print) message(msg)
      
      gene <- unique(pathway$Gene)
      exc.snps <- NULL
      exc.snps2 <- NULL
      comment <- NULL
      for(g in gene){
        chr <- unique(pathway$Chr[pathway$Gene == g])
        if(!(chr %in% large.chr)){
          next
        }
        
        snps.in.gene <- pathway$SNP[pathway$Gene == g]
        snps.in.gene <- intersect(snps.in.gene, rs)
        snps.in.gene <- sort(snps.in.gene)
        if(length(snps.in.gene) <= 1){
          next
        }
        
        rg <- raw.geno[, snps.in.gene, drop = FALSE]
        suppressWarnings(cor2 <- cor(rg, use = "pairwise.complete.obs", method = "pearson")^2)
        cor2[is.na(cor2)] <- 0 # Sometimes two SNPs are approximately independent due to the location of missing, see rs4673651 and rs16847776 in 1000 Genomes EUR as an example
        diag(cor2) <- -1
        
        tmp <- apply(cor2, 2, function(x){max(x, na.rm = TRUE) > options$huge.gene.R2})
        if(!any(tmp)){
          next
        }
        rg <- rg[, tmp, drop = FALSE]
        cor2 <- cor2[tmp, tmp, drop = FALSE]
        snps.in.gene <- snps.in.gene[tmp]
        
        maf <- apply(rg, 2, function(x){m <- mean(x, na.rm = TRUE)/2; pmin(m, 1-m)})
        names(maf) <- snps.in.gene
        
        if(length(snps.in.gene) > options$huge.gene.size){
          tmp <- order(maf)
          cor2 <- cor2[tmp, tmp]
          cor2[lower.tri(cor2)] <- -1
          maf <- maf[tmp]
          snps.in.gene <- snps.in.gene[tmp]
          count <- apply(cor2, 1, function(x){any(x > options$huge.gene.R2, na.rm = TRUE)})
          if(!any(count)){
            next
          }
          rm.snps <- snps.in.gene[count]
          exc.snps2 <- c(exc.snps2, rm.snps)
        }else{
          while(1){
            if(nrow(cor2) == 1){
              break
            }
            max.r2 <- max(cor2, na.rm = TRUE)
            if(max.r2 <= options$huge.gene.R2){
              break
            }
            id <- which(cor2 == max.r2, arr.ind = TRUE)
            s <- unique(rownames(id))
            snp.lower.maf <- names(which.min(maf[s]))
            k <- which(colnames(cor2) == snp.lower.maf)
            exc.snps <- c(exc.snps, snp.lower.maf)
            tmp <- c(cor2[snp.lower.maf, ], cor2[, snp.lower.maf])
            tmp <- tmp[tmp > 0]
            cc <- paste(names(tmp)[which(tmp == max.r2)[1]], round(max.r2, 3), sep = "_")
            comment <- c(comment, cc)
            cor2 <- cor2[-k, -k, drop = FALSE]
          }
        }
        rm(cor2)
        gc()
      }
      
      if(!is.null(exc.snps)){
        del.snps <- data.frame(SNP = exc.snps, reason = "HUGE_CHR", comment = comment, stringsAsFactors = FALSE)
        id <- which(colnames(raw.geno) %in% exc.snps)
        raw.geno <- raw.geno[, -id, drop = FALSE]
        rs <- setdiff(rs, exc.snps)
        if(length(rs) == 0){
          msg <- "No SNPs were left due to HUGE_CHR"
          #stop(msg)
          return(list(deleted.snps = deleted.snps, deleted.genes = deleted.genes))
        }
        deleted.snps <- rbind(deleted.snps, del.snps)
        pathway <- pathway[!(pathway$SNP %in% exc.snps), ]
      }
      
      if(!is.null(exc.snps2)){
        exc.snps2 <- unique(exc.snps2)
        if(!is.null(exc.snps)){
          exc.snps2 <- setdiff(exc.snps2, exc.snps)
        }
        if(!is.null(exc.snps2)){
          del.snps2 <- data.frame(SNP = exc.snps2, reason = "HUGE_CHR2", comment = "", stringsAsFactors = FALSE)
          id <- which(colnames(raw.geno) %in% exc.snps2)
          raw.geno <- raw.geno[, -id, drop = FALSE]
          rs <- setdiff(rs, exc.snps2)
          if(length(rs) == 0){
            msg <- "No SNPs were left due to HUGE_CHR2"
            #stop(msg)
            return(list(deleted.snps = deleted.snps, deleted.genes = deleted.genes))
          }
          deleted.snps <- rbind(deleted.snps, del.snps2)
          pathway <- pathway[!(pathway$SNP %in% exc.snps2), ]
        }
      }
    }
    
    ######
    
    if(options$huge.chr.R2 < 1){
      msg <- paste("Removing SNPs in high LD within huge chromosomes:", date())
      if(options$print && print) message(msg)
      
      chr <- unique(pathway$Chr)
      exc.snps <- NULL
      comment <- NULL
      for(c in chr){
        if(!(c %in% large.chr)){
          next
        }
        
        snps.in.chr <- pathway$SNP[pathway$Chr == c]
        snps.in.chr <- intersect(snps.in.chr, rs)
        snps.in.chr <- sort(snps.in.chr)
        if(length(snps.in.chr) <= 1){
          next
        }
        rg <- raw.geno[, snps.in.chr, drop = FALSE]
        suppressWarnings(cor2 <- cor(rg, use = "pairwise.complete.obs", method = "pearson")^2)
        cor2[is.na(cor2)] <- 0
        diag(cor2) <- -1
        
        tmp <- apply(cor2, 2, function(x){max(x, na.rm = TRUE) > options$huge.chr.R2})
        if(!any(tmp)){
          next
        }
        rg <- rg[, tmp, drop = FALSE]
        cor2 <- cor2[tmp, tmp, drop = FALSE]
        snps.in.chr <- snps.in.chr[tmp]
        
        maf <- apply(rg, 2, function(x){m <- mean(x, na.rm = TRUE)/2; pmin(m, 1-m)})
        names(maf) <- snps.in.chr
        while(1){
          if(nrow(cor2) == 1){
            break
          }
          max.r2 <- max(cor2, na.rm = TRUE)
          if(max.r2 <= options$huge.chr.R2){
            break
          }
          id <- which(!is.na(cor2) & (cor2 == max.r2), arr.ind = TRUE)
          s <- unique(rownames(id))
          snp.lower.maf <- names(which.min(maf[s]))
          k <- which(colnames(cor2) == snp.lower.maf)
          exc.snps <- c(exc.snps, snp.lower.maf)
          tmp <- c(cor2[snp.lower.maf, ], cor2[, snp.lower.maf])
          tmp <- tmp[tmp > 0]
          cc <- paste(names(tmp)[which.max(tmp)], round(max.r2, 3), sep = "_")
          comment <- c(comment, cc)
          cor2 <- cor2[-k, -k, drop = FALSE]
        }
        rm(cor2)
        gc()
      }
      
      if(!is.null(exc.snps)){
        del.snps <- data.frame(SNP = exc.snps, reason = "HUGE_CHR3", comment = comment, stringsAsFactors = FALSE)
        id <- which(colnames(raw.geno) %in% exc.snps)
        raw.geno <- raw.geno[, -id, drop = FALSE]
        rs <- setdiff(rs, exc.snps)
        if(length(rs) == 0){
          msg <- "No SNPs were left due to HUGE_CHR3"
          #stop(msg)
          return(list(deleted.snps = deleted.snps, deleted.genes = deleted.genes))
        }
        deleted.snps <- rbind(deleted.snps, del.snps)
        pathway <- pathway[!(pathway$SNP %in% exc.snps), ]
      }
    }
    
  }
  
  #########
  
  if(options$gene.miss.rate < 1){
    msg <- paste("Removing genes with high missing rate:", date())
    if(options$print && print) message(msg)
    
    gene <- unique(pathway$Gene)
    exc.snps <- NULL
    comment <- NULL
    for(g in gene){
      snps.in.gene <- pathway$SNP[pathway$Gene == g]
      snps.in.gene <- intersect(snps.in.gene, rs)
      snps.in.gene <- sort(snps.in.gene)
      cc <- !complete.cases(raw.geno[, snps.in.gene, drop = FALSE])
      gmr <- mean(cc)
      if(gmr > options$gene.miss.rate){
        exc.snps <- c(exc.snps, snps.in.gene)
        comment <- c(comment, gmr)
        pathway <- pathway[pathway$Gene != g, ]
        raw.geno <- raw.geno[, colnames(raw.geno) %in% pathway$SNP, drop = FALSE]
        rs <- colnames(raw.geno)
      }
    }
    
    if(!is.null(exc.snps)){
      del.snps <- data.frame(SNP = exc.snps, reason = "GENE_MISS_RATE", comment = comment, stringsAsFactors = FALSE)
      id <- which(colnames(raw.geno) %in% exc.snps)
      if(length(id) > 0){
        raw.geno <- raw.geno[, -id, drop = FALSE]
      }
      rs <- setdiff(rs, exc.snps)
      if(length(rs) == 0){
        msg <- "No SNPs were left due to GENE_MISS_RATE"
        #stop(msg)
        return(list(deleted.snps = deleted.snps, deleted.genes = deleted.genes))
      }
      deleted.snps <- rbind(deleted.snps, del.snps)
      pathway <- pathway[!(pathway$SNP %in% exc.snps), ]
    }
  }
  
  #########
  
  del.genes <- NULL
  reason <- NULL
  if(options$rm.gene.subset){
    msg <- paste("Removing genes which are subsets of other genes:", date())
    if(options$print && print) message(msg)
    
    chr <- unique(pathway$Chr)
    nc <- length(chr)
    for(i in 1:nc){
      pa <- pathway[pathway$Chr == chr[i], ]
      gene <- unique(pa$Gene)
      ng <- length(gene)
      if(ng == 1){
        next
      }
      
      ns <- rep(NA, ng)
      for(j in 1:ng){
        ns[j] <- sum(pa$Gene == gene[j])
      }
      
      gene <- gene[order(ns)]
      
      for(j in 1:(ng - 1)){
        g1 <- unique(pa$SNP[pa$Gene == gene[j]])
        for(k in (j+1):ng){
          g2 <- unique(pa$SNP[pa$Gene == gene[k]])
          if(length(g1) > length(g2)){ # impossible
            msg <- 'debug filter.raw.geno'
            stop(msg)
            next
          }
          
          if(all(g1 %in% g2)){
            del.genes <- c(del.genes, gene[j])
            reason <- c(reason, paste('Subset of ', gene[k], sep = ''))
            break
          }
        }
      }
    }
    
    if(!is.null(del.genes)){
      names(reason) <- del.genes
      deleted.genes <- data.frame(Gene = del.genes, reason = reason, stringsAsFactors = FALSE)
      pathway <- pathway[!(pathway$Gene %in% del.genes), ]
    }
  }
  
  #########
  
  list(deleted.snps = deleted.snps, deleted.genes = deleted.genes)
  
}


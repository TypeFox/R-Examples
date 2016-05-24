summary.minPtest <-
function(object,level=0.05,sign.SNP=FALSE, ...){
  if(length(which(object$p.adj.minp<=level))==0){
    stop("no permuation-based p-value is smaller or equal to the chosen level")
  }
  index.gene <- which(object$p.adj.minp<=level)
  gene.name <- rownames(object$p.adj.minp)[index.gene]
  p.list <- lapply(seq_along(gene.name), function(i){
    gene <- gene.name[i]
    snpnames <- object$SNPtoGene[which(object$SNPtoGene[,2]==gene.name[i]),1]
    if(sign.SNP){
      if(length(which(object$p.adj.psnp[snpnames,]<=level))==0){
        stop("no marginal p-value is smaller or equal to the chosen level for the SNPs located on one of the chosen genes", call.=FALSE)
      }
      p.adjust.temp <- object$p.adj.psnp[snpnames,]
      p.adjust.index <- which(p.adjust.temp<=level)
      snp.p.adjust <- sort(p.adjust.temp[p.adjust.index])
      p_value.temp <- object$psnp[snpnames,]
      if(length(snpnames)==1){
        names(snp.p.adjust) <- snpnames
        names(p_value.temp) <- snpnames
      }
      snp.p_value <- p_value.temp[names(snp.p.adjust)]
    }else{
      snp.p.adjust <- sort(object$p.adj.psnp[snpnames,])
      p_value.temp <- object$psnp[snpnames,]
      if(length(snpnames)==1){
        names(snp.p.adjust) <- snpnames
        names(p_value.temp) <- snpnames
      }
      snp.p_value <- p_value.temp[names(snp.p.adjust)]
    }
    snp.temp <- data.frame(names(snp.p.adjust), snp.p_value, snp.p.adjust)
    colnames(snp.temp) <- c("SNP", "snp.p_value","snp.p.adjust")
    rownames(snp.temp) <- NULL
    gene.temp <- data.frame(gene,object$minp[gene,],object$p.adj.minp[gene,])
    colnames(gene.temp) <- c("Gene","minP","gene.p.adjust")
    rownames(gene.temp) <- NULL
    list(gene.p.values = gene.temp, snp.p.values = snp.temp)
  })
  names(p.list) <- gene.name
  class(p.list) <- "summary.minPtest"
  return(p.list)
}

print.summary.minPtest <- function(x, ...){
  cat("\np-values:\n")
  p.table.temp <- lapply(x, function(i){
    test <- i
    list.1 <- do.call(cbind, i)
    if(dim(i$snp.p.values)[1]>1){
      list.1[2:(dim(list.1)[1]),1:3] <- as.factor(NA)
    }
    list.1
  })
  p.table <- do.call(rbind, p.table.temp)
  rownames(p.table) <- NULL
  colnames(p.table) <- c("Gene", "minP", "gene.p.adjust", "SNP", "snp.p.value", "snp.p.adjust")
  print(p.table, na.print="")
  invisible()
}


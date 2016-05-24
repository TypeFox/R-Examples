
norm.stat.test <- function(setup){
  
  start.time <- date()
  
  options <- setup$options
  norm.stat <- setup$norm.stat
  pathway <- setup$pathway
  allele.info <- setup$allele.info
  rm(setup)
  gc()
  
  options <- check.os(options)
  
  V <- norm.stat$V
  score0 <- norm.stat$score0
  
  ngrp <- length(V)
  
  group.id <- NULL
  gene.id <- NULL
  gene.pval <- NULL
  N.SNP <- NULL
  gene.name <- NULL
  chr <- NULL
  model <- NULL
  unadj.pvalue <- NULL
  for(g in 1:ngrp){
    msg <- paste("Permuting group ", names(V)[[g]], ": ", date(), sep = "")
    if(options$print) message(msg)
    
    U <- cov.svd(V[[g]], names(V)[[g]])
    sc <- score0[[g]]
    rs <- names(sc)
    group.setup <- create.group(pathway, rs)
    N.SNP <- c(N.SNP, group.setup$N.SNP)
    gene.name <- c(gene.name, group.setup$GeneInGroup)
    chr <- c(chr, rep(names(V)[[g]], length(group.setup$GeneInGroup)))
    gene.cutpoint.setup <- create.gene.cutpoint(pathway, rs, options)
    ngene <- length(group.setup$GeneInGroup)
    gene.id <- c(gene.id, 1:ngene)
    group.id <- c(group.id, rep(g, ngene))
    gpv <- artp2.chr(group.setup, gene.cutpoint.setup, U, sc, V[[g]], options, g)
    gene.pval <- c(gene.pval, gpv$gene.pval)
    model <- c(model, gpv$model)
  }
  
  msg <- paste0("Permutation completed: ", date())
  if(options$print) message(msg)
  
  pathway.cutpoint <- create.pathway.cutpoint(pathway, options)
  
  ppv <- artp2(group.id, gene.id, pathway.cutpoint, gene.name, options)
  
  pathway.pvalue <- ppv$pathway.pvalue
  most.sig.genes <- ppv$most.sig.genes
  arr.rank <- ppv$arr.rank
  
  gene.pvalue <- data.frame(Gene = gene.name, Group = as.integer(chr), N.SNP = N.SNP, Pvalue = gene.pval, stringsAsFactors = FALSE)
  if(!is.null(allele.info)){
    tmp <- merge(pathway[, c('SNP', 'Gene')], allele.info[, c('SNP', 'Chr')], by = 'SNP')[, c('Gene', 'Chr')]
    tmp <- tmp[!duplicated(tmp), ]
    gene.pvalue <- merge(gene.pvalue, tmp, by = 'Gene')
    gene.pvalue <- gene.pvalue[, c('Gene', 'Chr', 'Group', 'N.SNP', 'Pvalue')]
    if(all(gene.pvalue$Chr == gene.pvalue$Group)){
      gene.pvalue <- gene.pvalue[, c('Gene', 'Chr', 'N.SNP', 'Pvalue')]
    }
  }else{
    colnames(gene.pvalue) <- c('Gene', 'Chr', 'N.SNP', 'Pvalue')
  }
  
  gene.pvalue <- gene.pvalue[order(gene.pvalue$Pvalue), ]
  rownames(gene.pvalue) <- NULL
  
  ac1 <- sqrt(pathway.pvalue * (1 - pathway.pvalue) / options$nperm) / pathway.pvalue
  min.gene.pvalue <- min(gene.pvalue$Pvalue)
  ac2 <- ceiling(-log10(min.gene.pvalue))
  accurate <- ifelse(ac1 < .1 && ac2 <= log10(options$nperm), TRUE, FALSE)
  
  end.time <- date()
  test.timing <- as.integer(difftime(strptime(end.time, "%c"), strptime(start.time, "%c"), units = "secs"))
  
  list(gene.pvalue = gene.pvalue, pathway.pvalue = pathway.pvalue, most.sig.genes = most.sig.genes, 
       model = model, accurate = accurate, test.timing = test.timing, arr.rank = arr.rank)
  
}






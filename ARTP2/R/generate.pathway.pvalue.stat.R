
generate.pathway.pvalue.stat <- function(setup){
  
  start.time <- date()
  
  options <- setup$options
  norm.stat <- setup$norm.stat
  super.pathway <- setup$super.pathway
  pathway <- setup$pathway
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
  gene.stat.file <- NULL
  for(g in 1:ngrp){
    msg <- paste("Permuting group ", names(V)[[g]], ": ", date(), sep = "")
    if(options$print) message(msg)
    
    U <- cov.svd(V[[g]], names(V)[[g]])
    sc <- score0[[g]]
    rs <- names(sc)
    group.setup <- create.group(super.pathway, rs)
    N.SNP <- c(N.SNP, group.setup$N.SNP)
    gene.name <- c(gene.name, group.setup$GeneInGroup)
    chr <- c(chr, rep(names(V)[[g]], length(group.setup$GeneInGroup)))
    gene.cutpoint.setup <- create.gene.cutpoint(super.pathway, rs, options)
    ngene <- length(group.setup$GeneInGroup)
    gene.id <- c(gene.id, 1:ngene)
    group.id <- c(group.id, rep(g, ngene))
    gpv <- artp2.chr(group.setup, gene.cutpoint.setup, U, sc, V[[g]], options, g)
    gene.pval <- c(gene.pval, gpv$gene.pval)
    model <- c(model, gpv$model)
    
    gene.stat.file <- c(gene.stat.file, paste0(options$out.dir, "/", options$id.str, ".CID.", g - 1, ".GID.", 0:(ngene-1), '.bin'))
  }
  
  names(gene.stat.file) <- gene.name
  names(gene.id) <- gene.name
  names(group.id) <- gene.name
  names(gene.name) <- gene.name
  names(gene.pval) <- gene.name
  names(chr) <- gene.name
  names(N.SNP) <- gene.name
  
  msg <- paste0("Permutation completed: ", date())
  if(options$print) message(msg)
  
  npath <- length(pathway)
  pathway.pval.stat <- NULL
  gene.pvalue <- list()
  for(i in 1:npath){
    msg <- paste0("Computing p-value of pathway ", i, ": ", date())
    if(options$print) message(msg)
    
    pathway.cutpoint <- create.pathway.cutpoint(pathway[[i]], options)
    gn <- unique(pathway[[i]][, 'Gene'])
    pps <- artp2.select.genes(group.id[gn], gene.id[gn], pathway.cutpoint, gene.name[gn], options)
    if(is.null(pathway.pval.stat)){
      pathway.pval.stat <- pps
    }else{
      pathway.pval.stat <- cbind(pathway.pval.stat, pps)
    }
    
    gene.pvalue[[i]] <- data.frame(Gene = gene.name[gn], Chr = as.integer(chr[gn]), N.SNP = N.SNP[gn], Pvalue = gene.pval[gn], stringsAsFactors = FALSE)
    gene.pvalue[[i]] <- gene.pvalue[[i]][order(gene.pvalue[[i]]$Pvalue), ]
    rownames(gene.pvalue[[i]]) <- 1:nrow(gene.pvalue[[i]])
  }
  
  colnames(pathway.pval.stat) <- paste0('Pathway.', 1:npath)
  
  unlink(gene.stat.file)
  
  list(pathway.pval.stat = pathway.pval.stat, gene.pvalue = gene.pvalue)
  
}






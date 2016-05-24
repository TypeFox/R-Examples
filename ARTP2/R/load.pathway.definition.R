
load.pathway.definition <- function(pathway, options){
  
  msg <- paste("Loading definition of pathway:", date())
  if(options$print) message(msg)
  
  if(is.character(pathway)){
    tmp <- try(pd <- read.table(pathway, header = TRUE, as.is = TRUE), silent = TRUE)
    if(error.try(tmp)){
      msg <- paste0("Cannot load ", pathway)
      stop(msg)
    }else{
      if(nrow(pd) == 0){
        msg <- paste0("File below is empty: \n", pathway)
        stop(msg)
      }
      pathway <- pd
      rm(pd)
      gc()
    }
  }else{
    if(is.matrix(pathway)){
      pathway <- as.data.frame(pathway)
    }
  }
  
  header <- c("SNP", "Gene", "Chr")
  tmp <- (header %in% colnames(pathway))
  if(!all(tmp)){
    msg <- paste("Columns below were not found in pathway definition:\n", paste(header[!tmp], collapse = " "))
    stop(msg)
  }
  
  pathway <- pathway[, header]
  
  if(!is.null(options$selected.snps)){
    pathway <- pathway[pathway$SNP %in% options$selected.snps, , drop = FALSE]
    if(nrow(pathway) == 0){
      msg <- "No SNP is left if only use SNPs specified in options$selected.snps"
      stop(msg)
    }
  }
  
  if(!is.null(options$excluded.genes)){
    pathway <- pathway[!(pathway$Gene %in% options$excluded.genes), , drop = FALSE]
    if(nrow(pathway) == 0){
      msg <- "No SNP is left after removing genes specified by the users"
      stop(msg)
    }
  }
  
  pathway <- pathway[!duplicated(pathway), ]
  pathway <- pathway[order(pathway$Chr, pathway$Gene, pathway$SNP), ]
  
  tmp <- table(pathway$Gene, pathway$Chr)
  id <- apply(tmp, 1, function(x){sum(x > 0) > 1})
  if(any(id)){
    dup.genes <- rownames(tmp)[id]
    msg <- paste(c('The follow gene(s) are included in more than one chromosome:\n', dup.genes), collapse = ' ', sep = '')
    stop(msg)
  }
  
  rownames(pathway) <- NULL
  
  pathway
  
}

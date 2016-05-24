
load.reference.geno <- function(reference, snps.in.pathway, options){
  
  snps.in.pathway <- unique(snps.in.pathway)
  
  rt <- reference.type(reference)
  
  if(rt == 'ref.files'){
    
    msg <- paste("Loading genotypes from PLINK files:", date())
    if(options$print) message(msg)
    
    sel.subs <- options$selected.subs
    exc.subs <- options$excluded.subs
    
    if(is.null(sel.subs) && !is.null(exc.subs)){
      col.class <- rep("NULL", 6)
      col.class[2] <- "character"
      sel.subs <- read.table(reference$fam[1], header = FALSE, as.is = TRUE, colClasses = col.class)[, 1]
      sel.subs <- setdiff(sel.subs, exc.subs)
      exc.subs <- NULL
    }
    
    geno <- NULL
    for(i in 1:nrow(reference)){
      g <- read.bed(reference$bed[i], reference$bim[i], reference$fam[i], snps.in.pathway, sel.subs)
      if(!is.null(g)){
        if(is.null(geno)){
          geno <- g
        }else{
          geno <- cbind(geno, g)
        }
      }
    }
    
  }else{
    
    msg <- paste("Loading genotypes from reference:", date())
    if(options$print) message(msg)
    
    snps.in.pathway <- intersect(snps.in.pathway, colnames(reference))
    reference <- reference[, snps.in.pathway, drop = FALSE]
    
    foo1 <- function(u){
      if(any(is.na(u))){
        return(c(NA, NA))
      }
      id <- which(u %in% c('/', '-', ' ', '\\', '_'))
      if(length(id) > 0){
        u <- u[-id]
      }
      toupper(u)
    }
    
    foo2 <- function(s){
      s <- base::strsplit(toupper(s), '')
      tmp <- unique(unlist(s))
      tmp <- setdiff(tmp, c('/', '-', ' ', '\\', '_', NA))
      ea <- tmp[2]
      s <- sapply(s, foo1)
      apply(s, 2, function(u){ifelse(any(is.na(u)), NA, sum(u == ea))})
    }
    
    geno <- apply(reference, 2, foo2)
    
    colnames(geno) <- colnames(reference)
    
  }
  
  id <- which(apply(geno, 2, sd, na.rm = TRUE) == 0)
  if(length(id) > 0){
    geno <- geno[, -id, drop = FALSE]
  }
  
  geno
  
}

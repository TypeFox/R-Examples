
read.bed <- function(bed, bim, fam, sel.snps = NULL, sel.subs = NULL, encode012 = TRUE){
  
  col.class <- c("NULL", "character", "NULL", "NULL", "character", "character")
  bim.file <- read.table(bim, header = FALSE, as.is = TRUE, colClasses = col.class)
  colnames(bim.file) <- c('SNP', 'RefAllele', 'EffectAllele')
  nsnp <- nrow(bim.file)
  
  if(is.null(sel.snps)){
    sel.snps <- bim.file[, 1]
  }else{
    sel.snps <- as.character(sel.snps)
    if(any(duplicated(sel.snps))){
      msg <- 'Duplicated SNPs are detected and removed from sel.snps'
      warning(msg)
    }
    sel.snps <- unique(sel.snps)
    sel.snps <- intersect(bim.file[, 1], sel.snps)
  }
  
  sel.snp.id <- which(bim.file[, 1] %in% sel.snps)
  
  nsel <- length(sel.snp.id)
  if(nsel == 0){
    return(NULL)
  }
  
  bim.file <- bim.file[sel.snp.id, , drop = FALSE]
  sel.snps <- bim.file[, 1]
  
  col.class <- rep("NULL", 6)
  col.class[2] <- "character"
  sid <- read.table(fam, header = FALSE, as.is = TRUE, colClasses = col.class)[, 1]
  if(any(duplicated(sid))){
    msg <- paste0('Duplicated subjects exist in fam file: \n', fam)
    warning(msg)
  }
  nsub <- length(sid)
  
  geno <- rep(-1, nsub * nsel)
  
  tmp <- .C("ReadBED", as.character(bed), as.integer(nsub), 
            as.integer(nsnp), as.integer(nsel), as.integer(sel.snp.id), 
            geno = as.integer(geno), PACKAGE = "ARTP2")
  
  geno <- as.data.frame(matrix(tmp$geno, nrow = nsub, byrow = FALSE))
  rownames(geno) <- sid
  colnames(geno) <- sel.snps
  
  if(!is.null(sel.subs)){
    if(any(duplicated(sel.subs))){
      msg <- 'Duplicated subjects exist in sel.subs and are returnsed as duplicated lines'
      warning(msg)
    }
    sel.subs <- as.character(sel.subs)
    id <- which(sel.subs %in% rownames(geno))
    if(length(id) == 0){
      msg <- paste("No subjects were left in \n", bed)
      stop(msg)
    }
    sel.subs <- sel.subs[id]
    geno <- geno[sel.subs, , drop = FALSE]
  }
  
  geno[geno == -1] <- NA
  
  if(encode012){
    return(geno)
  }
  
  geno <- geno[, bim.file$SNP, drop = FALSE]
  
  for(i in 1:ncol(geno)){
    rs <- bim.file$SNP[i]
    g <- geno[, rs]
    ra <- bim.file$RefAllele[i]
    ea <- bim.file$EffectAllele[i]
    code <- paste0(c(ra, ra, ea), c(ra, ea, ea))
    geno[, rs] <- ifelse(is.na(g), NA, code[g + 1])
  }
  
  geno
  
}

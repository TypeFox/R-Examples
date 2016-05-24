
load.reference.allele <- function(reference, pathway, options){
  
  rt <- reference.type(reference)
  
  snps.in.pathway <- unique(pathway$SNP)
  
  if(rt == 'ref.files'){
    
    msg <- paste("Loading allele information from PLINK files:", date())
    if(options$print) message(msg)
    
    if("matrix" %in% class(reference)){
      reference <- as.data.frame(reference)
    }
    
    col.class <- rep("NULL", 6)
    col.class[c(2, 5, 6)] <- "character"
    col.class[c(1, 4)] <- 'integer'
    bim.files <- reference$bim
    nfiles <- length(bim.files)
    allele.info <- NULL
    for(i in 1:nfiles){
      tmp <- try(bim <- read.table(bim.files[i], header = FALSE, as.is = TRUE, colClasses = col.class), silent = TRUE)
      if(error.try(tmp)){
        msg <- paste0('Cannot load ', bim.files[i])
        stop(msg)
      }
      
      colnames(bim) <- c("Chr", "SNP", "Pos", "RefAllele", "EffectAllele")
      bim$Reference.ID <- i
      bim <- bim[bim$SNP %in% snps.in.pathway, ]
      bim$RefAllele <- toupper(bim$RefAllele)
      bim$EffectAllele <- toupper(bim$EffectAllele)
      allele.info <- rbind(allele.info, bim)
      
      rm(bim)
      gc()
    }
    
    if(is.null(allele.info)){
      msg <- "No SNPs are found in bim files"
      stop(msg)
    }
    
  }else{
    
    msg <- paste("Loading allele information from reference:", date())
    if(options$print) message(msg)
    
    snps.in.pathway <- intersect(snps.in.pathway, colnames(reference))
    if(length(snps.in.pathway) == 0){
      msg <- 'No SNPs are found in reference genotypes'
      stop(msg)
    }
    
    reference <- reference[, snps.in.pathway, drop = FALSE]
    
    foo1 <- function(s){
      s <- unique(unlist(base::strsplit(toupper(s), '')))
      s <- setdiff(s, c('/', '-', ' ', '\\', '_', NA))
      length(s)
    }
    
    id <- which(apply(reference, 2, foo1) == 2)
    if(length(id) == 0){
      msg <- 'All SNPs in reference genotypes are excluded'
      stop(msg)
    }
    
    reference <- reference[, id, drop = FALSE]
    
    foo2 <- function(s){
      s <- unique(unlist(base::strsplit(toupper(s), '')))
      s <- setdiff(s, c('/', '-', ' ', '\\', '_', NA))
      s
    }
    
    tmp <- apply(reference, 2, foo2)
    SNP <- colnames(reference)
    RefAllele <- as.vector(tmp[1, ])
    EffectAllele <- as.vector(tmp[2, ])
    allele.info <- data.frame(SNP = SNP, Pos = NA, RefAllele = RefAllele, EffectAllele = EffectAllele, stringsAsFactors = FALSE)
    path <- pathway[!duplicated(pathway$SNP), c('SNP', 'Chr')]
    allele.info <- merge(allele.info, path, by = 'SNP')
    allele.info <- allele.info[, c("Chr", "SNP", "Pos", "RefAllele", "EffectAllele")]
    
  }
  
  allele.info
  
}

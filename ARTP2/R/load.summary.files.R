
load.summary.files <- function(summary.files, lambda, sel.snps){
  
  msg <- paste("Loading summary files:", date())
  message(msg)
  
  header <- c('SNP', 'RefAllele', 'EffectAllele', 'BETA') # columns that must be provided by users
  opt.header <- c('P', 'SE')
  
  complete.header <- c(header, opt.header, 'Direction')
  
  nfiles <- length(summary.files)
  stat <- list()
  lam <- NULL
  
  fid <- 0
  for(i in 1:nfiles){
    st <- read.table(summary.files[i], header = TRUE, as.is = TRUE, nrows = 1e4)
    colnames(st) <- convert.header(colnames(st), complete.header)
    tmp <- (header %in% colnames(st))
    if(!all(tmp)){
      msg <- paste0("Columns below were not found in ", summary.files[i], ":\n", paste(header[!tmp], collapse = " "))
      stop(msg)
    }
    
    col.class <- sapply(st, class)
    col.id <- which(colnames(st) %in% complete.header)
    col.class[-col.id] <- "NULL"
    col.class[c('SNP', 'RefAllele', 'EffectAllele')] <- 'character'
    st <- read.table(summary.files[i], header = TRUE, as.is = TRUE, colClasses = col.class)
    colnames(st) <- convert.header(colnames(st), complete.header)
    if(!is.null(sel.snps)){
      st <- st[st$SNP %in% sel.snps, ]
    }
    
    if(nrow(st) == 0){
      next
    }
    
    if(!any(opt.header %in% colnames(st))){
      msg <- paste0("Neither SE nor P is not provided in ", summary.files[i])
      stop(msg)
    }
    
    
    if(!('P' %in% colnames(st))){
      st$P <- NA
    }
    
    if(!('SE' %in% colnames(st))){
      st$SE <- NA
    }
    
    if(!('Direction' %in% colnames(st))){
      msg <- paste0('Direction is absent in ', summary.files[i], '. Function meta() assumed equal sample sizes for all SNPs in that study. Invalidation of this assumption can lead to false positive if summary data of this study is used in pathway analysis')
      warning(msg)
      st$Direction <- ifelse(st$BETA == 0, '0', ifelse(st$BETA > 0, '+', '-'))
    }
    
    nc <- unique(nchar(st$Direction))
    if(length(nc) != 1){
      msg <- paste0('String lengths of Direction are unequal in ', summary.files[i])
      stop(msg)
    }
    
    st <- st[, complete.header]
    
    dup <- duplicated(st$SNP)
    if(any(dup)){
      dup.snps <- unique(st$SNP[dup])
      msg <- paste("SNPs below are duplicated: ", paste(dup.snps, collapse = " "))
      stop(msg)
    }
    
    id.no.SE.P <- which(is.na(st$SE) & is.na(st$P))
    if(length(id.no.SE.P) > 0){
      msg <- paste("For SNPs below, neither SE nor P is not provided in", summary.files[i], ":\n", paste(st$SNP[id.no.SE.P], collapse = " "))
      stop(msg)
    }
    
    st$RefAllele <- toupper(st$RefAllele)
    st$EffectAllele <- toupper(st$EffectAllele)
    
    id.no.SE <- which(is.na(st$SE))
    id.no.P <- which(is.na(st$P))
    
    if(length(id.no.SE) > 0){
      z2 <- qchisq(st$P[id.no.SE], df = 1, lower.tail = FALSE)
      st$SE[id.no.SE] <- abs(st$BETA[id.no.SE]/sqrt(z2))
    }
    
    if(length(id.no.P) > 0){
      st$P[id.no.P] <- pchisq((st$BETA[id.no.P]/st$SE[id.no.P])^2, df = 1, lower.tail = FALSE)
    }
    
    fid <- fid + 1
    lam <- c(lam, lambda[i])
    rownames(st) <- st$SNP
    st <- st[complete.cases(st), ]
    
    stat[[fid]] <- st
    rm(st)
    gc()
    
  }
  
  if(length(stat) == 0){
    msg <- "No SNPs to be included in analysis"
    stop(msg)
  }
  
  lambda <- lam
  
  list(stat = stat, lambda = lambda)
  
}


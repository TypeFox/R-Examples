map.SNPs.to.plinkfiles <- function(plink.files, pathway){
  
  bim <- plink.files$bim
  cc <- c('NULL', 'character', 'NULL', 'NULL', 'NULL', 'NULL')
  
  nf <- length(bim)
  SNP <- NULL
  iplink <- NULL
  for(i in 1:nf){
    re <- try(b <- read.table(bim[i], header = FALSE, as.is = TRUE, colClasses = cc)[, 1], silent = TRUE)
    if("try-error" %in% class(re)){
      msg <- paste0('Cannot load ', b)
      stop(msg)
    }
    
    b <- intersect(b, pathway$SNP)
    if(length(b) == 0){
      next
    }
    SNP <- c(SNP, b)
    iplink <- c(iplink, rep(i, length(b)))
    
  }
  
  if(length(SNP) == 0){
    msg <- 'No SNP is available in geno.files'
    stop(msg)
  }
  
  sf <- data.frame(SNP, iplink, stringsAsFactors = FALSE)
  sf <- sf[!duplicated(sf$SNP), ]
  rownames(sf) <- NULL
  
  sf
  
}


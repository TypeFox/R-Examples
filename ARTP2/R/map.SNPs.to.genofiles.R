
# read in the first row of each geno.files to obtain a list of available SNPs. Map each SNP to one geno.files
map.SNPs.to.genofiles <- function(geno.files, pathway){
  
  SNP <- NULL
  file <- NULL
  col <- NULL
  ncol <- NULL
  for(f in geno.files){
    rs <- scan(f, what='character',nlines = 1, quiet = TRUE)
    rs <- intersect(rs, pathway$SNP)
    nrs <- length(rs)
    
    if(nrs == 0){
      next
    }
    
    SNP <- c(SNP, rs)
    file <- c(file, rep(f, nrs))
    col <- c(col, 1:nrs)
    ncol <- c(ncol, rep(nrs, nrs))
  }
  
  if(length(SNP) == 0){
    msg <- 'No SNP is available in geno.files'
    stop(msg)
  }
  
  sf <- data.frame(SNP, file, col, ncol, stringsAsFactors = FALSE)
  sf <- sf[!duplicated(sf$SNP), ]
  rownames(sf) <- NULL
  
  sf
  
}


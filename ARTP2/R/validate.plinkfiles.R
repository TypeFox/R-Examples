
validate.plinkfiles <- function(geno.files){
  
  if(!is.data.frame(geno.files) && !is.matrix(geno.files)){
    msg <- 'geno.files should be a data frame'
    stop(msg)
  }
  
  if(is.matrix(geno.files)){
    geno.files <- as.data.frame(geno.files)
  }
  
  if(!all(c('fam', 'bim', 'bed') %in% colnames(geno.files))){
    msg <- 'geno.files should be a data frame with columns \'fam\', \'bim\', \'bed\''
    stop(msg)
  }
  
  geno.files <- geno.files[, c('fam', 'bim', 'bed')]
  
  fam <- geno.files$fam
  cc <- rep('NULL', 6)
  cc[1] <- 'character'
  for(f in fam){
    re <- try(tmp <- read.table(f, header = FALSE, as.is = TRUE, colClasses = cc, nrows = 1e3))
    if("try-error" %in% class(re)){
      msg <- paste0('Cannot load ', f)
      stop(msg)
    }
    rm(tmp)
    gc()
  }
  
  bim <- geno.files$bim
  cc <- rep('NULL', 6)
  cc[2] <- 'character'
  for(b in bim){
    re <- try(tmp <- read.table(b, header = FALSE, as.is = TRUE, colClasses = cc, nrows = 1e3))
    if("try-error" %in% class(re)){
      msg <- paste0('Cannot load ', b)
      stop(msg)
    }
    rm(tmp)
    gc()
  }
  
  
  geno.files
  
}

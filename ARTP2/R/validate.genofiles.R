
# check if geno.files can be read or not
validate.genofiles <- function(geno.files){
  
  miss.files <- NULL
  for(f in geno.files){
    re <- try(tmp <- read.table(f, header = TRUE, as.is = TRUE, nrows = 1), silent = TRUE)
    if("try-error" %in% class(re)){
      miss.files <- c(miss.files, f)
    }
  }
  
  if(setequal(miss.files, geno.files)){
    msg <- 'All genotype files are missing'
    stop(msg)
  }
  
  if(!is.null(miss.files)){
    msg <- paste(c('The following genotype files are not found: ', miss.files), sep = '', collapse = '\n')
    warning(msg)
  }
  
  geno.files <- setdiff(geno.files, miss.files)
  geno.files
  
}


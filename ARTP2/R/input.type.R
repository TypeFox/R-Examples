
input.type <- function(geno.files){
  
  if(is.null(geno.files)){
    return('data.frame')
  }
  
  if('character' %in% class(geno.files)){
    return('geno.files')
  }
  
  if(!('data.frame' %in% class(geno.files))){
    msg <- 'Invalid geno.files'
    stop(msg)
  }
  
  if(any(c('fam', 'bim', 'bed') %in% colnames(geno.files))){
    if(all(c('fam', 'bim', 'bed') %in% colnames(geno.files))){
      return('plink.files')
    }
    
    msg <- 'Invalid geno.files. Is your genotype data stored in plink files?'
    stop(msg)
  }
  
}


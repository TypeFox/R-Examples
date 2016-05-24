
reformat.reference.path <- function(reference){
  
  reference <- as.data.frame(reference)
  
  if(reference.type(reference) == 'ref.geno'){
    return(reference)
  }
  
  if("bed" %in% colnames(reference)){
    reference$bed <- as.character(reference$bed)
  }
  
  if("bim" %in% colnames(reference)){
    reference$bim <- as.character(reference$bim)
  }
  
  if("fam" %in% colnames(reference)){
    reference$fam <- as.character(reference$fam)
  }
  
  reference
  
}


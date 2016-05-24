
reference.type <- function(reference){
  
  header <- c("bed", "bim", "fam")
  tmp <- (header %in% colnames(reference))
  if(is.data.frame(reference) && !any(tmp)){
    return('ref.geno')
  }else{
    return('ref.files')
  }
  
}


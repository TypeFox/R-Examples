
parsePARline <- function(s){
  
  if(length(s) > 1)s <- paste(s, collapse= " ")
  
  sp <- strsplit(s, "=")[[1]]
  if(length(sp) != 2)stop("Fatal error in parsePAR.")
  
  parname <- str_trim(sp[1])
  parval <- str_trim(sp[2])
  
  # Single numeric
  tp <- trynumeric(parval)
  if(is.numeric(tp))return(tp)
  
  # Not returned yet - try splitting
  spl <- strsplit(parval, "[[:space:]]")[[1]]
  
  # Single character
  if(length(spl) == 1){
    v <- gsub("'","", parval)
    return(v) 
  }
  
  # Multiple numeric or character
  if(!is.numeric(trynumeric(spl))){
    parval <- gsub("'","", spl)
  } else {
    parval <- trynumeric(spl)
  }
  
  return(parval)
}

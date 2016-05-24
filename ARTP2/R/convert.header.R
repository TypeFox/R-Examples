
convert.header <- function(header1, header2){
  
  names(header1) <- toupper(header1)
  names(header2) <- toupper(header2)
  
  hd1 <- names(header1)
  hd2 <- names(header2)
  if(any(duplicated(hd1))){
    return(header1)
  }
  
  id <- which(hd1 %in% hd2)
  if(length(id) == 0){
    return(header1)
  }
  
  hd <- hd1[id]
  header1[id] <- header2[hd]
  header1
  
}

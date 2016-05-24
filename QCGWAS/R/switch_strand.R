switch_strand <-
function(input, strand_col = FALSE ) {
  output <- input[ , (1:3)[c(TRUE,TRUE,strand_col)] ]
  
  output[ ,1] <- ifelse(input[ ,1] == "A", "T", output[ ,1])
  output[ ,1] <- ifelse(input[ ,1] == "T", "A", output[ ,1])
  output[ ,1] <- ifelse(input[ ,1] == "C", "G", output[ ,1])
  output[ ,1] <- ifelse(input[ ,1] == "G", "C", output[ ,1])
  
  output[ ,2] <- ifelse(input[ ,2] == "A", "T", output[ ,2])
  output[ ,2] <- ifelse(input[ ,2] == "T", "A", output[ ,2])
  output[ ,2] <- ifelse(input[ ,2] == "C", "G", output[ ,2])
  output[ ,2] <- ifelse(input[ ,2] == "G", "C", output[ ,2])
  
  if(strand_col) {
    output[ ,3] <- ifelse(input[ ,3] == "+", "-", output[ ,3])
    output[ ,3] <- ifelse(input[ ,3] == "-", "+", output[ ,3])
  }
  return(output)
}

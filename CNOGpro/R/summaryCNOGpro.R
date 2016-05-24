summaryCNOGpro <-
function(experiment){
  cat("CNOGpro copy number analysis experiment of", experiment$Name, "\n")
  cat("Windowlength:", experiment$windowlength,"\n\n")
  cat("Summary statistics of chromosome:", experiment$accession,"\n")
  cat("Length: ", experiment$chrlength,"\tNumber of genetic elements: ", nrow(experiment$genes),"\tMean coverage: ",experiment$mean,"\tVariance: ",experiment$variance,"\n")
  if(experiment$is_GC_normalized) cat("Data have been normalized by GC content\n") 
  else cat("Data not yet normalized\n")  
}

printCNOGpro <-
function(experiment){
  cat("CNOGpro results for ", experiment$Name,". Chromosome: ",experiment$accession,"\n",sep="")
  if(!is.null(experiment$HMMtable)){print(experiment$HMMtable)}
  else if(!is.null(experiment$genes)){print(experiment$genes)}
  else{cat("You have not yet performed copy number analysis.\n")}


}

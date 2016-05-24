exclude.seqsite <-
function(sequences, site, verbose=TRUE){  
  sequences <- as.character(sequences)
  res.exc <- sequences[na.action(na.omit(sequences[vwhichPDict(PDict(site), DNAStringSet(sequences)) == 1]))]
  if(verbose == TRUE){
  cat(length(sequences) - length(res.exc), 
      " out of ", length(sequences), " sequences with the restriction site detected and removed \n",
      length(res.exc), 
      " sequence remaining \n", sep="")
  }
  return(res.exc)  
}

prob = function(combinations, freqs, rare, threshold, theta = 0){
  bIsAll = any(grepl("all", combinations, ignore.case = TRUE))
  
  if(bIsAll & length(combinations) > 1){
    warning("All option used, all other combinations will be ignored")
    combinations = "all";
  }
  
  .Call("DNAtools_Prob", combinations, freqs, rare, threshold, theta)
}
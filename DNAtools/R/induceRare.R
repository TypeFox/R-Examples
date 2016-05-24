induceRare = function(db, freqs, threshold = 0.01, rareCode = 99){
  
  nLoci = length(freqs)
  rareSet = sapply(freqs, function(x){which(x <= threshold)})
  
  for(loc in 1:nLoci){
    j1 = 2*loc
    j2 = j1 + 1
    
    ## order alleles
    idx = db[ ,j1] > db[ ,j2]
    if(any(idx)){
      tmp = db[idx, j1]
      db[idx, j1] = db[idx, j2]
      db[idx, j2] = tmp
    }
    
    i1 = db[, j1] %in% rareSet[[loc]]
    
    if(any(i1))
      db[i1,j1] = 99
    
    i2 = db[, j2] %in% rareSet[[loc]]
  
    if(any(i2))
      db[i2,j2] = 99
  }
  return(db)
}
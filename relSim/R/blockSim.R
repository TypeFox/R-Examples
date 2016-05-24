blockSim = function(N, Freqs, rel = "UN", ibsthresh = NULL, kithresh = NULL,
                    code = 1, falseNeg = TRUE, BlockSize = N/10){
  
  rel = toupper(rel)
  if(!grepl("(UN|FS|PC)", rel)){
    stop("rel must be one of 'UN', 'FS' or 'PC'")
  }
  
  nBlocks = N / BlockSize
  
  if(is.null(ibsthresh) & is.null(kithresh))
    stop("You must specify one or both of ibsthresh or kithresh")
  
  nResults = 0
  if(is.null(ibsthresh) & !is.null(kithresh)){
    nResults = length(kithresh)
    ibsthresh = rep(0, nResults) ## dummy vals
  }else if(is.null(kithresh) & !is.null(ibsthresh)){
    nResults = length(ibsthresh)
    kithresh = rep(0, nResults)
  }else{
    if(length(ibsthresh) != length(kithresh)){
      stop("ibsthresh and kithresh must be the same length")
    }else{
      nResults = length(ibsthresh)
    }
  }
  
  if(nResults == 0)
    stop("Nothing to count")
  
  nTotal = rep(0, nResults)
  
  pb = txtProgressBar(min = 0, max = nBlocks, style = 3)
  
  if(rel == "UN"){
    for(block in 1:nBlocks){
      Prof1 = randomProfiles(Freqs$freqs, BlockSize)
      Prof2 = randomProfiles(Freqs$freqs, BlockSize)
      
      count = blockStatCounts(Prof1, Prof2, BlockSize, Freqs$freqs, 
                              code, falseNeg, ibsthresh, kithresh, nResults)
      nTotal = nTotal + count
      setTxtProgressBar(pb, block)
    }
  }else if(rel == "FS"){
    for(block in 1:nBlocks){
      Prof1 = randomProfiles(Freqs$freqs, BlockSize)
      Prof2 = randomSibs(Prof1, Freqs$freqs, BlockSize)
      
      count = blockStatCounts(Prof1, Prof2, BlockSize, Freqs$freqs, 
                              code, falseNeg, ibsthresh, kithresh, nResults)
      
      nTotal = nTotal + count
      setTxtProgressBar(pb, block)
    }
  }else if(rel == "PC"){
    for(block in 1:nBlocks){
      Prof1 = randomProfiles(Freqs$freqs, BlockSize)
      Prof2 = randomChildren(Prof1, Freqs$freqs, BlockSize)
      
      count = blockStatCounts(Prof1, Prof2, BlockSize, Freqs$freqs, 
                              code, falseNeg, ibsthresh, kithresh, nResults)
      
      nTotal = nTotal + count
      setTxtProgressBar(pb, block)
    }
  }
  
  return(list(nTotal = nTotal, N = N, p = nTotal / N))
}

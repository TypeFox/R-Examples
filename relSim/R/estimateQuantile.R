estimateQuantile = function(p, Freqs, n, N, rel = "UN", w = 3,
                            BlockSize = N/10, resFile = NULL, ...){
  rel = toupper(rel)
  if(!grepl("(UN|FS|PC)", rel)){
    stop("rel must be one of 'UN', 'FS', or 'PC'")
  }
  
  v = p*(1-p)/n
  pL = max(0, p - w*sqrt(v))
  pU = min(p + w*sqrt(v), 1)
  
  if(n < (1/10^floor(log10(pU)))){
    stop("The sample size is insufficient to estimate the quantiles")
  }
  
  nLoci = length(Freqs$freqs)
  
  if(rel == "UN"){
    unrel = NULL
    if(is.null(resFile)){
      cat("Taking initial sample\n")
      unrel = sim(n, Freqs, ...)
      cat("\n")
    }else{
      cat(paste("Reading results from", resFile, "\n"))
      unrel = readResults(fileName = resFile)
    }
    
    qSib = quantile(unrel$sib, c(pL, pU), type = 8)
    qPC = quantile(unrel$pc, c(pL, pU), type = 8)
    
    if(qSib[2] - qSib[1] <= 0 | qPC[2] - qPC[1] <= 0){
      cat(sprintf("CI for sibs is [%6.1E, %6.1E]\n", qSib[1], qSib[2]))
      cat(sprintf("CI for PC is [%6.1E, %6.1E]\n", qPC[1], qPC[2]))
      stop("The sample size is insufficient")
    }
    
    cat("Taking secondary sample\n")
    cat(sprintf("CI for sibs is [%6.1E, %6.1E]\n", qSib[1], qSib[2]))
    cat(sprintf("CI for PC is [%6.1E, %6.1E]\n", qPC[1], qPC[2]))
    
    keep = vector(mode = "list", length = 2)
    names(keep) = c("pc","sib")
    
    nBlocks = N/BlockSize
    
    if(floor(nBlocks) - nBlocks < 0){
      stop("BlockSize must evenly divide N")
    }
    
    code = 4 ## lrPC and lrSib
    nSib = 0
    nPC = 0
    
    pb = txtProgressBar(min = 0, max = nBlocks, style = 3)
    
    for(b in 1:nBlocks){
      prof1 = randomProfiles(Freqs$freqs, BlockSize)
      prof2 = randomProfiles(Freqs$freqs, BlockSize) 
      
      blockRes = blockStats(prof1, prof2, BlockSize, Freqs$freqs, code)
      
      sibKeep = blockRes$sib[blockRes$sib >= qSib[1] &
                               blockRes$sib <= qSib[2]]
      
      pcKeep = blockRes$pc[blockRes$pc >= qPC[1] &
                             blockRes$pc <= qPC[2]]
      
      keep$pc = c(keep$pc, pcKeep)
      keep$sib = c(keep$sib, sibKeep)
      
      nSib = nSib + sum(blockRes$sib < qSib[1])
      nPC = nPC + sum(blockRes$pc < qPC[1])
      setTxtProgressBar(pb, b)
    }
    
    i = which.min(abs(p - 1:N/(N+1)))[1] ## [1] just to get rid of dupes
    qS = sort(keep$sib)[i - nSib]
    qP = sort(keep$pc)[i - nPC]
    
    return(list(pc = qP, sib = qS))
  }else if(rel == "FS"){
    sibs = NULL
    
    if(is.null(resFile)){
      cat("Taking initial sample\n")
      sibs = sim(n, Freqs, "FS",  ...)
    }else{
      cat(paste("Reading results from", resFile, "\n"))
      sibs = readResults(fileName = resFile)
    }
    
    qSib = quantile(sibs$sib, c(pL, pU), type = 8)
    
    if(qSib[2] - qSib[1] <= 0){
      stop("The sample size is insufficient")
    }
    
    cat("Taking secondary sample\n")
    cat(sprintf("CI for sibs is [%6.1E, %6.1E]\n", qSib[1], qSib[2]))
    
    keep = NULL
    nBlocks = N / BlockSize
    
    if(floor(nBlocks) - nBlocks < 0){
      stop("BlockSize must evenly divide N")
    }
    
    code = 1 ## lrSib
    nSib = 0
    
    pb = txtProgressBar(min = 0, max = nBlocks, style = 3)
    
    for(b in 1:nBlocks){
      sib1 = randomProfiles(Freqs$freqs, BlockSize)
      sib2 = randomSibs(sib1, Freqs$freqs, BlockSize)
      
      blockRes =  blockStats(sib1, sib2, BlockSize, Freqs$freqs, code)
      
      sibKeep = blockRes$sib[blockRes$sib >= qSib[1] &
                               blockRes$sib <= qSib[2]]
      
      keep = c(keep, sibKeep)
      nSib = nSib + sum(blockRes$sib < qSib[1])
      setTxtProgressBar(pb, b)
    }
    
    i = which.min(abs(p - 1:N/(N+1)))[1] ## [1] just to get rid of dupes
    qS = sort(keep)[i - nSib]
    
    return(list(sib = qS))
    
  }else if(rel == "PC"){
    pc = NULL
    
    if(is.null(resFile)){
      cat("Taking initial sample\n")
      pc = sim(n, Freqs, "PC",  ...)
    }else{
      cat(paste("Reading results from", resFile, "\n"))
      pc = readResults(fileName = resFile)
    }
    
    qPC = quantile(pc$pc, c(pL, pU), type = 8)
    
    if(qPC[2] - qPC[1] <= 0){
      stop("The sample size is insufficient")
    }
    
    cat("Taking secondary sample\n")
    cat(sprintf("CI for PC is [%6.1E, %6.1E]\n", qPC[1], qPC[2]))
    
    keep = NULL
    nBlocks = N / BlockSize
    
    if(floor(nBlocks) - nBlocks < 0){
      stop("BlockSize must evenly divide N")
    }
    
    code = 2 ## lrPC
    nSib = 0
    nPC = 0

    pb = txtProgressBar(min = 0, max = nBlocks, style = 3)
    
    for(b in 1:nBlocks){
      parent = randomProfiles(Freqs$freqs, BlockSize)
      child = randomProfiles(Freqs$freqs, BlockSize)
      
      
      
      blockRes = blockStats(parent, child, BlockSize, Freqs$freqs, code)
      
      pcKeep = blockRes$pc[blockRes$pc >= qPC[1] &
                             blockRes$pc <= qPC[2]]
      
      keep = c(keep, pcKeep)
      nPC = nPC + sum(blockRes$pc < qPC[1])
      setTxtProgressBar(pb, b)
    }
    
    i = which.min(abs(p - 1:N/(N+1)))[1] ## [1] just to get rid of dupes
    qP = sort(keep)[i - nPC]
    
    return(list(pc = qP))
  }
}



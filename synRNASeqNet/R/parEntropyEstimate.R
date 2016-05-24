parEntropyEstimate <-
function(idx, method = method, unit = unit,
                               priorHyperParam = priorHyperParam, 
                               shrinkageTarget = shrinkageTarget, boot = boot){
  getIndex <- function(idx){
    rrow <- ceiling((-1 + sqrt(8*idx + 1))/2)
    ccol <- idx - rrow*(rrow - 1)/2
    ans <- c(rrow, ccol)
    return(ans)
  }
  
  idx <- getIndex(idx)
  if(boot == F){
    if(idx[1] == idx[2]) cellCounts <- table(counts[idx[1], ]) else
      cellCounts <- table(counts[idx[1], ], counts[idx[2], ])
  }else{
    cellCounts <- table(counts[idx[1], ], counts[idx[2], sample(1:ncol(counts))])
  }
  
  if(method == "ML") ans <- entropyML(cellCounts, unit = unit) else
    if(method == "MM") ans <- entropyMM(cellCounts, unit = unit) else
      if(method == "Bayes"){
        #priorHyperParam <- match.arg(priorHyperParam)
        ans <- entropyBayes(cellCounts, unit = unit,
                            priorHyperParam = priorHyperParam)
      } else
        if(method == "CS") ans <- entropyCS(cellCounts, unit = unit) else
          if(method == "Shrink") ans <- entropyShrink(cellCounts, unit = unit,
                                                      shrinkageTarget = shrinkageTarget) else
                                                        stop("Unknown Entropy Estimate Method")
  
  return(ans)
}

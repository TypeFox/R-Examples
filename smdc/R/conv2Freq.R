conv2Freq <-
function(tmpMatrix,wordClass,breaks) {
  # 度数分布表に変換する
  freqDist <- matrix(0,nrow=length(breaks)-1,ncol=ncol(tmpMatrix))
  for (tmp in rownames(tmpMatrix)) {
    cat <- wordClass[tmp]
    if (!is.na(cat)) {
      freqDist[cat,] <- freqDist[cat,] + tmpMatrix[tmp,]
    }
  }
  colnames(freqDist) <- colnames(tmpMatrix)
  if (!is.null(names(breaks))) {
    rownames(freqDist) <- names(breaks)[2:length(breaks)]    
  }
  return(freqDist)
}

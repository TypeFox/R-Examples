simDic <-
function(docMatrix1,docMatrix2,scoreDict,breaks=seq(-1,1,length=11),norm=FALSE,method='cosine',scoreFunc=mean) {
  # ライブラリ読込み
  library('proxy')
  
  # 出現する単語をベクトル化する
  words <- unique(rbind(matrix(rownames(docMatrix1)),matrix(rownames(docMatrix2))))
  words <- words[order(words)] #ソートしておく
  
  # 単語ごとにスコアを算出する
  wordScores <- rep(NA,length(words))
  for (i in 1:length(words)) {
    cond <- (scoreDict[,1] == words[i])
    value <- scoreDict[cond,2]
    if (length(value)!=0) {
      wordScores[i] <- scoreFunc(value,na.rm=TRUE)
    }
  }
  
  # 単語をスコアごとのクラスに割り振る
  names(breaks) <- cut(breaks,breaks)
  wordClass <- cut(wordScores,breaks)
  names(wordClass) <- words
  
  # 入力行列を度数分布表に変換する
  docFreq1 <- conv2Freq(docMatrix1,wordClass,breaks)
  docFreq2 <- conv2Freq(docMatrix2,wordClass,breaks)
  
  # docMatrix1とdocMatrix2の列名が同一にならないようにしておく
  colnames(docFreq1) <- paste("r_",colnames(docMatrix1),sep="")
  colnames(docFreq2) <- paste("c_",colnames(docMatrix2),sep="")
  
  # 類似度を算出する
  sim <- as.matrix(simil(t(cbind(docFreq1,docFreq2)),method=method))[colnames(docFreq1),colnames(docFreq2)]
  rownames(sim) <- colnames(docMatrix1)
  colnames(sim) <- colnames(docMatrix2)
  
  # 正規化
  if (norm) {
    sim <- normalize(sim)
  }
  
  return(sim)
}

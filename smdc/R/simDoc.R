simDoc <-
function(docMatrix1,docMatrix2,norm=FALSE,method='cosine') {
  # ライブラリ読込み
  library('proxy')
  
  # 行列の行数を一致させる
  exDocMatrix <- uniform(docMatrix1,docMatrix2)
  exDocMatrix1 <- exDocMatrix[[1]]
  exDocMatrix2 <- exDocMatrix[[2]]

  # docMatrix1とdocMatrix2の列名が同一にならないようにしておく
  colnames(exDocMatrix1) <- paste("r_",colnames(docMatrix1),sep="")
  colnames(exDocMatrix2) <- paste("c_",colnames(docMatrix2),sep="")
  
  # 類似度行列を計算する
  sim <- as.matrix(simil(t(cbind(exDocMatrix1,exDocMatrix2)),method=method))[colnames(exDocMatrix1),colnames(exDocMatrix2)]
  rownames(sim) <- colnames(docMatrix1)
  colnames(sim) <- colnames(docMatrix2)
  
  # 正規化
  if (norm) {
    sim <- normalize(sim)
  }
  
  return(sim)
}

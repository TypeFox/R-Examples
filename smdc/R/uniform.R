uniform <-
function(matrix1,matrix2) {
  # 出現する特徴量をベクトル化する
  words <- unique(rbind(matrix(rownames(matrix1)),matrix(rownames(matrix2))))
  words <- words[order(words)] #ソートしておく
  
  # 行列の行数を揃える
  exMatrix1 <- matrix(0,nrow=length(words),ncol=ncol(matrix1))
  exMatrix2 <- matrix(0,nrow=length(words),ncol=ncol(matrix2))
  rownames(exMatrix1) <- words
  rownames(exMatrix2) <- words
  colnames(exMatrix1) <- colnames(matrix1)
  colnames(exMatrix2) <- colnames(matrix2)
  for (word in rownames(matrix1)) {
    exMatrix1[word,] <- matrix1[word,]
  }
  for (word in rownames(matrix2)) {
    exMatrix2[word,] <- matrix2[word,]
  }
  
  return(list(exMatrix1,exMatrix2))
}

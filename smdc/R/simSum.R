simSum <-
function(sim) {
  # ドキュメントごとに、最高スコアと類似度のもっとも高いクラスを算出する
  results <- rep(0,ncol(sim))
  names(results) <- colnames(sim)
  scores <- rep(0,ncol(sim))
  for (i in 1:ncol(sim)) {
    scores[i] <- max(sim[,i])
    results[i] <- rownames(sim)[which.max(sim[,i])]
  }
  
  # クラスごとに、類似度の高い順にドキュメントを集計する
  summary <- as.list(NULL,length=nrow(sim))
  for (i in 1:nrow(sim)) {
    cond <- results==rownames(sim)[i]
    summary[[i]] <- names(which(cond[order(-scores)]))
  }
  names(summary) <- rownames(sim)
  
  return(summary)
}

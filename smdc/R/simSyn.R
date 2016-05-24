simSyn <-
function(sims,weight) {  
  # 入力ベクトルの長さをチェック
  len = length(sims)
  if (len!=length(weight)) {
    stop(message="different lengths between sims and weight")
  }
  
  # 類似度を合成する
  sim <- weight[1] * sims[[1]]
  for (i in 2:len) {
    sim <- sim + weight[i] * sims[[i]]
  }
  
  return(sim)
}

normalize <-
function(sim) {
  # 各行ごとに平均0, 標準偏差1に正規化する
  meanVec <- apply(sim,1,mean,na.rm = TRUE)
  sdVec <- apply(sim,1,sd,na.rm = TRUE)
  sim <- t(scale(t(sim),meanVec,sdVec))
  return(sim)
}

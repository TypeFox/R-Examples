makeWeightsAIC <-
function(aics){
  min.i<-which(aics==min(aics))
  deltas<-aics-aics[min.i][1]
  w<-exp(-deltas/2)/sum(exp(-deltas/2))
  return(w)
}

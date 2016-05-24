portvol.Bayes <- function(tickers, weights = rep(1,length(tickers)),
                          start, end, data, sim.size = 1000){
  cctr.sim <- cctr.Bayes(tickers, weights, start, end, data, sim.size)
  portvol.sim <- rowSums(cctr.sim)
  return(portvol.sim)
}
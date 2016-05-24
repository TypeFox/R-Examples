cctr.Bayes <- function(tickers, weights = rep(1,length(tickers)),
                       start, end, data, sim.size = 1000){
  w <- weights/sum(weights)
  cctr.sim <- w*mctr.Bayes(tickers, w, start, end, data, sim.size)
  return(cctr.sim)
}
#########################################
#                                       # 
# RANK HISTOGRAM FOR ENSEMBLE FORECASTS #
#                                       #
#########################################
Rankhist <- function(ens, obs, reduce.bins=1) {


  # preprocess
  l <- Preprocess(ens=ens, obs=obs)
  ens <- l[["ens"]]
  obs <- l[["obs"]]

  stopifnot(nrow(ens) == length(obs))
  N <- dim(ens)[1]
  K <- dim(ens)[2]

  if ((K+1) %% reduce.bins != 0) {
    stop("number of histogram bins is not a multiple of reduce.bins")
  }

  ranks <- apply(cbind(obs, ens), 1, rank, ties.method="random")[1, ]
  rank.hist <- hist(ranks, breaks=seq(0.5, K+1.5, 
                    reduce.bins), plot=FALSE)[["counts"]]
  return(rank.hist)
}


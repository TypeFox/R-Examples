copynumberBootstrap <-
function(experiment, replicates, quantiles,windowlength){
  if (experiment$is_GC_normalized) {observations <- experiment$CorrReadsprWindow}
  else {observations <- experiment$ReadsprWindow}
  overallmean <- mean(observations)
  means <- vector("numeric",length=nrow(experiment$genes))
  lowerbound <- vector("numeric", length=nrow(experiment$genes))
  upperbound <- vector("numeric", length=nrow(experiment$genes))
  pb <- txtProgressBar(min=0,max=nrow(experiment$genes),style=3)
  for (generow in 1:nrow(experiment$genes)){
    from <- ceiling(experiment$genes$Left[generow]/windowlength)
    to <- ceiling(experiment$genes$Right[generow]/windowlength)
    theseobs <- observations[from:to]
    means[generow] <- mean(theseobs)
    setofmeans <- replicate(n=replicates,expr=mean(sample(x=theseobs,size=length(theseobs),replace=T)))
    lowerbound[generow] <- quantile(setofmeans,probs=quantiles[1])
    upperbound[generow] <- quantile(setofmeans,probs=quantiles[2])
    rm(setofmeans)
    setTxtProgressBar(pb,generow)
  }
  close(pb)
  experiment$genes$CN_boot <- means/overallmean
  experiment$genes$Lowerbound <- lowerbound/overallmean
  experiment$genes$Upperbound <- upperbound/overallmean
  
  return(experiment)
}

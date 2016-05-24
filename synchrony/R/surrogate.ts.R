surrogate.ts <- function (ts, distr.ts=NULL, trans.ts=NULL, nbreaks=10) {
  if (NCOL(ts)==1) {
    ts=cbind(1:NROW(ts), ts)
  }
  surr.ts=numeric(length=NROW(ts))*NA
  if (is.null(trans.ts)) {  
    ## Determine transition probabilities
    distr.ts <- cut(ts[,2], quantile(ts[,2], seq(0, 1, len = (nbreaks+1))), 
               include.lowest = TRUE, labels=FALSE)
    trans.ts=matrix(nrow=nbreaks, ncol=nbreaks, 0)
    for (i in 1:(NROW(ts)-1)) {
      trans.ts[distr.ts[i], distr.ts[i+1]]=trans.ts[distr.ts[i], distr.ts[i+1]]+1      
    }
    trans.ts=trans.ts/rowSums(trans.ts)
  }
 
  ## Construct surrogate time series
  surr.ts=numeric(length=NROW(ts))*NA
  v0.ind.ts=sample(1:NROW(ts), size=1)
  v0.ts=ts[v0.ind.ts, 2]
  b0.ts=distr.ts[v0.ind.ts]
  surr.ts[1]=v0.ts
  for (i in 2:(NROW(ts))) {
    b1.ts=sample(1:nbreaks, size=1, prob=trans.ts[b0.ts,])
    v1.ind.ts=sample(which(distr.ts==b1.ts), size=1)
    v1.ts=ts[v1.ind.ts, 2]
    surr.ts[i]=v1.ts
    v0.ind.ts=v1.ind.ts
    b0.ts=b1.ts
  }
  return(list(surr.ts=cbind(ts[,1], surr.ts), trans=trans.ts, distr=distr.ts))
}
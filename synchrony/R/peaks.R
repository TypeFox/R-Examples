peaks <- function  (t1, t2, nrands = 0, type = 1, quiet = FALSE) {
  
  if (NCOL(t1)==1 | NCOL(t2)==1) {
    t1=cbind(1:NROW(t1), t1)
    t2=cbind(1:NROW(t2), t2)    
  }
  observed.peaks=peaks.aux (t1, t2)
  
  if (nrands == 0) {
    results = observed.peaks    
  }
  else {
    randomized.peaks=numeric(length=nrands+1)
    if (!quiet)
      prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    for (n in 1:nrands) {
      if (type==1) {
        t1.tmp=cbind(t1[,1], sample(t1[,2]))
        t2.tmp=cbind(t2[,1], sample(t2[,2]))
      }
      else {
        lags=sample(1:NROW(t1), size=2)
        rands=mlag(cbind(t1[,2], t2[,2]), lags)
        t1.tmp=cbind(t1[,1], rands[,1])
        t2.tmp=cbind(t2[,1], rands[,2])
      }
      randomized.peaks[n]=peaks.aux(t1.tmp, t2.tmp)$obs
      if (!quiet)
        setTxtProgressBar(prog.bar, n)
    }
    randomized.peaks[n+1]=observed.peaks$obs
    pval=sum(randomized.peaks >= observed.peaks$obs)/(nrands+1)
    results=list(pval=pval, rands=randomized.peaks, obs=observed.peaks$obs, 
                 locations=observed.peaks$locations, index=observed.peaks$index)
  }  
  class(results)="synchrony"
  return (results)
}

peaks.aux <- function (t1, t2) {
  f1=find.minmax(t1)
  f2=find.minmax(t2)
  common.mins=f1$mins$index %in% f2$mins$index
  common.maxs=f1$maxs$index %in% f2$maxs$index
  peaks=(sum(common.mins)+sum(common.maxs))/sum(max(NROW(f1$mins), 
                                                    NROW(f2$mins)) +
                                                  max(NROW(f1$maxs), 
                                                      NROW(f2$maxs)))
  index=sort(c(f1$mins$index[common.mins], 
               f1$maxs$index[common.maxs]))
  return (list(obs=peaks, index=index))
}

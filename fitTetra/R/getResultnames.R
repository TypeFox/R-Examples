getResultnames <-
function(ng) {
  #these are the names of the allmodeldata component of the return value of fitTetra
  c("marker","markername","m","model","nsamp","nsel","npar","iter","dip","LL","AIC","BIC",
     "minsepar","selcrit","meanP","P80","P90","P95","P975","P99",
     paste(rep("muact",ng),0:(ng-1),sep=""), #actual means of the samples in each peak on arcsine-sqrt transformed scale
     paste(rep("sdact",ng),0:(ng-1),sep=""), #actual sd's of the samples in each peak on arcsine-sqrt transformed scale
     paste(rep("Pact",ng),0:(ng-1),sep=""), #actual freqs of the samples in each peak
     paste(rep("mutrans",ng),0:(ng-1),sep=""), #mu's on arcsine-sqrt transformed scale
     paste(rep("sdtrans",ng),0:(ng-1),sep=""), #sd's on arcsine-sqrt transformed scale
     paste(rep("P",ng),0:(ng-1),sep=""),   #mixture proportions
     paste(rep("mu",ng),0:(ng-1),sep=""),  #backtransformed mu's
     paste(rep("sd",ng),0:(ng-1),sep=""), #backtransformed sd's
     "message"
    )
}

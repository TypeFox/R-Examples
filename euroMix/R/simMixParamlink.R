simMixParamlink <-
function(y,alleles)
lapply(y$markerdata, function(m) alleles[sort(unique(m[m>0]))])

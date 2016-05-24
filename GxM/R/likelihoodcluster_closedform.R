likelihoodcluster_closedform <-
function(param, paramnames, zeroset, listdata, fMP, cl) {
  d1 = length(param);
  d2 = length(zeroset);
  paramlist = vector('list',length=(d1+d2));
  names(paramlist) = c(paramnames, zeroset);
  for (i in 1:d1)  paramlist[[i]] = param[i];
  if(d2>0) { for (j in 1:d2)  paramlist[[d1+j]] = 0; }
  fmp = parLapply(cl=cl, X=listdata, fun=fMP, param=paramlist);
  return(sum(unlist(fmp)));
}

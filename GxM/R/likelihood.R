likelihood <-
function(param, paramnames, zeroset, listdata, fMP, PointsW) {
  d1 = length(param);
  d2 = length(zeroset);
  paramlist = vector('list',length=(d1+d2));
  names(paramlist) = c(paramnames, zeroset);
  for (i in 1:d1)  paramlist[[i]] = param[i];
  if(d2>0) { for (j in 1:d2)  paramlist[[d1+j]] = 0; }
  fmp = lapply(X=listdata, FUN=fMP, param=paramlist, PointsW=PointsW);
  return(sum(unlist(fmp)));
}

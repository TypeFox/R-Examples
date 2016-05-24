rtruncated=function(n,lo,hi,pf,qf,...)
qf(pf(lo,...)+runif(n)*(pf(hi,...)-pf(lo,...)),...)
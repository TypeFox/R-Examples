Ims<-
function(n, k, N, pbar, cm = 1,type = c("exact", "napprox","ewmaSK","ewma2"),lam=1) 
n * cm + (N - n) * (1 - OC(pbar,n,k,type,lam))

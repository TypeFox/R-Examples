## speclagreg.K.crossval = function(X,Y,freq,lags=-3:3){  
##   r = c()
##   n = dim(X)[1]
##   for (K in 1:d){
##     tr = 1:floor(0.8*n)
##     ts = (length(tr)+1):n
##     A.est = speclagreg(X[tr,], as.matrix(Y[tr,]), lags=lags, K=K)
##     Ye = linproc(X,A.est)
##     E = sum( (Ye[ts,] - Y[ts,])^2 )
##     r = c(r,E)
##   }
##   Kb = freq
##   Kb[] = which.min(r)
##   Kb
## }  

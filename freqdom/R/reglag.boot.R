## @export
## reglag.boot = function(X,Y,A,rep=20,Kconst=1,K=NULL,freq=NULL,p=10,q=10,weights="Bartlett",plot=FALSE)
## {
##   CX = lagged.cov(X,lag=0)
##   CY = lagged.cov(Y,lag=0)
  
##   n = dim(X)[1]
##   dX = dim(X)[2]
##   dY = dim(Y)[2]
  
##   lags = A$lags

##   res = c()
##   cat(paste("Bootstraping,",rep,"repetitions.\n"))
  
##   for (i in 1:rep){
##     cat(".")
##     simX = rmvnorm(n,rep(0,dX),CX)
##     simY = rmvnorm(n,rep(0,dY),CY)
##     model = speclagreg(simX,simY,lags=lags,Kconst=Kconst,K=K,freq=freq,p=p,q=q,weights=weights)
##     res = rbind(res,timedom.norms(model)$norms)
##   }
  
##   cat("\n")
##   colnames(res) = lags
##   res = apply(res,2,function(X) { quantile(X,c(0.75,0.9,0.95)) })
  
##   Ans = timedom.norms(A)
##   if (plot){
##     ylim = c(0,max(res,Ans$norms))
##     plot(lags,Ans$norms,ylim=ylim,ylab="norm",xlab="lag")
##     matlines(colnames(res),t(res),type='l')
##   }
  
##   cat("Suggested lags")
##   suglags = A$lags[Ans$norms > res[3,]]
##   print(suglags)
##   list(lags=lags,quntiles=res,suglags=suglags)
## }

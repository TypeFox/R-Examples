compute.lr.stat <-
function(U, y, mu=NULL, sigma=NULL, exact=TRUE, verbose=FALSE, tol=10^-8){
  M = ncol(U)
  
  max.ll.obj = lapply (1:ncol(y), function(s){maximise.lr (y=y[,s], U=U, mu=mu, sigma=sigma, theta=NULL, exact=exact, verbose=verbose, tol=tol)})
  max.ll.0.obj = lapply (1:ncol(y), function(s){maximise.lr (y=y[,s], U=U, mu=mu, sigma=sigma, theta=0, exact=exact, verbose=verbose, tol=tol)})
  
  max.ll = sapply (max.ll.obj, function(x){x$ll})
  max.ll.0 = sapply (max.ll.0.obj, function(x){x$ll})
  
  mu.hat = sapply (max.ll.obj, function(x){x$mu})
  sigma.hat = sapply (max.ll.obj, function(x){x$sigma})
  
  list(ts=2*(max.ll - max.ll.0), mu.hat=mu.hat, sigma.hat=sigma.hat)
}

# mygllm.r: interface to modified gllm

mygllm <- function(y,s,X,maxit=1000, tol=0.00001, E=rep(1,length(s)))
{
   out <- .C("mygllm", as.integer(y), as.integer(s-1), 
             X=as.double(cbind(X,rep(0,nrow(X)))),
             as.integer(maxit), as.double(tol), estimate=as.double(E),
             length(s),length(y),ncol(X), PACKAGE="RecordLinkage")
   return(out$estimate)
}

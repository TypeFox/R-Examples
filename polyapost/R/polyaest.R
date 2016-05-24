polyaest<-function(A1, A2, b2, initsol, rep, ysamp, burnin)
{
  if(! is.matrix(A1)) stop("A1 not matrix")
  if(! is.matrix(A2)) stop("A2 not matrix")
  if(nrow(A2)!=length(b2)) {stop(" the no. of rows of the constraint matrix does not match the length of the rhs vector ")}
  P<-nullspace(A1)%*%t(nullspace(A1))
  chain<-means(P, ncol(P), A2, nrow(A2), b2, initsol, rep, ysamp)
  chain<-chain[burnin:rep]
  estimate<-mean(chain)
  quantiles<-quantile(chain, c(.025,.975))
  return(list(chain, estimate, quantiles))
}


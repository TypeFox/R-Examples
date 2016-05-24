# pfMLLik
# particle filter unbiased estimate of marginal likelihood (logged)

pfMLLik <- function(n,simx0,t0,stepFun,dataLik,data)
{
  times=c(t0,as.numeric(rownames(data)))
  deltas=diff(times)
  return(function(...){
    xmat=simx0(n,t0,...)
    ll=0
    for (i in 1:length(deltas)) {
      xmat=t(apply(xmat,1,stepFun,t0=times[i],deltat=deltas[i],...))
      w=apply(xmat,1,dataLik,t=times[i+1],y=data[i,],log=FALSE,...)
      if (max(w)<1e-20) {
        warning("Particle filter bombed")
        return(-1e99)
      }
      ll=ll+log(mean(w))
      rows=sample(1:n,n,replace=TRUE,prob=w)
      xmat=xmat[rows,]
    }
    ll
  })
}




# eof


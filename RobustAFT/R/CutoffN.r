CutoffN <- function(res,delta,cu=2.575829,zmax=10,gridsize=2000)
{
# computes cutoff on residual scale
  xs <- seq(from=cu,to=zmax,length=gridsize)
  fs <- apply(as.matrix(xs),1,RappN,res=res,delta=delta)
  alpha <- min(min(fs[2,],na.rm=TRUE),1)
  xu <- max(xs[fs[1,] < alpha],na.rm=TRUE)
  list(alpha=alpha,tu=xu)
}
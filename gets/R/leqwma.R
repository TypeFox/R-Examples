leqwma <-
function(x, length=5, lag=1, start=1, p=2, as.vector=FALSE)
{
  eqwma(x, length=length, lag=lag, start=start, p=p,
    log=TRUE, abs=TRUE, as.vector=as.vector)
}

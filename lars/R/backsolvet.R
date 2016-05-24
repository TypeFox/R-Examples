backsolvet <-
function(r, x, k=ncol(r))
{
  backsolve(r,x,k,transpose=TRUE)
}


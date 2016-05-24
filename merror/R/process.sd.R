"process.sd" <- function(x)
{
  # Log transform first, then exp back
  # Should be equivalent to:
  # sqrt(prod(var(x)[upper.tri(var(x))])^(2/(dim(x)[2]*(dim(x)[2]-1))))
  
  sqrt(exp(sum(log(var(x)[upper.tri(var(x))]))*(2/(dim(x)[2]*(dim(x)[2]-1)))))
}

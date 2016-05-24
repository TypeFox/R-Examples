
rfmc <- function(n, P, pi0)
{
  v = vector("numeric",n)
  r = length(pi0)
  v[1] = sample(r,1,prob=pi0)
  for (i in 2:n) {
    v[i] = sample(r,1,prob=P[v[i-1],])
  }
  ts(v)
}

# eof


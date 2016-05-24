`ddst.uniform.Nk` <-
function(x, base = ddst.base.legendre, Dmax = 10) {
  n = length(x)
  maxN = max(min(Dmax, n-2, 20),1)
  coord = numeric(maxN)
  for (j in 1:maxN) 
   coord[j] = ddst.phi(x, j, base)
  coord= cumsum(coord^2*n)
}


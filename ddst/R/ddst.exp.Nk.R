`ddst.exp.Nk` <-
function(x, base = ddst.base.legendre, Dmax = 5, n=length(x)) {
er = mean(x)
maxN = max(min(Dmax, n-2, 20),1)
coord = numeric(maxN)
u = numeric(maxN)
for (j in 1:maxN) { 
u[j] = ddst.phi(pexp(x,1/er), j, base)
coord[j] = t(u[1:j]) %*% MMexp[[j]] %*% u[1:j] * n
}
coord
}


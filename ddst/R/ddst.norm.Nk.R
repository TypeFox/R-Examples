`ddst.norm.Nk` <-
function(x, base = ddst.base.legendre, Dmax = 5, n=length(x)) {
er1 = mean(x)
sx  = sort(x)
H   = qnorm((1:n - 3/8)/(n+1/4))
er2 = mean((sx[-1] - sx[-n])/(H[-1] - H[-n]))
pp   = (x-er1)/er2
tmpp = c(mean(pp), (mean(pp^2) - 1)/2)
maxN = max(min(Dmax, length(x)-2, 20),1)
u = numeric(maxN)
for (j in 1:maxN) 
u[j] = ddst.phi(pnorm(x,er1,er2), j, base)
coord = numeric(Dmax)
for (k in 1:Dmax) {
korekta = u[1:k] - t(MMnorm12[[k]]) %*% tmpp
coord[k] = t(korekta) %*% MMnorm[[k]] %*% korekta * n
}
coord
}


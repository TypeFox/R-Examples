`ddst.extr.Nk` <-
function(x, base = ddst.base.legendre, Dmax = 5, n=length(x)) {
sx  = sort(x)
er2 = sum(sx*(1:n - n:1))/(n*(n-1)*log(2))
er1 = mean(x)+0.5772156649*er2
maxN = max(min(Dmax, length(x)-2, 20),1)

u = NULL
for (j in 1:maxN) 
u[j] = ddst.phi(1-pgumbel(-x,-er1,er2), j, base)
coord = NULL
gg1=mean(1-exp((x-er1)/er2))
gg2=mean((1-exp((x-er1)/er2))*(x-er1)/er2+1)
for (k in 1:Dmax) {
korekta = u[1:k] + t(MMextr12[[k]]) %*% sM22 %*% c(gg1,gg2)
coord[k] = t(korekta) %*% MMextr[[k]] %*% korekta * n
}
coord
}


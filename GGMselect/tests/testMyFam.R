library(GGMselect)
p=30
n=30
eta=0.13
dmax=3
iG = 4
iS = 20
set.seed(iG)
Gr <- simulateGraph(p,eta)
set.seed(iS*(pi/3.1415)**iG)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr$C)

K=2.5
mon1GRestALL <- selectFast(X, dmax, K, family="LA")
resFam <- selectMyFam(X, list(mon1GRestALL$LA$G), K)

print(all.equal(mon1GRestALL$LA$G, resFam$G))
print(all.equal(mon1GRestALL$LA$Neighb[,1, drop=F], resFam$Neighb))
print(all.equal(mon1GRestALL$LA$crit.min, resFam$crit.min))

K=c(2.5, 1)
zz=selectFast(X, dmax, K, family="LA")
resFam <- selectMyFam(X, list(zz$LA$G[,,1], zz$LA$G[,,2]), K=K)

print(all.equal(zz$LA$G, resFam$G))
print(all.equal(as.vector(zz$LA$Neighb[,1,1]), as.vector(resFam$Neighb[,,1])))
print(all.equal(zz$LA$crit.min, resFam$crit.min))

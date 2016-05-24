cat("Test of pen=Inf\n")
library(GGMselect)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++
p = 30
n = 20
dmax = 17
eta=0.05
iG = 1
set.seed(iG)
Gr <- simulateGraph(p,eta)
iS = 1
set.seed(iS*(pi/3.1415)**iG)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr$C)

# LONG !QEGraphEst <- selectQE(X, dmax=dmax, K=2,   verbose=T)
GraphEst <- selectFast(X, dmax=dmax, K=2,
                       verbose=T,
                       family=c("LA", "C01","EW" ))

#  print("Neighb")
#  print(GraphEst$LA$Neighb)
  print("crit.min");  print(GraphEst$LA$crit.min)
  print("crit.min");  print(GraphEst$EW$crit.min)
  print("crit.min");  print(GraphEst$C01$crit.min)
  print("crit.min");  print(GraphEst$C01.LA$crit.min)
  print("crit.min");  print(GraphEst$C01.LA.EW$crit.min)

library(GGMselect)
itest=2
# +++++++++++++++++++++++++++++++++++++++++++++++++
p=30
n=30
eta=0.13
dmax=3
iG = 7
iS = 9
set.seed(iG)
Gr <- simulateGraph(p,eta)
set.seed(iS*(pi/3.1415)**iG)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr$C)

K=2.5

# ptm <- proc.time()

print(selectQE(X, dmax, K, verbose=TRUE))

# cat("Elapsed time ",  proc.time()-ptm,"\n")
cat ("End of test ", itest, "\n")



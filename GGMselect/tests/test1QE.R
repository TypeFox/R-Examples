library(GGMselect)
itest=1
# +++++++++++++++++++++++++++++++++++++++++++++++++
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

# ptm <- proc.time()

print( selectQE(X, dmax, K, verbose=TRUE))

# cat("Elapsed time ",  proc.time()-ptm,"\n")

# Test of special options
dmax=c(rep(2,25), rep(3, 5)) # p-vector of maximal degree
max.mod = 10 # Force the stepwise
K= c(2.5, 3)
print(selectQE(X, dmax, K, verbose=TRUE, max.nG=max.mod))

cat ("End of test ", itest, "\n")



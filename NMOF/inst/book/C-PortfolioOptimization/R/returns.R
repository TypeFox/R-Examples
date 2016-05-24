# returns.R -- version 2010-09-11
## generate artificial price data: 
#  R = returns, P = prices
ns <- 100L  # number of scenarios
na <- 10L   # number of assets
R  <- 1 + array(rnorm(ns * na) * 0.01, dim = c(ns, na))
P  <- rbind(100,R)
P  <- apply(P, 2, cumprod)
matplot(P, type = "l")  # ... or ts.plot(P)

## discrete returns
# compute returns: rets should be equal to R
rets1 <- P[2L:nrow(P), ] / P[1L:(nrow(P) - 1L), ] 
# ... or
rets2 <- diff(P) / P[1L:(nrow(P)-1L), ] + 1 
max(abs(rets1 - R))  # not 'exactly' equal
max(abs(rets2 - R))  # not 'exactly' equal
max(abs(rets1 - rets1))  # 'exactly' equal

## log-returns
rets3 <- diff(log(P))
# ... almost like discrete returns
plot(as.vector(rets1) - as.vector(rets3) - 1L, cex = 0.5)
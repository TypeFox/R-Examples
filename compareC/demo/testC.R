
set.seed(2014)

### empirical type I error rate

rej.perc = numeric(2000)
for (i in 1:2000) {
nn = 100
lifetimes <- rexp(nn, rate = exp(1))
censtimes <- rexp(nn, rate = 0.1)
x1 <- rnorm(nn, lifetimes)
x2 <- rnorm(nn, lifetimes)
ztimes <- pmin(lifetimes, censtimes)
status <- as.numeric(censtimes > lifetimes)

rej.perc[i] = compareC(ztimes, status, x1, x2)$pval
}

mean(rej.perc < 0.05)   ### 0.0465

compareC(ztimes, status, x1, x2)   ### syntax to compare two Cs

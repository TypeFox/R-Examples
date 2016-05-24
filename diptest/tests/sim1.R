#### Very small scale simulation to make the point
#### --> See ../stuff/  for much more
library(diptest)

P.p <- c(1, 5, 10, 25)/100
(P.p <- c(P.p, 1/2, rev(1 - P.p)))

N.sim <- 9999
set.seed(94)
.p0 <- proc.time()
dU100 <- replicate(N.sim, dip(runif(100)))
cat('Time elapsed: ', (p1 <- proc.time()) - .p0,'\n'); .p0 <- p1
## Lynne (2003: P IV, 1.6 GHz): ~7   s
## 2010 (AMD Phenom II X4 925):  1.3 s

100 * round(q100 <- quantile(dU100, p = P.p), 4)

plot(density(sqrt(100) * dU100), lwd = 2, col=2,
     main = expression("Dip  distribution" ~~
         list(sqrt(n)* D[n], ~ n == 100)))
abline(h=0, col="dark gray", lty=3)

round(1e4 * quantile(dU100, p = seq(0,1, by = 0.01), names = FALSE))

##--- an extreme unimodal case -- i.e. very small dip():
set.seed(60); x <- rexp(301,1)^3
hist(x)
(dt.x <- dip.test(x))
(dt2 <- dip.test(x, simulate = TRUE))
(dt3 <- dip.test(x, simulate = TRUE, B = 10000))
stopifnot(dt.x$p.value == 1,## <- gave	NA  earlier
	  dt2$p.value == 1,
	  dt3$p.value == 1)


cat('Time elapsed: ', proc.time() - .p0,'\n') # "stats"

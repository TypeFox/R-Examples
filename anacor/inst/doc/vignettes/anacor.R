### R code from vignette source 'anacor.Rnw'

###################################################
### code chunk number 1: anacor.Rnw:357-361
###################################################
library("anacor")
data("tocher")
res <- anacor(tocher, scaling = c("standard", "centroid"))
res


###################################################
### code chunk number 2: anacor.Rnw:363-365 (eval = FALSE)
###################################################
## plot(res, plot.type = "jointplot", xlim = c(-2,1.5), ylim = c(-2, 1.5), asp = 1)
## plot(res, plot.type = "graphplot", xlim = c(-2,1.5), ylim = c(-2, 1.5), wlines = 5, asp = 1)


###################################################
### code chunk number 3: anacor.Rnw:381-382
###################################################
chisq.test(tocher)


###################################################
### code chunk number 4: anacor.Rnw:390-394
###################################################
data(bitterling)
res1 <- anacor(bitterling, ndim = 2, scaling = c("Benzecri", "Benzecri"))
res2 <- anacor(bitterling, ndim = 5, scaling = c("Benzecri", "Benzecri"))
res2


###################################################
### code chunk number 5: anacor.Rnw:396-398 (eval = FALSE)
###################################################
## plot(res1, plot.type = "benzplot", main = "Benzecri Distances (2D)")
## plot(res2, plot.type = "benzplot", main = "Benzecri Distances (5D)")


###################################################
### code chunk number 6: anacor.Rnw:416-419 (eval = FALSE)
###################################################
## data(glass)
## res <- anacor(glass)
## plot(res, plot.type = "regplot", xlab = "fathers occupation", ylab = "sons occupation", asp = 1)


###################################################
### code chunk number 7: anacor.Rnw:442-447
###################################################
data(maxwell)
res <- anacor(maxwell$table, row.covariates = maxwell$row.covariates, scaling = c("Goodman", "Goodman"))
res
plot(res, plot.type = "orddiag", asp = 1)
plot(res, plot.type = "transplot", legpos = "topright")



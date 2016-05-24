### R code from vignette source 'mmod-demo.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: start
###################################################
library(mmod)
data(nancycats)
nancycats


###################################################
### code chunk number 2: diffstat
###################################################
diff_stats(nancycats)


###################################################
### code chunk number 3: fig1plot
###################################################
nc.diff_stats <- diff_stats(nancycats, phi_st=TRUE)
with(nc.diff_stats, pairs(per.locus[,3:6], upper.panel=panel.smooth))


###################################################
### code chunk number 4: fig1
###################################################
nc.diff_stats <- diff_stats(nancycats, phi_st=TRUE)
with(nc.diff_stats, pairs(per.locus[,3:6], upper.panel=panel.smooth))


###################################################
### code chunk number 5: bs
###################################################
bs <- chao_bootstrap(nancycats, nreps=10)
bs.D <- summarise_bootstrap(bs, D_Jost)
bs.D


###################################################
### code chunk number 6: show_coords
###################################################
head(nancycats@other$xy, 4)
nc.pop_dists <- dist(nancycats@other$xy, method="euclidean")


###################################################
### code chunk number 7: pw
###################################################
nc.pw_D <- pairwise_D(nancycats, linearized=TRUE)


###################################################
### code chunk number 8: mantel
###################################################
mantel.rtest(nc.pw_D, log(nc.pop_dists), 999)


###################################################
### code chunk number 9: fig2plot
###################################################
fit <- lm(as.vector(nc.pw_D) ~ as.vector(nc.pop_dists))
plot(as.vector(nc.pop_dists), as.vector(nc.pw_D),
     ylab="pairwise D", xlab="physical distance")
abline(fit)


###################################################
### code chunk number 10: fig2
###################################################
fit <- lm(as.vector(nc.pw_D) ~ as.vector(nc.pop_dists))
plot(as.vector(nc.pop_dists), as.vector(nc.pw_D),
     ylab="pairwise D", xlab="physical distance")
abline(fit)



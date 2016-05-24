###############################################################################
## MAS 5.0 vs. robloxbioc for Uni-RNA samples
###############################################################################

## load MAQC-I data
library(MAQCsubsetAFX)
data(refA)
data(refB)
data(refC)
data(refD)

## Minimum Kolmogorov distance
## takes about 50 min for each reference set (Core i5 520M with 8 GByte RAM)
library(RobLoxBioC)
system.time(minKD.A <- KolmogorovMinDist(refA, Norm()))
system.time(minKD.B <- KolmogorovMinDist(refB, Norm()))
system.time(minKD.C <- KolmogorovMinDist(refC, Norm()))
system.time(minKD.D <- KolmogorovMinDist(refD, Norm()))

## load the results for random normal samples from R-forge
con <- url("http://robast.r-forge.r-project.org/data/minKD_norm.RData")
load(file = con)
close(con)

uni.n <- rep(c(11, 16), 4)

#######################################
## Figure 4 in Kohl and Deigner (2010)
#######################################
resA <- split(as.vector(minKD.A$dist), as.vector(minKD.A$n))[c(4,8)]
resB <- split(as.vector(minKD.B$dist), as.vector(minKD.B$n))[c(4,8)]
resC <- split(as.vector(minKD.C$dist), as.vector(minKD.C$n))[c(4,8)]
resD <- split(as.vector(minKD.D$dist), as.vector(minKD.D$n))[c(4,8)]
#setEPS(height = 6, width = 9)
#postscript(file = "Figure4.eps")
par(mar = c(4, 4, 3, 1))
boxplot(c(resA, resB, resC, resD), main = "Minimum Kolmogorov distance", 
        ylim = c(0, 0.49), at = c(1, 2, 4, 5, 7, 8, 10, 11), xlim = c(0.5, 14.5),
        ylab = "minimum Kolmogorov distance", xlab = "sample size",
        pch = 20)
boxplot(minKD.norm[,c(7,12)], at = c(13, 14), add = TRUE, pch = 20)
lines(c(1,2), 1/(2*c(11,16)), lwd = 2)
lines(c(4,5), 1/(2*c(11,16)), lwd = 2)
lines(c(7,8), 1/(2*c(11,16)), lwd = 2)
lines(c(10,11), 1/(2*c(11,16)), lwd = 2)
lines(c(13,14), 1/(2*c(11,16)), lwd = 2)
abline(h = c(0.1, 0.15), lty = 2, lwd = 1.5)
abline(v = c(3, 6, 9, 12), lwd = 1.5)
text(c(1.5, 4.5, 7.5, 10.5, 13.5), rep(0.48, 5), labels = c("refA", "refB", "refC", "refD", "normal"), font = 2)
legend("bottomleft", legend = "minimal possible distance", lty = 1, 
       bg = "white", cex = 0.8)
dev.off()

## Reference set B
pdf(file = "minKD_refB.pdf")
par(mfrow = c(1, 2))
boxplot(res, main = "Reference set B", ylim = c(0, 0.45), 
        ylab = "minimum Kolmogorov distance", xlab = "sample size")
lines(1:length(uni.n), 1/(2*uni.n), col = "orange", lwd = 2)
legend("topright", legend = "minimal possible distance", fill = "orange")
abline(h = c(0.1, 0.15), lty = 2, lwd = 1.5)

boxplot(minKD.norm[,c(7,12)], main = "Normal samples", ylim = c(0, 0.45), 
        ylab = "minimum Kolmogorov distance", xlab = "sample size")
lines(1:length(uni.n), 1/(2*uni.n), col = "orange", lwd = 2)
legend("topright", legend = "minimal possible distance", fill = "orange")
abline(h = c(0.1, 0.15), lty = 2, lwd = 1.5)
dev.off()

## Reference set C
pdf(file = "minKD_refC.pdf")
par(mfrow = c(1, 2))
boxplot(res, main = "Reference set C", ylim = c(0, 0.45), 
        ylab = "minimum Kolmogorov distance", xlab = "sample size")
lines(1:length(uni.n), 1/(2*uni.n), col = "orange", lwd = 2)
legend("topright", legend = "minimal possible distance", fill = "orange")
abline(h = c(0.1, 0.15), lty = 2, lwd = 1.5)

boxplot(minKD.norm[,c(7,12)], main = "Normal samples", ylim = c(0, 0.45), 
        ylab = "minimum Kolmogorov distance", xlab = "sample size")
lines(1:length(uni.n), 1/(2*uni.n), col = "orange", lwd = 2)
legend("topright", legend = "minimal possible distance", fill = "orange")
abline(h = c(0.1, 0.15), lty = 2, lwd = 1.5)
dev.off()

## Reference set D
pdf(file = "minKD_refD.pdf")
par(mfrow = c(1, 2))
boxplot(res, main = "Reference set D", ylim = c(0, 0.45), 
        ylab = "minimum Kolmogorov distance", xlab = "sample size")
lines(1:length(uni.n), 1/(2*uni.n), col = "orange", lwd = 2)
legend("topright", legend = "minimal possible distance", fill = "orange")
abline(h = c(0.1, 0.15), lty = 2, lwd = 1.5)

boxplot(minKD.norm[,c(7,12)], main = "Normal samples", ylim = c(0, 0.45), 
        ylab = "minimum Kolmogorov distance", xlab = "sample size")
lines(1:length(uni.n), 1/(2*uni.n), col = "orange", lwd = 2)
legend("topright", legend = "minimal possible distance", fill = "orange")
abline(h = c(0.1, 0.15), lty = 2, lwd = 1.5)
dev.off()


## MAS 5.0
## takes about 4.5 minutes for each reference set (Core i5 520M with 8 GByte RAM)
system.time(res.mas5.A <- mas5(refA))
system.time(res.mas5.B <- mas5(refB))
system.time(res.mas5.C <- mas5(refC))
system.time(res.mas5.D <- mas5(refD))

## MAS rmx
## takes about 45 seconds for each reference set (Core i5 520M with 8 GByte RAM)
library(RobLoxBioC)
system.time(res.rmx.A <- robloxbioc(refA, normalize = TRUE, add.constant = 0))
system.time(res.rmx.B <- robloxbioc(refB, normalize = TRUE, add.constant = 0))
system.time(res.rmx.C <- robloxbioc(refC, normalize = TRUE, add.constant = 0))
system.time(res.rmx.D <- robloxbioc(refD, normalize = TRUE, add.constant = 0))

## Spearman correlations
cor.mas5.A <- cor(exprs(res.mas5.A), method = "spearman")
cor.rmx.A <- cor(exprs(res.rmx.A), method = "spearman")
(diff.A <- cor.rmx.A-cor.mas5.A)
(rel.A <- cor.rmx.A/cor.mas5.A)

cor.mas5.B <- cor(exprs(res.mas5.B), method = "spearman")
cor.rmx.B <- cor(exprs(res.rmx.B), method = "spearman")
(diff.B <- cor.rmx.B-cor.mas5.B)
(rel.B <- cor.rmx.B/cor.mas5.B)

cor.mas5.C <- cor(exprs(res.mas5.C), method = "spearman")
cor.rmx.C <- cor(exprs(res.rmx.C), method = "spearman")
(diff.C <- cor.rmx.C-cor.mas5.C)
(rel.C <- cor.rmx.C/cor.mas5.C)

cor.mas5.D <- cor(exprs(res.mas5.D), method = "spearman")
cor.rmx.D <- cor(exprs(res.rmx.D), method = "spearman")
(diff.D <- cor.rmx.D-cor.mas5.D)
(rel.D <- cor.rmx.D/cor.mas5.D)

100*(range(rel.A[col(rel.A) > row(rel.A)])-1)
100*(range(rel.B[col(rel.B) > row(rel.B)])-1)
100*(range(rel.C[col(rel.C) > row(rel.C)])-1)
100*(range(rel.D[col(rel.D) > row(rel.D)])-1)
range(diff.A[col(diff.A) > row(diff.A)])
range(diff.B[col(diff.B) > row(diff.B)])
range(diff.C[col(diff.C) > row(diff.C)])
range(diff.D[col(diff.D) > row(diff.D)])

## Pearson correlations of log2-transformed data
cor.mas5.A1 <- cor(log2(exprs(res.mas5.A)))
cor.rmx.A1 <- cor(log2(exprs(res.rmx.A)))
(diff.A1 <- cor.rmx.A1-cor.mas5.A1)
(rel.A1 <- cor.rmx.A1/cor.mas5.A1)

cor.mas5.B1 <- cor(log2(exprs(res.mas5.B)))
cor.rmx.B1 <- cor(log2(exprs(res.rmx.B)))
(diff.B1 <- cor.rmx.B1-cor.mas5.B1)
(rel.B1 <- cor.rmx.B1/cor.mas5.B1)

cor.mas5.C1 <- cor(log2(exprs(res.mas5.C)))
cor.rmx.C1 <- cor(log2(exprs(res.rmx.C)))
(diff.C1 <- cor.rmx.C1-cor.mas5.C1)
(rel.C1 <- cor.rmx.C1/cor.mas5.C1)

cor.mas5.D1 <- cor(log2(exprs(res.mas5.D)))
cor.rmx.D1 <- cor(log2(exprs(res.rmx.D)))
(diff.D1 <- cor.rmx.D1-cor.mas5.D1)
(rel.D1 <- cor.rmx.D1/cor.mas5.D1)

100*(range(rel.A1[col(rel.A1) > row(rel.A1)])-1)
100*(range(rel.B1[col(rel.B1) > row(rel.B1)])-1)
100*(range(rel.C1[col(rel.C1) > row(rel.C1)])-1)
100*(range(rel.D1[col(rel.D1) > row(rel.D1)])-1)
range(diff.A1[col(diff.A1) > row(diff.A1)])
range(diff.B1[col(diff.B1) > row(diff.B1)])
range(diff.C1[col(diff.C1) > row(diff.C1)])
range(diff.D1[col(diff.D1) > row(diff.D1)])


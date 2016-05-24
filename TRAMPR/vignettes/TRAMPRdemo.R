### R code from vignette source 'TRAMPRdemo.Rnw'

###################################################
### code chunk number 1: TRAMPRdemo.Rnw:38-39
###################################################
library(TRAMPR)


###################################################
### code chunk number 2: TRAMPRdemo.Rnw:68-73
###################################################
files <- c("demo_samples_abi.txt",
           "demo_samples_abi_template_full.csv",
           "demo_samples_abi_info_full.csv",
           "demo_samples_abi_soilcore.csv")
file.copy(system.file(files, package="TRAMPR"), ".")


###################################################
### code chunk number 3: TRAMPRdemo.Rnw:85-86
###################################################
load.abi.create.template("demo_samples_abi.txt")


###################################################
### code chunk number 4: TRAMPRdemo.Rnw:100-102
###################################################
file.rename("demo_samples_abi_template_full.csv",
            "demo_samples_abi_template.csv")


###################################################
### code chunk number 5: TRAMPRdemo.Rnw:113-114
###################################################
load.abi.create.info("demo_samples_abi.txt")


###################################################
### code chunk number 6: TRAMPRdemo.Rnw:120-122
###################################################
file.rename("demo_samples_abi_info_full.csv",
            "demo_samples_abi_info.csv")


###################################################
### code chunk number 7: TRAMPRdemo.Rnw:131-132
###################################################
soilcore <- read.csv("demo_samples_abi_soilcore.csv")


###################################################
### code chunk number 8: TRAMPRdemo.Rnw:142-143
###################################################
primer.translate <- list(ITS1F="B", ITS4=c("G", "Y"))


###################################################
### code chunk number 9: TRAMPRdemo.Rnw:152-155
###################################################
demo.samples <- load.abi("demo_samples_abi.txt",
                         primer.translate=primer.translate,
                         soilcore=soilcore)


###################################################
### code chunk number 10: TRAMPRdemo.Rnw:165-166
###################################################
data(demo.knowns)


###################################################
### code chunk number 11: TRAMPRdemo.Rnw:181-182 (eval = FALSE)
###################################################
## plot(demo.knowns)


###################################################
### code chunk number 12: TRAMPRdemo.Rnw:187-188
###################################################
plot(demo.knowns)


###################################################
### code chunk number 13: TRAMPRdemo.Rnw:207-208
###################################################
fit <- TRAMP(demo.samples, demo.knowns)


###################################################
### code chunk number 14: TRAMPRdemo.Rnw:214-215 (eval = FALSE)
###################################################
## plot(fit, 565)


###################################################
### code chunk number 15: TRAMPRdemo.Rnw:220-221
###################################################
plot(fit, 565)


###################################################
### code chunk number 16: TRAMPRdemo.Rnw:235-238 (eval = FALSE)
###################################################
## pdf("TRAMP_fits.pdf")
## plot(fit)
## dev.off()


###################################################
### code chunk number 17: TRAMPRdemo.Rnw:244-246
###################################################
m <- summary(fit)
m[1:5, 1:5]


###################################################
### code chunk number 18: TRAMPRdemo.Rnw:254-255
###################################################
knowns.unk <- build.knowns(demo.samples)


###################################################
### code chunk number 19: TRAMPRdemo.Rnw:266-268 (eval = FALSE)
###################################################
## knowns.sporocarps <- build.knowns(demo.samples, restrict=TRUE,
##                                   min.ratio=0)


###################################################
### code chunk number 20: TRAMPRdemo.Rnw:275-276
###################################################
fit <- combine(fit, knowns.unk)


###################################################
### code chunk number 21: TRAMPRdemo.Rnw:281-282
###################################################
demo.knowns <- combine(demo.knowns, knowns.unk)


###################################################
### code chunk number 22: TRAMPRdemo.Rnw:289-290
###################################################
identical(demo.knowns, fit$knowns)


###################################################
### code chunk number 23: TRAMPRdemo.Rnw:298-299 (eval = FALSE)
###################################################
## demo.knowns <- fit$knowns


###################################################
### code chunk number 24: TRAMPRdemo.Rnw:310-312 (eval = FALSE)
###################################################
## samples <- setdiff(labels(demo.samples), labels(demo.knowns)) 
## plot(fit, samples)


###################################################
### code chunk number 25: TRAMPRdemo.Rnw:318-319 (eval = FALSE)
###################################################
## plot(fit, 221)


###################################################
### code chunk number 26: TRAMPRdemo.Rnw:324-325 (eval = FALSE)
###################################################
## plot(demo.samples, 221)


###################################################
### code chunk number 27: TRAMPRdemo.Rnw:334-335
###################################################
fit <- add.known(fit, 221, prompt=FALSE)


###################################################
### code chunk number 28: TRAMPRdemo.Rnw:340-341
###################################################
names(which(summary(fit)[,"221"]))


###################################################
### code chunk number 29: TRAMPRdemo.Rnw:349-350 (eval = FALSE)
###################################################
## plot(fit, 586)


###################################################
### code chunk number 30: TRAMPRdemo.Rnw:355-356
###################################################
plot(fit, 586)


###################################################
### code chunk number 31: TRAMPRdemo.Rnw:372-373 (eval = FALSE)
###################################################
## plot(fit, 586, grouped=TRUE)


###################################################
### code chunk number 32: TRAMPRdemo.Rnw:378-379
###################################################
plot(fit, 586, grouped=TRUE)


###################################################
### code chunk number 33: TRAMPRdemo.Rnw:401-402 (eval = FALSE)
###################################################
## fit <- update(fit)


###################################################
### code chunk number 34: TRAMPRdemo.Rnw:416-418
###################################################
m <- summary(fit, grouped=TRUE)
m <- m[, colSums(m) > 0]


###################################################
### code chunk number 35: matrix (eval = FALSE)
###################################################
## cores <- fit$samples$info$soilcore.fk
## m.bycore <- aggregate(m, by=list(cores=cores), any)
## rownames(m.bycore) <- m.bycore$cores
## m.bycore <- m.bycore[-1]
## m.bycore


###################################################
### code chunk number 36: TRAMPRdemo.Rnw:438-439
###################################################
cores <- fit$samples$info$soilcore.fk
m.bycore <- aggregate(m, by=list(cores=cores), any)
rownames(m.bycore) <- m.bycore$cores
m.bycore <- m.bycore[-1]
m.bycore


###################################################
### code chunk number 37: TRAMPRdemo.Rnw:451-453
###################################################
core.veg <- soilcore$vegetation[match(rownames(m.bycore),
                                      soilcore$soilcore.pk)]


###################################################
### code chunk number 38: rankabundance (eval = FALSE)
###################################################
## sp.freq <- sort(colSums(m.bycore), decreasing=TRUE)
## plot(sp.freq, xlab="Species rank", ylab="Species frequency", type="o",
##      pch=19)


###################################################
### code chunk number 39: TRAMPRdemo.Rnw:467-468
###################################################
tapply(rowSums(m.bycore), core.veg, mean)


###################################################
### code chunk number 40: TRAMPRdemo.Rnw:473-474
###################################################
sp.freq <- sort(colSums(m.bycore), decreasing=TRUE)
plot(sp.freq, xlab="Species rank", ylab="Species frequency", type="o",
     pch=19)


###################################################
### code chunk number 41: pca (eval = FALSE)
###################################################
## pca <- prcomp(m.bycore)
## 
## col <- c("blue", "red")[as.integer(core.veg)]
## pch <- c(1, 3)[as.integer(core.veg)]
## 
## plot(pca$x[,1], pca$x[,2], col=col, pch=pch, xlab="PC1", ylab="PC2")
## legend("top", levels(core.veg), col=c("blue", "red"), pch=c(1, 3))


###################################################
### code chunk number 42: TRAMPRdemo.Rnw:502-503
###################################################
pca <- prcomp(m.bycore)

col <- c("blue", "red")[as.integer(core.veg)]
pch <- c(1, 3)[as.integer(core.veg)]

plot(pca$x[,1], pca$x[,2], col=col, pch=pch, xlab="PC1", ylab="PC2")
legend("top", levels(core.veg), col=c("blue", "red"), pch=c(1, 3))


###################################################
### code chunk number 43: TRAMPRdemo.Rnw:520-521 (eval = FALSE)
###################################################
## write.TRAMPknowns(fit$knowns, "my_knowns")


###################################################
### code chunk number 44: TRAMPRdemo.Rnw:534-535 (eval = FALSE)
###################################################
## my.knowns <- read.TRAMPknowns("my_knowns")


###################################################
### code chunk number 45: TRAMPRdemo.Rnw:541-543 (eval = FALSE)
###################################################
## my.knowns <- fit$knowns
## save(my.knowns, file="my_knowns.Rdata")


###################################################
### code chunk number 46: TRAMPRdemo.Rnw:554-555 (eval = FALSE)
###################################################
## load("my_knowns.Rdata")



### R code from vignette source 'mefa.Rnw'

###################################################
### code chunk number 1: mefa.Rnw:56-57
###################################################
options(prompt = "R> ", continue = "+   ", useFancyQuotes = FALSE, width = 76)


###################################################
### code chunk number 2: mefa.Rnw:108-111
###################################################
library("mefa")
data("dol.count")
head(dol.count, 16)


###################################################
### code chunk number 3: mefa.Rnw:137-140
###################################################
x1 <- stcs(dol.count)
str(x1)
unique(x1$count)


###################################################
### code chunk number 4: mefa.Rnw:145-149
###################################################
x2 <- stcs(dol.count, expand = TRUE)
str(x2)
sum(x2$count)
unique(x2$count)


###################################################
### code chunk number 5: mefa.Rnw:158-161
###################################################
x3 <- stcs(dol.count, drop.zero = TRUE)
str(x3)
unique(x3$count)


###################################################
### code chunk number 6: mefa.Rnw:170-173
###################################################
m1 <- mefa(x1)
m1
m1$xtab["LT1", ]


###################################################
### code chunk number 7: mefa.Rnw:178-180
###################################################
str(m1$xtab)
str(m1$segm)


###################################################
### code chunk number 8: mefa.Rnw:185-186
###################################################
mefa(x1, nested = TRUE)


###################################################
### code chunk number 9: mefa.Rnw:191-192
###################################################
mefa(x1, segment = FALSE)


###################################################
### code chunk number 10: mefa.Rnw:199-201 (eval = FALSE)
###################################################
## data("dol.samp")
## str(dol.samp)


###################################################
### code chunk number 11: mefa.Rnw:204-206
###################################################
data("dol.samp")
str(dol.samp, strict.width = "cut")


###################################################
### code chunk number 12: mefa.Rnw:211-213 (eval = FALSE)
###################################################
## data("dol.taxa")
## str(dol.taxa)


###################################################
### code chunk number 13: mefa.Rnw:216-218
###################################################
data("dol.taxa")
str(dol.taxa, strict.width = "cut")


###################################################
### code chunk number 14: mefa.Rnw:223-230 (eval = FALSE)
###################################################
## mefa(x1, samp = dol.samp)
## mefa(x1, taxa = dol.taxa)
## m2 <- mefa(x1, samp = dol.samp, taxa = dol.taxa)
## m2
## str(m2$xtab)
## str(m2$samp)
## str(m2$taxa)


###################################################
### code chunk number 15: mefa.Rnw:233-240
###################################################
mefa(x1, samp = dol.samp)
mefa(x1, taxa = dol.taxa)
m2 <- mefa(x1, samp = dol.samp, taxa = dol.taxa)
m2
str(m2$xtab)
str(m2$samp, strict.width = "cut")
str(m2$taxa, strict.width = "cut")


###################################################
### code chunk number 16: mefa.Rnw:247-253 (eval = FALSE)
###################################################
## m2.sub <- mefa(x1, dol.samp[-c(1:5), ], dol.taxa[-c(1:80), ],
##   xtab.fixed = FALSE)
## m2.sub
## str(m2.sub$xtab)
## str(m2.sub$samp)
## str(m2.sub$taxa)


###################################################
### code chunk number 17: mefa.Rnw:256-262
###################################################
m2.sub <- mefa(x1, dol.samp[-c(1:5), ], dol.taxa[-c(1:80), ],
  xtab.fixed = FALSE)
m2.sub
str(m2.sub$xtab)
str(m2.sub$samp, strict.width = "cut")
str(m2.sub$taxa, strict.width = "cut")


###################################################
### code chunk number 18: mefastr-raw (eval = FALSE)
###################################################
## mefalogo()


###################################################
### code chunk number 19: mefastr
###################################################
opar <- par(mar = rep(0, 4))
mefalogo()
par(opar)


###################################################
### code chunk number 20: mefa.Rnw:299-301
###################################################
dim(m2)
dimnames(m2)


###################################################
### code chunk number 21: mefa.Rnw:306-307
###################################################
summary(m2)


###################################################
### code chunk number 22: mefa.Rnw:312-313
###################################################
names(summary(m2))


###################################################
### code chunk number 23: mefa.Rnw:318-320
###################################################
summary(m2)$s.rich
summary(m2)$mfill


###################################################
### code chunk number 24: plot-raw (eval = FALSE)
###################################################
## plot(m2, 1, main = "A")
## plot(m2, 4, type = "rank", trafo = "log", main = "B")


###################################################
### code chunk number 25: plot
###################################################
opar <- par(mfrow = c(1, 2))
plot(m2, 1, main = "A")
plot(m2, 4, type = "rank", trafo = "log", main = "B")
par(opar)


###################################################
### code chunk number 26: boxplot-raw (eval = FALSE)
###################################################
## boxplot(m2, 2, main = "A")
## boxplot(m2, 3, main = "B")


###################################################
### code chunk number 27: boxplot
###################################################
opar <- par(mfrow = c(1, 2))
boxplot(m2, 2, main = "A")
boxplot(m2, 3, main = "B")
par(opar)


###################################################
### code chunk number 28: mefa.Rnw:366-368
###################################################
molten <- melt(m2, "method")
str(molten)


###################################################
### code chunk number 29: mefa.Rnw:373-375
###################################################
m3 <- mefa(molten, dol.samp, dol.taxa)
m3


###################################################
### code chunk number 30: image-raw (eval = FALSE)
###################################################
## image(m3, trafo = "log", sub = "all segments", main="A")
## for (i in 1:2) image(m3, segm = i, trafo = "log",
##   sub = dimnames(m3)$segm[i], main = LETTERS[i + 1])


###################################################
### code chunk number 31: image
###################################################
opar <- par(mfrow = c(1, 3))
image(m3, trafo = "log", sub = "all segments", main="A")
for (i in 1:2) image(m3, segm = i, trafo = "log",
  sub = dimnames(m3)$segm[i], main = LETTERS[i + 1])
par(opar)


###################################################
### code chunk number 32: mefa.Rnw:410-414
###################################################
ex1 <- m2[1:20, 11:15, "fresh"]
dim(ex1)
dim(ex1$samp)
dim(ex1$taxa)


###################################################
### code chunk number 33: mefa.Rnw:421-425
###################################################
ex2 <- m2[m2$samp$method == "time"]
levels(ex2$samp$method)
ex3 <- m2[m2$samp$method == "time", drop = TRUE]
levels(ex3$samp$method)


###################################################
### code chunk number 34: mefa.Rnw:430-435
###################################################
size.5 <- as.factor(is.na(m3$taxa$size) | m3$taxa$size < 5)
levels(size.5) <- c("large", "small")
m4 <- aggregate(m3, "microhab", size.5)
t(m4$xtab)
lapply(m4$segm, t)


###################################################
### code chunk number 35: mefa.Rnw:487-491
###################################################
set.seed(1234)
m5 <- m2[ , sample(1:dim(m2)[2], 10)]
report(m5, "report.tex", tex = TRUE, segment = TRUE,
  taxa.name = 1, author.name = 2, drop.redundant = 1)


###################################################
### code chunk number 36: show-report
###################################################
writeLines(readLines("report.tex"))
invisible(file.remove("report.tex"))


###################################################
### code chunk number 37: mefa.Rnw:525-527
###################################################
mod.amin <- glm(m2$xtab[, "amin"] ~ ., data = m2$samp, family = poisson)
summary(mod.amin)


###################################################
### code chunk number 38: mefa.Rnw:539-542
###################################################
library("MASS")
mod.abu <- glm.nb(summary(m2)$s.abu ~ .^2, data = m2$samp)
summary(mod.abu)


###################################################
### code chunk number 39: mefa.Rnw:549-552
###################################################
prop.fr <- cbind(summary(m2[ , , "fresh"])$s.abu, summary(m2)$s.abu)
mod.fr <- glm(prop.fr ~ .^2, data = m2$samp, family = binomial)
summary(mod.fr)


###################################################
### code chunk number 40: mefa.Rnw:565-566
###################################################
library("vegan")


###################################################
### code chunk number 41: mefa.Rnw:569-572
###################################################
m6 <- as.mefa(m2, drop.zero = TRUE)
m6.ado <- adonis(m6$xtab ~ .^2, data = m6$samp, permutations = 100)
m6.ado


###################################################
### code chunk number 42: cca
###################################################
m2.cca <- cca(m2$segm[["fresh"]] ~ ., data=m2$samp, 
    subset=rowSums(m2$segm[["fresh"]]) > 0)
plot(m2.cca)


###################################################
### code chunk number 43: mefa.Rnw:599-605
###################################################
m.list <- list()
n1 <- rep(c("time", "quadrat"), each = 2)
n2 <- rep(c("fresh", "broken"), 2)
n3 <- paste(n1, n2, sep=".")
for (i in 1:4) m.list[[n3[i]]] <-
  aggregate(m2[m2$samp$method == n1[i], , n2[i]], "microhab")


###################################################
### code chunk number 44: clust-raw (eval = FALSE)
###################################################
## for (i in 1:4) {
##   tmp <- hclust(dist(m.list[[i]]$xtab), "ward")
##   plot(tmp, main = LETTERS[i], sub = names(m.list)[i], xlab = "")
## }


###################################################
### code chunk number 45: clust
###################################################
opar <- par(mfrow = c(2, 2))
for (i in 1:4) {
  tmp <- hclust(dist(m.list[[i]]$xtab), "ward")
  plot(tmp, main = LETTERS[i], sub = names(m.list)[i], xlab = "")
}
par(opar)



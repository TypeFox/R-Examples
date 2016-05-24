### R code from vignette source 'archetypes.Rnw'

###################################################
### code chunk number 1: archetypes.Rnw:70-71
###################################################
options(width=80, prompt='R> ')


###################################################
### code chunk number 2: archetypes.Rnw:298-299
###################################################
library(archetypes)


###################################################
### code chunk number 3: archetypes.Rnw:322-324 (eval = FALSE)
###################################################
## data("toy")
## plot(toy)


###################################################
### code chunk number 4: archetypes.Rnw:329-334
###################################################
data("toy")

par(mar = c(2,2,0,0) + 0.1, ps = 9)
plot(toy, xlab = "", ylab = "", xlim = c(0,20), ylim = c(0,20),
     pch = 19, col = gray(0.7), cex = 0.6)


###################################################
### code chunk number 5: archetypes.Rnw:343-345
###################################################
set.seed(1986)
a <- archetypes(toy, 3, verbose = TRUE)


###################################################
### code chunk number 6: archetypes.Rnw:356-357
###################################################
a


###################################################
### code chunk number 7: archetypes.Rnw:360-361
###################################################
parameters(a)


###################################################
### code chunk number 8: archetypes.Rnw:367-369 (eval = FALSE)
###################################################
## xyplot(a, toy, chull = chull(toy))
## xyplot(a, toy, adata.show = TRUE)


###################################################
### code chunk number 9: archetypes.Rnw:374-379
###################################################
par(mfrow = c(1,2), mar = c(2,2,0,0)+0.1, ps = 9)
xyplot(a, toy, chull = chull(toy),
       xlab = "", ylab = "", xlim = c(0, 20), ylim = c(0, 20), cex = 0.6)
xyplot(a, toy, adata.show = TRUE,
       xlab = "", ylab = "", xlim = c(0, 20), ylim = c(0, 20), cex = 0.6)


###################################################
### code chunk number 10: archetypes.Rnw:399-400 (eval = FALSE)
###################################################
## movieplot(a, toy)


###################################################
### code chunk number 11: archetypes.Rnw:405-411
###################################################
par(mfrow = c(2, 4), mar = c(0, 0, 0, 0) + 0.1, ps = 9)
movieplot(a, toy, xlim = c(0, 20), ylim = c(0, 20), cex = 0.6, axes = FALSE,
          postfn = function(iter) {
            box()
            text(1, 19, paste(iter + 1, ".", sep = ""), cex = 1)
          })


###################################################
### code chunk number 12: archetypes.Rnw:430-432
###################################################
set.seed(1986)
a4 <- stepArchetypes(data = toy, k = 3, verbose = FALSE, nrep = 4)


###################################################
### code chunk number 13: archetypes.Rnw:435-436
###################################################
a4


###################################################
### code chunk number 14: archetypes.Rnw:440-441
###################################################
summary(a4)


###################################################
### code chunk number 15: archetypes.Rnw:447-448 (eval = FALSE)
###################################################
## xyplot(a4, toy)


###################################################
### code chunk number 16: archetypes.Rnw:453-456
###################################################
par(mar = c(2, 2, 0, 0) + 0.1, ps = 9)
xyplot(a4, toy, cex = 0.6, xlim = c(0, 20), ylim = c(0, 20),
       xlab = "", ylab = "")


###################################################
### code chunk number 17: archetypes.Rnw:463-464 (eval = FALSE)
###################################################
## bestModel(a4)


###################################################
### code chunk number 18: archetypes.Rnw:466-467
###################################################
print(bestModel(a4), full = FALSE)


###################################################
### code chunk number 19: archetypes.Rnw:480-488
###################################################
file <- "as.RData"
if ( file.exists(file) ) {
  load(file = file)
} else {
  set.seed(1986)
  as <- stepArchetypes(data = toy, k = 1:10, verbose = FALSE, nrep = 4)
  save(as, file = file)
}


###################################################
### code chunk number 20: archetypes.Rnw:518-519
###################################################
rss(as)


###################################################
### code chunk number 21: archetypes.Rnw:522-523
###################################################
t(sapply(as, function(a) sapply(a, '[[', 'iters')))


###################################################
### code chunk number 22: archetypes.Rnw:534-535 (eval = FALSE)
###################################################
## screeplot(as)


###################################################
### code chunk number 23: archetypes.Rnw:540-545
###################################################
par(mar = c(3, 4, 0.1, 0) + 0.1, ps = 9)
screeplot(as, cex = 0.6, ylim = c(0, 0.08), axes = FALSE)
mtext("Archetypes", side = 1, line = 2)
axis(2, las = 2)
box()


###################################################
### code chunk number 24: archetypes.Rnw:553-556 (eval = FALSE)
###################################################
## a7 <- bestModel(as[[7]])
## xyplot(a7, toy, chull = chull(toy))
## xyplot(a7, toy, adata.show = TRUE)


###################################################
### code chunk number 25: archetypes.Rnw:561-568
###################################################
a7 <- bestModel(as[[7]])

par(mfrow = c(1, 2), mar = c(2, 2, 0, 0) + 0.1, ps = 9)
xyplot(a7, toy, chull = chull(toy),
       xlim = c(0, 20), ylim = c(0, 20), cex = 0.6, xlab = "", ylab = "")
xyplot(a7, toy, adata.show = TRUE,
       xlim = c(0, 20), ylim = c(0, 20), cex = 0.6, xlab = "", ylab = "")


###################################################
### code chunk number 26: archetypes.Rnw:582-593
###################################################
file <- "gas.RData"
if ( file.exists(file) ) {
  load(file = file)
} else {
  set.seed(1986)
  gas <- stepArchetypes(data = toy, k = 1:10,
                        family = archetypesFamily("original",
                        zalphasfn = archetypes:::ginv.zalphasfn),
                        verbose = FALSE, nrep = 4)
  save(gas, file = file)
}


###################################################
### code chunk number 27: archetypes.Rnw:632-633
###################################################
rss(gas)


###################################################
### code chunk number 28: archetypes.Rnw:637-638 (eval = FALSE)
###################################################
## movieplot(gas[[9]][[3]], toy)


###################################################
### code chunk number 29: archetypes.Rnw:643-649
###################################################
par(mfrow = c(1, 4), mar = c(0, 0, 0, 0) + 0.1, ps = 9)
movieplot(gas[[9]][[3]], toy, xlim = c(0, 20), ylim = c(0, 20), cex = 0.6,
          axes = FALSE, postfn = function(iter) {
            box()
            text(1, 19, paste(iter + 1, ".", sep = ""), cex = 1)
          })


###################################################
### code chunk number 30: archetypes.Rnw:667-670 (eval = FALSE)
###################################################
## ga7 <- bestModel(gas[[7]])
## xyplot(ga7, toy, chull = chull(toy))
## xyplot(ga7, toy, adata.show = TRUE)


###################################################
### code chunk number 31: archetypes.Rnw:675-682
###################################################
ga7 <- bestModel(gas[[7]])

par(mfrow = c(1, 2), mar = c(2, 2, 0, 0) + 0.1, ps = 9)
xyplot(ga7, toy, chull = chull(toy),
       xlim = c(0, 20), ylim = c(0, 20), cex = 0.6, xlab = "", ylab = "")
xyplot(ga7, toy, adata.show = TRUE,
       xlim = c(0, 20), ylim = c(0, 20), cex = 0.6, xlab = "", ylab = "")


###################################################
### code chunk number 32: archetypes.Rnw:692-693
###################################################
apply(coef(ga7, 'alphas'), 2, range)


###################################################
### code chunk number 33: archetypes.Rnw:795-797
###################################################
data("skel")
skel2 <- subset(skel, select = -Gender)


###################################################
### code chunk number 34: archetypes.Rnw:816-817 (eval = FALSE)
###################################################
## jd()


###################################################
### code chunk number 35: archetypes.Rnw:822-824
###################################################
par(mar = c(1, 4, 0, 0) + 0.1, ps = 9)
jd()


###################################################
### code chunk number 36: archetypes.Rnw:832-833 (eval = FALSE)
###################################################
## pcplot(skel2)


###################################################
### code chunk number 37: archetypes.Rnw:838-840
###################################################
datacol <- rgb(178, 178, 178, maxColorValue = 255,
               alpha = round(255 * 0.2))


###################################################
### code chunk number 38: archetypes.Rnw:842-844 (eval = FALSE)
###################################################
## par(mar = c(5, 0.4, 0, 0.4) + 0.1, ps = 9)
## pcplot(skel2, las = 2, col = datacol)


###################################################
### code chunk number 39: archetypes.Rnw:846-852
###################################################
png(filename = "body-pcplot-raw.png", units = "px",
    width = 590, height = 430, pointsize = 12)
par(mar = c(5.5, 0.4, 0, 0.4) + 0.1)
pcplot(skel2, las = 2, col = datacol)
graphics.off()
cat("\\includegraphics{body-pcplot-raw.png}\n")


###################################################
### code chunk number 40: archetypes.Rnw:861-869
###################################################
file <- "bas.RData"
if ( file.exists(file) ) {
  load(file = file)
} else {
  set.seed(1981)
  as <- stepArchetypes(skel2, k = 1:15, verbose = FALSE)
  save(as, file = file)
}


###################################################
### code chunk number 41: archetypes.Rnw:884-885 (eval = FALSE)
###################################################
## screeplot(as)


###################################################
### code chunk number 42: archetypes.Rnw:890-895
###################################################
par(mar = c(3, 4, 0.4, 0) + 0.1, ps = 9)
screeplot(as, cex = 0.6, axes = FALSE)
mtext("Archetypes", side = 1, line = 2)
axis(2, las = 2)
box()


###################################################
### code chunk number 43: archetypes.Rnw:903-904
###################################################
a3 <- bestModel(as[[3]])


###################################################
### code chunk number 44: archetypes.Rnw:907-908
###################################################
t(parameters(a3))


###################################################
### code chunk number 45: archetypes.Rnw:911-912 (eval = FALSE)
###################################################
## barplot(a3, skel2, percentiles = TRUE)


###################################################
### code chunk number 46: archetypes.Rnw:917-921
###################################################
par(mar = c(5, 4, 0.4, 0) + 0.1, ps = 9)
barplot(a3, skel2, percentiles = TRUE,
        below.compressed.height = 0.4,
        below.compressed.srt = 90)


###################################################
### code chunk number 47: archetypes.Rnw:933-934 (eval = FALSE)
###################################################
## pcplot(a3, skel2, data.col = as.numeric(skel$Gender))


###################################################
### code chunk number 48: archetypes.Rnw:939-943
###################################################
datacol <- c(rgb(0, 205, 0, maxColorValue = 255,
                 alpha = round(255 * 0.2)),
             rgb(0, 0, 255, maxColorValue = 255,
                 alpha=round(255 * 0.2)))


###################################################
### code chunk number 49: archetypes.Rnw:945-947 (eval = FALSE)
###################################################
## par(mar = c(5, 0.4, 0, 0.4) + 0.1, ps = 9)
## pcplot(a3, skel2, las = 2, data.col = datacol[skel$Gender])


###################################################
### code chunk number 50: archetypes.Rnw:949-955
###################################################
png(filename = "body-pcplot-gender.png", units = "px",
    width = 590, height = 430, pointsize = 12)
par(mar = c(5.5, 0.4, 0, 0.4) + 0.1)
pcplot(a3, skel2, las = 2, data.col = datacol[skel$Gender])
graphics.off()
cat("\\includegraphics{body-pcplot-gender.png}\n")


###################################################
### code chunk number 51: archetypes.Rnw:968-969 (eval = FALSE)
###################################################
## ternaryplot(coef(a3, 'alphas'), col = as.numeric(skel$Gender))


###################################################
### code chunk number 52: archetypes.Rnw:974-979
###################################################
library("vcd")
ternaryplot(coef(a3, 'alphas'), dimnames = 1:3, cex = 0.3,
            col = datacol[skel$Gender], main = NULL,
            labels = "none", grid = FALSE)
grid.text("1", x = 3, y = 3)


###################################################
### code chunk number 53: archetypes.Rnw:995-996 (eval = FALSE)
###################################################
## skeletonplot(parameters(a3))


###################################################
### code chunk number 54: archetypes.Rnw:1001-1003
###################################################
par(mar = c(1, 4, 0, 0) + 0.1, ps = 9)
skeletonplot(parameters(a3), skel.height = 190)


###################################################
### code chunk number 55: archetypes.Rnw:1053-1054
###################################################
sessionInfo()



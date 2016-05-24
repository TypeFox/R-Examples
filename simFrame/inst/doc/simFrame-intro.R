### R code from vignette source 'simFrame-intro.Rnw'

###################################################
### code chunk number 1: simFrame-intro.Rnw:95-98
###################################################
options(width=75, prompt="R> ")
library("simFrame")
set.seed(1234)


###################################################
### code chunk number 2: simFrame-intro.Rnw:237-238
###################################################
showMethods("setNA")


###################################################
### code chunk number 3: simFrame-intro.Rnw:401-405
###################################################
nc <- NAControl(NArate = 0.05)
getNArate(nc)
setNArate(nc, c(0.01, 0.03, 0.05, 0.07, 0.09))
getNArate(nc)


###################################################
### code chunk number 4: simFrame-intro.Rnw:507-510
###################################################
library("mvtnorm")
dc <- DataControl(size = 10, distribution = rmvnorm, dots = 
    list(mean = rep(0, 2), sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)))


###################################################
### code chunk number 5: simFrame-intro.Rnw:517-519
###################################################
foo <- generate(dc)
foo


###################################################
### code chunk number 6: simFrame-intro.Rnw:722-727
###################################################
data("eusilcP")
set <- setup(eusilcP, size = 10, k = 2)
summary(set)
set
draw(eusilcP[, c("id", "eqIncome")], set, i = 1)


###################################################
### code chunk number 7: simFrame-intro.Rnw:827-829
###################################################
cc <- DARContControl(target = "V2", epsilon = 0.2, 
    fun = function(x) x * 100)


###################################################
### code chunk number 8: simFrame-intro.Rnw:845-847
###################################################
bar <- contaminate(foo, cc)
bar


###################################################
### code chunk number 9: simFrame-intro.Rnw:963-964
###################################################
nc <- NAControl(NArate = 0.3)


###################################################
### code chunk number 10: simFrame-intro.Rnw:972-973
###################################################
setNA(bar, nc)


###################################################
### code chunk number 11: simFrame-intro.Rnw:1120-1129
###################################################
data("eusilcP")
sc <- SampleControl(size = 500, k = 50)
cc <- DARContControl(target = "eqIncome", epsilon = 0.02, 
    fun = function(x) x * 25)
sim <- function(x) {
    c(mean = mean(x$eqIncome), trimmed = mean(x$eqIncome, trim = 0.02))
}
set.seed(12345)
results <- runSimulation(eusilcP, sc, contControl = cc, fun = sim)


###################################################
### code chunk number 12: simFrame-intro.Rnw:1149-1153
###################################################
head(results)
aggregate(results)
tv <- mean(eusilcP$eqIncome)
tv


###################################################
### code chunk number 13: simFrame-intro.Rnw:1191-1192
###################################################
print(plot(results, true = tv))


###################################################
### code chunk number 14: simFrame-intro.Rnw:1194-1195
###################################################
print(simDensityplot(results, true = tv))


###################################################
### code chunk number 15: simFrame-intro.Rnw:1312-1315
###################################################
library("laeken")
data("eusilcP")
set.seed(12345)


###################################################
### code chunk number 16: simFrame-intro.Rnw:1322-1324
###################################################
set <- setup(eusilcP, design = "region", grouping = "hid", 
    size = c(75, 250, 250, 125, 200, 225, 125, 150, 100), k = 100)


###################################################
### code chunk number 17: simFrame-intro.Rnw:1337-1339
###################################################
cc <- DCARContControl(target = "eqIncome", epsilon = 0.005, 
    grouping = "hid", dots = list(mean = 500000, sd = 10000))


###################################################
### code chunk number 18: simFrame-intro.Rnw:1346-1356
###################################################
sim <- function(x, k) {
    g <- gini(x$eqIncome, x$.weight)$value
    eqIncHill <- fitPareto(x$eqIncome, k = k, 
        method = "thetaHill", groups = x$hid)
    gHill <- gini(eqIncHill, x$.weight)$value
    eqIncPDC <- fitPareto(x$eqIncome, k = k, 
        method = "thetaPDC", groups = x$hid)
    gPDC <- gini(eqIncPDC, x$.weight)$value
    c(standard = g, Hill = gHill, PDC = gPDC)
}


###################################################
### code chunk number 19: simFrame-intro.Rnw:1363-1365
###################################################
results <- runSimulation(eusilcP, set, contControl = cc, 
    design = "gender", fun = sim, k = 125)


###################################################
### code chunk number 20: simFrame-intro.Rnw:1372-1374
###################################################
head(results)
aggregate(results)


###################################################
### code chunk number 21: simFrame-intro.Rnw:1381-1382
###################################################
tv <- simSapply(eusilcP, "gender", function(x) gini(x$eqIncome)$value)


###################################################
### code chunk number 22: simFrame-intro.Rnw:1388-1389
###################################################
print(plot(results, true = tv, xlab = "Gini coefficient"))


###################################################
### code chunk number 23: simFrame-intro.Rnw:1454-1457
###################################################
library("robCompositions")
library("mvtnorm")
set.seed(12345)


###################################################
### code chunk number 24: simFrame-intro.Rnw:1486-1490
###################################################
crnorm <- function(n, mean, sigma) isomLRinv(rmvnorm(n, mean, sigma))
sigma <- matrix(c(1, -0.5, 1.4, -0.5, 1, -0.6, 1.4, -0.6, 2), 3, 3)
dc <- DataControl(size = 150, distribution = crnorm, 
    dots = list(mean = c(0, 2, 3), sigma = sigma))


###################################################
### code chunk number 25: simFrame-intro.Rnw:1497-1498
###################################################
nc <- NAControl(NArate = 0.05)


###################################################
### code chunk number 26: simFrame-intro.Rnw:1505-1512
###################################################
sim <- function(x, orig) {
    i <- apply(x, 1, function(x) any(is.na(x)))
    ni <- length(which(i))
    xKNNa <- impKNNa(x)$xImp
    xLS <- impCoda(x, method = "lm")$xImp
    c(knn = aDist(xKNNa, orig)/ni, LS = aDist(xLS, orig)/ni)
}


###################################################
### code chunk number 27: simFrame-intro.Rnw:1517-1518
###################################################
results <- runSimulation(dc, nrep = 50, NAControl = nc, fun = sim)


###################################################
### code chunk number 28: simFrame-intro.Rnw:1524-1526
###################################################
head(results)
aggregate(results)


###################################################
### code chunk number 29: simFrame-intro.Rnw:1533-1534
###################################################
print(plot(results, xlab = "Relative Aitchison distance"))


###################################################
### code chunk number 30: simFrame-intro.Rnw:1536-1538
###################################################
alpha <- if(names(dev.cur()) == "pdf") 0.6 else 1
print(simDensityplot(results, alpha = alpha, xlab = "Relative Aitchison distance"))


###################################################
### code chunk number 31: simFrame-intro.Rnw:1578-1579
###################################################
cl <- makeCluster(2, type="PSOCK")


###################################################
### code chunk number 32: simFrame-intro.Rnw:1585-1590
###################################################
clusterEvalQ(cl, {
        library("simFrame")
        library("robCompositions")
        library("mvtnorm")
    })


###################################################
### code chunk number 33: simFrame-intro.Rnw:1595-1596
###################################################
clusterSetRNGStream(cl, iseed=12345)


###################################################
### code chunk number 34: simFrame-intro.Rnw:1604-1616
###################################################
crnorm <- function(n, mean, sigma) isomLRinv(rmvnorm(n, mean, sigma))
sigma <- matrix(c(1, -0.5, 1.4, -0.5, 1, -0.6, 1.4, -0.6, 2), 3, 3)
dc <- DataControl(size = 150, distribution = crnorm, 
    dots = list(mean = c(0, 2, 3), sigma = sigma))
nc <- NAControl(NArate = c(0.01, 0.03, 0.05, 0.07, 0.09))
sim <- function(x, orig) {
    i <- apply(x, 1, function(x) any(is.na(x)))
    ni <- length(which(i))
    xKNNa <- impKNNa(x)$xImp
    xLS <- impCoda(x, method = "lm")$xImp
    c(knn = aDist(xKNNa, orig)/ni, LS = aDist(xLS, orig)/ni)
}


###################################################
### code chunk number 35: simFrame-intro.Rnw:1624-1625
###################################################
clusterExport(cl, c("crnorm", "sigma", "dc", "nc", "sim"))


###################################################
### code chunk number 36: simFrame-intro.Rnw:1630-1631
###################################################
results <- clusterRunSimulation(cl, dc, nrep = 50, NAControl = nc, fun = sim)


###################################################
### code chunk number 37: simFrame-intro.Rnw:1637-1638
###################################################
stopCluster(cl)


###################################################
### code chunk number 38: simFrame-intro.Rnw:1645-1647
###################################################
head(results)
aggregate(results)


###################################################
### code chunk number 39: simFrame-intro.Rnw:1653-1654
###################################################
print(plot(results, ylab = "Relative Aitchison distance"))


###################################################
### code chunk number 40: simFrame-intro.Rnw:1656-1659
###################################################
alpha <- if(names(dev.cur()) == "pdf") 0.6 else 1
print(simDensityplot(results, NArate=0.07, 
        alpha = alpha, xlab = "Relative Aitchison distance"))



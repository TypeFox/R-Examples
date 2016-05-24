### R code from vignette source 'simFrame-eusilc.Rnw'

###################################################
### code chunk number 1: simFrame-eusilc.Rnw:70-71
###################################################
options(width=72, prompt="R> ")


###################################################
### code chunk number 2: simFrame-eusilc.Rnw:131-134
###################################################
library("simFrame")
library("laeken")
data("eusilcP")


###################################################
### code chunk number 3: simFrame-eusilc.Rnw:141-152
###################################################
sim <- function(x, k) {
    x <- x[!is.na(x$eqIncome), ]
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
### code chunk number 4: simFrame-eusilc.Rnw:167-170
###################################################
set.seed(12345)
sc <- SampleControl(grouping = "hid", size = 1500, k = 100)
results <- runSimulation(eusilcP, sc, fun = sim, k = 175)


###################################################
### code chunk number 5: simFrame-eusilc.Rnw:179-183
###################################################
head(results)
aggregate(results)
tv <- gini(eusilcP$eqIncome)$value
plot(results, true = tv)


###################################################
### code chunk number 6: simFrame-eusilc.Rnw:198-199
###################################################
print(plot(results, true = tv))


###################################################
### code chunk number 7: simFrame-eusilc.Rnw:214-218
###################################################
set.seed(12345)
sc <- SampleControl(design = "region", grouping = "hid", 
    size = c(75, 250, 250, 125, 200, 225, 125, 150, 100), k = 100)
results <- runSimulation(eusilcP, sc, fun = sim, k = 175)


###################################################
### code chunk number 8: simFrame-eusilc.Rnw:222-226
###################################################
head(results)
aggregate(results)
tv <- gini(eusilcP$eqIncome)$value
plot(results, true = tv)


###################################################
### code chunk number 9: simFrame-eusilc.Rnw:232-233
###################################################
print(plot(results, true = tv))


###################################################
### code chunk number 10: simFrame-eusilc.Rnw:259-266
###################################################
set.seed(12345)
sc <- SampleControl(design = "region", grouping = "hid", 
    size = c(75, 250, 250, 125, 200, 225, 125, 150, 100), k = 100)
cc <- DCARContControl(target = "eqIncome", epsilon = 0.005, 
    grouping = "hid", dots = list(mean = 500000, sd = 10000))
results <- runSimulation(eusilcP, sc, 
    contControl = cc, fun = sim, k = 175)


###################################################
### code chunk number 11: simFrame-eusilc.Rnw:271-275
###################################################
head(results)
aggregate(results)
tv <- gini(eusilcP$eqIncome)$value
plot(results, true = tv)


###################################################
### code chunk number 12: simFrame-eusilc.Rnw:281-282
###################################################
print(plot(results, true = tv))


###################################################
### code chunk number 13: simFrame-eusilc.Rnw:319-326
###################################################
set.seed(12345)
sc <- SampleControl(design = "region", grouping = "hid", 
    size = c(75, 250, 250, 125, 200, 225, 125, 150, 100), k = 100)
cc <- DCARContControl(target = "eqIncome", epsilon = 0.005, 
    grouping = "hid", dots = list(mean = 500000, sd = 10000))
results <- runSimulation(eusilcP, sc, contControl = cc, 
    design = "gender", fun = sim, k = 125)


###################################################
### code chunk number 14: simFrame-eusilc.Rnw:335-339
###################################################
head(results)
aggregate(results)
tv <- simSapply(eusilcP, "gender", function(x) gini(x$eqIncome)$value)
plot(results, true = tv)


###################################################
### code chunk number 15: simFrame-eusilc.Rnw:345-346
###################################################
print(plot(results, true = tv))


###################################################
### code chunk number 16: simFrame-eusilc.Rnw:369-377
###################################################
set.seed(12345)
sc <- SampleControl(design = "region", grouping = "hid", 
    size = c(75, 250, 250, 125, 200, 225, 125, 150, 100), k = 100)
cc <- DCARContControl(target = "eqIncome", 
    epsilon = c(0, 0.0025, 0.005, 0.0075, 0.01),
    dots = list(mean = 500000, sd = 10000))
results <- runSimulation(eusilcP, sc, contControl = cc, 
    design = "gender", fun = sim, k = 125)


###################################################
### code chunk number 17: simFrame-eusilc.Rnw:382-386
###################################################
head(results)
aggregate(results)
tv <- simSapply(eusilcP, "gender", function(x) gini(x$eqIncome)$value)
plot(results, true = tv)


###################################################
### code chunk number 18: simFrame-eusilc.Rnw:392-393
###################################################
print(plot(results, true = tv))


###################################################
### code chunk number 19: simFrame-eusilc.Rnw:441-449
###################################################
set.seed(12345)
sc <- SampleControl(design = "region", grouping = "hid", 
    size = c(75, 250, 250, 125, 200, 225, 125, 150, 100), k = 50)
cc <- DCARContControl(target = "eqIncome", 
    epsilon = c(0, 0.005, 0.01), dots = list(mean = 500000, sd = 10000))
nc <- NAControl(target = "eqIncome", NArate = c(0, 0.05))
results <- runSimulation(eusilcP, sc, contControl = cc, 
    NAControl = nc, design = "gender", fun = sim, k = 125)


###################################################
### code chunk number 20: simFrame-eusilc.Rnw:456-460
###################################################
head(results)
aggregate(results)
tv <- simSapply(eusilcP, "gender", function(x) gini(x$eqIncome)$value)
plot(results, true = tv)


###################################################
### code chunk number 21: simFrame-eusilc.Rnw:467-468
###################################################
print(plot(results, true = tv))


###################################################
### code chunk number 22: simFrame-eusilc.Rnw:503-520
###################################################
cl <- makeCluster(2, type="PSOCK")
clusterEvalQ(cl, {
        library("simFrame")
        library("laeken")
        data("eusilcP")
    })
clusterSetRNGStream(cl, iseed=12345)
sc <- SampleControl(design = "region", grouping = "hid", 
    size = c(75, 250, 250, 125, 200, 225, 125, 150, 100), k = 50)
cc <- DCARContControl(target = "eqIncome", 
    epsilon = c(0, 0.005, 0.01), dots = list(mean = 500000, sd = 10000))
nc <- NAControl(target = "eqIncome", NArate = c(0, 0.05))
clusterExport(cl, c("sc", "cc", "nc", "sim"))
results <- clusterRunSimulation(cl, eusilcP, sc, 
    contControl = cc, NAControl = nc, design = "gender", 
    fun = sim, k = 125)
stopCluster(cl)


###################################################
### code chunk number 23: simFrame-eusilc.Rnw:524-528
###################################################
head(results)
aggregate(results)
tv <- simSapply(eusilcP, "gender", function(x) gini(x$eqIncome)$value)
plot(results, true = tv)


###################################################
### code chunk number 24: simFrame-eusilc.Rnw:535-536
###################################################
print(plot(results, true = tv))



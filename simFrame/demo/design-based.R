# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

## initializations
library("simFrame")
library("laeken")
data("eusilcP")
set.seed(12345)

## set up samples
set <- setup(eusilcP, design = "region", grouping = "hid", 
    size = c(75, 250, 250, 125, 200, 225, 125, 150, 100), k = 100)

## define contamination
cc <- DCARContControl(target = "eqIncome", epsilon = 0.005, 
    grouping = "hid", dots = list(mean = 500000, sd = 10000))

## define function for simulation runs
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

## run simulation
results <- runSimulation(eusilcP, set, contControl = cc, 
    design = "gender", fun = sim, k = 125)

## inspect results
head(results)
aggregate(results)

## compute true values
tv <- simSapply(eusilcP, "gender", function(x) gini(x$eqIncome)$value)

## plot results
plot(results, true = tv, xlab = "Gini coefficient")

# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

## initializations
library("simFrame")
data("eusilcP")

## control objects for sampling and contamination
sc <- SampleControl(size = 500, k = 50)
cc <- DARContControl(target = "eqIncome", epsilon = 0.02, 
    fun = function(x) x * 25)

## define function for simulation runs
sim <- function(x) {
    c(mean = mean(x$eqIncome), trimmed = mean(x$eqIncome, trim = 0.02))
}

## set seed and run simulation
set.seed(12345)
results <- runSimulation(eusilcP, sc, contControl = cc, fun = sim)

## inspect results
head(results)
aggregate(results)

## compute true values
tv <- mean(eusilcP$eqIncome)
tv

## plot results
plot(results, true = tv)
simDensityplot(results, true = tv)

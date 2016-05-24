## =============================================================================
## Import data from the grofit package (Kahm et al. 2010, J. Stat. Software,
##   doi:10.18637/jss.v033.i07).
##
## Convert its specific data structure into a generic data base form
##   with package reshape2 (Wickham, 2007, doi:10.18637/jss.v021.i12)
##   and fit several growth models with package growthrates.
##
## Author: Thomas Petzoldt, TU Dresden
## License: GPL >= 2, https://www.gnu.org/licenses/
## Please cite our work when using this package.
## =============================================================================

library("grofit")
library("reshape2")

## read data from package grofit
data(grofit.data)
data(grofit.time)

## convert data structure to generic format
m.data <- melt(grofit.data, id.vars = 1:3, measure.vars = 4:ncol(grofit.data))
m.time <- melt(grofit.time, id.vars = NULL)

## rename columns
names(m.data)[1:3] <- c("experiment", "info", "conc")
names(m.time)      <- c("variable2", "time")

## convert to "data base" format
database <- cbind(m.time, m.data)[c("experiment", "conc", "time", "value", "info")]

## plot the data
xyplot(value ~ time|experiment+conc, data = database)

## determine growth rate from maximum of 1st derivative of smoothing splines
sfit <- all_splines(value ~ time | experiment + conc, data = database, spar = .9)

## plot results in log scale
par(mfrow = c(3, 3))
plot(sfit, log = "y")

## plot results in linear scale
par(mfrow = c(3, 3))
plot(sfit)

## fit parametric model
p <- c(y0 = 0.01, mumax = 0.2, K = 0.2)

pfit <- all_growthmodels(value ~ grow_logistic(time, parms) | experiment + conc,
                         p = p, database)
par(mfrow = c(3, 3))
plot(pfit)

## easy linear fit with default parameters, "experiment" is redundant
efit <- all_easylinear(value ~ time | experiment + conc, data = database)
par(mfrow = c(3,3))
plot(efit, log = "y")

## easy linear fit with adapted settings
efit2 <- all_easylinear(value ~ time | conc, data = database, h = 10, quota = 0.95)
mfrow = c(3,3)
plot(efit2, log = "y")



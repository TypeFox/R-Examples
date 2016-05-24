pkgname <- "pastecs"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('pastecs')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("AutoD2")
### * AutoD2

flush(stderr()); flush(stdout())

### Name: AutoD2
### Title: AutoD2, CrossD2 or CenterD2 analysis of a multiple time-series
### Aliases: AutoD2 CrossD2 CenterD2
### Keywords: ts multivariate htest

### ** Examples

data(marphy)
marphy.ts <- as.ts(as.matrix(marphy[, 1:3]))
AutoD2(marphy.ts)
marphy.ts2 <- as.ts(as.matrix(marphy[, c(1, 4, 3)]))
CrossD2(marphy.ts, marphy.ts2)
# This is not identical to:
CrossD2(marphy.ts2, marphy.ts)
marphy.d2 <- CenterD2(marphy.ts, window=16)
lines(c(17, 17), c(-1, 15), col=4, lty=2)
lines(c(25, 25), c(-1, 15), col=4, lty=2)
lines(c(30, 30), c(-1, 15), col=4, lty=2)
lines(c(41, 41), c(-1, 15), col=4, lty=2)
lines(c(46, 46), c(-1, 15), col=4, lty=2)
text(c(8.5, 21, 27.5, 35, 43.5, 57), 11, labels=c("Peripheral Zone", "D1",
        "C", "Front", "D2", "Central Zone")) # Labels
time(marphy.ts)[marphy.d2$D2 > marphy.d2$chisq]



cleanEx()
nameEx("GetUnitText")
### * GetUnitText

flush(stderr()); flush(stdout())

### Name: GetUnitText
### Title: Format a nice time units for labels in graphs
### Aliases: GetUnitText
### Keywords: ts

### ** Examples

timeser <- ts(1:24, frequency=12)           # 12 observations per year
if (exists("is.R") && is.function(is.R) && is.R()) {  # We are in R
attr(timeser, "units") <- "years"           # time in years for 'ts' object
} else {                                              # We are in Splus
attr(attr(timeser, "tspar"), "units") <- "years" # Idem for Splus 'rts' object
}
GetUnitText(timeser)                        # formats unit (1/12 year=months)



cleanEx()
nameEx("abund")
### * abund

flush(stderr()); flush(stdout())

### Name: abund
### Title: Sort variables by abundance
### Aliases: abund extract.abund identify.abund lines.abund plot.abund
###   print.abund print.summary.abund summary.abund
### Keywords: multivariate

### ** Examples

data(bnr)
bnr.abd <- abund(bnr)
summary(bnr.abd)
plot(bnr.abd, dpos=c(105, 100))
bnr.abd$n <- 26
# To identify a point on the graph, use: bnr.abd$n <- identify(bnr.abd)
lines(bnr.abd)
bnr2 <- extract(bnr.abd)
names(bnr2)



cleanEx()
nameEx("buysbal")
### * buysbal

flush(stderr()); flush(stdout())

### Name: buysbal
### Title: Buys-Ballot table
### Aliases: buysbal
### Keywords: ts

### ** Examples

data(releve)
buysbal(releve$Day, releve$Melosul, frequency=4, units="days",
        datemin="21/03/1989", dateformat="d/m/Y")
buysbal(releve$Day, releve$Melosul, frequency=4, units="days",
        datemin="21/03/1989", dateformat="d/m/Y", count=TRUE)



cleanEx()
nameEx("daystoyears")
### * daystoyears

flush(stderr()); flush(stdout())

### Name: daystoyears
### Title: Convert time units from "days" to "years" or back
### Aliases: daystoyears yearstodays
### Keywords: ts

### ** Examples

# A vector with a "days" time-scale (25 values every 30 days)
A <- (1:25)*30
# Convert it to a "years" time-scale, using 23/05/2001 (d/m/Y) as first value
B <- daystoyears(A, datemin="23/05/2001", dateformat="d/m/Y")
B
# Convert it back to "days" time-scale
yearstodays(B, xmin=30)

# Here is an example showing the shift introduced, and its consequence:
C <- daystoyears(unclass(as.Date(c("1970-1-1","1971-1-1","1972-1-1","1973-1-1"),
	format = "%Y-%m-%d")))
C



cleanEx()
nameEx("decaverage")
### * decaverage

flush(stderr()); flush(stdout())

### Name: decaverage
### Title: Time series decomposition using a moving average
### Aliases: decaverage
### Keywords: ts smooth

### ** Examples

data(marbio)
ClausoB.ts <- ts(log(marbio$ClausocalanusB + 1))
ClausoB.dec <- decaverage(ClausoB.ts, order=2, times=10, sides=2, ends="fill")
plot(ClausoB.dec, col=c(1, 3, 2), xlab="stations")
# A stacked graph is more representative in this case
plot(ClausoB.dec, col=c(1, 3), xlab="stations", stack=FALSE, resid=FALSE,
        lpos=c(53, 4.3))



cleanEx()
nameEx("deccensus")
### * deccensus

flush(stderr()); flush(stdout())

### Name: deccensus
### Title: Time decomposition using the CENSUS II method
### Aliases: deccensus
### Keywords: ts smooth

### ** Examples

data(releve)
# Get regulated time series with a 'years' time-scale
rel.regy <- regul(releve$Day, releve[3:8], xmin=6, n=87, units="daystoyears",
    frequency=24, tol=2.2, methods="linear", datemin="21/03/1989", dateformat="d/m/Y")
rel.ts <- tseries(rel.regy)
# We must have complete cycles to allow using deccensus()
start(rel.ts)
end(rel.ts)
rel.ts2 <- window(rel.ts, end=c(1992,5))
rel.dec2 <- deccensus(rel.ts2[, "Melosul"], trend=TRUE)
plot(rel.dec2, col=c(1,4,3,2))



cleanEx()
nameEx("decdiff")
### * decdiff

flush(stderr()); flush(stdout())

### Name: decdiff
### Title: Time series decomposition using differences (trend elimination)
### Aliases: decdiff
### Keywords: ts smooth

### ** Examples

data(marbio)
ClausoB.ts <- ts(log(marbio$ClausocalanusB + 1))
ClausoB.dec <- decdiff(ClausoB.ts, lag=1, order=2, ends="fill")
plot(ClausoB.dec, col=c(1, 4, 2), xlab="stations")



cleanEx()
nameEx("decevf")
### * decevf

flush(stderr()); flush(stdout())

### Name: decevf
### Title: Time series decomposition using eigenvector filtering (EVF)
### Aliases: decevf
### Keywords: ts smooth

### ** Examples

data(releve)
melo.regy <- regul(releve$Day, releve$Melosul, xmin=9, n=87,
        units="daystoyears", frequency=24, tol=2.2, methods="linear",
        datemin="21/03/1989", dateformat="d/m/Y")
melo.ts <- tseries(melo.regy)
acf(melo.ts)
# Autocorrelation is not significant after 0.16
# That corresponds to a lag of 0.16*24=4 (frequency=24)
melo.evf <- decevf(melo.ts, lag=4, axes=1)
plot(melo.evf, col=c(1, 4, 2))
# A superposed graph is better in the present case
plot(melo.evf, col=c(1, 4), xlab="stations", stack=FALSE, resid=FALSE,
        lpos=c(0, 60000))



cleanEx()
nameEx("decloess")
### * decloess

flush(stderr()); flush(stdout())

### Name: decloess
### Title: Time series decomposition by the LOESS method
### Aliases: decloess
### Keywords: ts smooth

### ** Examples

data(releve)
melo.regy <- regul(releve$Day, releve$Melosul, xmin=9, n=87,
        units="daystoyears", frequency=24, tol=2.2, methods="linear",
        datemin="21/03/1989", dateformat="d/m/Y")
melo.ts <- tseries(melo.regy)
melo.dec <- decloess(melo.ts, s.window="periodic")
plot(melo.dec, col=1:3)



cleanEx()
nameEx("decmedian")
### * decmedian

flush(stderr()); flush(stdout())

### Name: decmedian
### Title: Time series decomposition using a running median
### Aliases: decmedian
### Keywords: ts smooth

### ** Examples

data(marbio)
ClausoB.ts <- ts(log(marbio$ClausocalanusB + 1))
ClausoB.dec <- decmedian(ClausoB.ts, order=2, times=10, ends="fill")
plot(ClausoB.dec, col=c(1, 4, 2), xlab="stations")
# This is a transect across a frontal zone:
plot(ClausoB.dec, col=c(0, 2), xlab="stations", stack=FALSE, resid=FALSE)
lines(c(17, 17), c(0, 10), col=4, lty=2)
lines(c(25, 25), c(0, 10), col=4, lty=2)
lines(c(30, 30), c(0, 10), col=4, lty=2)
lines(c(41, 41), c(0, 10), col=4, lty=2)
lines(c(46, 46), c(0, 10), col=4, lty=2)
text(c(8.5, 21, 27.5, 35, 43.5, 57), 8.7, labels=c("Peripheral Zone", "D1",
        "C", "Front", "D2", "Central Zone"))



cleanEx()
nameEx("decreg")
### * decreg

flush(stderr()); flush(stdout())

### Name: decreg
### Title: Time series decomposition using a regression model
### Aliases: decreg
### Keywords: ts smooth

### ** Examples

data(marphy)
density <- ts(marphy[, "Density"])
plot(density)
Time <- time(density)

# Linear model to represent trend
density.lin <- lm(density ~ Time)
summary(density.lin)
xreg <- predict(density.lin)
lines(xreg, col=3)
density.dec <- decreg(density, xreg)
plot(density.dec, col=c(1, 3, 2), xlab="stations")

# Order 2 polynomial to represent trend
density.poly <- lm(density ~ Time + I(Time^2))
summary(density.poly)
xreg2 <- predict(density.poly)
plot(density)
lines(xreg2, col=3)
density.dec2 <- decreg(density, xreg2)
plot(density.dec2, col=c(1, 3, 2), xlab="stations")

# Fit a sinusoidal model on seasonal (artificial) data
tser <- ts(sin((1:100)/12*pi)+rnorm(100, sd=0.3), start=c(1998, 4),
        frequency=24)
Time <- time(tser)
tser.sin <- lm(tser ~ I(cos(2*pi*Time)) + I(sin(2*pi*Time)))
summary(tser.sin)
tser.reg <- predict(tser.sin)
tser.dec <- decreg(tser, tser.reg)
plot(tser.dec, col=c(1, 4), xlab="stations", stack=FALSE, resid=FALSE,
        lpos=c(0, 4))
plot(tser.dec, col=c(1, 4, 2), xlab="stations")

# One can also use nonlinear models (see 'nls')
# or autoregressive models (see 'ar' and others in 'ts' library)



cleanEx()
nameEx("disjoin")
### * disjoin

flush(stderr()); flush(stdout())

### Name: disjoin
### Title: Complete disjoined coded data (binary coding)
### Aliases: disjoin
### Keywords: manip

### ** Examples

# Artificial data with 1/5 of zeros
Z <- c(abs(rnorm(8000)), rep(0, 2000))
# Let the program chose cuts
table(cut(Z, breaks=5))
# Create one class for zeros, and 4 classes for the other observations
Z2 <- Z[Z != 0]
cuts <- c(-1e-10, 1e-10, quantile(Z2, 1:5/5, na.rm=TRUE))
cuts
table(cut(Z, breaks=cuts))
# Binary coding of these data
disjoin(cut(Z, breaks=cuts))[1:10, ]



cleanEx()
nameEx("disto")
### * disto

flush(stderr()); flush(stdout())

### Name: disto
### Title: Compute and plot a distogram
### Aliases: disto
### Keywords: multivariate ts

### ** Examples

data(bnr)
disto(bnr)



cleanEx()
nameEx("escouf")
### * escouf

flush(stderr()); flush(stdout())

### Name: escouf
### Title: Choose variables using the Escoufier's equivalent vectors method
### Aliases: escouf extract.escouf identify.escouf lines.escouf plot.escouf
###   print.escouf print.summary.escouf summary.escouf
### Keywords: multivariate

### ** Examples

data(marbio)
marbio.esc <- escouf(marbio)
summary(marbio.esc)
plot(marbio.esc)
# The x-axis has short labels. For more info., enter: 
marbio.esc$vr
# Define a level at which to extract most significant variables
marbio.esc$level <- 0.90
# Show it on the graph
lines(marbio.esc)
# This can also be done interactively on the plot using:
# marbio.esc$level <- identify(marbio.esc)
# Finally, extract most significant variables
marbio2 <- extract(marbio.esc)
names(marbio2)



cleanEx()
nameEx("first")
### * first

flush(stderr()); flush(stdout())

### Name: first
### Title: Get the first element of a vector
### Aliases: first
### Keywords: manip

### ** Examples

a <- c(NA, 1, 2, NA, 3, 4, NA)
first(a)
first(a, na.rm=TRUE)



cleanEx()
nameEx("is.tseries")
### * is.tseries

flush(stderr()); flush(stdout())

### Name: is.tseries
### Title: Is this object a time series?
### Aliases: is.tseries
### Keywords: ts

### ** Examples

tser <- ts(sin((1:100)/6*pi)+rnorm(100, sd=0.5), start=c(1998, 4), frequency=12)
is.tseries(tser)      # TRUE
not.a.ts <- c(1,2,3)
is.tseries(not.a.ts)  # FALSE



cleanEx()
nameEx("last")
### * last

flush(stderr()); flush(stdout())

### Name: last
### Title: Get the last element of a vector
### Aliases: last
### Keywords: manip

### ** Examples

a <- c(NA, 1, 2, NA, 3, 4, NA)
last(a)
last(a, na.rm=TRUE)



cleanEx()
nameEx("local.trend")
### * local.trend

flush(stderr()); flush(stdout())

### Name: local.trend
### Title: Calculate local trends using cumsum
### Aliases: local.trend identify.local.trend
### Keywords: ts

### ** Examples

data(bnr)
# Calculate and plot cumsum for the 8th series
bnr8.lt <- local.trend(bnr[,8])
# To identify local trends, use:
# identify(bnr8.lt)
# and click points between which you want to compute local linear trends...



cleanEx()
nameEx("match.tol")
### * match.tol

flush(stderr()); flush(stdout())

### Name: match.tol
### Title: Determine matching observation with a tolerance in time-scale
### Aliases: match.tol
### Keywords: ts manip

### ** Examples

X <- 1:5
Table <- c(1, 3.1, 3.8, 4.4, 5.1, 6)
match.tol(X, Table)
match.tol(X, Table, tol=0.2)
match.tol(X, Table, tol=0.55)



cleanEx()
nameEx("pennington")
### * pennington

flush(stderr()); flush(stdout())

### Name: pennington
### Title: Calculate Pennington statistics
### Aliases: pennington
### Keywords: univar

### ** Examples

data(marbio)
pennington(marbio[, "Copepodits2"])
pennington(marbio[, "Copepodits2"], calc="mean", na.rm=TRUE)



cleanEx()
nameEx("pgleissberg")
### * pgleissberg

flush(stderr()); flush(stdout())

### Name: pgleissberg
### Title: Gleissberg distribution probability
### Aliases: pgleissberg
### Keywords: distribution

### ** Examples

# Until n=50, the exact probability is returned
pgleissberg(20, 10, lower.tail=TRUE, two.tailed=FALSE)
# For higher n values, it is approximated by a normal distribution
pgleissberg(60, 33, lower.tail=TRUE, two.tailed=FALSE)



cleanEx()
nameEx("regarea")
### * regarea

flush(stderr()); flush(stdout())

### Name: regarea
### Title: Regulate a series using the area method
### Aliases: regarea
### Keywords: ts manip chron smooth

### ** Examples

data(releve)
# The 'Melosul' series is regulated with a 25-days window
reg <- regarea(releve$Day, releve$Melosul, window=25)
# Then with a 50-days window
reg2 <- regarea(releve$Day, releve$Melosul, window=50)
# The initial and both regulated series are shown on the graph for comparison
plot(releve$Day, releve$Melosul, type="l")
lines(reg$x, reg$y, col=2)
lines(reg2$x, reg2$y, col=4)



cleanEx()
nameEx("regconst")
### * regconst

flush(stderr()); flush(stdout())

### Name: regconst
### Title: Regulate a series using the constant value method
### Aliases: regconst
### Keywords: ts manip chron smooth

### ** Examples

data(releve)
reg <- regconst(releve$Day, releve$Melosul)
plot(releve$Day, releve$Melosul, type="l")
lines(reg$x, reg$y, col=2)



cleanEx()
nameEx("reglin")
### * reglin

flush(stderr()); flush(stdout())

### Name: reglin
### Title: Regulation of a series using a linear interpolation
### Aliases: reglin
### Keywords: ts manip chron smooth

### ** Examples

data(releve)
reg <- reglin(releve$Day, releve$Melosul)
plot(releve$Day, releve$Melosul, type="l")
lines(reg$x, reg$y, col=2)



cleanEx()
nameEx("regspline")
### * regspline

flush(stderr()); flush(stdout())

### Name: regspline
### Title: Regulation of a time series using splines
### Aliases: regspline
### Keywords: ts manip chron smooth

### ** Examples

data(releve)
reg <- regspline(releve$Day, releve$Melosul)
plot(releve$Day, releve$Melosul, type="l")
lines(reg$x, reg$y, col=2)



cleanEx()
nameEx("regul")
### * regul

flush(stderr()); flush(stdout())

### Name: regul
### Title: Regulation of one or several time series using various methods
### Aliases: regul extract.regul hist.regul identify.regul lines.regul
###   plot.regul print.regul print.specs.regul print.summary.regul
###   specs.regul summary.regul
### Keywords: ts manip chron smooth

### ** Examples

data(releve)
# The series in this data frame are very irregularly sampled in time:
releve$Day
length(releve$Day)
intervals <- releve$Day[2:61]-releve$Day[1:60]
intervals
range(intervals)
mean(intervals)
# The series must be regulated to be converted in a 'rts' or 'ts object
rel.reg <- regul(releve$Day, releve[3:8], xmin=9, n=63, deltat=21,
        tol=1.05, methods=c("s","c","l","a","s","a"), window=21)
rel.reg
plot(rel.reg, 5)
specs(rel.reg)
# Now we can extract one or several regular time series
melo.ts <- extract(rel.reg, series="Melosul")
is.tseries(melo.ts)

# One can convert time-scale from "days" to "years" during regulation
# This is most useful for analyzing seasonal cycles in a second step
melo.regy <- regul(releve$Day, releve$Melosul, xmin=6, n=87,
        units="daystoyears", frequency=24, tol=2.2, methods="linear",
        datemin="21/03/1989", dateformat="d/m/Y")
melo.regy
plot(melo.regy, main="Regulation of Melosul")
# In this case, we have only one series in 'melo.regy'
# We can use also 'tseries()' instead of 'extract()'
melo.tsy <- tseries(melo.regy)
is.tseries(melo.tsy)



cleanEx()
nameEx("regul.adj")
### * regul.adj

flush(stderr()); flush(stdout())

### Name: regul.adj
### Title: Adjust regulation parameters
### Aliases: regul.adj
### Keywords: ts chron

### ** Examples

# This example follows the example for regul.screen()
# where we determined that xmin=9, deltat=21, n=63, with tol=1.05
# is a good choice to regulate the irregular time series in 'releve' 
data(releve)
regul.adj(releve$Day, xmin=9, deltat=21)
# The histogram indicates that it is not useful to increase tol
# more than 1.05, because few observations will be added
# except if we increase it to 5-7, but this value could be
# considered to be too large in comparison with deltat=22
# On the other hand, with tol <= 1, the number of matching
# observations will be almost divided by two!



cleanEx()
nameEx("regul.screen")
### * regul.screen

flush(stderr()); flush(stdout())

### Name: regul.screen
### Title: Test various regulation parameters
### Aliases: regul.screen
### Keywords: ts chron

### ** Examples

data(releve)
# This series is very irregular, and it is difficult
# to choose the best regular time-scale
releve$Day
length(releve$Day)
intervals <- releve$Day[2:61]-releve$Day[1:60]
intervals
range(intervals)
mean(intervals)
# A combination of xmin=1, deltat=22 and n=61 seems correct
# But is it the best one?
regul.screen(releve$Day, xmin=0:11, deltat=16:27, tol=1.05)
# Now we can tell that xmin=9, deltat=21, n=63, with tol=1.05
# is a much better choice! 



cleanEx()
nameEx("stat.desc")
### * stat.desc

flush(stderr()); flush(stdout())

### Name: stat.desc
### Title: Descriptive statistics on a data frame or time series
### Aliases: stat.desc
### Keywords: ts

### ** Examples

data(marbio)
stat.desc(marbio[,13:16], basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)



cleanEx()
nameEx("stat.pen")
### * stat.pen

flush(stderr()); flush(stdout())

### Name: stat.pen
### Title: Pennington statistics on a data frame or time series
### Aliases: stat.pen
### Keywords: ts

### ** Examples

data(marbio)
stat.pen(marbio[,c(4, 14:16)], basic=TRUE, desc=TRUE)



cleanEx()
nameEx("stat.slide")
### * stat.slide

flush(stderr()); flush(stdout())

### Name: stat.slide
### Title: Sliding statistics
### Aliases: stat.slide lines.stat.slide plot.stat.slide print.stat.slide
### Keywords: ts

### ** Examples

data(marbio)
# Sliding statistics with fixed-length blocs
statsl <- stat.slide(1:68, marbio[, "ClausocalanusA"], xmin=0, n=7, deltat=10)
statsl
plot(statsl, stat="mean", leg=TRUE, lpos=c(55, 2500), xlab="Station",
        ylab="ClausocalanusA")

# More information on the series, with predefined blocs
statsl2 <- stat.slide(1:68, marbio[, "ClausocalanusA"],
        xcut=c(0, 17, 25, 30, 41, 46, 70), basic=TRUE, desc=TRUE, norm=TRUE,
        pen=TRUE, p=0.95)
statsl2
plot(statsl2, stat="median", xlab="Stations", ylab="Counts",
        main="Clausocalanus A")              # Median
lines(statsl2, stat="min")                   # Minimum
lines(statsl2, stat="max")                   # Maximum
lines(c(17, 17), c(-50, 2600), col=4, lty=2) # Cuts
lines(c(25, 25), c(-50, 2600), col=4, lty=2)
lines(c(30, 30), c(-50, 2600), col=4, lty=2)
lines(c(41, 41), c(-50, 2600), col=4, lty=2)
lines(c(46, 46), c(-50, 2600), col=4, lty=2)
text(c(8.5, 21, 27.5, 35, 43.5, 57), 2300, labels=c("Peripheral Zone", "D1",
        "C", "Front", "D2", "Central Zone")) # Labels
legend(0, 1900, c("series", "median", "range"), col=1:3, lty=1)
# Get cuts back from the object
statsl2$xcut



cleanEx()
nameEx("trend.test")
### * trend.test

flush(stderr()); flush(stdout())

### Name: trend.test
### Title: Test if an increasing or decreasing trend exists in a time
###   series
### Aliases: trend.test
### Keywords: ts htest

### ** Examples

data(marbio)
trend.test(marbio[, 8])
# Run a bootstrap test on the same series
marbio8.trend.test <- trend.test(marbio[, 8], R=99)
# R=999 is a better value... but it is very slow!
marbio8.trend.test  
plot(marbio8.trend.test)
marbio8.trend.test$p.value



cleanEx()
nameEx("tsd")
### * tsd

flush(stderr()); flush(stdout())

### Name: tsd
### Title: Decomposition of one or several regular time series using
###   various methods
### Aliases: tsd extract.tsd plot.tsd print.specs.tsd print.summary.tsd
###   print.tsd specs.tsd summary.tsd
### Keywords: ts smooth loess nonparametric

### ** Examples

data(releve)
# Regulate the series and extract them as a time series object
rel.regy <- regul(releve$Day, releve[3:8], xmin=6, n=87, units="daystoyears",
        frequency=24, tol=2.2, methods="linear", datemin="21/03/1989",
        dateformat="d/m/Y")
rel.ts <- tseries(rel.regy)

# Decompose all series in the set with the "loess" method
rel.dec <- tsd(rel.ts, method="loess", s.window=13, trend=FALSE)
rel.dec
plot(rel.dec, series=5, col=1:3)    # An plot series 5

# Extract "deseasoned" components
rel.des <- extract(rel.dec, series=3:6, components="deseasoned")
rel.des[1:10,]

# Further decompose these components with a moving average
rel.des.dec <- tsd(rel.des, method="average", order=2, times=10)
plot(rel.des.dec, series=3, col=c(2, 4, 6))
# In this case, a superposed graph is more appropriate:
plot(rel.des.dec, series=3, col=c(2,4), stack=FALSE, resid=FALSE,
        labels=c("without season cycle", "trend"), lpos=c(0, 55000))
# Extract residuals from the latter decomposition
rel.res2 <- extract(rel.des.dec, components="residuals")



cleanEx()
nameEx("tseries")
### * tseries

flush(stderr()); flush(stdout())

### Name: tseries
### Title: Convert a 'regul' or a 'tsd' object into a time series
### Aliases: tseries
### Keywords: ts manip

### ** Examples

data(releve)
rel.regy <- regul(releve$Day, releve[3:8], xmin=6, n=87, units="daystoyears",
        frequency=24, tol=2.2, methods="linear", datemin="21/03/1989",
        dateformat="d/m/Y")
# This object is not a time series
is.tseries(rel.regy)     # FALSE
# Extract all time series contained in the 'regul' object
rel.ts <- tseries(rel.regy)
# Now this is a time series
is.tseries(rel.ts)       # TRUE



cleanEx()
nameEx("turnogram")
### * turnogram

flush(stderr()); flush(stdout())

### Name: turnogram
### Title: Calculate and plot a turnogram for a regular time series
### Aliases: turnogram extract.turnogram identify.turnogram plot.turnogram
###   print.summary.turnogram print.turnogram summary.turnogram
### Keywords: ts htest

### ** Examples

data(bnr)
# Let's transform series 4 into a time series (supposing it is regular)
bnr4 <- as.ts(bnr[, 4])
plot(bnr4, type="l", main="bnr4: raw data", xlab="Time")
# A simple turnogram is calculated
bnr4.turno <- turnogram(bnr4)
summary(bnr4.turno)
# A complete turnogram confirms that "level=3" is a good value: 
turnogram(bnr4, complete=TRUE)
# Data with maximum info. are extracted (thus taking 1 every 3 observations)
bnr4.interv3 <- extract(bnr4.turno)
plot(bnr4, type="l", lty=2, xlab="Time")
lines(bnr4.interv3, col=2)
title("Original bnr4 (dotted) versus max. info. curve (plain)")
# Choose another level (for instance, 6) and extract the corresponding series
bnr4.turno$level <- 6
bnr4.interv6 <- extract(bnr4.turno)
# plot both extracted series on top of the original one
plot(bnr4, type="l", lty=2, xlab="Time")
lines(bnr4.interv3, col=2)
lines(bnr4.interv6, col=3)
legend(70, 580, c("original", "interval=3", "interval=6"), col=1:3, lty=c(2, 1, 1))
# It is hard to tell on the graph which series contains more information
# The turnogram shows us that it is the "interval=3" one!



cleanEx()
nameEx("turnpoints")
### * turnpoints

flush(stderr()); flush(stdout())

### Name: turnpoints
### Title: Analyze turning points (peaks or pits)
### Aliases: turnpoints extract.turnpoints lines.turnpoints plot.turnpoints
###   print.summary.turnpoints print.turnpoints summary.turnpoints
### Keywords: ts

### ** Examples

data(marbio)
plot(marbio[, "Nauplii"], type="l")
# Calculate turning points for this series
Nauplii.tp <- turnpoints(marbio[, "Nauplii"])
summary(Nauplii.tp)
plot(Nauplii.tp)
# Add envelope and median line to original data
plot(marbio[, "Nauplii"], type="l")
lines(Nauplii.tp)
# Note that lines() applies to the graph of original dataset!!!
title("Raw data, envelope maxi., mini. and median line")



cleanEx()
nameEx("vario")
### * vario

flush(stderr()); flush(stdout())

### Name: vario
### Title: Compute and plot a semi-variogram
### Aliases: vario
### Keywords: ts

### ** Examples

data(bnr)
vario(bnr[, 4])



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

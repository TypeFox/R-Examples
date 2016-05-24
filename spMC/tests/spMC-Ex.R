pkgname <- "spMC"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('spMC')
set.seed(0)

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("All-spMC")
flush(stderr()); flush(stdout())

### Name: setCores
### Title: Set the number of CPU cores for HPC
### Aliases: setCores
### Keywords: programming

### ** Examples

#Set 2 CPU cores for parallel computation
setCores(2)

### Name: ACM
### Title: ACM Data
### Aliases: ACM
### Keywords: datasets

### ** Examples

data(ACM)
str(ACM)
summary(ACM)

### Name: boxplot.lengths
### Title: Stratum Lengths Boxplot
### Aliases: boxplot.lengths
### Keywords: spatial distribution hplot

### ** Examples

direction <- c(0,0,1)

# Compute the appertaining directional line for each location
loc.id <- which_lines(ACM[, 1:3], direction)

# Estimate stratum lengths
gl <- getlen(ACM$MAT3, ACM[, 1:3], loc.id, direction)

# Make the boxplot of the object gl
par(mfrow = c(1, 1))
boxplot(gl)
graphics::par(get("par.postscript", pos = 'CheckExEnv'))

### Name: sim_ck
### Title: Conditional Simulation Based on Indicator Cokriging
### Aliases: sim_ck
### Keywords: spatial distribution

### ** Examples

data(ACM)
# Estimate transition rates matrices and 
# proportions for the categorical variable MAT5
x <- multi_tpfit(ACM$MAT5, ACM[, 1:3])

# Generate the simulation grid
mygrid <- list()
mygrid$X <- seq(min(ACM$X), max(ACM$X), length = 3)
mygrid$Y <- seq(min(ACM$Y), max(ACM$Y), length = 3)
mygrid$Z <- -180 * 0:2 - 1
mygrid <- as.matrix(expand.grid(mygrid$X, mygrid$Y, mygrid$Z))

# Simulate the random field through
# Simple Indicator Cokriging algorithm and
# optimize by Simulated Annealing
myANSimP <- sim_ck(x, ACM$MAT5, ACM[, 1:3], mygrid)

# Simulate the random field through
# Ordinary Indicator Cokriging algorithm and
# optimize by Genetic Algorithm
myGASimS <- sim_ck(x, ACM$MAT5, ACM[, 1:3], mygrid)

### Name: density.lengths
### Title: Empirical Densities Estimation of Stratum Lengths
### Aliases: density.lengths
### Keywords: spatial distribution

### ** Examples

# Compute the empirical densities of stratum lengths
dgl <- density(gl)

### Name: embed_MC
### Tested in multi_tpfit() during the boxplot.lengths() example

### Name: getlen
### Tested during the boxplot.lengths() example

### Name: hist.lengths
### Title: Histograms of Stratum Lengths for Each Observed Category
### Aliases: hist.lengths
### Keywords: spatial distribution hplot

### ** Examples

# Plot the histograms
hist(gl)

### Name: sim_ik
### Title: Conditional Simulation Based on Indicator Kriging
### Aliases: sim_ik
### Keywords: spatial distribution

### ** Examples

# Simulate the random field through
# Simple Indicator Kriging algorithm and
# optimize by Simulated Annealing
myANSimP <- sim_ik(x, ACM$MAT5, ACM[, 1:3], mygrid, ordinary = FALSE)
myANSimP <- quench(x, ACM$MAT5, ACM[, 1:3], myANSimP, optype = "param", max.it = 2)
myANSimF <- quench(x, ACM$MAT5, ACM[, 1:3], myANSimP, optype = "fullprobs", max.it = 2)

# Simulate the random field through
# Ordinary Indicator Kriging algorithm and
# optimize by Genetic Algorithm
myGASimS <- sim_ik(x, ACM$MAT5, ACM[, 1:3], mygrid)
myGASimS <- quench(x, ACM$MAT5, ACM[, 1:3], myGASimS, GA = TRUE, optype = "semiprobs", max.it = 2)
myGASimC <- quench(x, ACM$MAT5, ACM[, 1:3], myGASimS, GA = TRUE, optype = "coordprobs", max.it = 2)

### Name: tpfit_ils
### It will be tested during the multi_tpfit_ils() example

### Name: image.multi_tpfit
### Title: Images with Multidimensional Transiograms
### Aliases: image.multi_tpfit
### Keywords: spatial distribution hplot

### ** Examples

# Set short names for categories 3 and 4
names(x$prop)[3:4] <- c("Clay and Sand", "Gravel and Sand")

# Plot 2-D theoretical sections of
# a multidimensional transiogram
image(x, 3, max.dist=c(20,10,5), which.dire=2:3,
      mar = .7, col=rev(heat.colors(500)),
      breaks=0:500/500, nlevels = 5)

# 3D-Plot for a 2-D theoretical sections of
# a multidimensional transiogram
persp(x, 3, max.dist=c(20,10,5), which.dire=2:3,
      mar = .7, col=rev(heat.colors(500)))

### Name: image.pemt (is.pemt)
### Title: Images with Pseudo-empirical Multidimensional Transiograms
### Aliases: imgMultiTransiogram
### Keywords: spatial distribution hplot

### ** Examples

# Compute a 2-D section of a pseudo-empirical
# multidimensional transiogram
psEmpTr <- pemt(ACM$MAT3, ACM[, 1:3], 2,
                max.dist = c(20, 10, 5), 
                which.dire=c(1, 3), 
                mle = TRUE)

# Contour plots of 2-D pseudo-empirical 
# sections of multidimensional transiograms
contour(psEmpTr, mar = .7)

# Plot 2-D pseudo-empirical sections of
# a multidimensional transiogram
image(psEmpTr, col = rev(heat.colors(500)), 
      breaks = 0:500 / 500, mar = .7)

# 3D-Plot for a 2-D pseudo-empirical sections of
# a multidimensional transiogram
persp(psEmpTr, col = rainbow(500), mar = .7,
      theta = 15, phi = 45)

# Test the object psEmpTr
is.pemt(psEmpTr)

### Name: is.lengths
### Title: Object test for lengths class
### Aliases: is.lengths
### Keywords: spatial attribute

### ** Examples

# Test the object gl
is.lengths(gl)

### Name: is.multi_tpfit
### Title: Object test for multi_tpfit class
### Aliases: is.multi_tpfit
### Keywords: spatial attribute

### ** Examples

# Test the object x
is.multi_tpfit(x)

### Name: is.multi_transiogram
### Title: Object test for multi_transiogram class
### Aliases: is.multi_transiogram
### Keywords: spatial attribute

### ** Examples

# Generate the matrix of 
# multidimensional lags
lags <- expand.grid(X=-1:1, Y=-1:1, Z=-1:1)
lags <- as.matrix(lags)

# Compute transition probabilities 
# from the multidimensional MC model
TrPr <- predict(x, lags)

# Test the object TrPr
is.multi_transiogram(TrPr)

### Name: is.tpfit
### Title: Object test for tpfit class
### Aliases: is.tpfit
### Keywords: spatial attribute

### ** Examples

# Estimate the parameters of a 
# one-dimensional MC model
MoPa <- tpfit(ACM$PERM, ACM[, 1:3], c(0, 0, 1))

# Test the object MoPa
is.tpfit(MoPa)

### Name: is.transiogram
### Title: Object test for transiogram class
### Aliases: is.transiogram
### Keywords: spatial attribute

### ** Examples

# Compute theoretical transition probabilities 
# from the one-dimensional MC model
TTPr <- predict(MoPa, lags = 0:2/2)

# Compute empirical transition probabilities 
ETPr <- transiogram(ACM$MAT5, ACM[, 1:3], c(0, 0, 1), 200, 20, reverse = TRUE)

# Test the objects TTPr and ETPr
is.transiogram(TTPr)
is.transiogram(ETPr)

### Name: sim_mcs
### Title: Multinomial Categorical Simulation
### Aliases: sim_mcs
### Keywords: spatial distribution

### ** Examples

# Simulate the random field
myMCSim <- sim_mcs(x, ACM$MAT5, ACM[, 1:3], mygrid)

### Name: metpfit
### It will be tested during the multi_tpfit_me() example

### Name: mixplot
### Title: Plot of Multiple One-dimensional Transiograms
### Aliases: mixplot
### Keywords: spatial distribution hplot

### ** Examples

# Estimate empirical transition 
# probabilities by points
ETr <- transiogram(ACM$MAT3, ACM[, 1:3], c(0, 0, 1), 100)

# Estimate the transition rate matrix
RTm <- tpfit(ACM$MAT3, ACM[, 1:3], c(0, 0, 1))

# Compute transition probabilities 
# from the one-dimensional MC model
TPr <- predict(RTm, lags = ETr$lags)

# Plot empirical vs. theoretical transition probabilities
mixplot(list(ETr, TPr), type = c("p", "l"), pch = "+", col = c(3, 1))

### Name: mlen
### Tested in multi_tpfit() during the boxplot.lengths() example

### Name: multi_tpfit_ils
### Title: Iterated Least Squares Method for Multidimensional Model
###   Parameters Estimation
### Aliases: multi_tpfit_ils
### Keywords: spatial distribution models

### ** Examples

# Estimate the parameters of a 
# multidimensional MC model
multi_tpfit_ils(ACM$MAT3, ACM[, 1:3], 200)

### Name: multi_tpfit_me
### Title: Maximum Entropy Method for Multidimensional Model Parameters
###   Estimation
### Aliases: multi_tpfit_me
### Keywords: spatial distribution models

### ** Examples

# Estimate transition rates matrices and
# proportions for the categorical variable MAT3
multi_tpfit_me(ACM$MAT3, ACM[, 1:3], mle = "mdn")

### Name: multi_tpfit
### Tested during the sim_ck() example

### Name: sim_path
### Title: Conditional Simulation Based on Path Algorithms
### Aliases: sim_path
### Keywords: spatial distribution

### ** Examples

# Simulate the random field through
# the fixed path algorithm
myFixPathSim <- sim_path(x, ACM$MAT5, ACM[, 1:3], mygrid,
                         radius = 50, fixed = TRUE)

# Simulate the random field through
# the random path algorithm
myRndPathSim <- sim_path(x, ACM$MAT5, ACM[, 1:3], mygrid, radius = 50)

### Name: plot.density.lengths
### Title: Plot Empirical Densities Estimates of Stratum Lengths
### Aliases: plot.density.lengths
### Keywords: spatial distribution hplot

### ** Examples

# Plot the empirical densities of stratum log-lengths
plot(dgl)

### Name: plot.hist.lengths
### Title: Plot Histograms of Stratum Lengths
### Aliases: plot.hist.lengths
### Keywords: spatial distribution hplot

### ** Examples

# Compute the histograms
hgl <- hist(gl, plot = FALSE)

# Plot the histograms
plot(hgl, col = "#efffef")

### Name: plot.lengths
### Title: Plot Stratum Lengths
### Aliases: plot.lengths
### Keywords: spatial distribution hplot

### ** Examples

# Plot the object gl
par(mfrow = c(1,1))
plot(gl)
graphics::par(get("par.postscript", pos = 'CheckExEnv'))

### Name: plot.transiogram
### Title: Plot One-dimensional Transiograms
### Aliases: plot.transiogram
### Keywords: spatial distribution hplot

### ** Examples

# Plot empirical transition probabilities with confidence interval
plot(ETr, type = "l", ci = .95)

# Plot theoretical transition probabilities
plot(TPr, type = "l")

### Name: predict.multi_tpfit
### Tested during the is.multi_transiogram() example

### Name: predict.tpfit
### Tested during the is.transiogram() example

### Name: print.density.lengths
### Title: Printing Empirical Densities Estimates of Stratum Lengths
### Aliases: print.density.lengths
### Keywords: spatial distribution

### ** Examples

# Print the empirical densities of stratum lengths
print(dgl)

### Name: print.lengths
### Title: Printing Stratum Lengths for Each Observed Category
### Aliases: print.lengths
### Keywords: spatial

### ** Examples

# Print stratum lengths
print(gl)

### Name: print.multi_tpfit
### Title: Printing Model Parameters for Multidimensional Continuos Lag
###   Spatial MC
### Aliases: print.multi_tpfit
### Keywords: spatial

### ** Examples

# Print results
print(x)

### Name: print.multi_transiogram
### Title: Printing Theoretical Multidimensional Transiograms
### Aliases: print.multi_transiogram
### Keywords: spatial

### ** Examples

# Print results
print(TrPr)

### Name: print.summary.lengths
### Title: Printing Stratum Lengths Summary for Each Observed Category
### Aliases: print.summary.lengths
### Keywords: spatial distribution

### ** Examples

# Summarize the stratum lengths
sgl <- summary(gl)

# Print the summary of stratum lengths
print(sgl)

### Name: print.tpfit
### Title: Printing Model Parameters for One-dimensional Continuos Lag
###   Spatial MC
### Aliases: print.tpfit
### Keywords: spatial

### ** Examples

# Print results
print(MoPa)

### Name: print.transiogram
### Title: Printing Theoretical or Empirical One-dimensional Transiograms
### Aliases: print.transiogram
### Keywords: spatial

### ** Examples

# Print results
print(TTPr)
print(ETPr)

rownames(ETPr$Tmat) <- c("Clay", "Gravel", "Sand and Clay", "Sand and Gravel", "Sand")
colnames(ETPr$Tmat) <- rownames(ETPr$Tmat)
plot(ETPr, type = "l")


### Name: summary.lengths
### Tested during the print.summary.lengths() example

### Name: tpfit
### Tested in multi_tpfit() during the boxplot.lengths() example

### Name: transiogram
### Tested during the is.transiogram() example

### Name: which_lines
### Tested during the boxplot.lengths() example

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

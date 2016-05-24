###################################################
### Chap02Start
###################################################
library(mistat)
data(YARNSTRG)
data(CYCLT)  


###################################################
### LibraryDataPrint
###################################################
# This is a comment
install.packages("mistat",          # Install mistat package 
                 dependencies=TRUE) # and its dependencies
                 #
library(mistat)  # A command to make our datasets
                 # and functions available
                 #
data(CYCLT)      # Load specified data set
                 # CYCLT is a vector of values
                 #
help(CYCLT)      # Read the help page about CYCLT
                 #
CYCLT            # Print CYCLT to Console


###################################################
### Seed123
###################################################
set.seed(123)


###################################################
### rbinom 50
###################################################
X <-                  # Assign to object X
  rbinom(n = 50,      # 50 pairs of binomial variates 
         size = 1,    # set the number of trials
         prob = 0.5)  # set the probability of success
                      #
X                     # Equivalent to command print(X)
                      #
ls()                  # List the available objects
                      #
rm(X)                 # Remove object X


###################################################
### PlotSteelrods50
###################################################
data(STEELROD)                  # STEELROD is a vector
#
plot(STEELROD,                  # Plot vector STEELROD 
     ylab = "Steel rod Length", # set y axis title 
     xlab = "Index")            # set x axis title


###################################################
### Steelrods50
###################################################
data(STEELROD)                  # STEELROD is a vector
                                #
plot(STEELROD,                  # Plot vector STEELROD 
     ylab = "Steel rod Length", # set y axis title 
     xlab = "Index")            # set x axis title


###################################################
### STEELRODshift
###################################################
X <- STEELROD
X[26:length(X)] <- X[26:length(X)] - 3
plot(X, ylab = "Length of steel rod", xlab = "Index" )
lines(x = c(1,  25), y = rep(mean(X[1:25]), 2))
lines(x = c(26, 50), y = rep(mean(X[26:50]), 2))
rm(X)


###################################################
### Chap02.Rnw:246-247
###################################################
set.seed(123)


###################################################
### RandomVariationAroundTrend
###################################################
X <-  seq(from=1,           # Assign to X a sequence from 1
          to=50,            # to 50
          by=1)             # increment of sequence
#
# Equivalent to:
# X <- 1:50                 # Integer sequence from 1 to 50
#
X <- sin(                   # Reassign X with sine of
  X*(2*pi)/50)              # X*2*pi/50
#
X <- X + rnorm(n=length(X), # Add to X a random normal 
               mean=0,      # component with mean 0 
               sd=0.05)     # and standard deviation 0.05
# 
plot(X,                     
     ylab="Values")         
#
abline(h=0,                 # Add a horizontal line at y=0 
       lty="dashed",        # set line type
       col="lightgray")     # set line color


###################################################
### PlotRandomVariationAroundTrend
###################################################
set.seed(123)
X <-  seq(from=1,           # Assign to X a sequence from 1
          to=50,            # to 50
          by=1)             # increment of sequence
#
# Equivalent to:
# X <- 1:50                 # Integer sequence from 1 to 50
#
X <- sin(                   # Reassign X with sine of
  X*(2*pi)/50)              # X*2*pi/50
#
X <- X + rnorm(n=length(X), # Add to X a random normal 
               mean=0,      # component with mean 0 
               sd=0.05)     # and standard deviation 0.05
# 
plot(X,                     
     ylab="Values")         
#
abline(h=0,                 # Add a horizontal line at y=0 
       lty="dashed",        # set line type
       col="lightgray")     # set line color
rm(X)


###################################################
### PlotSampleOfMeasurementsOnThreeInstruments
###################################################
set.seed(123)
X <- rep(c(5, 2, 5), each=10)
X <- X + c(rnorm(10, sd=0.5), rnorm(10, sd=0.2), rnorm(10, sd=0.1)) 
plot(X, 
     ylim=c(0,8),
     ylab="Weight")
lines(c(1,10), c(5,5))
lines(c(11,20), c(2,2))
lines(c(21,30), c(5,5))
text(5.5, 7, "A")
text(15.5, 4, "B")
text(25.5, 7, "C")
rm(X)


###################################################
### PlotBarDiagram
###################################################
X <- matrix(c(5, 13, 45, 25, 12), nrow=1)
barplot(X, 
        width=1, 
        space=4, 
        col="grey50",
        names.arg=c(10, 20, 30, 40, 50), 
        ylab="Frequency")


###################################################
### PlotCumulativeFrequencyDistribution
###################################################
plot(c(10, 20, 30, 40, 50), cumsum(X), 
     type="s", 
     xlab="", ylab="Cumulative Frequency", 
     ylim=c(0,100))
rm(X)


###################################################
### BarDiagramBlemishes
###################################################
data(BLEMISHES)  # BLEMISHES is a matrix-like structure           


###################################################
### BLEMISHEShead
###################################################
BLEMISHES[1:3, ] # Return rows 1 to 3, all columns


###################################################
### BLEMISHEShead
###################################################
head(BLEMISHES,  # Equivalently head returns part of an object 
     n=3)        # set number of elements


###################################################
### BarDiagramBlemishes
###################################################
X <- factor(BLEMISHES$count, # Encode count vector as a factor
  levels=0:5)                # specify labels for the levels
                             #
X <- table(X                 # Reassign to X a frequency table
          )                  # of encoded variable count
                             #
X <- prop.table(X)           # Reassign with proportional table
                             #
barplot(X,                   # Plot a Bar diagram
        width=1,             # set width of bars
        space=4,             # set space between bars
        col="grey50",       
        ylab="Proportional Frequency")


###################################################
### PlotBarDiagramBlemishes
###################################################
X <- factor(BLEMISHES$count, # Encode count vector as a factor
            levels=0:5)                # specify labels for the levels
#
X <- table(X                 # Reassign to X a frequency table
)                  # of encoded variable count
#
X <- prop.table(X)           # Reassign with proportional table
#
barplot(X,                   # Plot a Bar diagram
        width=1,             # set width of bars
        space=4,             # set space between bars
        col="grey50",       
        ylab="Proportional Frequency")


###################################################
### PlotCumulativeFrequancyDistributionBlemishes
###################################################
plot(names(X), cumsum(X), 
     type="s", 
     xlab="", ylab="Cumulative Frequency")
rm(X)


###################################################
### YARNSTRGhist02
###################################################
hist(YARNSTRG,  # Plot an histogram of the given data values
     breaks=6,  # set the number of cells for the histogram 
     main="",   # set the main title to void
     xlab = "Log yarn strength")


###################################################
### YARNSTRGecdf
###################################################
plot.ecdf(YARNSTRG,  # Plot empirical cumulative distribution 
          pch=NA,    # set no symbol in steps
          main="",   
          xlab="Log Yarn Strength")


###################################################
### YARNSTRGhist01
###################################################
data(YARNSTRG)  # YARNSTRG is a vector of values
                #
hist(YARNSTRG)  # Plot an histogram of the given data values


###################################################
### YARNSTRGhist02code
###################################################
hist(YARNSTRG,  # Plot an histogram of the given data values
     breaks=6,  # set the number of cells for the histogram 
     main="",   # set the main title to void
     xlab = "Log yarn strength")
                     #      
plot.ecdf(YARNSTRG,  # Plot empirical cumulative distribution 
          pch=NA,    # set no symbol in steps
          main="",   
          xlab="Log Yarn Strength")


###################################################
### PlotEcdfSampleStatistics
###################################################
plot.ecdf(YARNSTRG, 
          main="")
abline(h=c(0.25,0.5,0.75), 
       lty=2, 
       col="gray70")
arrows(sort(YARNSTRG)[25], 0.25, 
       sort(YARNSTRG)[25], 0)
arrows(sort(YARNSTRG)[50], 0.50, 
       sort(YARNSTRG)[50], 0)
arrows(sort(YARNSTRG)[75], 0.75, 
       sort(YARNSTRG)[75], 0)
text(sort(YARNSTRG)[25] + 0.2, 0.09, "Q1")
text(sort(YARNSTRG)[50] + 0.2, 0.09, "Q2")
text(sort(YARNSTRG)[75] + 0.2, 0.09, "Q3")


###################################################
### CYCLTsummary
###################################################
quantile(CYCLT,                # Sample quantiles
         probs = seq(from=0,   # set given probabilities
                     to=1,     #
                     by=0.25), #
         type=6                # set algorithm to be used:
)                              # 6 for MINITAB like one
#
mean(CYCLT,                    # Arithmetic mean
     trim=0.0,                 # set fraction of observations
     # to be trimmed from each end
     na.rm=TRUE                # set whether NA values should
)                              # be stripped out
#
# summary(CYCLT)               # As above but uses R default 
# type algorithm


###################################################
### BLEMISHESSkeKurt
###################################################
library(e1071)             # Load package e1071
                           #
skewness(BLEMISHES$count)  # Computes the skewness
                           #
kurtosis(BLEMISHES$count)  # Computes the kurtosis


###################################################
### PlotDistributionSkewness
###################################################
X <- seq(0, 1, length.out=200)
plot(X, dbeta(X, shape1=2, shape2=8), 
     type="l", 
     xlab=expression(x),
     ylab=expression(y),
     ylim=c(0, 4.5))
lines(X, dbeta(X, shape1=10, shape2=10))
lines(X, dbeta(X, shape1=8, shape2=2))
text(0.1, 4.2, labels="Positive \nSkewness")
text(0.5, 4.2, labels="Symmetric")
text(0.9, 4.2, labels="Negative \nSkewness")


###################################################
### PlotDistributionSteepness
###################################################
X <- seq(0, 1, length.out=200)
plot(X*6-3, dbeta(X, shape1=8, shape2=8)/6, 
     type="l", 
     xlab=expression(x),
     ylab=expression(y),
     ylim=c(0, 0.55))
lines(X*6-3, dbeta(X, shape1=2.5, shape2=2.5)/6, lty="dotdash")
X <- seq(-3, 3, length.out=200)
lines(X, dnorm(X), lty="dashed")
text(1.0, 0.52, labels="Steep")
text(1.5, 0.40, labels="Normal")
text(2.5, 0.28, labels="Flat")
rm(X)


###################################################
### BoxPlotLogYarnStrength
###################################################
boxplot(YARNSTRG,             # Produce box-and-whisker plot
        ylab="Log Yarn Strength")


###################################################
### PlotBoxPlotLogYarnStrength
###################################################
boxplot(YARNSTRG,             # Produce box-and-whisker plot
        ylab="Log Yarn Strength")


###################################################
### PlotQQPlotLogYarnStrength
###################################################
qqplot(qunif(seq(0,1,length.out=100)), YARNSTRG,
       xlab="Fraction of Sample",
       ylab="Sample Quantiles")


###################################################
### OELECTrobust
###################################################
File <- paste(                # Compose a string with
  path.package("mistat"),     # mistat package path and
  "/csvFiles/OELECT.csv",     # /DatFiles/OELECT.csv 
  sep = "")                   # separate the terms with ""
#
Oelect <- read.csv(file=File) # Read a csv file and assign to
# Oelect
#
rm(File)                      # 
#
mySummary <- function(        # Define a new function
  x, trim=0, type=6 )         # with arguments x, trim and type                          
  #
{                             # The new function does:
  #
  qq <-quantile(x, type=type) # Calculate quantiles
  #
  qq <- c(qq[1L:3L],          # Concatenate quantiles and mean   
          mean(x, trim=trim), 
          qq[4L:5L])
  #
  names(qq) <- c("Min.", "1st Qu.",  # Assign names to values
                 "Median", "Mean", 
                 "3rd Qu.", "Max.")
  #
  qq <- signif(qq)            # Round to significant values
  #
  return(qq)                  # Return results
  #
}                             # Function end 
#
mySummary(Oelect$OELECT)      # Apply mySummary to Oelect data
#                         
mySummary(Oelect$OELECT,      # Apply mySummary to Oelect data 
          trim=0.05)          # set trim to 5%
#
sd(Oelect$OELECT)             # Computes the standard deviation


###################################################
### OELECTrobust2
###################################################
OutVolt <- sort(Oelect$OELECT)  # Assign a vector 
# of sorted values
#
OutVolt[99] <- 2289.86          # Assign a specific value
# at position 99
#
mySummary(OutVolt)    # Apply function mySummary
#
mySummary(OutVolt,    # Apply mySummary with
          trim=0.05)  # trim = 5%
#
sd(OutVolt)           # Computes the standard deviation


###################################################
### OELECTrobust3
###################################################
IQR(OutVolt)/1.349  # Robust estimate of S


###################################################
### Chap02End
###################################################
rm(Oelect, OutVolt, mySummary)
detach(package:mistat)

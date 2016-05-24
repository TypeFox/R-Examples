### R code from vignette source 'CDF_Analysis.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
# Load the spsurvey package
library(spsurvey)



###################################################
### code chunk number 2: data
###################################################
# Load the data set and determine the number of rows in the data frame
data(FL_lakes)
nr <- nrow(FL_lakes)



###################################################
### code chunk number 3: data
###################################################
# Display the initial six lines in the data file
head(FL_lakes)



###################################################
### code chunk number 4: data
###################################################
# Display the final six lines in the data file
tail(FL_lakes)



###################################################
### code chunk number 5: Quanteval
###################################################
# Use the summary function to summarize the data structure of the dissolved
# oxygen variable
cat("\nSummarize the data structure of the dissolved oxygen variable:\n")
summary(FL_lakes$Oxygen)



###################################################
### code chunk number 6: Quanteval
###################################################
#
# Conduct an analysis of lake condition variables
#

# Create the sites data frame, which identifies sites to use in the analysis
# Note that only sampled sites are used
sites <- data.frame(siteID=FL_lakes$siteID,
                    Use=FL_lakes$Status == "Sampled")



###################################################
### code chunk number 7: Quanteval
###################################################
# Create the subpop data frame, which defines populations and subpopulations for
# which estimates are desired
subpop <- data.frame(siteID=FL_lakes$siteID,
                     Basin=FL_lakes$Basin)



###################################################
### code chunk number 8: Quanteval
###################################################
# Create the design data frame, which identifies the stratum code, weight,
#    x-coordinate, and y-coordinate for each site ID
design <- data.frame(siteID=FL_lakes$siteID,
                     wgt=FL_lakes$wgt,
                     xcoord=FL_lakes$xcoord,
                     ycoord=FL_lakes$ycoord)



###################################################
### code chunk number 9: Quanteval
###################################################
# Create the data.cont data frame, which specifies the variables to use in the
# analysis
data.cont <- data.frame(siteID=FL_lakes$siteID,
                        DissolvedOxygen=FL_lakes$Oxygen)



###################################################
### code chunk number 10: Conditionevalpop
###################################################
#
# Conduct an analysis of the dissolved oxygen variables correcting for
# population size
#

# Assign frame size values
framesize <- c("NWFWMD-1"=451, "NWFWMD-2"=394, "SFWMD-9"=834, "SJRWMD-1"=1216,
               "SRWMD-1"=1400, "SWFWMD-4"=851)



###################################################
### code chunk number 11: Quanteval
###################################################
# Calculate CDF estimates for the quantitative variables
CDF_Estimates <- cont.analysis(sites, subpop, design, data.cont,
   popsize=list(Basin=as.list(framesize)))



###################################################
### code chunk number 12: Quanteval
###################################################
# Write CDF estimates as a csv file
write.table(CDF_Estimates$CDF, file="CDF_Estimates.csv", sep=",",
            row.names=FALSE)



###################################################
### code chunk number 13: Quanteval
###################################################
cont.cdfplot("CDF_Estimates.pdf", CDF_Estimates$CDF)



###################################################
### code chunk number 14: Quanteval
###################################################
# Test for statistical difference between CDFs for basins
CDF_Tests <- cont.cdftest(sites, subpop, design, data.cont,
   popsize=list(Basin=as.list(framesize)))



###################################################
### code chunk number 15: Quanteval
###################################################
# Print results of the statistical tests for difference between CDFs from
# basins for dissolved oxygen
print(CDF_Tests, digits=3)



###################################################
### code chunk number 16: figure2
###################################################
# Display basins that have significantly different CDFs
n1 <- length(levels(CDF_Tests$Subpopulation_1))
n2 <- length(levels(CDF_Tests$Subpopulation_2))
plot(1:n2, 1:n1, type="n", xlab="Second Basin", ylab="First Basin", xaxt="n",
     yaxt="n")
count=1
for(i in 1:n1) {
   for(j in i:n2) {
      text(j, i, ifelse(CDF_Tests$p_Value[count] < 0.01, "X", " "))
      count <- count+1
   }
}
axis(side=1, at=1:n2, labels=levels(CDF_Tests$Subpopulation_2), cex.axis=0.75)
axis(side=2, at=1:n1, labels=levels(CDF_Tests$Subpopulation_1), cex.axis=0.75)
title("Significantly Different CDFs")
abline(1, 1, col="red", lwd=2)



###################################################
### code chunk number 17: Quanteval
###################################################
# Write CDF test results as a csv file
write.table(CDF_Tests, file="CDF_Tests.csv", sep=",", row.names=FALSE)




### R code from vignette source 'Linear_Analysis.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
# Load the spsurvey package
library(spsurvey)



###################################################
### code chunk number 2: data
###################################################
# Load the data set and determine the number of rows in the data frame
data(IN_streams)
nr <- nrow(IN_streams)



###################################################
### code chunk number 3: data
###################################################
# Display the initial six lines in the data file
head(IN_streams)



###################################################
### code chunk number 4: data
###################################################
# Display the final six lines in the data file
tail(IN_streams)



###################################################
### code chunk number 5: Statuseval
###################################################
# Use the table and addmargins functions to create a table displaying the count
# for each code of the status variable
cat("\nA table displaying the number of values for each level of the status
variable follows:\n")
addmargins(table(IN_streams$Status))



###################################################
### code chunk number 6: Statuseval
###################################################
# Use the table and addmargins functions to create a table displaying the count
# for each code of the TNT variable
cat("\nA table displaying the number of values for each level of the TNT
variable follows:\n")
addmargins(table(IN_streams$TNT))



###################################################
### code chunk number 7: Statuseval
###################################################
#
# Conduct an analysis of site status evaluation variables
#

# Create the sites data frame, which identifies sites to use in the analysis
# Note that all sites will be used to estimate number of streams in each category
sites <- data.frame(siteID=IN_streams$siteID,
                    Use=rep(TRUE, nr))



###################################################
### code chunk number 8: Statuseval
###################################################
# Create the subpop data frame, which defines populations and subpopulations for
# which estimates are desired
subpop <- data.frame(siteID=IN_streams$siteID,
                     Upper_Wabash=rep("Upper Wabash", nr), 
							       Strahler_Order=IN_streams$Strahler_Cat)



###################################################
### code chunk number 9: Statuseval
###################################################
# Create the design data frame, which identifies the stratum code, weight,
#    x-coordinate, and y-coordinate for each site ID
design <- data.frame(siteID=IN_streams$siteID,
                     wgt=IN_streams$wgt,
                     xcoord=IN_streams$xcoord,
                     ycoord=IN_streams$ycoord)



###################################################
### code chunk number 10: Statuseval
###################################################
# Create the data.cat data frame, which specifies the variables to use in the
# analysis
data.cat <- data.frame(siteID=IN_streams$siteID,
                       Status=IN_streams$Status,
                       Target_NonTarget=IN_streams$TNT)



###################################################
### code chunk number 11: Statuseval
###################################################
# Calculate extent estimates for the site status evaluation variables
Extent_Estimates <- cat.analysis(sites, subpop, design, data.cat)



###################################################
### code chunk number 12: Statuseval
###################################################
# Print the extent estimates for all basins combined
print(Extent_Estimates[c(1:8, 32:34),])



###################################################
### code chunk number 13: Statuseval
###################################################
# Write results as a comma-separated value (csv) file
write.csv(Extent_Estimates, file="Extent_Estimates.csv", sep=",",
            row.names=FALSE)



###################################################
### code chunk number 14: Conditioneval
###################################################
# Use the table and addmargins functions to create a table displaying the count
# for each code of the IBI status variable
cat("\nA table displaying the number of values for each level of the IBI status
variable follows:\n")
addmargins(table(IN_streams$IBI_Status))



###################################################
### code chunk number 15: Conditioneval
###################################################
# Use the table and addmargins functions to create a table displaying the count
# for each code of the QHEI status variable
cat("\nA table displaying the number of values for each level of the QHEI status
variable follows:\n")
addmargins(table(IN_streams$QHEI_Status))



###################################################
### code chunk number 16: Conditioneval
###################################################
#
# Conduct an analysis of stream condition variables
#

# Create the sites data frame
# Note that only sampled sites are used
sites <- data.frame(siteID=IN_streams$siteID,
                    Use=IN_streams$Status == "Sampled")

# Note that the existing subpop and design data frames can be reused



###################################################
### code chunk number 17: Conditioneval
###################################################
# Create the data.cat data frame, which specifies the variables to use in the
# analysis
data.cat <- data.frame(siteID=IN_streams$siteID,
                       IBI_Status=IN_streams$IBI_Status,
                       QHEI_Status=IN_streams$QHEI_Status)



###################################################
### code chunk number 18: Conditioneval
###################################################
# Calculate estimates for the categorical variables
Condition_Estimates <- cat.analysis(sites, subpop, design, data.cat)



###################################################
### code chunk number 19: Conditioneval
###################################################
# Print the condition estimates for all basins combined
print(Condition_Estimates[c(1:3, 16:18),])



###################################################
### code chunk number 20: Conditioneval
###################################################
# Write results as a csv file
write.csv(Condition_Estimates, file="Condition_Estimates.csv")



###################################################
### code chunk number 21: Conditionevalpop
###################################################
#
# Conduct an analysis of stream condition variables correcting for population
# size
#

# Note that the existing sites, subpop, design, and data.cont data frames can be
# reused

# Assign frame size values
framesize <- c("1st"=4514.450, "2nd"=1443.260, "3rd"=740.146, "4th"=660.294)



###################################################
### code chunk number 22: Conditionevalpop
###################################################
# Calculate estimates for the stream condition variables
Condition_Estimates_popsize <- cat.analysis(sites, subpop, design, data.cat,
   popsize=list(Upper_Wabash=sum(framesize),
                Strahler_Order=as.list(framesize)))



###################################################
### code chunk number 23: Conditionevalpop
###################################################
# Print the stream condition estimates for all sites combined
print(Condition_Estimates_popsize[c(1:3, 16:18),])



###################################################
### code chunk number 24: Conditionevalpop
###################################################
# Write results as a csv file
write.csv(Condition_Estimates_popsize, file="Condition_Estimates_popsize.csv")



###################################################
### code chunk number 25: Quanteval
###################################################
# Use the summary function to summarize the data structure of the IBI score
# variable
cat("\nSummarize the data structure of the IBI score variable:\n")
summary(IN_streams$IBI_Score)



###################################################
### code chunk number 26: Quanteval
###################################################
# Use the summary function to summarize the data structure of the QHEI score
# variable
cat("\nSummarize the data structure of the QHEI score variable:\n")
summary(IN_streams$QHEI_Score)



###################################################
### code chunk number 27: Quanteval
###################################################
#
# Conduct an analysis of quantitative variables
#

# Note that the existing sites, subpop, and design data frames can be reused

# Create the data.cont data frame, which specifies the variables to use in the
# analysis
data.cont <- data.frame(siteID=IN_streams$siteID,
                        IBI_Score=IN_streams$IBI_Score,
                        QHEI_Score=IN_streams$QHEI_Score)



###################################################
### code chunk number 28: Quanteval
###################################################
# Calculate CDF estimates for the quantitative variables
CDF_Estimates <- cont.analysis(sites, subpop, design, data.cont,
   popsize=list(Upper_Wabash=sum(framesize),
                Strahler_Order=as.list(framesize)))



###################################################
### code chunk number 29: Quanteval
###################################################
# Write CDF estimates as a csv file
write.csv(CDF_Estimates$CDF, file="CDF_Estimates.csv")



###################################################
### code chunk number 30: Quanteval
###################################################
cont.cdfplot("CDF_Estimates.pdf", CDF_Estimates$CDF)



###################################################
### code chunk number 31: Quanteval
###################################################
# Print the percentile estimates for IBI score for all sites combined
print(CDF_Estimates$Pct[1:10,])



###################################################
### code chunk number 32: Quanteval
###################################################
# Write percentile estimates as a csv file
write.csv(CDF_Estimates$Pct, file="Percentile_Estimates.csv")



###################################################
### code chunk number 33: Quanteval
###################################################
# Test for statistical difference between CDFs for Strahler order categories
CDF_Tests <- cont.cdftest(sites, subpop[,c(1,3)], design, data.cont,
   popsize=list(Strahler_Order=as.list(framesize)))



###################################################
### code chunk number 34: Quanteval
###################################################
# Print results of the statistical tests for difference between CDFs from
# Strahler order categories for IBI score
print(CDF_Tests, digits=2)



###################################################
### code chunk number 35: Quanteval
###################################################
# Write CDF test results as a csv file
write.csv(CDF_Tests, file="CDF_Tests.csv")




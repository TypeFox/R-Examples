### R code from vignette source 'Finite_Analysis.Rnw'

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
### code chunk number 5: Statuseval
###################################################
# Use the table and addmargins functions to create a table displaying the count
# for each code of the status variable
cat("\nA table displaying the number of values for each level of the status
variable follows:\n")
addmargins(table(FL_lakes$Status))



###################################################
### code chunk number 6: Statuseval
###################################################
# Use the table and addmargins functions to create a table displaying the count
# for each code of the TNT variable
cat("\nA table displaying the number of values for each level of the TNT
variable follows:\n")
addmargins(table(FL_lakes$TNT))



###################################################
### code chunk number 7: Statuseval
###################################################
#
# Conduct an analysis of site status evaluation variables
#

# Create the sites data frame, which identifies sites to use in the analysis
# Note that all sites will be used to estimate number of lakes in each category
sites <- data.frame(siteID=FL_lakes$siteID,
                    Use=rep(TRUE, nr))



###################################################
### code chunk number 8: Statuseval
###################################################
# Create the subpop data frame, which defines populations and subpopulations for
# which estimates are desired
subpop <- data.frame(siteID=FL_lakes$siteID,
                     CombinedBasins=rep("All Basins", nr), 
							       Basin=FL_lakes$Basin)



###################################################
### code chunk number 9: Statuseval
###################################################
# Create the design data frame, which identifies the stratum code, weight,
#    x-coordinate, and y-coordinate for each site ID
design <- data.frame(siteID=FL_lakes$siteID,
                     wgt=FL_lakes$wgt,
                     xcoord=FL_lakes$xcoord,
                     ycoord=FL_lakes$ycoord)



###################################################
### code chunk number 10: Statuseval
###################################################
# Create the data.cat data frame, which specifies the variables to use in the
# analysis
data.cat <- data.frame(siteID=FL_lakes$siteID,
                       Status=FL_lakes$Status,
                       Target_NonTarget=FL_lakes$TNT)



###################################################
### code chunk number 11: Statuseval
###################################################
# Calculate extent estimates for the site status evaluation variables
Extent_Estimates <- cat.analysis(sites, subpop, design, data.cat)



###################################################
### code chunk number 12: Statuseval
###################################################
# Print the extent estimates for all basins combined
print(Extent_Estimates[c(1:7, 45:47),])



###################################################
### code chunk number 13: Statuseval
###################################################
# Write results as a comma-separated value (csv) file
write.csv(Extent_Estimates, file="Extent_Estimates.csv")



###################################################
### code chunk number 14: Conditioneval
###################################################
# Use the table and addmargins functions to create a table displaying the count
# for each code of the pH category variable
cat("\nA table displaying the number of values for each level of the pH category
variable follows:\n")
addmargins(table(FL_lakes$pH_Cat))



###################################################
### code chunk number 15: Conditioneval
###################################################
# Use the table and addmargins functions to create a table displaying the count
# for each code of the fecal coliform category variable
cat("\nA table displaying the number of values for each level of the fecal
coliform category variable follows:\n")
addmargins(table(FL_lakes$Coliform_Cat))



###################################################
### code chunk number 16: Conditioneval
###################################################
#
# Conduct an analysis of lake condition variables
#

# Create the sites data frame
# Note that only sampled sites are used
sites <- data.frame(siteID=FL_lakes$siteID,
                    Use=FL_lakes$Status == "Sampled")

# Note that the existing subpop and design data frames can be reused



###################################################
### code chunk number 17: Conditioneval
###################################################
# Create the data.cat data frame, which specifies the variables to use in the
# analysis
data.cat <- data.frame(siteID=FL_lakes$siteID,
                       pHCat=FL_lakes$pH_Cat,
                       ColiformCat=FL_lakes$Coliform_Cat)



###################################################
### code chunk number 18: Conditioneval
###################################################
# Calculate estimates for the categorical variables
Condition_Estimates <- cat.analysis(sites, subpop, design, data.cat)



###################################################
### code chunk number 19: Conditioneval
###################################################
# Print the condition estimates for all basins combined
print(Condition_Estimates[c(1:4, 28:32),])



###################################################
### code chunk number 20: Conditioneval
###################################################
# Write results as a csv file
write.csv(Condition_Estimates, file="Condition_Estimates.csv")



###################################################
### code chunk number 21: Conditionevalpop
###################################################
#
# Conduct an analysis of lake condition variables correcting for population size
#

# Note that the existing sites, subpop, design, and data.cont data frames can be
# reused

# Assign frame size values
framesize <- c("NWFWMD-1"=451, "NWFWMD-2"=394, "SFWMD-9"=834, "SJRWMD-1"=1216,
               "SRWMD-1"=1400, "SWFWMD-4"=851)



###################################################
### code chunk number 22: Conditionevalpop
###################################################
# Calculate estimates for the lake condition variables
Condition_Estimates_popsize <- cat.analysis(sites, subpop, design, data.cat,
   popsize=list(CombinedBasins=sum(framesize),
                Basin=as.list(framesize)))



###################################################
### code chunk number 23: Conditionevalpop
###################################################
# Print the lake condition estimates for all basins combined
print(Condition_Estimates_popsize[c(1:4, 28:32),])



###################################################
### code chunk number 24: Conditionevalpop
###################################################
# Write results as a csv file
write.csv(Condition_Estimates_popsize, file="Condition_Estimates_popsize.csv")



###################################################
### code chunk number 25: Quanteval
###################################################
# Use the summary function to summarize the data structure of the dissolved
# oxygen variable
cat("\nSummarize the data structure of the dissolved oxygen variable:\n")
summary(FL_lakes$Oxygen)



###################################################
### code chunk number 26: Quanteval
###################################################
# Use the summary function to summarize the data structure of the turbidity
# variable
cat("\nSummarize the data structure of the turbidity variable:\n")
summary(FL_lakes$Turbidity)



###################################################
### code chunk number 27: Quanteval
###################################################
#
# Conduct an analysis of quantitative variables
#

# Note that the existing sites, subpop, and design data frames can be reused

# Create the data.cont data frame, which specifies the variables to use in the
# analysis
data.cont <- data.frame(siteID=FL_lakes$siteID,
                        DissolvedOxygen=FL_lakes$Oxygen,
                        Turbidity=FL_lakes$Turbidity)



###################################################
### code chunk number 28: Quanteval
###################################################
# Calculate CDF estimates for the quantitative variables
CDF_Estimates <- cont.analysis(sites, subpop, design, data.cont,
   popsize=list(CombinedBasins=sum(framesize),
                Basin=as.list(framesize)))



###################################################
### code chunk number 29: Quanteval
###################################################
# Write CDF estimates as a csv file
write.csv(CDF_Estimates$CDF, file="CDF_Estimates.csv")



###################################################
### code chunk number 30: Quanteval
###################################################
cont.cdfplot("CDF_Estimates.pdf", CDF_Estimates$CDF, logx=c("","x"))



###################################################
### code chunk number 31: Quanteval
###################################################
# Print the percentile estimates for dissolved oxygen for all basins combined
print(CDF_Estimates$Pct[1:10,])



###################################################
### code chunk number 32: Quanteval
###################################################
# Write percentile estimates as a csv file
write.csv(CDF_Estimates$Pct, file="Percentile_Estimates.csv")



###################################################
### code chunk number 33: Quanteval
###################################################
# Test for statistical difference between CDFs for basins
CDF_Tests <- cont.cdftest(sites, subpop[,c(1,3)], design, data.cont,
   popsize=list(Basin=as.list(framesize)))



###################################################
### code chunk number 34: Quanteval
###################################################
# Print results of the statistical tests for difference between CDFs from
# basins for dissolved oxygen
print(CDF_Tests, digits=3)



###################################################
### code chunk number 35: Quanteval
###################################################
# Write CDF test results as a csv file
write.csv(CDF_Tests, file="CDF_Tests.csv")




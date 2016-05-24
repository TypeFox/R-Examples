### R code from vignette source 'Area_Analysis.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
# Load the spsurvey package
library(spsurvey)



###################################################
### code chunk number 2: data
###################################################
# Load the data set and determine the number of rows in the data frame
data(SC_estuaries)
nr <- nrow(SC_estuaries)



###################################################
### code chunk number 3: data
###################################################
# Display the initial six lines in the data file
head(SC_estuaries)



###################################################
### code chunk number 4: data
###################################################
# Display the final six lines in the data file
tail(SC_estuaries)



###################################################
### code chunk number 5: Statuseval
###################################################
# Use the table and addmargins functions to create a table displaying the count
# for each code of the status variable
cat("\nA table displaying the number of values for each level of the status
variable follows:\n")
addmargins(table(SC_estuaries$Status))



###################################################
### code chunk number 6: Statuseval
###################################################
#
# Conduct an analysis of site status evaluation variables
#

# Create the sites data frame, which identifies sites to use in the analysis
# Note that all sites will be used to estimate number of estuaries in each category
sites <- data.frame(siteID=SC_estuaries$siteID,
                    Use=rep(TRUE, nr))



###################################################
### code chunk number 7: Statuseval
###################################################
# Create the subpop data frame, which defines populations and subpopulations for
# which estimates are desired
subpop <- data.frame(siteID=SC_estuaries$siteID,
                     All_Estuaries=rep("All Estuaries", nr), 
							       Estuary_Type=SC_estuaries$Stratum)



###################################################
### code chunk number 8: Statuseval
###################################################
# Create the design data frame, which identifies the stratum code, weight,
#    x-coordinate, and y-coordinate for each site ID
design <- data.frame(siteID=SC_estuaries$siteID,
                     wgt=SC_estuaries$wgt,
                     xcoord=SC_estuaries$xcoord,
                     ycoord=SC_estuaries$ycoord)



###################################################
### code chunk number 9: Statuseval
###################################################
# Create the data.cat data frame, which specifies the variables to use in the
# analysis
data.cat <- data.frame(siteID=SC_estuaries$siteID,
                       Status=SC_estuaries$Status)



###################################################
### code chunk number 10: Statuseval
###################################################
# Calculate extent estimates for the site status evaluation variables
Extent_Estimates <- cat.analysis(sites, subpop, design, data.cat)



###################################################
### code chunk number 11: Statuseval
###################################################
# Print the extent estimates
print(Extent_Estimates)



###################################################
### code chunk number 12: Statuseval
###################################################
# Write results as a comma-separated value (csv) file
write.csv(Extent_Estimates, file="Extent_Estimates.csv")



###################################################
### code chunk number 13: Conditioneval
###################################################
# Use the table and addmargins functions to create a table displaying the count
# for each code of the IBI status variable
cat("\nA table displaying the number of values for each level of the IBI status
variable follows:\n")
addmargins(table(SC_estuaries$IBI_status))



###################################################
### code chunk number 14: Conditioneval
###################################################
# Use the table and addmargins functions to create a table displaying the count
# for each code of the WQ status variable
cat("\nA table displaying the number of values for each level of the WQ status variable follows:\n")
addmargins(table(SC_estuaries$WQ_status))



###################################################
### code chunk number 15: Conditioneval
###################################################
#
# Conduct an analysis of estuary condition variables
#

# Create the sites data frame
# Note that only sampled sites are used
sites <- data.frame(siteID=SC_estuaries$siteID,
                    Use=SC_estuaries$Status == "Sampled")

# Note that the existing subpop and design data frames can be reused



###################################################
### code chunk number 16: Conditioneval
###################################################
# Create the data.cat data frame, which specifies the variables to use in the
# analysis
data.cat <- data.frame(siteID=SC_estuaries$siteID,
                       IBI_Status=SC_estuaries$IBI_status,
                       WQ_Status=SC_estuaries$WQ_status)



###################################################
### code chunk number 17: Conditioneval
###################################################
# Calculate estimates for the categorical variables
Condition_Estimates <- cat.analysis(sites, subpop, design, data.cat)



###################################################
### code chunk number 18: Conditioneval
###################################################
# Print the condition estimates for all basins combined
print(Condition_Estimates[c(1:4, 13:16),])



###################################################
### code chunk number 19: Conditioneval
###################################################
# Write results as a csv file
write.csv(Condition_Estimates, file="Condition_Estimates.csv")



###################################################
### code chunk number 20: Conditionevalpop
###################################################
#
# Conduct an analysis of estuary condition variables correcting for population
# size
#

# Note that the existing sites, subpop, design, and data.cont data frames can be
# reused

# Assign frame size values
framesize <- c("Open Water"=628.509298, "Tidal Creek"=105.829522)



###################################################
### code chunk number 21: Conditionevalpop
###################################################
# Calculate estimates for the estuary condition variables
Condition_Estimates_popsize <- cat.analysis(sites, subpop, design, data.cat,
   popsize=list(All_Estuaries=sum(framesize),
                Estuary_Type=as.list(framesize)))



###################################################
### code chunk number 22: Conditionevalpop
###################################################
# Print the estuary condition estimates for all sites combined
print(Condition_Estimates_popsize[c(1:4, 13:16),])



###################################################
### code chunk number 23: Conditionevalpop
###################################################
# Write results as a csv file
write.csv(Condition_Estimates_popsize, file="Condition_Estimates_popsize.csv")



###################################################
### code chunk number 24: Quanteval
###################################################
# Use the summary function to summarize the data structure of the IBI score
# variable
cat("\nSummarize the data structure of the IBI score variable:\n")
summary(SC_estuaries$IBI_score)



###################################################
### code chunk number 25: Quanteval
###################################################
# Use the summary function to summarize the data structure of the WQ score
# variable
cat("\nSummarize the data structure of the WQ score variable:\n")
summary(SC_estuaries$WQ_score)



###################################################
### code chunk number 26: Quanteval
###################################################
#
# Conduct an analysis of quantitativen variables
#

# Note that the existing sites, subpop, and design data frames can be reused

# Create the data.cont data frame, which specifies the variables to use in the
# analysis
data.cont <- data.frame(siteID=SC_estuaries$siteID,
                        IBI_Score=SC_estuaries$IBI_score,
                        WQ_Score=SC_estuaries$WQ_score)



###################################################
### code chunk number 27: Quanteval
###################################################
# Calculate CDF estimates for the quantitative variables
CDF_Estimates <- cont.analysis(sites, subpop, design, data.cont,
   popsize=list(All_Estuaries=sum(framesize),
                Estuary_Type=as.list(framesize)))



###################################################
### code chunk number 28: Quanteval
###################################################
# Write CDF estimates as a csv file
write.csv(CDF_Estimates$CDF, file="CDF_Estimates.csv")



###################################################
### code chunk number 29: Quanteval
###################################################
cont.cdfplot("CDF_Estimates.pdf", CDF_Estimates$CDF)



###################################################
### code chunk number 30: Quanteval
###################################################
# Print the percentile estimates for IBI score for all sites combined
print(CDF_Estimates$Pct[1:10,])



###################################################
### code chunk number 31: Quanteval
###################################################
# Write percentile estimates as a csv file
write.csv(CDF_Estimates$Pct, file="Percentile_Estimates.csv")



###################################################
### code chunk number 32: Quanteval
###################################################
# Test for statistical difference between CDFs for strata
CDF_Tests <- cont.cdftest(sites, subpop[,c(1,3)], design, data.cont,
   popsize=list(Estuary_Type=as.list(framesize)))



###################################################
### code chunk number 33: Quanteval
###################################################
# Print results of the statistical tests for difference between strata CDFs for
# IBI score and WQ score
print(CDF_Tests, digits=3)



###################################################
### code chunk number 34: Quanteval
###################################################
# Write CDF test results as a csv file
write.csv(CDF_Tests, file="CDF_Tests.csv")




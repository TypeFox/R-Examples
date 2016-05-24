### R code from vignette source 'Risk_Analysis.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
# Load the spsurvey package
library(spsurvey)



###################################################
### code chunk number 2: data
###################################################
# Load the data set and determine the number of rows in the data frame
data(NLA_2007)
nr <- nrow(NLA_2007)



###################################################
### code chunk number 3: data
###################################################
# Display the initial six lines in the data file
head(NLA_2007)



###################################################
### code chunk number 4: data
###################################################
# Display the final six lines in the data file
tail(NLA_2007)



###################################################
### code chunk number 5: Riskeval
###################################################
#
# Conduct a relative risk analysis
#

# Create the sites data frame, which identifies sites to use in the analysis
sites <- data.frame(siteID=NLA_2007$siteID,
                    Use=rep(TRUE, nr))



###################################################
### code chunk number 6: Riskeval
###################################################
# Create the subpop data frame, which defines populations and subpopulations for
# which estimates are desired
subpop <- data.frame(siteID=NLA_2007$siteID,
                     Western_US=rep("Western_US", nr),
                     Lake_Origin=NLA_2007$Lake_Origin)



###################################################
### code chunk number 7: Riskeval
###################################################
# Create the design data frame, which identifies the stratum code, weight,
#    x-coordinate, and y-coordinate for each site ID
design <- data.frame(siteID=NLA_2007$siteID,
                     wgt=NLA_2007$wgt,
                     xcoord=NLA_2007$xcoord,
                     ycoord=NLA_2007$ycoord)



###################################################
### code chunk number 8: Riskeval
###################################################
# Create the data.risk data frame, which specifies the variables to use in the
# analysis
data.risk <- data.frame(siteID=NLA_2007$siteID,
                        Chlorophyll_a=NLA_2007$Chla_cond,
                        MacroInvert_OE=NLA_2007$OE5_cond,
                        Total_Nitrogen=NLA_2007$NTL_cond,
                        Total_Phosphorus=NLA_2007$PTL_cond,
                        Turbidity=NLA_2007$Turbidity_cond)



###################################################
### code chunk number 9: Riskeval
###################################################
#
# Conduct a relative risk analysis
#

# Assign response and stressor variable names
resp_vars <- c("Chlorophyll_a", "MacroInvert_OE")
stress_vars <- c("Total_Nitrogen", "Total_Phosphorus", "Turbidity")



###################################################
### code chunk number 10: Riskeval
###################################################
# Calculate relative risk estimates
RelRisk_Estimates <- relrisk.analysis(sites, subpop, design, data.risk,
   response.var= rep(resp_vars, each=length(stress_vars)),
   stressor.var=rep(stress_vars, length(resp_vars)))



###################################################
### code chunk number 11: Riskeval
###################################################
# Print the relative risk estimates
print(RelRisk_Estimates)



###################################################
### code chunk number 12: Riskeval
###################################################
# Write results as a comma-separated value (csv) file
write.csv(RelRisk_Estimates, file="RelRisk_Estimates.csv")



###################################################
### code chunk number 13: Riskeval
###################################################
# Calculate attributable risk estimates
AttRisk_Estimates <- attrisk.analysis(sites, subpop, design, data.risk,
   response.var= rep(resp_vars, each=length(stress_vars)),
   stressor.var=rep(stress_vars, length(resp_vars)))



###################################################
### code chunk number 14: Riskeval
###################################################
# Print the attributable risk estimates
print(AttRisk_Estimates)



###################################################
### code chunk number 15: Riskeval
###################################################
# Write results as a csv file
write.csv(AttRisk_Estimates, file="AttRisk_Estimates.csv")




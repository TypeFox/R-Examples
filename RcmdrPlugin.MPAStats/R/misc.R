# Last modified: 2012-06-15 by Andrew Heiss
#--------------------------------------------

# Run `str()` on the active dataset
checkStructure <- function() {
  doItAndPrint(paste("str(", ActiveDataSet(), ")", sep = ""))
}


# Run `DataframeSummary()` on the active dataset
# TODO: Create a dialog that allows the user to toggle and specify confidence levels
summarizeDataframe <- function() {
  doItAndPrint(paste("DataframeSummary(", ActiveDataSet(), ")", sep = ""))
}


# Calculate the factor change coefficients (e^b) for the 
# current logistic regression model
factorChange <- function() {
  .activeModel <- ActiveModel()
  doItAndPrint(paste("exp(coef(", .activeModel, "))", sep = ""))
}
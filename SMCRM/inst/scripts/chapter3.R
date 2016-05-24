# Reproduction of Analyses of Chapter 3
# 
# Author: Tobias Verbeke <tobias.verbeke@openanalytics.eu>
###############################################################################

library(SMCRM)

# load data
data(customerAcquisition)

# response probability
responseModel <- glm(acquisition ~ acq_expense + acq_expense_sq + industry + revenue + employees,
    family = binomial, data = customerAcquisition)
summary(responseModel)

# number of acquired customers

# initial order quantity

# duration - time

# firm performance

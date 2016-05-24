########################################################################
# This is a longer version of the example to demostrate more #functions.
########################################################################

rm(list = ls())

require(EditImputeCont)

##################
# Input data formats
##################

data(NestedEx) ; makeDemoTxt(NestedEx)     # make demo textfiles from NestedEx dataset
 
D_obs1 = read.table("Demo_D_obs.txt", header = TRUE) 
  # The original data file has a header with variable names 
  # to be used in ratio edit description.
	
Ratio1 = read.table("Demo_Ratio_edit.txt", header = FALSE) 
Range1 = read.table("Demo_Range_edit.txt", header = FALSE) 
Balance1 = read.table("Demo_Balance_edit.txt", header = FALSE) 

### Alternatively, you can use the syntax of R pacakge 'editrules' ### 
#
# Ratio1 <- editmatrix(c(
#  "X1 <= 1096.63*X5", "X1 <= 2980.96*X7", "X1 <= 148.41*X8", "X1 <= 7.39*X9",
#  "X5 <= 0.37*X1", "X5 <= 54.60*X7", "X5 <= 2.72*X8", "X5 <= 0.14*X9",
#  "X7 <= 0.14*X1", "X7 <= 1.65*X5", "X7 <= 7.39*X8", "X7 <= 0.05*X9",
#  "X8 <= 1.65*X1", "X8 <= 54.60*X5", "X8 <= 403.43*X7", "X8 <= 1.65*X9",
#  "X9 <= 20.09*X1", "X9 <= 403.43*X5", "X9 <= 13359.73*X7", "X9 <= 148.41*X8"
# ))
# Range1 <- editmatrix(c(
#  "X1 >= 2", "X2 <= 1.2e+06", "X11 >= 0.002", "X11 <= 1.2e+04"
# ))
# Balance1 <- editmatrix(c(
#  "X1 == X2+X3+X4", "X5 == X6 + 0.4*X10 + 0.6*X11", "X7 == 0.4*X10 + 0.6*X11"
# ))

##################
# Read the data
##################

data1 = readData(Y.original=D_obs1, ratio=Ratio1, range=Range1, 
balance=Balance1, eps.bal=0.6)

print(data1$Edit.editmatrix)
plot(data1$Edit.editmatrix)	  ## function of 'editrules' package

##################
# Make and initialize the model
##################

model1 = createModel(data1, K=50)

##################
# Run MCMC
##################

# run 1 iteration
model1$Iterate()

# run additional 10 iterations
model1$Run(10)
scatterPlot(model1, data1, xvar=1) 
  # compare Y.input and Y.edited at iter = 11
EI_data1 = model1$Y.edited    # store the first edit-imputed dataset

# run additional 10 iterations
model1$Run(10)
scatterPlot(model1, data1) 
  # compare Y.input and Y.edited at iter = 21
EI_data2 = model1$Y.edited    # store the second edit-imputed dataset

# You can draw multiple edit-imputed data by multiple EI functions
# Get 5 edit-imputed data from MCMC with 500 iterations after 50 burn-in
result1 = multipleEI(model1, 50, 5, 10)
dim(result1)
# [1]   5 1000   11
# 5 Edit-imputed datasets of 1000 records with 11 variables
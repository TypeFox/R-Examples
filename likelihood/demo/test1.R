##
## TEST 1 - calculation of negative log likelihood
## Model: Radius = a + b * DBH
## Dataset: included crown_rad dataset
##

library(likelihood)

# Set up our dataset
data(crown_rad)
dataset <- crown_rad

# Create our model function
model <- function (a, b, DBH) {a + b * DBH}

# Create a place to put our parameters and values,
# set values for a and b, and indicate
# that DBH comes from the column marked "DBH"
# in the dataset
par <- list(a = 1.12, b = 0.07)
var <- list(DBH = "DBH")

# We'll use the normal probability distribution function -
# add the options for it to our parameter list
# "x" value in PDF is observed value
var$x<-"Radius"

# Mean in normal PDF
var$mean<-"predicted"
var$sd<-0.815585

# Have it calculate log likelihood
var$log<-TRUE

result<-likeli(model, par, var, dataset, dnorm)
cat("This result should be -119.969: ", result)

# Try other values for a and b
par$a<-0
par$b<-0.115
var$sd<-0.962325
result<-likeli(model, par, var, dataset, dnorm)
cat("This result should be -139.382: ", result)

par$a<-2.36
par$b<-0
var$sd<-1.094228
result<-likeli(model, par, var, dataset, dnorm)
cat("This result should be -148.89: ", result)

# Do it through predicted_results
predicted <- predicted_results(model, par, var, dataset)
result <- sum(dnorm(x=dataset$Radius, mean=predicted, sd=var$sd,log=TRUE))
cat("This result should be -148.89: ", result)

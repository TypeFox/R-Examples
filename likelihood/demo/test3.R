##
## TEST 3 - calculation of negative log likelihood
## with variance of probability distribution function
## a function of DBH
## Model: Radius = a + b * DBH
## Variance: Variance = c * DBH (we'll take the square
## root to make standard deviation)
## Our parameters are a, b, and c
## Dataset: included crown_rad dataset
##

library(likelihood)

# Set up our dataset
data(crown_rad)
dataset <- crown_rad

# Create our model function
model <- function (a, b, DBH, species) {
  a[species] + b * DBH
}

# Create our variance function
varfun <- function(c, DBH) {sqrt(c * DBH)}

# Create a place to put our parameters and
# set initial values for a and b, and indicate
# that DBH comes from the column marked "DBH"
# in the dataset
par <- list(a = c(1.12, 1.0, 1.5), b = 0.07, c = 0.05)

# Create a place for all other values, starting with
# values from the dataset
var <- list(species="Species", DBH = "DBH")

# We'll use the normal probability distribution function -
# add the options for it to our parameter list
# "x" value in PDF is observed value
var$x<-"Radius"

# Mean in normal PDF
var$mean<-"predicted"
var$sd<-varfun

# Have it calculate log likelihood
var$log<-TRUE

result<-likeli(model, par, var, dataset, dnorm)

cat("This result should be -138.5: ", result)

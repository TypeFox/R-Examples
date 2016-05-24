##
## TEST 2 - calculation of negative log likelihood
## with multiple functions called from the model
## Model: Height = a + b * N where N = r * X
## and a = q + p * DBH
## Parameters are b, q, and r
## Dataset: included from_sortie dataset
##

library(likelihood)

# Set up our dataset
data(from_sortie)

# Create our model function
model <- function (a, b, N) {a + b * N}
# Put in a default which we will not provide in
# the parameter list
afun <- function(q, p=1.3, DBH) {q + p * DBH}

# Create a place to put our parameters and
# set initial values for b, q, and r
par<-list(b = 3.4, q = -6.5, r = 0.01)

# Create the N function
nfun <- function(r, X) {r * X}

# Set up our var list, with non-parameters information.
# Say that a and N are function results
var <- list( a = afun, N = nfun)

# Put in the parameters for sumneigh
var$target_data <- from_sortie
var$DBH <- "DBH"
var$X <- "X"

# We'll use the normal probability distribution function -
# add the options for it to our parameter list
# "x" value in PDF is observed value
var$x<-"Height"

# Mean in normal PDF
var$mean<-"predicted"
var$sd<-4.123106

# Have it calculate log likelihood
var$log<-TRUE

result<-likeli(model, par, var, from_sortie, dnorm)

cat("This result should be -429.678: ", result)

##
## Simulated annealing to find maximum likelihood estimates of a and b in
## the following model:
## Model: Radius = a + b * DBH
## Dataset: included crown_rad dataset
## This is a global, unbounded search.
##

library(likelihood)

# Set up our dataset
data(crown_rad)
dataset <- crown_rad

# Create our model function
modelfun <- function (a, b, DBH) {a + b * DBH}

# Create a place to put our parameters and
# set initial values for a and b, and indicate
# that DBH comes from the column marked "DBH"
# in the dataset
par <- list(a = 0, b = 0)

# Create a place to put all the other values
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

results<-anneal(model = modelfun, par = par, var = var, 
  source_data = dataset, pdf = dnorm, dep_var = "Radius", max_iter = 20000)

cat("Max likelihood should be about -119.7454.  Result: ", results$max_likeli, "\n")
cat("Best a should be about 1.117194.  Result: ", results$best_pars$a, "\n")
cat("Best b should be about 0.06755417.  Result: ", results$best_pars$b, "\n")
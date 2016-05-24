##
## ANNEALING 1: simulated annealing to maximize negative log
## likelihood for the following:
## Model: Radius = a + b * DBH
## Dataset: included crown_rad dataset
##

library(likelihood)

# Set up our dataset
data(crown_rad)
dataset = crown_rad

# Create our model function
modelfun = function (a, b, DBH) {a + b * DBH}

# Create a place to put our parameters and
# set initial values for a and b
par = list(a = 0, b = 0)

# Create a place to put all the other values, and indicate
# that DBH comes from the column marked "DBH"
# in the dataset
var = list(DBH = "DBH")

# Set bounds
par_lo=list(a = 0, b = 0)
par_hi=list(a = 50, b = 50)

# We'll use the normal probability distribution function 
# add the options for it to our parameter list
# "x" value in PDF is observed value
var$x="Radius"

# Mean in normal PDF: use the reserved word "predicted"
# to indicate that this is where the output from the scientific 
# model function goes
var$mean="predicted"
var$sd=0.815585

# Have it calculate log likelihood
var$log=TRUE

results=anneal(model = modelfun, par = par, var = var, 
  source_data = dataset, par_lo = par_lo, par_hi = par_hi,
  pdf = dnorm, dep_var = "Radius", max_iter = 20000)
  
cat("Max likelihood should be about negative 119.7454.  Result: ", results$max_likeli, "\n")
cat("Best a should be about 1.117194.  Result: ", results$best_pars$a, "\n")
cat("Best b should be about 0.06755417.  Result: ", results$best_pars$b, "\n")
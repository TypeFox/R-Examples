##
## Annealing run with variance of probability distribution
## function a function of DBH and sub-indexing on parameters
## Model: Radius = a + b * DBH
## Variance: Variance = c * DBH (we'll take the square
## root to make standard deviation)
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
par<-list(a = c(0.8, 1.0, 1.5), b = 0.07, c = 0.05)
var<-list(species="Species", DBH = "DBH")

par_lo <- list( a = c(0, 0, 0), b = 0, c = 0)
par_hi <- list( a = c(2, 2, 2), b = 10, c = 10)

# We'll use the normal probability distribution function -
# add the options for it to our parameter list
# "x" value in PDF is observed value
var$x<-"Radius"

# Mean in normal PDF
var$mean<-"predicted"
var$sd<-varfun

# Have it calculate log likelihood
var$log<-TRUE

results<-anneal(model=model, par=par, var=var, source_data=dataset, 
par_lo=par_lo, par_hi=par_hi, pdf=dnorm, dep_var="Radius",
max_iter=20000)

cat("Max likelihood should be about -110.163.  Result: ", results$max_likeli, "\n")
cat("Best a1 should be about 0.503.  Result: ", results$best_pars$a[1], "\n")
cat("Best a2 should be about 0.945.  Result: ", results$best_pars$a[2], "\n")
cat("Best a3 should be about 0.744.  Result: ", results$best_pars$a[3], "\n")
cat("Best b should be about 0.088.  Result: ", results$best_pars$b, "\n")
cat("Best c should be about 0.0375.  Result: ", results$best_pars$c, "\n")
# This file shows the approximator package used on a simple 1-d test
# case.  It generates some data randomly, then attempts to infer what
# the parameters used to generate that data are.  One can then compare
# the estimates with the true values.  


# load  the libraries:
library(approximator)
library(emulator)


# set seed:
set.seed(0)

# First a design matrix:
D1.1d <- matrix(1:6)



# And a subsets object
source("subsets_1d.R")

# and some basis functions:
"basis.1d" <- 
function (x) 
{
  out <- cbind(1,x)
  colnames(out) <- c("const","x")
  return(out)
}


# create a hyperparameter function:
source("hpafun_1d.R")

#...and call it, to creat hyperparameter object hpa.1d:
hpa.1d <- hpa.fun.1d(1:9)


# Now a function that creates data:
source("datamaker_1d.R")
z.1d <- generate.1d.observations(D1=D1.1d, subsets=subsets.1d, basis.fun=basis.1d, hpa=hpa.1d, betas = NULL, export.truth=FALSE)


# Now some checks.  First, look at H:
jj.H <- H.fun.app(D1=D1.1d, subsets=subsets.1d , basis=basis.1d , hpa=hpa.1d)
#  Look at jj.H and verify that it is right.




# Now optimize the hyperparameters:
a1 <- opt.1(D=D1.1d , z=z.1d , basis=basis.1d , subsets=subsets.1d , hpa=hpa.1d)
a2 <- opt.gt.1(level=2 , D=D1.1d , z=z.1d , basis=basis.1d , subsets=subsets.1d , hpa=hpa.1d)

# And use the second-level optimized hyperpareters (ie a2) to give the emulator mean:

jj.ans <- mdash.fun(3,D1=D1.1d,subsets=subsets.1d,hpa=a2,z=z.1d,basis=basis.1d)

# (preceding line gives a weird error, under MacOSX 10.5.6, with R-GUI
# 1.28 but no error when running R from the command line).




# And its variance:
jj.var <- c_fun(x=as.matrix(4),xdash=as.matrix(5),subsets=subsets.1d,hpa=hpa.1d)

#####################################
## Load stuff for testing purposes ##
#####################################

library(hdi)

data(riboflavin)

x <- riboflavin[,-1]
y <- riboflavin[,1]

dim(x)
##- [1]   71 4088
length(y)
##- [1] 71

doExtras <- interactive() # i.e., FALSE for  routine R CMD check

p. <- if(doExtras) 50 else 16 # smaller for speed
p.

x.use <- x[,1:p.]

######################
## Lasso projection ##
######################

## set seed because of cv
set.seed(3) ; fit.lasso  <- lasso.proj(x = x.use, y = y)
## Check standardization, i.e., equivariance :
set.seed(3) ; fit.lasso2 <- lasso.proj(x = 2 + 4 * x.use, y = y)

## verbose
ncores <- if(.Platform$OS.type == "windows") 1 else getOption("mc.cores", 2L)
  
set.seed(3) ; fit.tmp  <- lasso.proj(x = x.use, y = y, verbose = TRUE)
set.seed(3) ; fit.tmp2 <- lasso.proj(x = x.use, y = y,
                                     parallel = TRUE, ncores = ncores,
                                     verbose = TRUE)

## confidence intervals
ci.lasso  <- confint(fit.lasso,  level = 0.95)
ci.lasso2 <- confint(fit.lasso2, level = 0.95)

stopifnot(
    all.equal(fit.lasso$pval, fit.lasso2$pval)
   ,
    all.equal(c(0,0), range(fit.lasso$bhat / fit.lasso2$bhat - 4))
   ,
    all.equal(ci.lasso, ci.lasso2 * 4)
   ,
    TRUE)

if(!doExtras) {
    stopifnot(
        all.equal(as.vector(fit.lasso$bhat),
                  c(0.54650099, -0.64364814, 0.079821945, 0.26406221, -0.21405501,
                    -0.63576549, -0.095448048, 0.40801737, -1.2194818, -0.11113313,
                    0.3474404, 1.1425587, -0.54460967, 0.45298509, -0.31922868, 0.42184791),
                  tol = 4e-7)# 1e-8
        ,
        all.equal(unname(ci.lasso),
                  matrix(c(-0.186587, -1.88574, -1.03822, -0.274364, -0.808168, -1.30643,
                           -0.994648, -0.446016, -2.29657, -1.23525, -0.728214, 0.256838,
                           -1.297,    -0.416467, -1.16564, -0.593977, 1.27959, 0.598448,
                           1.19786,   0.802489, 0.380058, 0.0348979, 0.803752, 1.26205,
                           -0.142394, 1.01298,  1.42309, 2.02828, 0.207783, 1.32244,
                           0.527183, 1.43767),
                         16, 2), tol = 5e-5)
    )
}

##################################
## Code of example of help page ##
##################################

x <- matrix(rnorm(100*20), nrow = 100, ncol = 10)
y <- x[,1] + x[,2] + rnorm(100)
fit.lasso <- lasso.proj(x, y)
which(fit.lasso$pval.corr < 0.05) # typically: '1' and '2' and no other

## Group-wise testing of the first two coefficients
fit.lasso$groupTest(1:2)

## Hierarchical testing using distance matrix based on
## correlation matrix
out.clust <- fit.lasso$clusterGroupTest()
plot(out.clust)

## Fit the lasso projection method without doing the preparations
## for group testing (saves time and memory)
fit.lasso.faster <- lasso.proj(x, y, suppress.grouptesting = TRUE)

## Use the scaled lasso for the initial estimate
fit.lasso.scaled <- lasso.proj(x, y, betainit = "scaled lasso")
which(fit.lasso.scaled$pval.corr < 0.05)

## Use a robust estimate for the standard error
fit.lasso.robust <- lasso.proj(x, y, robust = TRUE)
which(fit.lasso.robust$pval.corr < 0.05)


## ----mint_create, eval=TRUE, results="hide"------------------------------
rm(list = ls()) # Clear workspace
library(MInt)

# Specify path to the design matrix
dFile <- system.file("extdata", "x.txt", package="MInt");

# Specify path to the response matrix
rFile <- system.file("extdata", "y.txt", package="MInt");


# Create the MInt model object
m <- mint(y=rFile, x=dFile, fmla = ~ feature1 + feature2)

## ----mint_estimate, eval=TRUE--------------------------------------------
m <- estimate(m)

## ----mint_compare, eval=TRUE---------------------------------------------
P_true <- as.matrix( read.csv(system.file("extdata", "P_true.txt", package="MInt"), header=FALSE) );

-cov2cor(P_true)
-cov2cor(m$param$P)


library( "miscTools" )

# Construct a simple OLS regression:
set.seed( 123 )
x1 <- runif(100)
x2 <- runif(100)
y <- 3 + 4*x1 + 5*x2 + rnorm(100)
m <- lm(y~x1+x2)  # estimate it
nObs(m)
nParam(m)

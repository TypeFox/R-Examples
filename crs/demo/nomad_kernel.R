require(crs)

set.seed(42)

## Example - simulated data

n <- 10000
num.eval <- 50
x1 <- runif(n)
x2 <- runif(n)
z <- rbinom(n,1,.5)
dgp <- cos(2*pi*x1)+sin(2*pi*x2)+z
z <- factor(z)
y <- dgp + rnorm(n,sd=.5)

## Estimate a model with specified degree and bandwidth

model.kernel <- crs(y~x1+x2+z,
                    degree=c(5,5),
                    segments=c(1,1),                    
                    lambda=c(0.1),
                    basis="additive",
                    complexity="degree-knots",
                    cv="none",
                    kernel=TRUE)

summary(model.kernel)

## Use initial value starting points in degree, segments, and lambda

model.kernel <- crs(y~x1+x2+z,
                    degree=c(5,5),
                    segments=c(1,1),
                    lambda=c(0.1),
                    basis="auto",
                    complexity="degree-knots",
                    cv="nomad",
                    kernel=TRUE)

summary(model.kernel)

## We could compare with exhaustive search, but that takes some time
## as numerical search is conducted for each degree/segment combination.

#model.kernel.multiple <- crs(y~x1+x2+z,
#                             basis="auto",
#                             complexity="degree-knots",
#                             cv="exhaustive",
#                             kernel=TRUE,
#                             nmulti=10)


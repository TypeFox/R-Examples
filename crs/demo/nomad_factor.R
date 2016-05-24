require(crs)

set.seed(42)

n <- 10000
num.eval <- 50

x1 <- runif(n,-5,5)
x2 <- runif(n,-5,5)

dgp <- sin(sqrt(x1^2+x2^2))/sqrt(x1^2+x2^2)

y <- dgp + rnorm(n,sd=.1)

## No initial points, generate by sample

model.nomad <- crs(y~x1+x2,
                   basis="auto",
                   cv="nomad",
                   complexity="degree-knots",
                   knots="uniform",
                   deriv=1,
                   cv.func="cv.aic")

summary(model.nomad)

## Use degree and segments as initial points

model.nomad <- crs(y~x1+x2,
                   basis="auto",
                   degree=c(3,3),
                   segments=c(2,4),
                   cv="nomad",
                   complexity="degree-knots",
                   knots="uniform",
                   deriv=1,
                   cv.func="cv.aic")

summary(model.nomad)

## Compare with exhaustive search (cv="exhaustive")

model <- crs(y~x1+x2,
             basis="auto",
             cv="exhaustive",
             complexity="degree-knots",
             knots="uniform",
             deriv=1,
             cv.func="cv.aic")

summary(model)

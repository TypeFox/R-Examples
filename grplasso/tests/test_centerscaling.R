library(grplasso)

set.seed(7198)

n <- 100
p <- 10

x <- matrix(runif(n * p, min = -2.5, max = 2.5), nrow = n, ncol = p)
y <- 4 * sin(x[,1]) + x[,2]^2 + rnorm(n)

x.new <- matrix(runif(n * p, min = -2.5, max = 2.5), nrow = n, ncol = p)

shift <- 10
scale <- 10

## Use Lasso shift, i.e. group-size = 1
fit <- grplasso(cbind(1,x), y, index = c(NA, 1:10),
                lambda = c(50, 10), center = TRUE, standardize = TRUE, 
                model = LinReg(), control = grpl.control(trace = 0))

## Rescale
x.resc <- shift + scale * x
fit.resc <- grplasso(cbind(1,x.resc), y, index = c(NA, 1:10),
                     lambda = c(50, 10), center = TRUE, standardize = TRUE, 
                     model = LinReg(), control = grpl.control(trace = 0))

## Compare estimators, without intercept
stopifnot(all.equal(coef(fit)[-1,], coef(fit.resc)[-1,] * scale))

## Check intercepts
mu.x <- apply(x, 2, mean)
int  <- mean(y) - apply(coef(fit)[-1,] * mu.x, 2, sum)

mu.x.resc <- apply(x.resc, 2, mean)
int.resc  <- mean(y) - apply(coef(fit.resc)[-1,] * mu.x.resc, 2, sum)

stopifnot(all.equal(int, coef(fit)[1,], tol = 10^-7))
stopifnot(all.equal(int.resc, coef(fit.resc)[1,], tol = 10^-7))
          
## Compare predictions
stopifnot(all.equal(predict(fit, newdata = cbind(1,x.new)),
                    predict(fit.resc,
                            newdata = cbind(1, shift + scale * x.new))))

## Check situation where we have an unpenalized covariate and re-scaling

index <- c(NA, 1:9, NA)

fit.a <- grplasso(cbind(1,x), y, index = index,
                  lambda = c(50, 10), center = TRUE, standardize = TRUE, 
                  model = LinReg(), control = grpl.control(trace = 0))

fit.b <- grplasso(cbind(1,x.resc), y, index = index,
                  lambda = c(50, 10), center = TRUE, standardize = TRUE, 
                  model = LinReg(), control = grpl.control(trace = 0))

stopifnot(all.equal(coef(fit.a)[10,], scale * coef(fit.b)[10,]))

########################################################################
##                                                                    ##
## Check whether every case is running, including function lambda.max ##
##                                                                    ##
########################################################################

## center = TRUE & unpenalized intercept

x.use <- cbind(1, x)
index <- c(NA, 1:10)

lambda.max <- lambdamax(x.use, y, index, model = LinReg())
lambda     <- lambda.max * c(1, 0.1)

fit1 <- grplasso(x.use, y, index, model = LinReg(), lambda = lambda,
                 center = TRUE)

## center = TRUE & penalized intercept

x.use <- cbind(1, x)
index <- c(99, 1:10)

lambda.max <- lambdamax(x.use, y, index, model = LinReg())
lambda     <- lambda.max * c(1, 0.1)

fit2 <- grplasso(x.use, y, index, model = LinReg(), lambda = lambda,
                 center = TRUE)

## center = TRUE & *no* intercept
x.use <- cbind(x)
index <- c(1:10)

lambda.max <- lambdamax(x.use, y, index, model = LinReg())
lambda     <- lambda.max * c(1, 0.1)

fit3 <- grplasso(x.use, y, index, model = LinReg(), lambda = lambda,
                 center = TRUE)

## center = FALSE & *no* intercept
x.use <- cbind(x)
index <- c(1:10)

lambda.max <- lambdamax(x.use, y, index, model = LinReg(),
                        center = FALSE)
lambda     <- lambda.max * c(1, 0.1)

fit4 <- grplasso(x.use, y, index, model = LinReg(), lambda = lambda,
                 center = FALSE)

## Check whether fit3 and fit4 are the same
stopifnot(all.equal(coef(fit3), coef(fit4)))


## center = FALSE & standardize = TRUE
x.use <- cbind(1,x)
index <- c(NA, 1:10)

lambda.max <- lambdamax(x.use, y, index, model = LinReg(),
                        center = FALSE)
lambda     <- lambda.max * c(1, 0.1)

fit5 <- grplasso(x.use, y, index, model = LinReg(), lambda = lambda,
                 center = FALSE)


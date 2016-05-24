###
# test functionality of multivariate mstop

#detach("package:gamboostLSS", unload = TRUE)
require(gamboostLSS)

### create some data first
set.seed(1907)
x1 <- rnorm(1000)
x2 <- rnorm(1000)
x3 <- rnorm(1000)
x4 <- rnorm(1000)
x5 <- rnorm(1000)
x6 <- rnorm(1000)
mu    <- exp(1.5 +1 * x1 +0.5 * x2 -0.5 * x3 -1 * x4)
sigma <- exp(-0.4 * x3 -0.2 * x4 +0.2 * x5 +0.4 * x6)
y <- numeric(1000)
for( i in 1:1000)
    y[i] <- rnbinom(1, size = sigma[i], mu = mu[i])
dat <- data.frame(x1, x2, x3, x4, x5, x6, y)


## check if model with different mstops is what we expect
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 10),
                     center = TRUE)
mstop(model)
model2 <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = list(mu = 10, sigma = 20)),
                     center = TRUE)
mstop(model2)

f1 <- fitted(model, parameter = "mu", type = "response")
f2 <- fitted(model, parameter = "sigma", type = "response")
model3 <- glmboost(y ~ ., family = NBinomialSigma(mu = f1, sigma = f2,
                          stabilization = "none"),
                   data = dat,
                   control = boost_control(mstop = 10),
                   center = TRUE)

tmp <- coef(model3) + coef(model)$sigma
stopifnot(max(abs(tmp - coef(model2)$sigma)) < sqrt(.Machine$double.eps))

# Plot
layout(matrix(c(1:4, 6, 5), byrow = TRUE, ncol = 2))
plot(model, xlim = c(0,20), ylim = range(sapply(coef(model2), range)))
plot(model2, xlim = c(0, 20), ylim = range(sapply(coef(model2), range)))
plot(model, xlim = c(0,20), ylim = range(sapply(coef(model2), range)),
     parameter = "sigma")
cp <- coef(model3, aggregate = "cumsum")
cp <- matrix(unlist(cp), nrow = length(cp), byrow = TRUE)
cp <- cp + coef(model)$sigma
cp <- cbind(coef(model)$sigma, cp)
cf <- cp[, ncol(cp)]
col <- hcl(h = 40, l = 50, c = abs(cf)/max(abs(cf)) * 490)
matlines(10:20, t(cp), type = "l", xlab = "Number of boosting iterations",
         ylab = "Coefficients", col = col)



### check subset method
ms <- list(mu = 10, sigma = 20)
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                      control = boost_control(mstop = ms, trace = TRUE),
                      center = TRUE)

model[c(20, 30)]   # check if two values can be specified
ms <- list(mu = 20, sigma = 30)
modela <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = ms, trace = TRUE),
                     center = TRUE)
stopifnot(max(abs(coef(model)[[1]] - coef(modela)[[1]]))
          < sqrt(.Machine$double.eps))
stopifnot(max(abs(coef(model)[[2]] - coef(modela)[[2]]))
          < sqrt(.Machine$double.eps))

model[40]          # check if one value can be specified
mstop(model)
modelb <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                       control = boost_control(mstop = 40, trace = TRUE),
                       center = TRUE)
stopifnot(all.equal(risk(model), risk(modelb)))

model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = 10, trace = TRUE),
                     center = TRUE)
model[20]
model2 <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                      control = boost_control(mstop = 20, trace = TRUE),
                      center = TRUE)
stopifnot(all.equal(risk(model), risk(model2)))

ms <- list(mu = 10, sigma = 20)
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = ms, trace = TRUE),
                     center = TRUE)
model[c(5,10)]
ms <- list(mu = 5, sigma = 10)
model2 <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                      control = boost_control(mstop = ms, trace = TRUE),
                      center = TRUE)
stopifnot(all.equal(risk(model), risk(model2)))


### check subset method where only bigger mstop-value is touched
# increase model
ms <- list(mu = 10, sigma = 20)
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = ms, trace = TRUE),
                     center = TRUE)
model[c(10,25)]
mstop(model)
ms <- list(mu = 10, sigma = 25)
model2 <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                      control = boost_control(mstop = ms, trace = TRUE),
                      center = TRUE)
stopifnot(all.equal(risk(model), risk(model2)))

# reduce model
ms <- list(mu = 10, sigma = 20)
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = ms, trace = TRUE),
                     center = TRUE)
model[c(10,15)]
mstop(model)
ms <- list(mu = 10, sigma = 15)
model2 <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                      control = boost_control(mstop = ms, trace = TRUE),
                      center = TRUE)
stopifnot(all.equal(risk(model), risk(model2)))

# reduce such that mu needs to be touced again
ms <- list(mu = 10, sigma = 20)
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                     control = boost_control(mstop = ms, trace = TRUE),
                     center = TRUE)
model[c(10,9)]
mstop(model)
ms <- list(mu = 10, sigma = 9)
model2 <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                      control = boost_control(mstop = ms, trace = TRUE),
                      center = TRUE)
stopifnot(all.equal(risk(model), risk(model2)))

### check multiple values of nu

nus <- list(mu = 0, sigma = 0.2)
model <- glmboostLSS(y ~ ., families = NBinomialLSS(), data = dat,
                      control = boost_control(mstop = 10, nu = nus, trace = TRUE),
                      center = TRUE)
stopifnot(all(coef(model)[[1]] == 0))
stopifnot(any(coef(model)[[2]] != 0))

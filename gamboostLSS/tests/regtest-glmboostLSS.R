###
# check glmboostLSS()

require("gamboostLSS")
require("gamlss")

set.seed(1907)
n <- 5000
x1  <- runif(n)
x2 <- runif(n)
mu <- 2 -1*x1 - 3*x2
sigma <- exp(-1*x1 + 3*x2)
df <- exp(1 + 3*x1 + 1*x2)
y <- rTF(n = n, mu = mu, sigma = sigma, nu = df)

### check subset method
model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                     control = boost_control(mstop = 10),
                     center = TRUE)
model2 <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                          control = boost_control(mstop = 20),
                          center = TRUE)
model[20]
stopifnot(all.equal(coef(model),coef(model2)))
stopifnot(length(coef(model2, aggregate = "none")[[1]][[1]]) ==
          length(coef(model, aggregate = "none")[[1]][[1]]))
stopifnot(length(coef(model2, aggregate = "none")[[1]][[1]]) == 20)

model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                     control = boost_control(mstop = 10),
                     center = TRUE)
model2[10]
stopifnot(all.equal(coef(model),coef(model2)))
stopifnot(length(coef(model2, aggregate = "none")[[1]][[1]]) ==
          length(coef(model, aggregate = "none")[[1]][[1]]))
stopifnot(length(coef(model2, aggregate = "none")[[1]][[1]]) == 10)

### check trace
model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(),
                     control = boost_control(mstop = 10, trace =TRUE),
                     center = TRUE)
model[100]

### check formula-interface with lists
set.seed(1907)
n <- 5000
x1  <- runif(n)
x2 <- runif(n)
mu <- 2 - 3*x2
sigma <- exp(-1*x1 + 3*x2)
df <- exp(1 + 3*x1)
y <- rTF(n = n, mu = mu, sigma = sigma, nu = df)
model <- glmboostLSS(list(mu = y ~ x2,
                          sigma = y ~ x1 + x2,
                          df = y ~ x1),
                     families = StudentTLSS(),
                     control = boost_control(mstop = 10, trace =TRUE),
                     center = TRUE)

stopifnot(all.equal(lapply(coef(model, which = ""), function(x) names(x)[-1]),
                    list(mu = "x2", sigma = c("x1", "x2"), df = "x1")))

model <- glmboostLSS(list(mu = y ~ x2,
                          df = y ~ x1,
                          sigma = y ~ x1 + x2),
                     families = StudentTLSS(),
                     control = boost_control(mstop = 10, trace =TRUE),
                     center = TRUE)

stopifnot(all.equal(lapply(coef(model, which = ""), function(x) names(x)[-1]),
                    list(mu = "x2", sigma = c("x1", "x2"), df = "x1")))

### check formula-interface with lists and different responses
### (not really sensible with the current families)
y2 <- y + 1
model2 <- glmboostLSS(list(mu = y ~ x2,
                           sigma = y ~ x1 + x2,
                           df = y2 ~ x1),
                      families = StudentTLSS(),
                      control = boost_control(mstop = 10, trace =TRUE),
                      center = TRUE)
sapply(model, function(comps) comps$offset)
sapply(model2, function(comps) comps$offset)
stopifnot((model2[[1]]$offset - model[[1]]$offset) < sqrt(.Machine$double.eps))
stopifnot((model2[[2]]$offset - model[[2]]$offset) < sqrt(.Machine$double.eps))
stopifnot(model2[[3]]$offset != model[[3]]$offset)

### even better check for offset-issue when different responses are used
set.seed(0804)
x1 <- runif(1000)
x2 <- runif(1000)
x3 <- runif(1000)
x4 <- runif(1000)
mu    <- 1.5 +1 * x1 +4 * x2
sigma <- exp(1 - 0.2 * x3 -0.4 * x4)
y <- rnorm(mean=mu, sd=1, n=length(mu))
y2 <- rnorm(mean=0, sd=sigma, n=length(sigma))
dat <- data.frame(x1, x2, x3, x4, y, y2)
model <- glmboostLSS(formula=list(mu=y~x1+x2, sigma=y2~x3+x4),
                     families=GaussianLSS(), data=dat)
## offset for mu must be equal to the mean of y per default
stopifnot((model$mu$offset - mean(y)) < sqrt(.Machine$double.eps))
## offset for sigma must be equal to the log of the standard deviation of y2 per
## default
stopifnot((model$sigma$offset - log(sd(y2))) < sqrt(.Machine$double.eps))


### check weights interface
set.seed(1907)
n <- 2500
x1  <- runif(n)
x2 <- runif(n)
mu <- 2 - 3*x2
sigma <- exp(-1*x1 + 3*x2)
df <- exp(1 + 3*x1)
y <- rTF(n = n, mu = mu, sigma = sigma, nu = df)
dat <- data.frame(x1, x2, y)
dat2 <- rbind(dat, dat) # data frame with duplicate entries

model <- glmboostLSS(list(mu = y ~ x2,
                          sigma = y ~ x1 + x2,
                          df = y ~ x1),
                     data = dat, weights = rep(2, nrow(dat)),
                     families = StudentTLSS(),
                     control = boost_control(mstop = 10, trace =TRUE),
                     center = TRUE)

model2 <- glmboostLSS(list(mu = y ~ x2,
                           sigma = y ~ x1 + x2,
                           df = y ~ x1),
                      data = dat2, families = StudentTLSS(),
                      control = boost_control(mstop = 10, trace =TRUE),
                      center = TRUE)

stopifnot(all.equal(coef(model), coef(model2)))

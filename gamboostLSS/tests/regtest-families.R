###
# check families

require("gamboostLSS")
require("gamlss")


### check families with only one offset specified (other to choose via optim)
set.seed(1907)
n <- 5000
x1  <- runif(n)
x2 <- runif(n)
mu <- 2 -1*x1 - 3*x2
sigma <- exp(-1*x1 + 3*x2)
df <- exp(1 + 3*x1 + 1*x2)
y <- rTF(n = n, mu = mu, sigma = sigma, nu = df)
data <- data.frame(y = y, x1 = x1, x2 = x2)
rm("y", "x1", "x2")

model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(mu = 0.5),
                     data = data,
                     control = boost_control(mstop = 10), center = TRUE)
coef(model)

model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(sigma = 1),
                     data = data,
                     control = boost_control(mstop = 10), center = TRUE)
coef(model)

model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(df = 4),
                     data = data,
                     control = boost_control(mstop = 10), center = TRUE)
coef(model)

model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(mu = 0.5, df = 1),
                     data = data,
                     control = boost_control(mstop = 10), center = TRUE)
coef(model)

### multivariate minimum for offset
loss <- function(df, sigma,y, f){
    -1 * (lgamma((df+1)/2) - log(sigma) - lgamma(1/2) - lgamma(df/2) - 0.5 *
          log(df) - (df+1)/2 * log(1 + (y-f)^2 / (df * sigma^2)))
}
riski <- function(x, y, w = rep(1, length(y))){
    f <- x[[1]]
    df <- exp(x[[2]])
    sigma <- exp(x[[3]])
    sum(w * loss(y = y, f = f, df = df, sigma = sigma))
}

res <- optim(fn = riski, par = c(0, 1, 1), y = data$y, w = rep(1, length(data$y)))$par
model <- glmboostLSS(y ~ x1 + x2, families = StudentTLSS(mu = res[1], sigma =
                                  exp(res[2]), df = exp(res[3])),
                     data = data,
                     control = boost_control(mstop = 10), center = TRUE)
model[100]
coef(model)

### check survival families

x1 <- runif(1000)
x2 <- runif(1000)
x3 <- runif(1000)
w <- rnorm(1000)

time <-  exp(3 + 1*x1 +2*x2  + exp(0.2 * x3) * w)
status <- rep(1, 1000)

## check if error messages are correct
try(glmboost(time ~ x1 + x2 + x3, family = Lognormal(),
             control = boost_control(trace = TRUE), center = TRUE))
try(glmboostLSS(time ~ x1 + x2 + x3, families = LogNormalLSS(),
                control = boost_control(trace = TRUE), center = TRUE))
require("survival")
try(glmboostLSS(list(mu = Surv(time, status) ~ x1 + x2 + x3,
                     sigma = time ~ x1 + x2 + x3), families = LogNormalLSS(),
                control = boost_control(trace = TRUE), center = TRUE))

## check results

# LogNormalLSS()
(m1 <- survreg(Surv(time, status) ~ x1 + x2 + x3, dist="lognormal"))
model <- glmboostLSS(Surv(time, status) ~ x1 + x2 + x3, families = LogNormalLSS(),
                     control = boost_control(trace = TRUE), center = TRUE)
stopifnot(sum(abs(coef(model, off2int = TRUE)[[1]] - c(3, 1, 2, 0)))
          < sum(abs(coef(m1) - c(3, 1, 2, 0))))
stopifnot(sum(abs(coef(model, off2int = TRUE)[[2]] - c(0, 0, 0, 0.2))) < 0.25)

# LogLogLSS()
etamu <- 3 + 1*x1 +2*x2
etasi <- exp(rep(0.2, 1000))
for (i in 1:1000)
    time[i] <- exp(rlogis(1, location = etamu[i], scale = etasi[i]))
status <- rep(1, 1000)
(m1 <- survreg(Surv(time, status) ~ x1 + x2 + x3, dist="loglogistic"))
model <- glmboostLSS(Surv(time, status) ~ x1 + x2 + x3, families = LogLogLSS(),
                     control = boost_control(trace = TRUE), center = TRUE)
model[350]
stopifnot(sum(abs(coef(model, off2int = TRUE, which ="")[[1]] - c(3, 1, 2, 0)))
          < sum(abs(coef(m1) - c(3, 1, 2, 0))))
stopifnot(sum(abs(coef(model, off2int = TRUE)[[2]] - c(0.2, 0, 0, 0))) < 0.25)

# WeibullLSS()
etamu <- 3 + 1*x1 +2*x2
etasi <- exp(rep(0.2, 1000))
status <- rep(1, 1000)
time <- rep(NA, 1000)
for (i in 1:1000)
    time[i] <- rweibull(1, shape = exp(- 0.2), scale = exp(etamu[i]))
(m1 <- survreg(Surv(time, status) ~ x1 + x2 + x3, dist="weibull"))
model <- glmboostLSS(Surv(time, status) ~ x1 + x2 + x3,
                     families = WeibullLSS(),
                     control = boost_control(trace = TRUE), center = TRUE)
model[300]
stopifnot(sum(abs(coef(model, off2int = TRUE, which ="")[[1]] - c(3, 1, 2, 0)))
          < sum(abs(coef(m1) - c(3, 1, 2, 0))))
stopifnot(sum(abs(coef(model, off2int = TRUE)[[2]] - c(0.2, 0, 0, 0))) < 0.4)


### Check that "families"-object contains a response function
NBinomialMu2 <- function(...){
    RET <- NBinomialMu(...)
    RET@response <- function(f) NA
    RET
}

NBinomialSigma2 <- function(...){
    RET <- NBinomialSigma(...)
    RET@response <- function(f) NA
    RET
}

NBinomialLSS2 <- function(mu = NULL, sigma = NULL){
    if ((!is.null(sigma) && sigma <= 0) || (!is.null(mu) && mu <= 0))
        stop(sQuote("sigma"), " and ", sQuote("mu"),
             " must be greater than zero")
    RET <- Families(mu = NBinomialMu2(mu = mu, sigma = sigma),
                        sigma = NBinomialSigma2(mu = mu, sigma = sigma))
    return(RET)
}

try(NBinomialLSS2())

#detach("package:gamboostLSS", unload = TRUE)


### Check that weighed median works well, which is needed if the negative
### gradient is stabilized using MAD
set.seed(122)

## some tests
x <- c(1, 2, 3)
stopifnot(weighted.median(x) == median(x))

## this doesn't work a.t.m.
w <- c(0, 1, 2)
stopifnot(weighted.median(x, w) == median(rep(x, w)))

x <- c(1, 2, 3, 4)
stopifnot(weighted.median(x) == median(x))

w <- c(0, 1, 2, 3)
stopifnot(weighted.median(x, w) == median(rep(x, w)))

w <- rep(0:1, each = 50)
x <- rnorm(100)
stopifnot(weighted.median(x, w) == median(rep(x, w)))

w <- rep(0:1, each = 51)
x <- rnorm(102)
stopifnot(weighted.median(x, w) == median(rep(x, w)))

## hope this is correct:
w <- runif(100, 0, 1)
x <- rnorm(100)
(wm <- weighted.median(x, w))

## check weighted.median with NAs.

# Missing in weights
tmp <- w[2]
w[2] <- NA
stopifnot(is.na(weighted.median(x, w)))
stopifnot(weighted.median(x, w, na.rm = TRUE) == wm)
## not always true but here it is

# Missing in x
w[2] <- tmp
x[5] <- NA
stopifnot(is.na(weighted.median(x, w)))
stopifnot(weighted.median(x, w, na.rm = TRUE) == wm)
## not always true but here it is

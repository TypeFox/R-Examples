### ===== actuar: An R Package for Actuarial Science =====
###
### Demo of the risk theory facilities provided by actuar
###
### AUTHORS: Christophe Dutang, Vincent Goulet <vincent.goulet@act.ulaval.ca>

require(actuar)
require(graphics)


### DISCRETIZATION OF CONTINUOUS DISTRIBUTIONS

## Upper and lower discretization of a Gamma(2, 1) distribution with a
## step (or span, or lag) of 0.5. The value of 'to' is chosen so as to
## cover most of the distribution.
x <- seq(0, qgamma(1 - 1E-6, 2, 1), by = 0.5)
xu <- tail(x, 1)
fu <- discretize(pgamma(x, 2, 1), method = "upper",
                 from = 0, to = xu, step = 0.5)
fl <- discretize(pgamma(x, 2, 1), method = "lower",
                 from = 0, to = xu, step = 0.5)
curve(pgamma(x, 2, 1), xlim = range(x), lwd = 2)
par(col = "blue")
plot(stepfun(head(x, -1), diffinv(fu)), pch = 19, add = TRUE)
par(col = "green")
plot(stepfun(x, diffinv(fl)), pch = 19, add = TRUE)
par(col = "black")

## Discretization with the rounding method, which has the true cdf
## pass through the midpoints of the intervals [x - step/2, x +
## step/2).
fr <- discretize(pgamma(x, 2, 1), method = "rounding",
                 from = 0, to = xu, step = 0.5)
curve(pgamma(x, 2, 1), xlim = range(x), lwd = 2)
par(col = "blue")
plot(stepfun(head(x, -1), diffinv(fr)), pch = 19, add = TRUE)
par(col = "black")

## Local matching of the first moment. This requires a function to
## compute the limited expected value of the true distribution in any
## point.
fb <- discretize(pgamma(x, 2, 1), method = "unbiased",
                 lev = levgamma(x, 2, 1),
                 from = 0, to = xu, step = 0.5)
curve(pgamma(x, 2, 1), xlim = range(x), lwd = 2)
par(col = "blue")
plot(stepfun(x, diffinv(fb)), pch = 19, add = TRUE)
par(col = "black")

all.equal(diff(pgamma(range(x), 2, 1)),
          sum(fb))                      # same total probability
all.equal(levgamma(xu, 2, 1) - xu * pgamma(xu, 2, 1, lower.tail = FALSE),
          drop(crossprod(x, fb)))       # same expected value

## Comparison of all four methods
fu <- discretize(plnorm(x), method = "upper", from = 0, to = 5)
fl <- discretize(plnorm(x), method = "lower", from = 0, to = 5)
fr <- discretize(plnorm(x), method = "rounding", from = 0, to = 5)
fb <- discretize(plnorm(x), method = "unbiased", from = 0, to = 5,
                 lev = levlnorm(x))
curve(plnorm(x), from = 0, to = 5, lwd = 2)
par(col = "blue")
plot(stepfun(0:4, diffinv(fu)), pch = 19, add = TRUE)
par(col = "red")
plot(stepfun(0:5, diffinv(fl)), pch = 19, add = TRUE)
par(col = "green")
plot(stepfun(0:4, diffinv(fr)), pch = 19, add = TRUE)
par(col = "magenta")
plot(stepfun(0:5, diffinv(fb)), pch = 19, add = TRUE)
legend("bottomright", legend = c("upper", "lower", "rounding", "unbiased"),
       col = c("blue", "red", "green", "magenta"), lty = 1, pch = 19,
       text.col = "black")
par(col = "black")


### CALCULATION OF THE AGGREGATE CLAIM AMOUNT DISTRIBUTION

## Calculation of the aggregate claim amount distribution using the
## recursive method (Panjer). Argument 'x.scale' is used to specify
## how much a value of 1 is really worth.
fx.b <- discretize(pgamma(x, 2, 1), from = 0, to = 22, step = 0.5,
                   method = "unbiased", lev = levgamma(x, 2, 1))
Fs.b <- aggregateDist("recursive", model.freq = "poisson",
                      model.sev = fx.b, lambda = 10, x.scale = 0.5)
summary(Fs.b)                           # summary method
knots(Fs.b)                             # support of Fs.b (knots)
Fs.b(knots(Fs.b))                       # evaluation at knots
plot(Fs.b, do.points = FALSE, verticals = TRUE, xlim = c(0, 60)) # graphic
mean(Fs.b)                              # empirical mean
quantile(Fs.b)                          # quantiles

## Convolutions (exact calculation). Requires a vector of
## probabilities for the frequency model. This method can quickly
## become impractical for a large expected number of claims.
pn <- dpois(0:qpois(1-1E-6, 10), 10)
Fs <- aggregateDist("convolution", model.freq = pn, model.sev = fx.b,
                    x.scale = 0.5)
summary(Fs)                             # summary method
knots(Fs)                               # support of Fs (knots)
Fs(knots(Fs))                           # evaluation at knots
plot(Fs, do.points = FALSE, verticals = TRUE, xlim = c(0, 60)) # graphic
mean(Fs)                                # empirical mean
quantile(Fs)                            # quantiles

## Normal approximation. Not hugely useful, but simple to implement...
Fs.n <- aggregateDist("normal", moments = c(20, 60))
summary(Fs.n)                           # summary method
plot(Fs.n, xlim = c(0, 60))             # graphic
mean(Fs.n)                              # true mean
quantile(Fs.n)                          # normal quantiles

## Normal Power II approximation. The approximation is valid for
## values above the expected value only.
Fs.np <- aggregateDist("npower", moments = c(20, 60, 0.516398))
summary(Fs.np)                          # summary method
plot(Fs.np, xlim = c(0, 60))            # truncated graphic

## Simulation method. Function 'simul' is used to simulate the data
## (see the 'simulation' demo for examples).
Fs.s <- aggregateDist("simulation",
                      model.freq = expression(y = rpois(10)),
                      model.sev = expression(y = rgamma(2, 1)),
                      nb.simul = 10000)
summary(Fs.s)                           # summary method
plot(Fs.s, do.points = FALSE, verticals = TRUE, xlim = c(0, 60)) # graphic
mean(Fs.s)                              # empirical mean
quantile(Fs.s)                          # quantiles

## Graphic comparing the cdfs obtained by a few methods.
fx.u <- discretize(pgamma(x, 2, 1), from = 0, to = 22, step = 0.5,
                   method = "upper")
Fs.u <- aggregateDist("recursive", model.freq = "poisson",
                      model.sev = fx.u, lambda = 10, x.scale = 0.5)
fx.l <- discretize(pgamma(x, 2, 1), from = 0, to = 22, step = 0.5,
                   method = "lower")
Fs.l <- aggregateDist("recursive", model.freq = "poisson",
                      model.sev = fx.l, lambda = 10, x.scale = 0.5)
par(col = "black")
plot(Fs.b, do.points = FALSE, verticals = TRUE, xlim = c(0, 60), sub = "")
par(col = "blue")
plot(Fs.u, do.points = FALSE, verticals = TRUE, add = TRUE, sub = "")
par(col = "red")
plot(Fs.l, do.points = FALSE, verticals = TRUE, add = TRUE, sub = "")
par(col = "green")
plot(Fs.s, do.points = FALSE, verticals = TRUE, add = TRUE, sub = "")
par(col = "magenta")
plot(Fs.n, add = TRUE, sub = "")
legend("bottomright",
       legend = c("recursive + unbiased", "recursive + upper",
                  "recursive + lower", "simulation",
                  "normal approximation"),
       col = c("black", "blue", "red", "green", "magenta"),
       lty = 1, text.col = "black", cex = 1.2)
par(col = "black")

## Table of quantiles for the same methods as graphic above.
x <- knots(Fs.l)
m <- which.min(x[round(Fs.l(x), 6) > 0])
M <- which.max(x[round(Fs.u(x), 6) < 1])
x <- x[round(seq.int(from = m, to = M, length = 30))]
round(cbind(x = x,
            Lower = Fs.l(x),
            Unbiased = Fs.b(x),
            Upper = Fs.u(x),
            Simulation = Fs.s(x),
            Normal = Fs.n(x)), 6)


### CALCULATION OF THE ADJUSTMENT COEFFICIENT

## No reinsurance, generalized Erlang claim amounts, inverse gamma
## interarrival times and independence. The adjustment coefficient is
## increasing with the safety loading.
mgf <- function(x) 1/(1 - x) * 2/(2 - x) * 3/(3 - x)
adjCoef(mgf, mgfinvgamma(x, 2, 6/11), 1.1, 1)
adjCoef(mgf, mgfinvgamma(x, 2, 6/11), 1.2, 1)
adjCoef(mgf, mgfinvgamma(x, 2, 6/11), 1.3, 1)

## More sophisticated example: comparison of the effect of dependence
## on the adjustment coefficient in the case of proportional
## reinsurance. Use a Clayton copula with exponential marginals.
rclayton <- function(alpha, n)
{
    val <- cbind(runif(n), runif(n))
    val[, 2] <- (val[, 1]^(-alpha) * (val[, 2]^(-alpha/(alpha + 1)) - 1) + 1)^(-1/alpha)
    val
}
u <- rclayton(2, 1000)             # variates with positive dependence
x <- qexp(u[, 1])                  # claim amounts
w <- qexp(u[, 2])                  # interarrival times

## Premium rate and Lundberg's functions of the retention rate. We
## assume a safety loading of 20% for the insurer and 30% for the
## reinsurer and premium calculated with the expected value principle.
p <- function(a) mean(x)/mean(w) * (1.2 - 1.3 + 1.3 * a)
h <- function(r, a) mean(exp(r * (a * x - p(a) * w)))
R1 <- adjCoef(h = h, upper = 1, reinsurance = "prop", from = 1/3, to = 1)
plot(R1)

## Repeat the above with independent claim amounts and interarrival
## times.
u <- rclayton(1, 1000)             # independent variates
x <- qexp(u[,1])                   # claim amounts
w <- qexp(u[,2])                   # interarrival times
R2 <- adjCoef(h = h, upper = 1, reinsurance = "prop", from = 1/3, to = 1)
plot(R2, add = TRUE, col = "green")
legend("bottomright", legend = c("dependence", "independence"),
       col = c("black", "green"), lty = 1)

## Similar example with excess-of-loss reinsurance.
## positive dependence
u <- rclayton(2, 1000)             # variates with positive dependence
x <- qexp(u[,1])                   # claim amounts
w <- qexp(u[,2])                   # interarrival times
p <- function(L)
    mean(x)/mean(w) * (1.2 - 1.3) + 1.3 * mean(pmin(L, x))/mean(w)
h <- function(r, L)
    mean(exp(r * (pmin(L, x) - p(L) * w)))
R3 <- adjCoef(h = h, upper = 1, reinsurance = "prop", from = 0, to = 10)
plot(R3)

u <- rclayton(1, 1000)             # independent variates
x <- qexp(u[,1])                   # claim amounts
w <- qexp(u[,2])                   # interarrival times
R4 <- adjCoef(h = h, upper = 1, reinsurance = "prop", from = 0, to = 10)
plot(R4, add = TRUE, col = "green")
legend("bottomright", legend = c("dependence", "independence"),
       col = c("black", "green"), lty = 1)


### CALCULATION OF RUIN PROBABILITIES

## Case with an explicit formula: exponential claims and interarrival
## times. Safety loading is always 20% and premiums are always
## calculated according to the expected value principle.
psi <- ruin(claims = "exponential",
            par.claims = list(rate = 1),
            wait = "exponential",
            par.wait = list(rate = 1),
            premium = 1.2)
psi(0:10)
plot(psi, from = 0, to = 10)

## Exponential claims and hyper-exponential interarrival times.
psi <- ruin(claims = "exponential",
            par.claims = list(rate = 2),
            wait   = "exponential",
            par.wait = list(rate = c(2, 3, 1)/2, w = c(2, 3, 1)/6),
            premium = 1.2)
psi(0:10)

## Hyper-exponential claims and interarrival times.
psi <- ruin(claims = "exponential",
            par.claims = list(rate = c(2, 3, 1)/2, w = c(2, 3, 1)/6),
            wait = "exponential",
            par.wait = list(rate = c(2, 3, 1)/4, w = c(2, 3, 1)/6),
            premium = 0.6)
psi(0:10)

## Exponential claims and Erlang interarrival times
psi <- ruin(claims = "exponential",
            par.claims = list(rate = 2),
            wait = "Erlang",
            par.wait = list(shape = 2, rate = 1),
            premium = 1.2)
psi(0:10)

## Erlang claims and interarrival times
psi <- ruin(claims = "Erlang",
            par.claims = list(shape = 2, rate = 2),
            wait = "Erlang",
            par.wait   = list(shape = 2, rate = 1),
            premium = 0.6)
psi(0:10)

## Mixture of Erlang for claims and Erlang interarrival times
psi <- ruin(claims = "Erlang",
            par.claims = list(shape = c(2, 4), rate = c(1, 3),
                              w = c(1, 2)/3),
            wait = "Erlang",
            par.wait = list(shape = 2, rate = 1),
            premium = 1.2)
psi(0:10)

## Generalized Erlang claims and mixture of two generalized Erlang
## interarrival times. These must be given as phase-type distributions
## to 'ruin'.
prob.c <- c(1, 0, 2, 0)/3
rate.c <- cbind(c(-1, 0, 0, 0), c(1, -3, 0, 0),
                c(0, 0, -2, 0), c(0, 0, 2, -3))
mean.c <- mphtype(1, prob.c, rate.c)
prob.w <- c(1, 0, 0)
rate.w <- cbind(c(-1, 0, 0), c(1, -2, 0), c(0, 2, -3))
mean.w <- mphtype(1, prob.w, rate.w)
psi <- ruin(claims = "phase-type",
            par.claims = list(prob = prob.c, rate = rate.c),
            wait = "phase-type",
            par.wait = list(prob = prob.w, rate = rate.w),
            premium = 1.2 * mean.c/mean.w)
psi(0:10)

par(op)

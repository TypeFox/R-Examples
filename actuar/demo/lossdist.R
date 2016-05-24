### ===== actuar: An R Package for Actuarial Science =====
###
### Demo of the loss distributions facilities provided by actuar
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

require(actuar)
require(graphics)


### A utility function to create graphs for probability laws
showgraphs <- function(fun, par, what = c("d", "p", "m", "lev"), xlim)
{
    dist <- switch(fun,
                   trbeta = "TRANSFORMED BETA DISTRIBUTION",
                   genpareto = "GENERALIZED PARETO DISTRIBUTION",
                   burr = "BURR DISTRIBUTION",
                   invburr = "INVERSE BURR DISTRIBUTION",
                   pareto = "PARETO DISTRIBUTION",
                   invpareto = "INVERSE PARETO DISTRIBUTION",
                   llogis = "LOGLOGISTIC DISTRIBUTION",
                   paralogis = "PARALOGISTIC DISTRIBUTION",
                   invparalogis = "INVERSE PARALOGISTIC DISTRIBUTION",
                   trgamma = "TRANSFORMED GAMMA DISTRIBUTION",
                   invtrgamma = "INVERSE TRANSFORMED GAMMA DISTRIBUTION",
                   invgamma = "INVERSE GAMMA DISTRIBUTION",
                   weibull = "WEIBULL DISTRIBUTION",
                   invweibull = "INVERSE WEIBULL DISTRIBUTION",
                   invexp = "INVERSE EXPONENTIAL DISTRIBUTION",
                   pareto1 = "SINGLE PARAMETER PARETO DISTRIBUTION",
                   lgamma = "LOGGAMMA DISTRIBUTION",
                   genbeta = "GENERALIZED BETA DISTRIBUTION",
                   phtype = "PHASE-TYPE DISTRIBUTION",
                   gamma = "GAMMA DISTRIBUTION",
                   exp = "EXPONENTIAL DISTRIBUTION",
                   chisq = "CHI-SQUARE DISTRIBUTION",
                   lnorm = "LOGNORMAL DISTRIBUTION",
                   invgauss = "INVERSE GAUSSIAN DISTRIBUTION",
                   norm = "NORMAL DISTRIBUTION",
                   beta = "BETA DISTRIBUTION",
                   unif = "UNIFORM DISTRIBUTION")

    if (missing(xlim))
    {
        qf <- match.fun(paste("q", fun, sep = ""))
        formals(qf)[names(par)] <- par
        xlim <- c(0, qf(0.999))
    }
    k <- seq.int(4)
    limit <- seq(0, xlim[2], len = 10)
    mfrow = c(ceiling(length(what) / 2), 2)
    op <- par(mfrow = mfrow, oma = c(0, 0, 2, 0))

    for (t in what)
    {
        f   <- match.fun(paste(t, fun, sep = ""))
        formals(f)[names(par)] <- par

        main <- switch(t,
                       "d" = "Probability Density Function",
                       "p" = "Cumulative Distribution Function",
                       "m" = "Raw Moments",
                       "lev" = "Limited Expected Value Function",
                       "mgf" = "Moment Generating Function")

        if (t == "m")
            plot(k, f(k), type = "l", col = 4, lwd = 2, main = main)
        else if (t == "lev")
            plot(limit, f(limit), type = "l", col = 4, lwd = 2, main = main)
        else if (t == "mgf")
            curve(f(x), xlim = c(0, 2), col = 4, lwd = 2, main = main)
        else
            curve(f(x), xlim = xlim, col = 4, lwd = 2, main = main)
        title(main = dist, outer = TRUE)
    }
    par(op)
}

###
### DATA SETS
###

## The package includes the individual dental claims and grouped
## dental claims data sets often referred to in Klugman, Panjer &
## Willmot (1998, 2004)
data(dental); dental
data(gdental); gdental


###
### PROBABILITY LAWS
###

## Illustration of the new probability laws functions provided by the
## package.

## TRANSFORMED BETA FAMILY

## Transformed beta distribution
showgraphs("trbeta", list(shape1 = 3, shape2 = 4, shape3 = 5, scale = 10))

## Generalized Pareto distribution
showgraphs("genpareto", list(shape1 = 10, shape2 = 4, scale = 10))

## Burr distribution
showgraphs("burr", list(shape1 = 3, shape2 = 4, scale = 10))

## Inverse Burr distribution
showgraphs("invburr", list(shape1 = 3, shape2 = 6, scale = 10))

## Pareto distribution
showgraphs("pareto", list(shape = 10, scale = 10))

## Inverse Pareto distribution
showgraphs("invpareto", list(shape = 4, scale = 1), what = c("d", "p"))

## Loglogistic distribution
showgraphs("llogis", list(shape = 6, scale = 10))

## Paralogistic distribution
showgraphs("paralogis", list(shape = 3, scale = 10))

## Inverse paralogistic distribution
showgraphs("invparalogis", list(shape = 6, scale = 10))

## TRANSFORMED GAMMA FAMILY

## Transformed gamma distribution
showgraphs("trgamma", list(shape1 = 3, shape2 = 1, scale = 10))

## Inverse transformed gamma distribution
showgraphs("invtrgamma", list(shape1 = 3, shape2 = 2, scale = 10))

## Inverse gamma distribution
showgraphs("invgamma", list(shape = 6, scale = 10))

## Weibull distribution ('mweibull' and 'levweibull')
showgraphs("weibull", list(shape = 1.5, scale = 10))

## Inverse Weibull distribution
showgraphs("invweibull", list(shape = 6, scale = 10))

## Inverse exponential distribution
showgraphs("invexp", list(rate = 1), what = c("d", "p"))

## OTHER DISTRIBUTIONS

## Single parameter Pareto distribution
showgraphs("pareto1", list(shape = 5, min = 10), xlim = c(0, 50))

## Loggamma distribution
showgraphs("lgamma", list(shapelog = 2, ratelog = 5))

## Generalized beta distribution
showgraphs("genbeta", list(shape1 = 1, shape2 = 2, shape3 = 3, scale = 2))

## Phase-type distribution
showgraphs("phtype", list(prob = c(0.5614, 0.4386), rates = matrix(c(-8.64, 0.101, 1.997, -1.095), 2, 2)), what = c("d", "p", "m", "mgf"), xlim = c(0.001, 5))

## DISTRIBUTIONS ALREADY IN R

## Gamma distribution
showgraphs("gamma", list(shape = 3, rate = 5), what = c("m", "lev", "mgf"))

## Chi-square distribution
showgraphs("chisq", list(df = 3), what = c("m", "lev", "mgf"))

## Exponential distribution
showgraphs("exp", list(rate = 5), what = c("m", "lev", "mgf"))

## Lognormal distribution
showgraphs("lnorm", list(meanlog = 1, sdlog = 1), what = c("m", "lev"))

## Inverse gaussian distribution (from package SuppDists)
showgraphs("invgauss", list(nu = 1, lambda = 10), what = c("m", "lev", "mgf"), xlim = c(0, 10))

## Normal distribution
showgraphs("norm", list(mean = 0, sd = 1), what = c("m", "mgf"))

## Beta distribution
showgraphs("beta", list(shape1 = 1, shape2 = 2), what = c("m", "lev"))

## Uniform distribution
showgraphs("unif", list(min = 0, max = 1), what = c("m", "lev", "mgf"))


###
### GROUPED DATA MANIPULATION
###

## Creation of grouped data objects
x <- grouped.data(groups = c(0, 25, 50, 100, 150, 250, 500),
                  line1 = c(30, 31, 57, 42, 65, 84),
                  line2 = c(26, 33, 31, 19, 16, 11))
x

## Extraction and replacement: only "[" and "[<-" are officially
## supported.
x[, 1]                                  # group boundaries
x[1]                                    # notice the difference
x[, -1]                                 # group frequencies
x[1:3,]                                 # first 3 groups
x[1, 2] <- 22; x                        # frequency replacement
x[1, 1] <- c(0, 20); x                  # boundary replacement

## A method for 'mean' exists for grouped data objects.
mean(x)

## Function 'hist' handles individual data only. We provide a method
## for grouped data.  Only the first frequencies column is considered.
hist(x[, -3])

## Function 'ogive' returns a function to compute the ogive of grouped
## data in any point, much like 'ecdf' does for individual
## data. Methods also exist to extract the group boundaries ('knots')
## and plot the ogive.
Fnt <- ogive(x)
summary(Fnt)
knots(Fnt)                              # group boundaries
Fnt(knots(Fnt))                         # ogive at group boundaries
plot(Fnt)                               # plot of the ogive

## A method of 'quantile' for grouped data objects computes linearly
## smoothed quantiles from a sample, that is the inverse of the ogive
## in various points.
quantile(x)
Fnt(ogive(x))


###
### EMPIRICAL MOMENTS CALCULATION
###

## Function 'emm' computes the k-th empirical moment of a sample,
## whether it is individual or grouped data.
emm(dental)                             # == mean(dental)
emm(gdental)                            # == mean(gdental)
emm(dental, order = 1:3)                # first three moments
emm(gdental, order = 1:3)               # idem

## Function 'elev' is similar to 'ecdf' and 'ogive' in that it returns
## a function to compute the empirical limited expected value (first
## limited moment) for any limit. There are methods for individual and
## grouped data.
lev <- elev(dental)
lev(knots(lev))                         # ELEV at data points
plot(lev, type = "o", pch = 19)         # plot of the ELEV function

lev <- elev(gdental)
lev(knots(lev))                         # ELEV at data points
plot(lev, type = "o", pch = 19)         # plot of the ELEV function


###
### MINIMUM DISTANCE ESTIMATION
###

## Maximum likelihood estimation (for individual data) is well covered
## by 'fitdistr' in package MASS. We provide function 'mde' to fit
## models using three distance minimization techniques: Cramer-von
## Mises (for individual and grouped data), chi-square and layer
## average severity (both grouped data only). Usage (and inner
## working) is very similar to 'fitdistr'.
mde(dental, pexp, start = list(rate = 1/200), measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200), measure = "CvM")
mde(gdental, pexp, start = list(rate = 1/200), measure = "chi-square")
mde(gdental, levexp, start = list(rate = 1/200), measure = "LAS")


###
### COVERAGE MODIFICATIONS
###

## Function 'coverage' is useful to obtain the probability density
## function (pdf) or cumulative distribution function (cdf) of a loss
## random variable under coverage modifications.
f <- coverage(dgamma, pgamma, deductible = 1, limit = 7)
curve(dgamma(x, 3), xlim = c(0, 10), ylim = c(0, 0.3))    # original
curve(f(x, 3), xlim = c(0.01, 5.99), col = 4, add = TRUE) # modified

x <- rgamma(1000, 3, 1)                 # sample of claim amounts
x <- pmin(x, 7)[x > 1] - 1              # deductible and limit

library(MASS)                           # for ML estimation
m <- mean(x)                            # empirical mean
v <- var(x)                             # empirical variance
(p <- fitdistr(x, f, start = list(shape = m^2/v, rate = m/v))$estimate ) # MLE
hist(x + 1, breaks = 0:10, prob = TRUE)  # histogram of observed data
curve(dgamma(x, p[1], p[2]), add = TRUE) # fit of underlying distribution

par(op)

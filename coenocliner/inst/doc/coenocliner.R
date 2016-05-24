## ----prelim, echo = FALSE, results = "hide"------------------------------
library("knitr")
opts_chunk$set(fig.lp = "fig:", size = "small", out.width = ".7\\linewidth",
               fig.align = "center", fig.show = "hold")

## ----load-library--------------------------------------------------------
library("coenocliner")

## ----show-matrix-grad-locs-----------------------------------------------
xy <- cbind(x = seq(from = 4, to = 7, length.out = 100),
            y = seq(from = 1, to = 100, length.out = 100))

## ----show-list-grad-locs-------------------------------------------------
xy <- list(x = seq(from = 4, to = 6, length.out = 100),
           y = seq(from = 1, to = 100, length.out = 100))

## ----show-gaussian-params------------------------------------------------
showParams("gaussian")

## ----param-example-------------------------------------------------------
opt <- c(4,5,6)
tol <- rep(0.25, 3)
h <- c(10,20,30)
parm <- cbind(opt = opt, tol = tol, h = h)     # matrix form
parl <- list(opt = opt, tol = tol, h = h)      # list form

## ----param-example-2-----------------------------------------------------
opty <- c(25, 50, 75)
tol <- c(5, 10, 20)
pars <- list(px = parm,
             py = cbind(opt = opty, tol = tol))

## ----count-models, echo = FALSE------------------------------------------
c("poisson", "negbin", "bernoulli", "binary", "binomial", "betabinomial", "ZIP", "ZINB", "ZIB", "ZIBB")

## ----example1-params-----------------------------------------------------
set.seed(2)
M <- 20                                    # number of species
ming <- 3.5                                # gradient minimum...
maxg <- 7                                  # ...and maximum
locs <- seq(ming, maxg, length = 100)      # gradient locations
opt  <- runif(M, min = ming, max = maxg)   # species optima
tol  <- rep(0.25, M)                       # species tolerances
h    <- ceiling(rlnorm(M, meanlog = 3))    # max abundances
pars <- cbind(opt = opt, tol = tol, h = h) # put in a matrix

## ----example1-expectations-----------------------------------------------
mu <- coenocline(locs, responseModel = "gaussian", params = pars,
                 expectation = TRUE)

## ----example1-head-expectations------------------------------------------
class(mu)
dim(mu)
mu

## ----example1-plot-expectations, fig.cap = "Gaussian species response curves along a hypothetical pH gradient", fig.height = 5----
plot(mu, lty = "solid", type = "l", xlab = "pH", ylab = "Abundance")

## ----example1-sim--------------------------------------------------------
simp <- coenocline(locs, responseModel = "gaussian", params = pars,
                   countModel = "poisson")

## ----example1-plot-simulations, fig.cap = "Simulated species abundances with Poisson errors from Gaussian response curves along a hypothetical pH gradient", fig.height = 5----
plot(simp, lty = "solid", type = "p", pch = 1:10, cex = 0.8,
     xlab = "pH", ylab = "Abundance")

## ----example1-nb-sim-----------------------------------------------------
simnb <- coenocline(locs, responseModel = "gaussian", params = pars,
                    countModel = "negbin", countParams = list(alpha = 0.5))

## ----example1-plot-nb-simulations, fig.cap = "Simulated species abundance with negative binomial errors from Gaussian response curves along a hypothetical pH gradient", fig.height = 5----
plot(simnb, lty = "solid", type = "p", pch = 1:10, cex = 0.8,
     xlab = "pH", ylab = "Abundance")

## ----example2-beta-pars--------------------------------------------------
A0    <- c(5,4,7,5,9,8) * 10               # max abundance
m     <- c(25,85,10,60,45,60)              # location on gradient of modal abundance
r     <- c(3,3,4,4,6,5) * 10               # species range of occurence on gradient
alpha <- c(0.1,1,2,4,1.5,1)                # shape parameter
gamma <- c(0.1,1,2,4,0.5,4)                # shape parameter
locs  <- 1:100                             # gradient locations
pars  <- list(m = m, r = r, alpha = alpha,
              gamma = gamma, A0 = A0)      # species parameters, in list form

## ----example2-beta-expectations------------------------------------------
mu <- coenocline(locs, responseModel = "beta", params = pars, expectation = TRUE)

## ----example2-head-------------------------------------------------------
mu

## ----example2-beta-plot-expectations, fig.cap = "Generalised beta function species response curves along a hypothetical environmental gradient recreating Figure 2 in Minchin (1987).", fig.height = 5----
plot(mu, lty = "solid", type = "l", xlab = "Gradient", ylab = "Abundance")

## ----example3-params-----------------------------------------------------
set.seed(10)
N <- 30                                           # number of samples
M <- 20                                           # number of species
## First gradient
ming1 <- 3.5                                      # 1st gradient minimum...
maxg1 <- 7                                        # ...and maximum
loc1 <- seq(ming1, maxg1, length = N)             # 1st gradient locations
opt1 <- runif(M, min = ming1, max = maxg1)        # species optima
tol1 <- rep(0.5, M)                               # species tolerances
h    <- ceiling(rlnorm(M, meanlog = 3))           # max abundances
par1 <- cbind(opt = opt1, tol = tol1, h = h)      # put in a matrix
## Second gradient
ming2 <- 1                                        # 2nd gradient minimum...
maxg2 <- 100                                      # ...and maximum
loc2 <- seq(ming2, maxg2, length = N)             # 2nd gradient locations
opt2 <- runif(M, min = ming2, max = maxg2)        # species optima
tol2 <- ceiling(runif(M, min = 5, max = 50))      # species tolerances
par2 <- cbind(opt = opt2, tol = tol2)             # put in a matrix
## Last steps...
pars <- list(px = par1, py = par2)                # put parameters into a list
locs <- expand.grid(x = loc1, y = loc2)           # put gradient locations together

## ----example3-expectations-----------------------------------------------
mu2d <- coenocline(locs, responseModel = "gaussian",
                   params = pars, extraParams = list(corr = 0.5),
                   expectation = TRUE)

## ----example3-head-locs--------------------------------------------------
head(locs)

## ----example3-persp-plots, fig.cap = "Bivariate Gaussian species responses for four selected species."----
layout(matrix(1:4, ncol = 2))
op <- par(mar = rep(1, 4))
persp(mu2d, species = c(2, 8, 13, 19), ticktype = "detailed",
      zlab = "Abundance")
par(op)
layout(1)

## ----example3-nb-sim-----------------------------------------------------
sim2d <- coenocline(locs, responseModel = "gaussian",
                    params = pars, extraParams = list(corr = 0.5),
                    countModel = "negbin", countParams = list(alpha = 1))

## ----example3-persp-plots2, fig.cap = "Simulated counts using negative binomial errors from bivariate Gaussian species responses for four selected species."----
layout(matrix(1:4, ncol = 2))
op <- par(mar = rep(1, 4))
persp(sim2d, species = c(2, 8, 13, 19), ticktype = "detailed",
      zlab = "Abundance")
par(op)
layout(1)

## ----sessionInfo, results = "asis", echo = FALSE-------------------------
toLatex(sessionInfo())


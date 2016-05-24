## ----opts, echo = FALSE, message = FALSE---------------------------------
library("knitr")
#knitr::opts_chunk$set(eval = FALSE)

## ----eval=TRUE, echo=FALSE, results="hide"-------------------------------
suppressMessages(require("growthrates"))
#require("growthrates")

## ----eval=FALSE----------------------------------------------------------
#  library("growthrates")

## ------------------------------------------------------------------------
grow_logistic_yshift <- function(time, parms) {
  with(as.list(parms), {
    y <- (K * y0) / (y0 + (K - y0) * exp(-mumax * time)) + y_shift
    as.matrix(data.frame(time = time, y = y, log_y = log(y)))
  })
}

## ---- fig.width=5, fig.height=4------------------------------------------
time <- 1:10
out <- grow_logistic_yshift(time, parms = list(y0 = 1, mumax = 0.5, K = 10, y_shift = 2))
plot(time, out[, "y"], type = "b")

## ---- fig.width=5, fig.height=4------------------------------------------
x <- seq(5, 100, 5)
y <- c(2.1, 2.3, 5, 4.7, 4.3, 6.9, 8.2, 11.5, 8.8, 10.2, 14.5, 12.5,
       13.6, 12.7, 14.2, 12.5, 13.8, 15.1, 12.7, 14.9)

fit <- fit_growthmodel(grow_logistic_yshift,
                       p = c(y0 = 1, mumax = 0.1, K = 10, K = 10, y_shift = 1),
                       time = x, y = y)
plot(fit)
summary(fit)

## ------------------------------------------------------------------------
ode_K_linear <- function (time, init, parms, ...) {
  with(as.list(c(parms, init)), {
    dy <- mumax * y * (1 - y/K)
    dK <- dK
    list(c(dy, dK), log_y = unname(log(y)))
  })
}

grow_K_linear <- function(time, parms, ...) {
  init    <- parms[c("y0", "K")]           # initial values
  names(init) <- c("y", "K")               # force names of state variables
  odeparms <- parms[c("mumax", "dK")]      # the parms of the ODE model
  out <- ode(init, time, ode_K_linear, parms = odeparms)
  out
}

## ------------------------------------------------------------------------
grow_K_linear <- growthmodel(grow_K_linear,
                             pnames = c("y0", "K", "mumax", "deltaK"))
head(grow_K_linear(time = 1:10, c(y0 = .1, K = 1, mumax = 0.1, dK = 0.5)))

## ---- fig.width=4, fig.height=3------------------------------------------
x <- seq(5, 100, 5)
y <- c(0.1, 2.2, 3.1, 1.5, 8.9, 8, 8.4, 9.8, 9.3, 10.6, 12, 13.6,
  13.1, 13.3, 11.6, 14.7, 12.6, 13.9, 16.9, 14.4)
fit <- fit_growthmodel(grow_K_linear,
                       p = c(y0 = 0.1, mumax = 0.2, K = 10, dK = .1), time = x, y = y)
plot(fit)
summary(fit)

## ---- eval=FALSE---------------------------------------------------------
#  ## The following example shows how to use compiled growth models
#  ## from inline code, by using the 'cOde' package of Daniel Kaschek
#  ## Note: This example needs the R development tools.
#  ##  - suitable compilers on Linux and Mac
#  ##  - Rtools on Windows from
#  
#  library("growthrates")
#  library("cOde")
#  
#  ## define a system of ODEs and compile it --------------------------------------
#  ode_K_linear <- funC(c(
#    y = "mumax * y * (1-y/K)",
#    K = "dK"
#  ))
#  
#  yini <- c(y = 1, K = 10)
#  parms = c(mumax = 0.1, dK = 0.05)
#  
#  ## run the model
#  out1 <- odeC(yini, times = 0:100, ode_K_linear, parms = parms)
#  
#  ## generate artificial test data with normal distributed noise
#  x <- seq(5, 100, 5)
#  y <- odeC(yini, x, ode_K_linear, parms)[, "y"] + rnorm(x)
#  
#  
#  ## create a "growthmodel" with interfaces compatible to package growthrates
#  ## It is essential to use consistent names for parameters and initial values!
#  
#  grow_K_linear <- function(time, parms, ...) {
#    init    <- parms[c("y0", "K")]  # initial values
#    names(init) <- c("y", "K")      # force names
#    out <- odeC(init, time, ode_K_linear, parms)
#    cbind(out, log_y = log(out[, "y"]))
#  }
#  
#  ## convert this to an object, (maybe needed by future extensions)
#  grow_K_linear <- growthmodel(grow_K_linear, pnames = c("y0", "mumax", "K", "dK"))
#  
#  ## Test the growthmodel
#  ## Columns with names 'time', 'y' and 'log_y' are mandatory.
#  head(grow_K_linear(time = x, c(y0 = 1, mumax = 0.1, K = 10, dK = 0.1)))
#  
#  
#  ## Fit the model ---------------------------------------------------------------
#  fit <- fit_growthmodel(grow_K_linear,
#                         p = c(y0 = 1, mumax = 0.1, K = 10, dK = 0.1), time = x, y = y)
#  plot(fit)
#  summary(fit)
#  
#  ## Unload DLL and cleanup ------------------------------------------------------
#  ## DLL creation should ideally be directed to a temporary directory.
#  dll <- paste(ode_K_linear, .Platform$dynlib.ext, sep = "")
#  dyn.unload(dll)
#  unlink(dll)
#  unlink(paste(ode_K_linear, ".c", sep = ""))
#  unlink(paste(ode_K_linear, ".o", sep = ""))


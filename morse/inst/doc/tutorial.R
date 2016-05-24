## ----include=FALSE, echo=FALSE-------------------------------------------

knitr::opts_chunk$set(fig.width = 7, 
                      fig.height = 4, 
                      cache = TRUE)


## ---- echo=FALSE, cache=TRUE, results='hide'-----------------------------
library(morse)
library(ggplot2)

## ---- cache=TRUE---------------------------------------------------------
data(cadmium2)
survDataCheck(cadmium2)

## ---- cache=TRUE---------------------------------------------------------
dat <- survData(cadmium2)
head(dat)

## ---- cache=TRUE---------------------------------------------------------
plot(dat, style = "ggplot", pool.replicate = FALSE)

## ---- cache=TRUE---------------------------------------------------------
plot(dat, concentration = 124, addlegend = TRUE,
     pool.replicate = FALSE, style = "ggplot")

## ---- cache=TRUE---------------------------------------------------------
plot(dat, target.time = 21, pool.replicate = FALSE, style = "ggplot",
     addlegend = TRUE)

## ---- cache=TRUE---------------------------------------------------------
plot(dat, concentration = 232, target.time = 21)

## ---- cache=TRUE---------------------------------------------------------
summary(dat)

## ---- results="hide", cache=TRUE-----------------------------------------
fit <- survFitTT(dat,
                 target.time = 21,
                 lcx = c(10, 20, 30, 40, 50))

## ---- cache=TRUE---------------------------------------------------------
summary(fit)

## ---- cache=TRUE---------------------------------------------------------
plot(fit, log.scale = TRUE, ci = TRUE, style = "ggplot",
     addlegend = TRUE)

## ---- results="hide", cache=TRUE-----------------------------------------
data("cadmium1")
wrong_fit <- survFitTT(survData(cadmium1),
                       target.time = 21,
                       lcx = c(10, 20, 30, 40, 50))
plot(wrong_fit, log.scale = TRUE, ci = TRUE, style = "ggplot",
     addlegend = TRUE)

## ----cache=TRUE, results="hide"------------------------------------------
ppc(fit, style = "ggplot")

## ---- cache=TRUE---------------------------------------------------------
# (1) load dataset
data(cadmium2)

# (2) check structure and integrity of the dataset
reproDataCheck(cadmium2)

# (3) create a `reproData` object
dat <- reproData(cadmium2)

# (4) represent the cumulated number of offspring as a function time
plot(dat, style = "ggplot", pool.replicate = FALSE)
plot(dat, target.time = 21, addlegend = TRUE, style = "ggplot",
     pool.replicate = FALSE)
plot(dat, concentration = 124, addlegend = TRUE, style = "ggplot",
     pool.replicate = FALSE)
plot(dat, concentration = 124, target.time = 21, style = "ggplot")

# (5) check information on the experimental design
summary(dat)

# (6) fit an exposure-response model at target-time
fit <- reproFitTT(dat, stoc.part = "bestfit",
                  target.time = 21,
                  ecx = c(10, 20, 30, 40, 50),
                  quiet = TRUE)
summary(fit)

## ---- cache=TRUE---------------------------------------------------------
plot(fit, log.scale = TRUE, ci = TRUE,
     style = "ggplot")

## ---- cache=TRUE---------------------------------------------------------
ppc(fit, style = "ggplot")

## ---- cache=TRUE---------------------------------------------------------
summary(fit)

## ---- cache=TRUE---------------------------------------------------------
dat <- reproData(cadmium2)
plot(as.survData(dat))

## ---- cache=TRUE, eval=FALSE---------------------------------------------
#  ppc(out, style = "ggplot")


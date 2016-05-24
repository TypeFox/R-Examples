## ----setup, include=FALSE, cache=FALSE-----------------------------------
library(knitr)
my_pdf = function(file, width, height) {
  pdf(file, width = width, height = height, pointsize = 8)
}

## ------------------------------------------------------------------------
library(MortalitySmooth)

## ----eval=FALSE----------------------------------------------------------
#  citation("MortalitySmooth")

## ------------------------------------------------------------------------
years <- 1950:2006
death <- selectHMDdata("Japan", "Deaths", "Females",
                       ages = 80, years = years)
exposure <- selectHMDdata("Japan", "Exposures", "Females",
                          ages = 80, years = years)
fitBIC <- Mort1Dsmooth(x=years, y=death,
                       offset=log(exposure))

## ------------------------------------------------------------------------
fitBIC$psi2 

## ------------------------------------------------------------------------
fitBICover <- Mort1Dsmooth(x=years, y=death,
                           offset=log(exposure),
                           overdispersion=TRUE)

## ------------------------------------------------------------------------
fitBIC$lambda;fitBICover$lambda

## ----fig.keep='none'-----------------------------------------------------
plot(fitBICover)
lines(years, fitBIC$logmortality, col=4, lwd=2, lty=2) 

## ----echo=FALSE----------------------------------------------------------
plot(fitBICover)
lines(years, fitBIC$logmortality, col=4, lwd=2, lty=2) 

## ------------------------------------------------------------------------
exposure.int <- exposure
exposure.int[20:40] <- NA

## ------------------------------------------------------------------------
w.int <- rep(1, length(years))
w.int[20:40] <- 0
fit1Dint <- Mort1Dsmooth(x=years, y=death,
                         offset=log(exposure.int), w=w.int)

## ----fig.keep='none'-----------------------------------------------------
plot(fit1Dint)

## ----echo=FALSE----------------------------------------------------------
plot(fit1Dint)

## ------------------------------------------------------------------------
years.ext <- seq(years[1], years[length(years)]+20)
exposure.ext <- c(exposure, rep(NA, 20))
death.ext <- c(death, rep(NA, 20))

## ------------------------------------------------------------------------
w.ext <- c(rep(1, length(years)), rep(0, 20))

## ----fig.keep='none'-----------------------------------------------------
fit1Dext <- Mort1Dsmooth(x=years.ext, y=death.ext,
                         offset=log(exposure.ext), w=w.ext,
                         method=3, lambda=fitBICover$lambda)
plot(fit1Dext)

## ----echo=FALSE----------------------------------------------------------
plot(fit1Dext)

## ------------------------------------------------------------------------
ages <- 50:100
years <- 1950:2006
death <- selectHMDdata("Sweden", "Deaths", "Females",
                       ages = ages, years = years) 
exposure <- selectHMDdata("Sweden", "Exposures", "Females",
                          ages = ages, years = years)
fitBIC <- Mort2Dsmooth(x=ages, y=years, Z=death,
                       offset=log(exposure)) 

## ------------------------------------------------------------------------
fitBIC$psi2

## ------------------------------------------------------------------------
fitBICover <- Mort2Dsmooth(x=ages, y=years, Z=death,
                           offset=log(exposure),
                           overdispersion=TRUE)

## ------------------------------------------------------------------------
fitBIC$lambdas
fitBICover$lambdas

## ----fig.keep='none'-----------------------------------------------------
new.exposure <- exposure
new.exposure[20:30, ] <- 0
W <- matrix(1, length(ages), length(years))
W[20:30, ] <- 0
fit2Dint <- Mort2Dsmooth(x=ages, y=years, Z=death,
                         offset=log(new.exposure), W=W)
plot(fit2Dint)

## ----echo=FALSE----------------------------------------------------------
plot(fit2Dint)

## ----fig.keep='none'-----------------------------------------------------
W <- matrix(1, length(ages), length(years))
set.seed(3)
zeros <- sample(x=1:prod(dim(W)),
                size=round(prod(dim(W))/100*90),
                replace=FALSE)
W[zeros] <- 0
fit2Dint10 <- Mort2Dsmooth(x=ages, y=years, Z=death,
                           offset=log(exposure), W=W)
plot(fit2Dint10)

## ----echo=FALSE----------------------------------------------------------
plot(fit2Dint10)

## ----fig.keep='none'-----------------------------------------------------
years.new <- 1950:2025
M <- matrix(0, nrow=length(ages),
            ncol=length(years.new)-length(years))
death.new <- cbind(death, M)
exposure.new <- cbind(exposure, M)
W <- matrix(1, length(ages), length(years))
W.new <- cbind(W, M)
fit2Dext <- Mort2Dsmooth(x=ages, y=years.new, Z=death.new,
                         offset=log(exposure.new), W=W.new)
plot(fit2Dext)

## ----echo=FALSE----------------------------------------------------------
plot(fit2Dext)


## ----setup, include=FALSE------------------------------------------------
library(knitr)
library(tabplot)

## ------------------------------------------------------------------------
require(ggplot2)
data(diamonds)
## add some NA's
is.na(diamonds$price) <- diamonds$cut == "Ideal"
is.na(diamonds$cut) <- (runif(nrow(diamonds)) > 0.8)

## ----message=FALSE-------------------------------------------------------
n <- nrow(diamonds)
N <- 200L * n

## convert to ff format (not enough memory otherwise)
require(ffbase)
diamondsff <- as.ffdf(diamonds)
nrow(diamondsff) <- N

# fill with identical data
for (i in chunk(from=1, to=N, by=n)){
  diamondsff[i,] <- diamonds
}

## ----message=FALSE-------------------------------------------------------
system.time(
	p <- tablePrepare(diamondsff)
)

## ----message=FALSE-------------------------------------------------------
system.time(
	tab <- tableplot(p, plot=FALSE)
)

## ----message=FALSE-------------------------------------------------------
system.time(
	tab <- tableplot(p, sample=TRUE, sampleBinSize=1e2, plot=FALSE)
)

## ----message=FALSE-------------------------------------------------------
system.time(
	tab <- tableplot(p, sample=TRUE, sampleBinSize=1e3, plot=FALSE)
)

## ----message=FALSE-------------------------------------------------------
system.time(
	tab <- tableplot(p, sample=TRUE, sampleBinSize=1e4, plot=FALSE)
)


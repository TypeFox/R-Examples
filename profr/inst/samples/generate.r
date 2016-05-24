# Code to generate sample Rprof files

a <- function() rnorm(1e7)
b <- function() rnorm(1e7)
c <- function() {
  rnorm(1e7)
  a()
}
d <- function() {
  a()
  b()
  c()
  a()
  return()
}

Rprof("nesting.rprof")
d()
Rprof(NULL)


data(diamonds, package="ggplot2")
library(reshape)

Rprof("reshape.rprof")
dm <- melt(diamonds)
cast(dm, cut + color + clarity ~ variable, mean)
Rprof(NULL)

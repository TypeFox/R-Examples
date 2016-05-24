## ----, echo = FALSE------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
library(lazyeval)
f1 <- function(x) lazy(x)
g1 <- function(y) f1(y)

g1(a + b)

## ------------------------------------------------------------------------
f2 <- function(x) lazy(x, .follow_symbols = FALSE)
g2 <- function(y) f2(y)

g2(a + b)

## ------------------------------------------------------------------------
a <- 10
b <- 1

lazy_eval(g1(a + b))
lazy_eval(g2(a + b))


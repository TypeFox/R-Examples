## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(fig.width=7, fig.height=5, comment="")
library(BrailleR)

## ----hist, fig.cap="A histogram of 1000 random values from a normal distribution"----
x=rnorm(1000)
VI(hist(x))

## ----BrailleRHistExample-------------------------------------------------
example(hist)


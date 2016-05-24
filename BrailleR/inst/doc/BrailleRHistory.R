## ----setup, results="hide"-----------------------------------------------
library(BrailleR)

## ----hist, fig.cap="A histogram of 1000 random values from a normal distribution"----
x=rnorm(1000)
VI(hist(x))


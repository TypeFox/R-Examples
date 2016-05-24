## ----include=FALSE-------------------------------------------------------
knitr::opts_chunk$set(highlight=FALSE)  # remove candy shop colors
# options(prompt="R> ", continue="+  ", width=75, useFancyQuotes=FALSE)
# opts_chunk$set(fig.path="figures/figure-", fig.align="center")
# render_sweave()           # use Sweave environments
# set_header(highlight="")  # do not use the Sweave.sty package

## ----message=FALSE-------------------------------------------------------
library("ccaPP")
data("diabetes")
x <- diabetes$x
y <- diabetes$y

## ------------------------------------------------------------------------
spearman <- maxCorGrid(x, y)
spearman

## ------------------------------------------------------------------------
spearman$a
spearman$b

## ------------------------------------------------------------------------
maxCorGrid(x, y, method = "kendall")
maxCorGrid(x, y, method = "M")
maxCorGrid(x, y, method = "pearson")

## ------------------------------------------------------------------------
maxCorGrid(x, y, consistent = TRUE)
maxCorGrid(x, y, method = "kendall", consistent = TRUE)

## ------------------------------------------------------------------------
permTest(x, y, R = 100, nCores = 2, seed = 2016)

## ------------------------------------------------------------------------
permTest(x, y, R = 100, method = "kendall", nCores = 2, seed = 2016)
permTest(x, y, R = 100, method = "M", nCores = 2, seed = 2016)
permTest(x, y, R = 100, method = "pearson", nCores = 2, seed = 2016)

## ------------------------------------------------------------------------
y[1, "GlucoseIntolerance"] <- 8.1

## ------------------------------------------------------------------------
permTest(x, y, R = 100, nCores = 2, seed = 2016)
permTest(x, y, R = 100, method = "kendall", nCores = 2, seed = 2016)
permTest(x, y, R = 100, method = "M", nCores = 2, seed = 2016)
permTest(x, y, R = 100, method = "pearson", nCores = 2, seed = 2016)


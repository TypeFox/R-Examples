## ----example, fig.align='center', fig.width=5, fig.height=4--------------
x <- c(rnorm(200, 4, 1), rnorm(200, 5, 2), rnorm(200, 6, 1.5))
groups <- c(rep("Cond1", 200), rep("Cond2", 200), rep("Cond3", 200))

library(sinaplot)
sinaplot(x, groups)

## ----blood, echo=FALSE, results='asis'-----------------------------------
knitr::kable(head(blood, 10))

## ----bloodDensity, fig.align='center', fig.height=6, fig.width=7---------
sinaplot(blood$value, blood$type)

## ----bloodNeighb, fig.align='center', fig.height=6, fig.width=7----------
sinaplot(blood$value, blood$type, method = "neighbourhood")

## ----bloodScaleOff, fig.align='center', fig.height=6, fig.width=7--------
sinaplot(blood$value, blood$type, method = "neighbourhood", 
         groupwiseScale = FALSE)

## ----sessionInfo---------------------------------------------------------
sessionInfo()


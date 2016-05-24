## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center")

## ----message=FALSE-------------------------------------------------------
library(nullabor)
library(ggplot2)
library(dplyr)

## ------------------------------------------------------------------------
d <- lineup(null_permute("mpg"), mtcars)
head(d)
# Position of actual data plot
attr(d, "pos")

## ----, fig.height = 10, fig.width = 11-----------------------------------
qplot(mpg, wt, data = d) + facet_wrap(~ .sample)

## ----, fig.height = 10, fig.width = 11-----------------------------------
d <- rorschach(null_permute("mpg"), mtcars, n = 20, p = 0)
qplot(mpg, wt, data = d) + facet_wrap(~ .n)

## ------------------------------------------------------------------------
head(null_dist("mpg", dist = "normal")(mtcars))

## ------------------------------------------------------------------------
head(null_permute("mpg")(mtcars))

## ------------------------------------------------------------------------
head(null_lm(wt~mpg, method = 'rotate')(mtcars))

## ------------------------------------------------------------------------
uni_dist(rnorm(100), rpois(100, 2))

## ------------------------------------------------------------------------
with(mtcars, reg_dist(data.frame(wt, mpg), data.frame(sample(wt), mpg)))

## ------------------------------------------------------------------------
with(mtcars, box_dist(data.frame(as.factor(am), mpg),  data.frame(as.factor(sample(am)), mpg)))

## ------------------------------------------------------------------------
with(mtcars, sep_dist(data.frame(wt, mpg,  as.numeric(as.factor(mtcars$cyl))), data.frame(sample(wt), mpg,  as.numeric(as.factor(mtcars$cyl))), nclust = 3))

## ------------------------------------------------------------------------
with(mtcars, bin_dist(data.frame(wt, mpg), data.frame(sample(wt), mpg), lineup.dat = NULL, X.bin = 5, Y.bin = 5))

## ------------------------------------------------------------------------
calc_mean_dist(lineup(null_permute('mpg'), mtcars, pos = 10), var = c('mpg', 'wt'), met = 'reg_dist', pos = 10)

## ------------------------------------------------------------------------
calc_diff(lineup(null_permute('mpg'), mtcars, pos = 10), var = c('mpg', 'wt'), met = 'reg_dist', dist.arg = NULL, pos = 10)

## ----, fig.height = 5, fig.width = 5.5-----------------------------------
opt.diff <- opt_bin_diff(lineup(null_permute('mpg'), mtcars, pos = 10), var = c('mpg', 'wt'), 2, 4, 2, 4, pos = 10, plot = TRUE)
opt.diff$p

## ----, fig.height = 10, fig.width = 11-----------------------------------
lineup.dat <- lineup(null_permute('mpg'), mtcars, pos = 10)
qplot(mpg, wt, data = lineup.dat, geom = 'point') + facet_wrap(~ .sample)

## ------------------------------------------------------------------------
#decrypt('...') 
#[1] 'True data in position 10' # Use pos = 10

## ----, message = FALSE---------------------------------------------------
dist.vals <- distmet(lineup.dat, var = c('mpg', 'wt'),'reg_dist', null_permute('mpg'), pos = 10, repl = 100, dist.arg = NULL) 

## ------------------------------------------------------------------------
head(dist.vals$lineup)
dist.vals$diff
head(dist.vals$closest)
head(dist.vals$null_values)
dist.vals$pos

## ----, message = FALSE---------------------------------------------------
dist.vals <- distmet(lineup.dat, var = c('mpg', 'wt'),'bin_dist', null_permute('mpg'), pos = 10, repl = 100, dist.arg = list(lineup.dat = lineup.dat, X.bin = 5, Y.bin = 5)) 

## ----, fig.height = 5, fig.width = 5.5-----------------------------------
distplot(dist.vals)


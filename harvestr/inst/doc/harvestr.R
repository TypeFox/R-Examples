### R code from vignette source 'harvestr.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library(harvestr)
library(plyr)
library(MCMCpack)
library(dostats)
print.data.frame <- function(x, ...){
  print(ascii(df, include.rownames = FALSE), type = 'rest')
}


###################################################
### code chunk number 2: ex1-setup
###################################################
library(harvestr)
library(plyr)
seeds <- gather(100, seed=12345)


###################################################
### code chunk number 3: ex1-generate_data
###################################################
datasets <- farm(seeds, {
  x <- rnorm(100)
  y <- rnorm(100, mean=x)
  data.frame(y,x)
})


###################################################
### code chunk number 4: r ex1-analysis
###################################################
analyses <-  harvest(datasets, lm)


###################################################
### code chunk number 5: ex1-summarize
###################################################
library(dostats)
coefs <- t(sapply(analyses, coef))
adply(coefs,2, dostats, mean, sd)


###################################################
### code chunk number 6: stochastic
###################################################
library(MCMCpack)
library(plyr)
posteriors <- harvest(datasets, MCMCregress, formula=y~x)
dataframes <- harvest(posteriors, as.data.frame)
X.samples  <- harvest(dataframes, `[[`, "x")
densities  <- harvest(X.samples, density)


###################################################
### code chunk number 7: harvestr.Rnw:123-125
###################################################
plot(densities[[1]])
l_ply(densities, lines)


###################################################
### code chunk number 8: caching_run1
###################################################
unlink("harvestr-cache", recursive=TRUE)  # reset cache
system.time({
    posteriors1 <- harvest(datasets, MCMCregress, formula=y~x, cache=TRUE)
})


###################################################
### code chunk number 9: cache_run2
###################################################
system.time({
    posteriors2 <- harvest(datasets, MCMCregress, formula=y~x, cache=TRUE)
})



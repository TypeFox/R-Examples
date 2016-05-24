library("psychomix")
set.seed(1)

### Rost
r <- simRaschmix(design = "rost2", extremes = FALSE)
re <- simRaschmix(design = "rost2", extremes = TRUE)

mr <- raschmix(r, k = 2, nrep = 1, scores = "saturated")
mrs <- raschmix(r, k = 1:2, nrep = 1, scores = "saturated")
mre <- raschmix(data = re, k = 2, nrep = 1, scores = "saturated")
mres <- raschmix(data = re, k = 1:3, nrep = 1, scores = "saturated")

mr
mrs
mre
mres

options(digits = 4)
parameters(mr)

## mrrefit <- refit(mr)
## summary(mrrefit)

## -------------------------------------------------

## ### DIFSim
## data("DIFSim", package = "psychotree")
## #data("DIFSim", package = "psychotools")
## DIFSim.na <- DIFSim
## DIFSim.na$resp[1,1] <- NA
## DIFSim.na$resp[2,] <- NA
## DIFSim.na$age[3] <- NA
## #m2 <- raschmix(DIFSim$resp, k = 2, nrep = 1, type = "rost")
## m1 <- raschmix(DIFSim$resp, k = 2, nrep = 1, scores = "saturated",
##                control = list(iter.max = 20))
## m2 <- raschmix(resp ~ 1, DIFSim.na, k = 1, nrep = 1, scores = "saturated")
## ## concomitant
## m3 <- raschmix(resp ~ age + gender, data = DIFSim, k = 2, nrep = 1,
##                scores = "saturated")
## m3mv <- raschmix(resp ~ age + gender, data = DIFSim, k = 2, nrep = 1,
##                scores = "meanvar")
## #m4 <- raschmix(resp ~ age + gender, data = DIFSim[-(1:10),], k = 2:3, nrep = 1,
## #               scores = "saturated")

## ## print and summary
## m3
## summary(m3)
## ## m4
## ## m4.1 <- getModel(m4, which = 1)
## ## m4.1
## ## summary(m4.1)

## ## logLik
## logLik(m3)

## ## parameters
## parameters(m3, which = "concomitant", component = 2:1)
## parameters(m3, which = "score")
## score.probs(m3mv)

## ## weights
## ## weights(m3)
## ## ## flexmix requires integer weights
## ## w <- sample(1:nrow(DIFSim), nrow(DIFSim))
## ## m3w <- raschmix(resp ~ age + gender, data = DIFSim, k = 2, nrep = 1,
## ##                scores = "saturated", weights = w)
## ## weights(m3w)

## ## refit

## ## plot
## ## plot(m3)
## ## histogram(m3)
## ## m.nident <- raschmix(data = cbind(0,DIFSim$resp[,1:5],1,1,DIFSim$resp[,-(1:5)]),
## ##                      scores = "saturated", k = 3, nrep = 1)
## ## plot(m.nident)
## ## plot(m.nident, pch = 19:21, cex = matrix(rep((1+1:23)/10, 3), ncol = 3))
## ## plot(m.nident, index = FALSE, component = 1:2)
## ## plot(m.nident, index = TRUE, component = 1:2)


## ---- eval=TRUE----------------------------------------------------------
library(hisse)

## ---- eval=TRUE----------------------------------------------------------
library(diversitree)
set.seed(4)
# Essentially we are setting up a model that models the evolution of two binary characters
# Thus, we are assuming the following state combinations 1=00, 2=10, 3=01, 4=11:
pars <- c(0.1,0.1,0.1,0.2, rep(0.03, 4), 0.01,0.01,0, 0.01, 0, 0.01, 0.01,0,0.01, 0,0.01,0.01)
phy <- tree.musse(pars, max.taxa=50, x0=1, include.extinct=FALSE)
sim.dat <- data.frame(names(phy$tip.state), phy$tip.state)
# Now we want to make the states associated with the second character hidden from us. So, 
# we remove states 3 and 4 and make them 1 and 2
sim.dat[sim.dat[,2]==3,2] = 1
sim.dat[sim.dat[,2]==4,2] = 2
# This next step simply forces the character to be binary:
sim.dat[,2] = sim.dat[,2] - 1

## ---- eval=TRUE----------------------------------------------------------
head(sim.dat)

## ---- eval=TRUE----------------------------------------------------------
turnover.anc = c(1,1,0,0)
eps.anc = c(1,1,0,0)

## ---- eval=TRUE----------------------------------------------------------
turnover.anc = c(1,2,0,0)

## ---- eval=TRUE----------------------------------------------------------
turnover.anc = c(1,2,3,4)

## ---- eval=TRUE----------------------------------------------------------
eps.anc = c(0,0,0,0)

## ---- eval=TRUE----------------------------------------------------------
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates

## ---- eval=TRUE----------------------------------------------------------
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual

## ---- eval=TRUE----------------------------------------------------------
trans.rates.nodual.equal16 = ParEqual(trans.rates.nodual, c(1,6))
trans.rates.nodual.equal16

## ---- eval=TRUE----------------------------------------------------------
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

## ---- eval=TRUE----------------------------------------------------------
trans.rates.nodual.allequal = trans.rates.nodual
trans.rates.nodual.allequal[!is.na(trans.rates.nodual.allequal) & !trans.rates.nodual.allequal == 0] = 1
trans.rates.nodual.allequal

## ---- eval=TRUE----------------------------------------------------------
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
trans.rates.bisse

## ---- eval=FALSE---------------------------------------------------------
#  pp = hisse(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc,
#             eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal)

## ---- eval=TRUE----------------------------------------------------------
turnover.anc = c(1,2,0,3)
eps.anc = c(1,2,0,3)

## ---- eval=TRUE----------------------------------------------------------
trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates.nodual.no0B <- ParDrop(trans.rates, c(2,3,5,7,8,9,10,12))
trans.rates.nodual.no0B

## ---- eval=FALSE---------------------------------------------------------
#  pp = hisse(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc,
#             eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal, output.type="net.div")

## ---- eval=TRUE----------------------------------------------------------
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)

## ---- eval=TRUE----------------------------------------------------------
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))

## ---- eval=TRUE----------------------------------------------------------
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

## ---- eval=TRUE----------------------------------------------------------
# Now we want three specific rates:
trans.rates.nodual.threerates <- trans.rates.nodual
# Set all transitions from 0->1 to be governed by a single rate:
to.change <- cbind(c(1,3), c(2,4))
trans.rates.nodual.threerates[to.change] = 1
# Now set all transitions from 1->0 to be governed by a single rate:
to.change <- cbind(c(2,4), c(1,3))
trans.rates.nodual.threerates[to.change] = 2
# Finally, set all transitions between the hidden state to be a single rate (essentially giving 
# you an estimate of the rate by which shifts in diversification occur:
to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.nodual.threerates[to.change] = 3
trans.rates.nodual.threerates

## ---- eval=FALSE---------------------------------------------------------
#  pp = hisse(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc,
#             eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal)

## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=5)

## ---- eval=TRUE----------------------------------------------------------
load("testrecon1.Rsave")
class(pp.recon)
pp.recon

## ---- eval=TRUE----------------------------------------------------------
plot.hisse.states(pp.recon, rate.param="net.div", show.tip.label=FALSE)

## ---- eval=TRUE----------------------------------------------------------
plot.hisse.states(pp.recon, rate.param="net.div", show.tip.label=FALSE, rate.range=c(0,0.072))

## ---- eval=TRUE----------------------------------------------------------
pp.recon$aic

## ---- eval=FALSE---------------------------------------------------------
#  pp.recon = MarginRecon(phy, sim.dat, f=c(1,1), hidden.states=TRUE, pars=pp$solution,
#                         aic=pp$aic, n.cores=2)

## ---- eval=TRUE----------------------------------------------------------
hisse.results.list = list()
load("testrecon1.Rsave")
hisse.results.list[[1]] = pp.recon
load("testrecon2.Rsave")
hisse.results.list[[2]] = pp.recon
load("testrecon3.Rsave")
hisse.results.list[[3]] = pp.recon
# Now supply the list the plotting function
plot.hisse.states(hisse.results.list, rate.param="net.div", show.tip.label=FALSE, rate.range=c(0,0.072))

## ---- eval=FALSE---------------------------------------------------------
#  # First, suck in all the files with .Rsave line ending in your working directory:
#  files = system("ls -1 | grep .Rsave", intern=TRUE)
#  # Create an empty list object
#  hisse.results.list = list()
#  # Now loop through all files, adding the embedded pp.recon object in each
#  for(i in sequence(length(files))){
#    load(files[i])
#    hisse.results.list[[i]] = pp.recon
#    rm(pp.recon)
#  }


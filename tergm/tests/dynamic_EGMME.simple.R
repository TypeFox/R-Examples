#  File tests/dynamic_EGMME.simple.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(tergm)
n<-40
do.plot <- TRUE
g0<-network.initialize(n,dir=FALSE)

#                    edges, mean.age
target.stats<-c(     n*1/2,       20)
logit<-function(p)log(p/(1-p))
coef.exact<-function(density,duration)
    list(form=-log(((1+exp(logit(1-1/duration)))/(density/(1-density)))-1),
         diss=logit(1-1/duration))


truth <- coef.exact(target.stats[1]/network.dyadcount(g0),
                    target.stats[2])

# Get a deliberately bad starting network.
g1<-san(g0~meandeg,target.stats=target.stats[1],verbose=TRUE)

# Fit the model with very poor starting values.
set.seed(1)
dynfit<-stergm(g1,formation=~edges,dissolution=~edges, targets=~edges+mean.age, estimate="EGMME",target.stats=target.stats[-3],verbose=TRUE,control=control.stergm(SA.plot.progress=do.plot,SA.restart.on.err=FALSE,init.form=c(-log(.95/.05)),init.diss=c(1)))

print(summary(dynfit))
mcmc.diagnostics(dynfit)

stopifnot(all.equal(unlist(truth),unlist(coef(dynfit)),tol=0.01,check.attributes=FALSE))
},"simple EGMME")

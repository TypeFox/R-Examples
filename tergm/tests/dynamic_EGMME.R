#  File tests/dynamic_EGMME.R in package tergm, part of the Statnet suite
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

#                    edges, degree(1), mean.age
target.stats<-c(      n*1/2,    n*0.6,       20)

# Get a deliberately bad starting network.
g1<-san(g0~meandeg+degree(1),target.stats=target.stats[-3],verbose=TRUE)

coef.form <- c(-6.57, 1.01)
coef.diss <- c(2.944439)

# Fit the model with very poor starting values.
set.seed(1)
dynfit<-stergm(g1,formation=~edges+degree(1),dissolution=~offset(edges), targets="formation", estimate="EGMME", offset.coef.diss=coef.diss,target.stats=target.stats[-3],verbose=TRUE,control=control.stergm(SA.plot.progress=do.plot,SA.restart.on.err=FALSE,init.form=c(-log(.95/.05),0)))

print(summary(dynfit))
mcmc.diagnostics(dynfit)

stopifnot(all.equal(c(coef.form,coef.diss),unlist(coef(dynfit)),tol=0.01,check.attributes=FALSE))

# All parameters free, edges, degree(1), and edge.ages as target.
set.seed(1)
dynfit2<-stergm(g1,formation=~edges+degree(1),dissolution=~edges, targets=~edges+degree(1)+mean.age, estimate="EGMME", target.stats=target.stats,control=control.stergm(SA.plot.progress=do.plot,SA.plot.stats=TRUE))

print(summary(dynfit2))
mcmc.diagnostics(dynfit2)

stopifnot(all.equal(c(coef.form,coef.diss),unlist(coef(dynfit2)),tol=0.01,check.attributes=FALSE))
}, "EGMME")

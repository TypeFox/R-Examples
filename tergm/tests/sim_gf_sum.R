#  File tests/sim_gf_sum.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
library(tergm)

logit<-function(p)log(p/(1-p))

coef.form.f<-function(coef.diss,density) -log(((1+exp(coef.diss))/(density/(1-density)))-1)

simtest <- function(S, edges, dur, n, dir=FALSE, bip=0){
  m <- if(bip) bip else n
  
  target.stats<-edges
  coef.diss<-logit(1-1/dur)
  
  dc <- matrix(rnorm(n*m),n,m)
  if(!dir && !bip) dc <- dc+t(dc)
  # Note that if nw1!=nw2, then
  # summary(nw1~edgecov(dc)+edgecov.ages(dc))!=summary(nw2~edgecov(dc)+edgecov.ages(dc))
  # with probability 1. Thus, they offer a sort of a checksum of a network.
  
  g0<-network.initialize(n, dir=dir, bipartite=bip)
  g0 %n% "dc" <- dc
  g0 %n% "time" <- 0
  g0 %n% "lasttoggle" <- -1-rgeom(network.dyadcount(g0),1/4)


  
  dyads <- network.dyadcount(g0)
  density<-edges/dyads
  
  coef.form<-coef.form.f(coef.diss,density)
  
  # Get a reasonably close starting network.
  g1<-san(g0~edges,target.stats=target.stats,verbose=TRUE)
  
  print(coef.form)
  print(coef.diss)
  
  # Simulate. Starting from an ordinary network:
  dynsim<-simulate(g1,formation=~edges,dissolution=~edges,coef.form=coef.form,coef.diss=coef.diss,time.slices=S,verbose=TRUE,monitor=~edgecov("dc")+edgecov.ages("dc"))
  
  # Resuming from a networkDynamic:
  dynsim2<-simulate(dynsim,formation=~edges,dissolution=~edges,coef.form=coef.form,coef.diss=coef.diss,time.slices=S,verbose=TRUE,monitor=~edgecov("dc")+edgecov.ages("dc"))
  
  # Resuming from a resumed networkDynamic:
  dynsim3<-simulate(dynsim2,formation=~edges,dissolution=~edges,coef.form=coef.form,coef.diss=coef.diss,time.slices=S,verbose=TRUE,monitor=~edgecov("dc")+edgecov.ages("dc"))

  sim.stats <- attr(dynsim3,"stats")
  
  # Replay the simulation using a networkDynamic:
  gf1.stats <- as.matrix(tergm.godfather(dynsim3~edgecov("dc")+edgecov.ages("dc"), start=0, end=S*3))
  
  # Replay the sim using networkDynamic with no explicit time params (should be using net.obs.period)
  gf1b.stats <- as.matrix(tergm.godfather(dynsim3~edgecov("dc")+edgecov.ages("dc")))
  
  # Replay the simulation using a network + list of changes:
  gf2.stats <- as.matrix(tergm.godfather(g1~edgecov("dc")+edgecov.ages("dc"), start=0, end=S*3, changes=attr(dynsim3,"changes")))
  
  # Evaluate summary statistics directly:
  summ.stats <- summary(dynsim3~edgecov("dc")+edgecov.ages("dc"), at=1:(S*3))
  
  # If they aren't all identical, we are in trouble.
  stopifnot(all.equal(sim.stats,gf1.stats),
            all.equal(sim.stats,gf1b.stats),
            all.equal(sim.stats,gf2.stats),
            all.equal(sim.stats,summ.stats))
}

set.seed(1)
cat("Undirected, unipartite:\n")
simtest(20, 10, 10, 10, FALSE, FALSE)
cat("Directed, unipartite:\n")
simtest(20, 10, 10, 10, TRUE, FALSE)
cat("Undirected, biipartite:\n")
simtest(20, 10, 10, 10, FALSE, 3)

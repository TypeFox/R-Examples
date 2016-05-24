#  File tests/simulate_networkDynamic.R in package tergm, part of the Statnet suite
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

# Construct a network with 20 nodes and 20 edges
n<-20
target.stats<-edges<-20
g0<-network.initialize(n,dir=TRUE)
g1<-san(g0~edges,target.stats=target.stats,verbose=TRUE)

S<-10

# To get an average duration of 10...
duration<-10
coef.diss<-logit(1-1/duration)

# To get an average of 20 edges...
dyads<-network.dyadcount(g1)
density<-edges/dyads
coef.form<-coef.form.f(coef.diss,density)



# Simulate a networkDynamic for five steps
dynsim<-simulate(g1,formation=~edges,dissolution=~edges,coef.form=coef.form,coef.diss=coef.diss,time.slices=S)

# ----- check net.obs.period encoding ----

# should have a net.obs.period object from 0 to 5
if (!all(unlist((dynsim%n%'net.obs.period')$observations)==c(0,1,1,11))){
  stop("simulate.network did not encode net.obs.period$observation correctly")
}

# "Resume" the simulation for five steps
dynsim2<-simulate(dynsim,time.slices=S)

if (!all(unlist((dynsim2%n%'net.obs.period')$observations)==c(0,1,1,11,11,21))){
  stop("simulate.networkDynamic did not encode net.obs.period$observation correctly")
}

# check a realistic example of starting a sim from scratch
mean.rel.dur <- 10
msm.sim <- network.initialize(1000,directed=F)
activate.vertices(msm.sim,-Inf,Inf)
set.network.attribute(msm.sim,'net.obs.period',list(observations=list(c(-1,0)),
                                                    mode="discrete", time.increment=1,time.unit="step"))
formation <- ~edges
dissolution <- ~offset(edges)
target.stats <- 400
coef.diss <- log(mean.rel.dur-1)
formation.with.stnet <- update.formula(formation,msm.startnet~.)
msm.startnet <- network.collapse(msm.sim,at=0)
msm.est <- ergm(formation.with.stnet,target.stats=target.stats)
coef.form <- msm.est$coef
coef.form[1] <- coef.form[1] - coef.diss
msm.edgelist <- as.edgelist(simulate(msm.est))
add.edges(msm.sim,msm.edgelist[,1],msm.edgelist[,2])
activate.edges(msm.sim, -Inf, Inf)

#first timestep 
msm.sim <- simulate(msm.sim,
                    formation=formation,
                    dissolution=~edges,
                    coef.form=coef.form,
                    coef.diss=coef.diss,
                    time.slices = 1,
                    monitor="all",
                    verbose=T
)
#second timestep
msm.sim <- simulate(msm.sim,
                    formation=formation,
                    dissolution=~edges,
                    coef.form=coef.form,
                    coef.diss=coef.diss,
                    time.slices = 1,
                    monitor="all",
                    verbose=T
)

if(!all(unlist((msm.sim%n%'net.obs.period')$observations)==c(-1,  0,  0,  1,  1,  2))){
  stop('net.obs.period (and simulation) was not constructed as expected when stopping and restarting simulate.networkdynamic')
}


# test a sim with vertex dynamics applied after the sim stage
mean.rel.dur <- 10
msm.sim <- network.initialize(1000,directed=F)
activate.vertices(msm.sim,-Inf,Inf)
set.network.attribute(msm.sim,'net.obs.period',list(observations=list(c(-1,0)),
                                                    mode="discrete", time.increment=1,time.unit="step"))
formation <- ~edges
dissolution <- ~offset(edges)
target.stats <- 400
coef.diss <- log(mean.rel.dur-1)
formation.with.stnet <- update.formula(formation,msm.startnet~.)
# simulate a set of edges to use as the starting point for the network
msm.startnet <- network.collapse(msm.sim,at=0)
msm.est <- ergm(formation.with.stnet,target.stats=target.stats)
coef.form <- msm.est$coef
coef.form[1] <- coef.form[1] - coef.diss
msm.edgelist <- as.edgelist(simulate(msm.est))
add.edges(msm.sim,msm.edgelist[,1],msm.edgelist[,2])
activate.edges(msm.sim, -Inf, Inf)

# simulate first timestep (0,1)
msm.sim <- simulate(msm.sim,
                    formation=formation,
                    dissolution=~edges,
                    coef.form=coef.form,
                    coef.diss=coef.diss,
                    time.slices = 1,
                    monitor="all",
                    verbose=T
)

# toggle off vertices for the step just simulated
msm.sim<-deactivate.vertices(msm.sim,v=sample(which(is.active(msm.sim,v=1:network.size(msm.sim),at=0)),size=10),onset=0,terminus=Inf,deactivate.edges=TRUE)
add.vertices.active(msm.sim,nv=5,onset=0,terminus=Inf)


# simulate second timestep (1,2)
msm.sim <- simulate(msm.sim,
                    formation=formation,
                    dissolution=~edges,
                    coef.form=coef.form,
                    coef.diss=coef.diss,
                    time.slices = 1,
                    monitor="all",
                    verbose=T
)

msm.sim<-deactivate.vertices(msm.sim,v=sample(which(is.active(msm.sim,v=1:network.size(msm.sim),at=1)),size=10),onset=1,terminus=Inf,deactivate.edges=TRUE)
add.vertices.active(msm.sim,nv=1,onset=1,terminus=Inf)

# check for correct activity
if(!all(sapply(-1:2,function(t){network.size.active(msm.sim,at=t)})==c(1000,995,986,986))){
  stop('vertex dynamics did not seem to generate correct activity in stergm sim')
}

if(!all(network.dynamic.check(msm.sim,complete=FALSE)$dyad.checks)){
  stop('vertex dynamics in stergm sim created inconsistent dyads')
}

# at this point msm.sim has observations as far as 2.  Try sampling at 2 and upating from (2,3)
msm.sim <- simulate(msm.sim,
                    formation=formation,
                    dissolution=~edges,
                    coef.form=coef.form,
                    coef.diss=coef.diss,
                    time.slices = 1, time.start=2,time.offset=0,
                    monitor="all",
                    verbose=T
)

if (max(unlist((msm.sim%n%'net.obs.period')$observations)) > 3){
  stop("simulate.networkDynamic updated net.obs.period inappropriately when time.start and time.offset was explicitly set" )
}
if (!all(get.change.times(msm.sim)==0:2)){
  stop("simulate.networkDynamic updated edge spells inappropriately when time.start and time.offset was explicitly set")
}

# try it again to make sure nothing got mangled
msm.sim <- simulate(msm.sim,
                    formation=formation,
                    dissolution=~edges,
                    coef.form=coef.form,
                    coef.diss=coef.diss,
                    time.slices = 1, time.start=3,time.offset=0,
                    monitor="all",
                    verbose=T
)
if (max(unlist((msm.sim%n%'net.obs.period')$observations)) > 4){
  stop("simulate.networkDynamic updated net.obs.period inappropriately when time.start and time.offset was explicitly set" )
}
if (!all(get.change.times(msm.sim)==0:3)){
  stop("simulate.networkDynamic updated edge spells inappropriately when time.start and time.offset was explicitly set")
}


# check that vertex ids in changes are correctly translated
dyn<-as.networkDynamic(network.initialize(4))
deactivate.vertices(dyn,v=1)
# define stergm that should toggle on all ties
changes<-simulate.networkDynamic(dyn,formation=~edges,dissolution=~edges,coef.form=1,coef.diss=0,output='changes')
# check if any changes involve vertex 1 (shouldn't because it is inactive)
if(any(changes[,2:3]==1)){
  stop("simulate.networkDynamic returned changes involving an inactive vertex" )
}


mean.rel.dur <- 10
msm.sim <- as.networkDynamic(network.initialize(1000,directed=F))
#set.network.attribute(msm.sim,'net.obs.period',list(observations=list(c(-1,0)),
#                                                    mode="discrete", time.increment=1,time.unit="step"))
formation <- ~edges
dissolution <- ~offset(edges)
target.stats <- 400
coef.diss <- log(mean.rel.dur-1)
formation.with.stnet <- update.formula(formation,msm.startnet~.)
# simulate a set of edges to use as the starting point for the network
msm.startnet <- network.collapse(msm.sim,at=0)
msm.est <- ergm(formation.with.stnet,target.stats=target.stats)
coef.form <- msm.est$coef
coef.form[1] <- coef.form[1] - coef.diss
msm.edgelist <- as.edgelist(simulate(msm.est))
add.edges(msm.sim,msm.edgelist[,1],msm.edgelist[,2])

msm.sim <- simulate(msm.sim,
                    formation=formation,
                    dissolution=~edges,
                    coef.form=coef.form,
                    coef.diss=coef.diss,
                    time.slices = 1, time.start=0,time.offset=0,
                    monitor="all",
                    verbose=T
)


if(length((msm.sim%n%'net.obs.period')$observations)>1){
  stop('simulate.networkDynamic mangled net.obs.period with time.offset=0 argument')
}


# --- check that adding and deleting vertices behaves properly ---
g0<-network.initialize(20,dir=TRUE)
g1<-san(g0~edges,target.stats=20,verbose=TRUE)
# Simulate a networkDynamic from static 
dynsim<-simulate(g1,formation=~edges,dissolution=~edges,
                 coef.form=-9.326322,coef.diss=2.197225,
                 time.slices=1,verbose=TRUE)
# add more vertices
add.vertices(dynsim,3)
# simulate another step
dynsim<-simulate(dynsim,formation=~edges,dissolution=~edges,
                 coef.form=-9.326322,coef.diss=2.197225,
                 time.slices=1,verbose=TRUE)
# expect to see a pid defined, but will just be numeric values
if (dynsim%n%'vertex.pid'!="tergm_pid"){
  stop("simulate.networkDynamic net not add a vertex.pid to the network as expected")
}
if(!all(is.numeric(dynsim%v%"tergm_pid"))){
  stop("simulate.networkDynamic did not create new numeric vertex pid as expected")
}


# add more (at this point there should be some non-numeric pids created)
add.vertices(dynsim,3)
if(all(is.numeric(dynsim%v%"tergm_pid"))){
  stop("add.vertices did not create new non-numeric vertex pids as expected")
}

# simulate
dynsim<-simulate(dynsim,formation=~edges,dissolution=~edges,
                 coef.form=-9.326322,coef.diss=2.197225,
                 time.slices=1,verbose=TRUE)
# delete some vertices
delete.vertices(dynsim,vid = 5:10)
dynsim<-simulate(dynsim,formation=~edges,dissolution=~edges,
                 coef.form=-9.326322,coef.diss=2.197225,
                 time.slices=1,verbose=TRUE)

# test with an externally defined PID  --------------

# Simulate a networkDynamic from static 
dynsim<-simulate(g1,formation=~edges,dissolution=~edges,
                 coef.form=-9.326322,coef.diss=2.197225,
                 time.slices=1,verbose=TRUE)
dynsim%n%'vertex.pid'<-'letters'
dynsim%v%'letters'<-LETTERS[1:20]
# simulate another step
dynsim<-simulate(dynsim,formation=~edges,dissolution=~edges,
                 coef.form=-9.326322,coef.diss=2.197225,
                 time.slices=1,verbose=TRUE)
# expect to see a pid defined, but will just be numeric values
if (dynsim%n%'vertex.pid'!="letters"){
  stop("simulate.networkDynamic vertex pid missing")
}
if(!all(is.na(dynsim%v%"tergm_pid"))){
  stop("simulate.networkDynamic created a tergm_pid when a vertex.pid already existed")
}



# ---- check summary.statistics.networkDynamic ---
# this was giving error in #958
my.nD <- network.initialize(100,directed=F)
activate.vertices(my.nD, onset=0, terminus = 10)
class(my.nD)
summary(my.nD~isolates, at=1)


# --- checking simulate with Inf coefficients ---
# was giving error to #995
n <- 50
nw <- network.initialize(n, directed = FALSE)
nw <- set.vertex.attribute(nw, "loc", rep(0:1, each = n/2))

fit <- ergm(nw ~ edges + offset(nodemix("loc", base=c(1, 3))),
            target.stats = 15,
            offset.coef = -Inf)
coef.form <- fit$coef
coef.form[1] <- coef.form[1] - log(40-1)

fit.sim <- simulate(fit)
diag.sim <- simulate(fit.sim,
                     formation = ~ edges + offset(nodemix("loc", base=c(1, 3))),
                     dissolution = ~offset(edges),
                     coef.form = coef.form,
                     coef.diss = log(40-1),
                     time.slices = 100,
                     monitor = "formation",
                     nsim = 1)

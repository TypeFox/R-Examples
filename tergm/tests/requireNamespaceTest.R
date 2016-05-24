library(statnet.common)
opttest({
library(network)

requireNamespace('ergm')
requireNamespace('tergm')
par(ask=FALSE)
n<-30
g0<-network.initialize(n,dir=FALSE)

#                     edges, degree(1), mean.age
target.stats<-c(      n*1/2,    n*0.6,        20)

dynfit<-tergm::stergm(g0,formation = ~edges+degree(1), dissolution = ~edges,
               targets = ~edges+degree(1)+mean.age,
               target.stats=target.stats, estimate="EGMME",
               control=tergm::control.stergm(SA.restart.on.err=FALSE))

par(ask=TRUE)
ergm::mcmc.diagnostics(dynfit)
summary(dynfit)

# use a summary formula to display number of isolates and edges
                                        # at discrete time points
my.nD <- simulate(dynfit, nsim=1, time.slices=100, output="networkDynamic")
ergm::summary.formula(my.nD~isolates+edges+mean.age, at=1:10)
ergm::summary.statistics.formula(my.nD~mean.age)
}, "requireNamespace")

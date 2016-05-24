# tests for temporal density functions
require(tsna)
require(testthat)
require(networkDynamicData)

# tests of argument processing
data(moodyContactSim)

expect_error(tEdgeDensity(0),regexp='first argument is a networkDynamic')
expect_error(tEdgeDensity(moodyContactSim,mode='foo'))
expect_error(tEdgeDensity(moodyContactSim,agg.unit='foo'))


# density of a network with no edges
expect_equal(tEdgeDensity(as.networkDynamic(network.initialize(3)),mode='event'),0)
# density of network size 0
expect_equal(tEdgeDensity(as.networkDynamic(network.initialize(0)),mode='event'),0)

# complete network, explicit timing
complete<-as.networkDynamic(network.initialize(5))
complete[,]<-1
activate.edges(complete,onset=0,terminus=100)
expect_equal(tEdgeDensity(complete,mode='event'),0.01)
expect_equal(tEdgeDensity(complete,agg.unit='dyad'),1)
expect_equal(tEdgeDensity(complete),1)

# density of a complete network with no explicit timing
complete<-as.networkDynamic(network.initialize(5))
complete[,]<-1
expect_equal(tEdgeDensity(complete,mode='event'),1)
expect_equal(tEdgeDensity(complete,agg.unit='dyad'),1)
expect_equal(tEdgeDensity(complete),1)

# density of a complete network with no explicit timing, active.default=FALSE
expect_equal(tEdgeDensity(complete,active.default=FALSE,agg.unit='dyad'),0)
expect_equal(tEdgeDensity(complete,active.default=FALSE,agg.unit='edge'),0)

# density of a complete network with no explicit timing
complete<-as.networkDynamic(network.initialize(5))
complete[,]<-1
expect_equal(tEdgeDensity(complete,mode='event'),1)
expect_equal(tEdgeDensity(complete,agg.unit='dyad'),1)
expect_equal(tEdgeDensity(complete),1)


# single edge, implicit timing
half<-as.networkDynamic(network.initialize(2))
half[1,2]<-1
expect_equal(tEdgeDensity(half,mode='event'),1)
expect_equal(tEdgeDensity(half,agg.unit='dyad'),0.5)
expect_equal(tEdgeDensity(half),1)

# single undirected edge, implicit timing
half<-as.networkDynamic(network.initialize(2,directed=FALSE))
half[1,2]<-1
expect_equal(tEdgeDensity(half,mode='event'),1)
expect_equal(tEdgeDensity(half,agg.unit='dyad'),1)
expect_equal(tEdgeDensity(half),1)

# two edges, each half active
half<-network.initialize(2)
add.edges.active(half,tail=1:2,head=2:1,onset=0:1,terminus=1:2)
expect_equal(tEdgeDensity(half,mode='event'),0.5)
expect_equal(tEdgeDensity(half,agg.unit='dyad'),0.5)
expect_equal(tEdgeDensity(half),0.5)

# two edges, each half active in a network of 4
half<-network.initialize(4)
add.edges.active(half,tail=1:2,head=2:1,onset=0:1,terminus=1:2)
expect_equal(tEdgeDensity(half,mode='event'),0.5)
expect_equal(tEdgeDensity(half,agg.unit='dyad'),2/24)
expect_equal(tEdgeDensity(half),0.5)

# single edge in range defined by net.obs
obs<-as.networkDynamic(network.initialize(2))
add.edges.active(obs,tail=1,head=2,onset=1,terminus=2)
obs%n%'net.obs.period'<-list(observations=list(c(0,3)),mode="discrete", time.increment=1,time.unit="step")
expect_equal(tEdgeDensity(obs,mode='event'),1/3)
expect_equal(tEdgeDensity(obs,agg.unit='dyad'),1/6)
expect_equal(tEdgeDensity(obs),1/3)

# single edge in network with loop
loop<-as.networkDynamic(network.initialize(2,loops=TRUE))
loop[1,2]<-1
expect_equal(tEdgeDensity(loop,mode='event'),1)
expect_equal(tEdgeDensity(loop,agg.unit='dyad'),0.25) # loop increases possible dyads to 4
expect_equal(tEdgeDensity(loop),1)


# test on example networks 
# have no valid baseline to compare, so just comparing to previous output
data(concurrencyComparisonNets)
expect_equal(tEdgeDensity(base,mode='event'),0.009809139,tolerance=0.0000001)
expect_equal(tEdgeDensity(base),0.1996973,tolerance=0.0000001)
expect_equal(tEdgeDensity(base,agg.unit='dyad'),0.000766006,tolerance=0.0000001)

expect_equal(tEdgeDensity(middle, mode='event'),0.00981934,tolerance=0.0000001)
expect_equal(tEdgeDensity(middle),0.1921979,tolerance=0.0000001)
expect_equal(tEdgeDensity(middle,agg.unit='dyad'),0.0007387788,tolerance=0.0000001)

expect_equal(tEdgeDensity(monog, mode='event'),0.009819324,tolerance=0.0000001)
expect_equal(tEdgeDensity(monog),0.1979761,tolerance=0.000001)
expect_equal(tEdgeDensity(monog,agg.unit='dyad'),0.0007617818,tolerance=0.0000001)




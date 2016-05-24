# some tests for the basic durations function 
require(tsna)
require(testthat)
require(networkDynamicData)

# ------ edgeDuration tests -------

data(moodyContactSim)
expect_equal(edgeDuration(moodyContactSim),c(32, 33, 32, 26, 30, 24, 32, 30, 26, 27, 27, 32, 31, 27, 26, 34, 44, 26))
# check with counts
expect_equal(edgeDuration(moodyContactSim,mode='counts'),rep(1,18))

# basic testing network
test<-network.initialize(5)
test[,]<-1
activate.edges(test,onset=0,terminus=2)
activate.edges(test,onset=5,terminus=6)

# spell level
expect_equal(edgeDuration(test,subject='spells'),c(2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1))
# edge level
expect_equal(edgeDuration(test,subject='edges'),rep(3,network.edgecount(test)))
expect_equal(edgeDuration(test,subject='dyads'),rep(3,network.edgecount(test)))


# test for dyads
test<-network.initialize(2)
add.edges(test,1,2)
add.edges(test,1,2)
add.edges(test,1,2)
expect_equal(edgeDuration(test,subject='edges',start=0,end=2),c(2,2,2))
expect_equal(edgeDuration(test,subject='dyads',start=0,end=2),6)
# test counts
expect_equal(edgeDuration(test,mode='count',subject='edges',start=0,end=2),c(1,1,1))
expect_equal(edgeDuration(test,mode='count',subject='dyads',start=0,end=2),3)


# ------ vertexDuration  tests ------

test<-network.initialize(3)
activate.vertices(test,onset = c(0,1),terminus=c(1,2),v=1:2)
vertexDuration(test) # this returns 1,1,Inf,  inf should be truncated by default

# mode option switch to counts?
expect_equal(vertexDuration(test,mode='counts'),c(1,1,1))

# does it trucate to net obs period appropriately?
test%n%'net.obs.period'<-list(observations=list(c(0,3)),mode="discrete", time.increment=1,time.unit="step")
expect_equal(vertexDuration(test),c(1,1,3))

# does it do spell level appropriately?
activate.vertices(test,v=1,onset=2,terminus=3)
expect_equal(vertexDuration(test,subject='spells'),c(1,1,1,3))
expect_equal(vertexDuration(test,subject='vertices'),c(2,1,3))

# test v argument
expect_equal(vertexDuration(test,v=2),1)

# does output for real example stay consistant?
data(windsurfers)
expect_equal(vertexDuration(windsurfers),c(27, 13, 19, 22, 5, 8, 5, 2, 1, 4, 5, 9, 7, 8, 7, 18, 9, 4, 2, 1, 5, 3, 6, 6, 4, 6, 7, 6, 14, 5, 4, 1, 1, 2, 4, 1, 8, 3, 14, 8, 8, 3, 3, 7,3, 3, 1, 8, 2, 3, 6, 4, 3, 1, 2, 3, 2, 4, 7, 7, 2, 3, 3, 2, 5, 6, 5, 3, 6, 7, 2, 1, 1, 14, 2, 2, 2, 5, 3, 2, 2, 2, 2, 4, 3, 1, 3, 2, 3, 2, 1, 2, 1, 1, 1))
expect_equal(vertexDuration(windsurfers,mode='counts'),c(4, 6, 8, 7, 4, 5, 3, 2, 1, 3, 4, 4, 6, 5, 6, 9, 7, 3, 2, 1, 4, 2, 3, 4, 3, 3, 5, 4, 5, 4, 3, 1, 1, 1, 3, 1, 5, 3, 8, 7, 6, 2, 3, 5, 2, 3, 1, 6, 2, 2, 4, 4, 3, 1, 1, 3, 1, 3, 6, 4, 1, 3, 3, 2, 4, 5, 3, 2, 5, 4, 2, 1, 1, 7, 2, 2, 2, 4, 3, 2, 2, 1, 2, 2, 2, 1, 3, 1, 3, 2, 1, 2, 1, 1, 1))


# ------- tiedDuration ------

# directed case (out)
test<-network.initialize(3)
add.edges.active(test,1,2,onset=0,terminus=2)
add.edges.active(test,1,3,onset=1,terminus=2)
expect_equal(tiedDuration(test),c(3,0,0))
# un directed
expect_equal(tiedDuration(test,neighborhood = 'combined'),c(3,2,1))
# reversed
expect_equal(tiedDuration(test,neighborhood = 'in'),c(0,2,1))

# counts
expect_equal(tiedDuration(test,mode='counts'),c(2,0,0))

# active default on non-dynamic netowrk
test2<-network.initialize(3)
add.edges(test2,1,2)
expect_equal(tiedDuration(test2,active.default = TRUE), c(1,0,0))
expect_equal(tiedDuration(test2,active.default = FALSE),c(0,0,0))

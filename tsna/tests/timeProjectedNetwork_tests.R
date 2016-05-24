# tests for the time-projected functions

require(tsna)
library(networkDynamicData)
require(testthat)

# trivial example network, directed case
test<-network.initialize(3)
add.edges.active(test,tail=1,head=2,onset=0,terminus=1)
add.edges.active(test,tail=2,head=3,onset=1,terminus=2)
testProj<-timeProjectedNetwork(test,start=0,end=2)
expect_equal(as.matrix(testProj),
              matrix(
              c(0, 1, 0, 1, 0, 0,
                0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 1,
                0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 1,
                0, 0, 0, 0, 0, 0),ncol=6,byrow=TRUE),check.attributes=FALSE)

# undirected case
test%n%'directed'<-FALSE
testProj<-timeProjectedNetwork(test,start=0,end=2)
expect_equal(as.matrix(testProj),
             matrix(
               c(0, 1, 0, 1, 0, 0,
                 1, 0, 0, 0, 1, 0,
                 0, 0, 0, 0, 0, 1,
                 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0, 1,
                 0, 0, 0, 0, 1, 0),ncol=6,byrow=TRUE),check.attributes=FALSE)


# test edge attribute
data("nd_test_nets")
test<-nd_test_nets[[27]]
activate.edge.attribute(test,'foo',onset=1,terminus=2,'A')
activate.edge.attribute(test,'foo',onset=2,terminus=3,'B')
activate.edge.attribute(test,'foo',onset=3,terminus=4,'C')
activate.edge.attribute(test,'foo',onset=4,terminus=5,'D')
activate.edge.attribute(test,'foo',onset=5,terminus=6,'E')
activate.edge.attribute(test,'foo',onset=6,terminus=7,'F')
activate.edge.attribute(test,'foo',onset=7,terminus=8,'G')
testProj<-timeProjectedNetwork(test,start=0,end=10)
expect_equal(get.edge.attribute(testProj,'foo'),c("A", "B", "C", "D", "E", "F", "G", NA , NA))
expect_equal(get.edge.attribute(testProj,'pid'),1:9)
# in the undirected case, edges should be copied twice
test%n%'directed'<-FALSE
testProj<-timeProjectedNetwork(test,start=0,end=10)
expect_equal(get.edge.attribute(testProj,'foo'),c("A", "A","B","B", "C","C","D", "D", "E","E","F", "F", "G","G",NA,NA, NA , NA))
expect_equal(get.edge.attribute(testProj,'pid'),c(1 ,1 ,2 ,2 ,3 ,3 ,4 ,4 ,5, 5, 6, 6 ,7, 7, 8, 8 ,9 ,9))

data(moodyContactSim)
changes<-get.change.times(moodyContactSim)
moodyProj<-timeProjectedNetwork(moodyContactSim,onsets=changes,termini=changes)
# make sure the listing didn't get mangled for NA
expect_equal(length(moodyProj%e%'na'),length(get.edge.attribute(moodyProj,'na',unlist=FALSE)))

data(harry_potter_support)
hpProj<-timeProjectedNetwork(harry_potter_support)
plot(hpProj,arrowhead.cex = 0,edge.col=ifelse(hpProj%e%'edge.type'=='within_slice','black','gray'),vertex.cex=0.7)

# check that specific slices copied correctly
# WHY DOES THIS FAIL?
#expect_equal(as.matrix(network.extract(harry_potter_support,at=5)),as.matrix(hpProj)[(64*4+1):(64*5),(64*4+1):(64*5)])

# check that vertex attributes copied
expect_equal((hpProj%v%'gender')[1:64],harry_potter_support%v%'gender')
expect_equal((hpProj%v%'gender')[65:128],harry_potter_support%v%'gender')

expect_equal(length(network.vertex.names(hpProj)),network.size(hpProj))

# check edge type added
expect_true("edge.type"%in%list.edge.attributes(hpProj))


moodyProj<-timeProjectedNetwork(moodyContactSim,time.increment=100)

# correct size of new network?
expect_equal(network.size(moodyContactSim)*(moodyContactSim%n%'net.obs.period')$observations[[1]][2]/100,network.size(moodyProj))

# create network from changes
changes<-get.change.times(moodyContactSim)
moodyProjChange<-timeProjectedNetwork(moodyContactSim,onsets=changes,termini =changes)

# test for some vertex inactivity
data(windsurfers)
windProj<-timeProjectedNetwork(windsurfers,start=0,end=5)

# gplot3d(windProj,edge.col=ifelse(proj%e%'edge.type'=='within_slice','black','gray'))

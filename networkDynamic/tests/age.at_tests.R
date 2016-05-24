# tests for age.at functions
library(networkDynamic)
library(testthat)

# ---------- edges.age.at tests ------
test<-network.initialize(5)
add.edges(test,tail = 1:4,head=2:5 )
activate.edges(test,onset=0:2,terminus=1:3,e=1:3)


expect_equal(edges.age.at(test,at=0),c(0,NA,NA,Inf))
expect_equal(edges.age.at(test,at=2.5),c(NA,NA,0.5,Inf))

# test active default
expect_equal(edges.age.at(test,at=0,active.default=FALSE),c(0,NA,NA,NA))

# test deleted edge
delete.edges(test,e=2)

expect_equal(edges.age.at(test,at=1),c(NA,NA,NA,Inf))

# test restricted edge set
expect_equal(edges.age.at(test,at=1,e=c(1,3,4)),c(NA,NA,Inf))


# network size 0
expect_equal(edges.age.at(network.initialize(0),at=0),logical(0))

# no edges in network
expect_equal(edges.age.at(network.initialize(3),at=0),logical(0))


# ------------------- test for dyads.age.at -------
test<-network.initialize(5)
add.edges(test,tail = 1:4,head=2:5 )
activate.edges(test,onset=0:2,terminus=1:3,e=1:3)

expect_equal(dyads.age.at(test,at=0,tail=1,head=2),0)
expect_equal(dyads.age.at(test,at=0,tail=2,head=1),NA) # reverse order
expect_equal(dyads.age.at(test,at=1,tail=1,head=2),NA)

# check for out of sync tail and head
expect_error(dyads.age.at(test,at=0,tail=1,head=2:3),regexp = 'must be the same length')

# check for multiple edges
expect_equal(dyads.age.at(test,at=1,tail=1:2,head=2:3),c(NA,0))
expect_equal(dyads.age.at(test,at=0,tail=1:2,head=2:3),c(0,NA))

# non existing edges
expect_equal(dyads.age.at(test,at=1,tail=5,head=1),NA)

# head or tail out of range
expect_equal(dyads.age.at(test,at=1,tail=5,head=6),NA) # no error returned

# deleted edges
delete.edges(test,e=2)
expect_equal(dyads.age.at(test,at=1,tail=1:2,head=2:3),c(NA,NA))

# active default test
expect_equal(dyads.age.at(test,at=0,tail=4,head=5),Inf)
expect_equal(dyads.age.at(test,at=0,tail=4,head=5,active.default=FALSE),NA)

# test edge format arge
expect_equal(dyads.age.at(test,at=0),c(0,NA,Inf))
out<-dyads.age.at(test,at=0,format.out='edgelist')
expect_equal(out[,'tails'],c(1,3,4))
expect_equal(out[,'heads'],c(2,4,5))
expect_equal(out[,'ages'],c(0,NA,Inf))

out<-dyads.age.at(test,at=0,format.out='matrix')
expect_equal(out,matrix(c(NA,  NA,  NA,  NA,  NA,   0,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA, Inf, NA),nrow=5,ncol=5))

# test symmetry for non directed case
test%n%'directed'<-FALSE
out<-dyads.age.at(test,at=0,format.out='matrix')
expect_equal(out,matrix(c(NA,   0,  NA,  NA,  NA,   0,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA,  NA, Inf,  NA,  NA,  NA, Inf,  NA),nrow=5,ncol=5))


# ----- test vertices.age.at -----
test<-network.initialize(5)
activate.vertices(test,v = 1:4,onset=0:3,terminus=1:4)
expect_equal(vertices.age.at(test,at=0),c(0,NA,NA,NA,Inf))
expect_equal(vertices.age.at(test,at=0,active.default=FALSE),c(0,NA,NA,NA,NA))
expect_equal(vertices.age.at(network.initialize(0),at=0),logical(0))
expect_equal(vertices.age.at(test,at=1,v=c(1,2,4,5)),c(NA,0,NA,Inf))



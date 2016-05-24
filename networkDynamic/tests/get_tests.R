########################################################################
# This file contains the testing suite for the "get" methods, i.e.:
#      get.edgeIDs.active
#      get.edges.active
#      get.neighborhood.active     
#      get.change.times
#########################################################################

require(networkDynamic)
require(testthat)

#---------------- GET.CHANGE.TIMES TESTS ------------------------
# Notes:
#  -- minimal level of testing done here. sorry -- alc 6/28/12
#-------------------------------------------------------------------

cat("testing get.change.times ... ")

# a case with no change times
data(flo)
net1 = network(flo)
t1 = get.change.times(net1)
a1 = (length(t1) == 0)

# a case with only edge change times
activate.edges(net1, onset=1:20, terminus=101:120)
t2 = get.change.times(net1, edge.activity=F)
t3 = get.change.times(net1)
a2 = (length(t2) == 0)
a3 = all(t3==c(1:20, 101:120))

# a case with only vertex change times
net2 = network(flo)
activate.vertices(net2, at=seq(2,32,2))
t4 = get.change.times(net2, vertex.activity=F)
t5 = get.change.times(net2)
a4 = (length(t4) == 0)
a5 = all(t5 == seq(2,32,2))

# a case with both types of change times
activate.edges(net2, at=60:99)
t6 = get.change.times(net2, vertex.activity = F)
t7 = get.change.times(net2, edge.activity = F)
t8 = get.change.times(net2)
a6 = all(t6 == 60:99)
a7 = all(t7 == seq(2,32,2))
a8 = all(t8 == c(seq(2,32,2), 60:99))

# a case with infinity
deactivate.edges(net1)
activate.edges(net2)
t9 = get.change.times(net1)
t10 = get.change.times(net2)
t11 = get.change.times(net1, ignore.inf=F)
t12 = get.change.times(net2, ignore.inf=F)
a9 = (length(t9) == 0)
a10 = all(t10 == seq(2,32,2))
a11 = length(t11)==0
a12 = all(t12 == c(-Inf,seq(2,32,2),Inf))
a.tests = paste("a", seq(1,12), sep="")
a.results= sapply(a.tests, function(x){eval(parse(text=x))})
if(any(!a.results)){
  bad.tests = paste("a", which(!a.results), sep="", collapse=" ")
  stop(paste("get.change.times is returning incorrect times in tests",
             bad.tests))
}


# test when called on a nD with no edge times

# test vertex attributes stuff
net<-network.initialize(5)
activate.vertex.attribute(net,"number","one",onset=0,terminus=2)
activate.vertex.attribute(net,"number","two",onset=3,terminus=4)
activate.vertex.attribute(net,"number","three",onset=5,terminus=6)

expect_equal(get.change.times(net),c(0,2,3,4,5,6))
expect_equal(get.change.times(net,vertex.attribute.activity=FALSE),numeric(0))

# test if only some vertices have attribute
net<-network.initialize(5)
activate.vertex.attribute(net,"number","one",onset=0,terminus=2,v=1:2)
expect_equal(get.change.times(net),c(0,2),info="checking if only some vertices have activity attribute")

# test edge attribute stuff
net<-network.initialize(5)
net[1,2]<-1
activate.edge.attribute(net,"number","one",onset=0,terminus=2)
activate.edge.attribute(net,"number","two",onset=3,terminus=4)
activate.edge.attribute(net,"number","three",onset=5,terminus=6)

expect_equal(get.change.times(net),c(0,2,3,4,5,6),info='check edge attribute activity activity')
expect_equal(get.change.times(net,edge.attribute.activity=FALSE),numeric(0),info='ignore check edge attribute activity')

# test if only some edges have attributes
net<-network.initialize(5)
net[1,2]<-1
net[2,3]<-1
activate.edge.attribute(net,"number","one",onset=0,terminus=2,e=2)
expect_equal(get.change.times(net),c(0,2),info='check if only some edges have activity attribute')

# test with deleted edge
delete.edges(net,eid=1)
expect_equal(get.change.times(net),c(0,2),info='check if some edges deleted')

# test network attribute stuff
net<-network.initialize(5)
activate.network.attribute(net,"number","one",onset=0,terminus=2)
activate.network.attribute(net,"number","two",onset=3,terminus=4)
activate.network.attribute(net,"number","three",onset=5,terminus=6)

expect_equal(get.change.times(net),c(0,2,3,4,5,6),info='check network attribute activity activity')
expect_equal(get.change.times(net,network.attribute.activity=FALSE),numeric(0),info='ignore check network attribute activity')

# check if null spell
net<-network.initialize(2)
deactivate.vertices(net,onset=-Inf,terminus=Inf)
expect_equal(get.change.times(net,ignore.inf=FALSE),numeric(0),info="check for null spell")

expect_equal(get.change.times(network.initialize(0)),numeric(0))

cat("ok\n")


#---------------- GET.EDGEIDS.ACTIVE TESTS ------------------------
# Notes:
#  --this function is basicallly 2 lines, which call functions that
#    are well tested (get.edgeIDs by time, is.active by me), thus
#    the testing here is quite minimal
#-------------------------------------------------------------------

cat("testing get.edgeIDs.active ... ")
anet <- network.initialize(100)
set.seed(10)
heads <- sample(100, 150, replace=TRUE)
set.seed(25)
tails <- sample(100, 150, replace=TRUE)
add.edges(anet, tails, heads)

anet.copy <-anet


# only checking a few simple cases
iov = as.matrix(anet, matrix.type="edgelist")[1:8,]
activate.edges(anet, c(-Inf, -Inf,  10, 10, 10),
                     c( Inf,   20, Inf, 20, 10), e=2:6)
activate.edges(anet, 10, 20, e=7)
activate.edges(anet, 20, 30, e=7)
activate.edges(anet, 30, 30, e=7)
activate.edges(anet, 40, 50, e=7)
deactivate.edges(anet, e=8)
spells = lapply(lapply(anet$mel[1:8],"[[","atl"),"[[","active")

# point queries that should return an empty vectors
b1 = length(get.edgeIDs.active(anet, iov[4,1], at=-Inf, active.default=F))==0    # at -Inf
b2 = length(get.edgeIDs.active(anet, iov[5,1], at=Inf, active.default=F))==0     # at Inf
b3 = length(get.edgeIDs.active(anet, iov[5,1], at=20, active.default=F))==0      # at b
b4 = length(get.edgeIDs.active(anet, iov[3,1], at=30, active.default=F))==0      # at H

# point queries that should return non-empty vectors
b5 = get.edgeIDs.active(anet, iov[3,1], at=-Inf, active.default=F)==3  # at -Inf
b6 = get.edgeIDs.active(anet, iov[4,1], at=Inf, active.default=F)==4   # at Inf
b7 = get.edgeIDs.active(anet, iov[5,1], at=10, active.default=F)==5    # at a
b8 = get.edgeIDs.active(anet, iov[7,1], at=30, active.default=F)==7    # at b

# interval queries that should return empty vectors
b9 = length(get.edgeIDs.active(anet, iov[3,1], Inf, Inf, active.default=F))==0   # over (Inf,Inf)
b10 = length(get.edgeIDs.active(anet, iov[8,1], -Inf, Inf, active.default=F))==0  # over (-Inf,Inf)
b11 = length(get.edgeIDs.active(anet, iov[3,1], 30, Inf, active.default=F))==0    # over (H,Inf)
b12 = length(get.edgeIDs.active(anet, iov[4,1], -Inf, 0, active.default=F))==0    # over (-Inf, L)
b13 = length(get.edgeIDs.active(anet, iov[5,1], 0, 5, active.default=F))==0       # over (L1,L2)
b14 = length(get.edgeIDs.active(anet, iov[5,1], 0, 15, active.default=F, rule="all"))==0   # over (L,M)  
b15 = length(get.edgeIDs.active(anet, iov[7,1], 0, 60, active.default=F, rule="all"))==0   # over (L,H)
b16 = length(get.edgeIDs.active(anet, iov[7,1], 45, 55, active.default=F, rule="all"))==0  # over (M, H)
b17 = length(get.edgeIDs.active(anet, iov[6,1], 15, 20, active.default=F))==0     # over (H1, H2)      

# interval queries that should return non-empty vectors
b18 = min(get.edgeIDs.active(anet, iov[1,1], -Inf, Inf, active.default=T))==1  # over null
b19 = get.edgeIDs.active(anet, iov[2,1], -Inf, Inf, active.default=F)==2  # over (-Inf,Inf)
b20 = get.edgeIDs.active(anet, iov[4,1], 10, Inf, active.default=F)==4    # over (a,Inf)
b21 = get.edgeIDs.active(anet, iov[3,1], -Inf, 20, active.default=F)==3   # over (-Inf, b)
b22 = get.edgeIDs.active(anet, iov[5,1], 10, 20, active.default=F)==5     # over (a,b)
b23 = get.edgeIDs.active(anet, iov[6,1], 10, 10, active.default=F)==6     # over (a,a)
b24 = get.edgeIDs.active(anet, iov[4,1], 0, 15, active.default=F)==4      # over (L,M)  
b25 = get.edgeIDs.active(anet, iov[7,1], 0, 60, active.default=F)==7      # over (L,H)
b26 = get.edgeIDs.active(anet, iov[5,1], 15, 18, active.default=F)==5     # over (M1, M2)
b27 = get.edgeIDs.active(anet, iov[7,1], 45, 60, active.default=F)==7     # over (M, H)

b.tests = paste("b", seq(1,27), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  stop(paste("get.edgeIDs.active is returning incorrect edge IDs in tests",
             bad.tests))
}


# some other tests that Zack has written                   
# Test case check that given onset all get.edges.active reproduces get.edges.active
a1 <-activate.edges(anet,onset=0,terminus=3)
v<-1:100
c1=all(sapply(v,function(x){identical(get.edgeIDs.active(a1,x,at=1),get.edgeIDs(a1,x))}))
# same as above, but check alter functionality, alters present
c2=all(sapply(v,function(x)
	{
	all(sapply(which(a1[1,]==1),function(y){identical(get.edgeIDs.active(a1,x,at=1,alter=y),get.edgeIDs(a1,x,alter=y))}))
	}))
# alters not present
c3=all(sapply(v,function(x)
	{
	all(sapply(which(a1[1,]!=1),function(y){identical(get.edgeIDs.active(a1,x,at=1,alter=y),get.edgeIDs(a1,x,alter=y))}))
	}))
# Check that only empty cases are identical 
c4=all(sapply(v,function(x){identical(get.edgeIDs.active(a1,x,at=-1),get.edgeIDs(a1,x))})
  ==sapply(v,function(x){length(get.edgeIDs(a1,x))==0}))

c.tests = paste("c", seq(1,4), sep="")
c.results= sapply(c.tests, function(x){eval(parse(text=x))})
if(any(!c.results)){
  bad.tests = paste("c", which(!b.results), sep="", collapse=" ")
  stop(paste("get.edgeIDs.active is returning incorrect edge IDs in tests",
             bad.tests))
}




#---------------- GET.EDGES.ACTIVE TESTS ------------------------
# Notes:
#  --this function is nearly identical to get.edgeIDs.active,
#    so we do even less testing here.
#-------------------------------------------------------------------

cat("testing get.edges.active ... ")

## Zack's tests
# within interval
b1=all(sapply(v,function(x){identical(get.edges(a1,x),get.edges.active(a1,v=x,at=1))}))
# outside interval, should only be same on empty cases
b2=all(sapply(v,function(x){identical(get.edges(a1,x),get.edges.active(a1,v=x,at=-1))})
  ==
  sapply(v,function(x){length(get.edges(a1,x))==0}))
## skye's tests
net <-network.initialize(5)
net[1,2]<-1;
net[2,3]<-1;
net <-activate.edges(net,onset=1,terminus=Inf,e=1)
net <-activate.edges(net,onset=2,terminus=3,e=2)
b3=length(get.edges.active(net,v=2,at=1))== 0
b4=length(get.edges.active(net,v=2,at=2, neighborhood="combined"))== 2

b.tests = paste("b", seq(1,4), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  stop(paste("get.edges.active is returning incorrect edges in tests",
             bad.tests))
}

get.edges.active(network.initialize(0),v=1,at=1)

cat("ok\n")



#---------------- GET.NEIGHBORHOOD.ACTIVE TESTS ------------------------
# Notes:
#  --this function largely relies on 'get.edges.active', so rather than
#    test an extensive set of inputs, I'm white-box testing each of
#    the branches in the main "if" loop
#-------------------------------------------------------------------

cat("testing get.neighborhood.active ... ")

# branch 1: undirected network
data(flo)
net1 <- network(flo, directed=FALSE)
el1 <- as.matrix(net1, matrix.type="edgelist")
activate.edges(net1, 10, 20, e=seq(2,20,2))
activate.edges(net1, 10, 10, e=seq(1,19,2))
n1 = get.neighborhood.active(net1, v=9, 10, 20, rule="all")
n2 = get.neighborhood.active(net1, v=9, at=10)
b1 = all(n1==c(2,3,13,16))
b2 = all(n2==c(1,2,3,13,14,16))

# branch 2: undirected network with loops
#  um, I can't get network to accept a network with self loops

# branch 3:  undirected network, no neighborhood
n3 = get.neighborhood.active(net1, v=12, at=10)
b3 = length(n3) == 0

# branch 4:  directed network, out neighborhood
net2 <- network(flo)
el2 <- as.matrix(net2, matrix.type="edgelist")
activate.edges(net2, -Inf, 20, e=seq(2,40,2))
activate.edges(net2, 20, 20, e=seq(1,39,2))
n4 = get.neighborhood.active(net2, v=15, 10, 20, type="out")
n5 = get.neighborhood.active(net2, v=15, at=20, type="out")
b4 = all(n4==c(5,11,13))
b5 = (n5==4)

# branch 5:  directed network, no out neighborhood
n6 = get.neighborhood.active(net2, v=1, 10, 20, type="out")
n7 = get.neighborhood.active(net2, v=15, 30, 40, type="out")
b6 = length(n6)==0
b7 = length(n7)==0

# branch 6:  directed network, in neighborhood
net3 <- network(flo)
activate.edges(net3, 10, Inf, e=seq(2,40,2))
activate.edges(net3, 10, 20, e=seq(1,39,2))
n8 = get.neighborhood.active(net3, v=9, 10, 20, type="in")
n9 = get.neighborhood.active(net3, v=9, 10, 30, type="in", rule="all")
b8 = all(n8==c(1,2,3,13,14,16))
b9 = all(n9==c(2,13,16))

# branch 7:  directed network, no in neighborhood
n10 = get.neighborhood.active(net3, v=12, 10, 20, type="in")
n11 = get.neighborhood.active(net3, v=9, -10, -5, type="in")
b10 = length(n10)==0
b11 = length(n11)==0

# branch 8:  directed network, combined neighborhood, having in and out ties
net4 <- network(flo)
activate.edges(net4, at=10, e=seq(2,40,2))
activate.edges(net4, at=20, e=seq(1,39,2))
n12 = get.neighborhood.active(net4, v=7, 10, 30, type="combined")
n13 = get.neighborhood.active(net4, v=7, at=20, type="combined")
b12 = all(n12==c(2,4,8,16))
b13 = all(n13==c(2,4,16))

# branch 9:  directed network, combined neighborhood, having only in ties
net5 <- network(flo)
net5[15,] = 0   # remove all out edges from node 15
el5 <- as.matrix(net5, matrix.type="edgelist")
activate.edges(net5, e=seq(2,36,2))
n14 = get.neighborhood.active(net5, v=15, -Inf, Inf, type="combined")
n15 = get.neighborhood.active(net5, v=15, -Inf, Inf, type="combined", active.default=F)
b14 = all(n14==c(4,5,11,13))
b15 = all(n15==c(4,11))

# branch 10:  directed network, combined neighborhood, having only out ties
net6 <- network(flo)
net6[,14] = 0   # remove all in edges to node 14
el6 <- as.matrix(net6, matrix.type="edgelist")
activate.edges(net6, e=23)
deactivate.edges(net6, e=25)
n16 = get.neighborhood.active(net6, v=14, at=6, type="combined")
b16 = (n16==9)

# branch 11:  directed network, combined neighborhood, but no in or out ties
n17 = get.neighborhood.active(net5, v=12, at=6, type="combined")
b17 = length(n17)==0

b.tests = paste("b", seq(1,17), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  stop(paste("get.neighborhood.active is returning wrong neighborhood in tests",
             bad.tests))
}


# Zack's tests (adjusted a touch)
# within interval
net7 <- network(flo)
activate.edges(net7, 1,2)
v = 1:network.size(net7)
c1=all(sapply(v,function(x){identical(sort(get.neighborhood(net7,x)),get.neighborhood.active(net7,v=x,at=1))}))
# outside interval
c2=all(sapply(v,function(x){length(get.neighborhood.active(net7,v=x,at=-1))==0}))
#c2=all(sapply(v,function(x){
#              nc<-get.neighborhood(a1,x)
#              nd<-get.neighborhood.active(a1,v=x,at=-1)
#              identical(nc[order(nc)],nd[order(nd)])})
#        ==
#       sapply(v,function(x){length(get.edges(a1,x))==0}))

c.tests = paste("c", seq(1,2), sep="")
c.results= sapply(c.tests, function(x){eval(parse(text=x))})
if(any(!c.results)){
  bad.tests = paste("c", which(!c.results), sep="", collapse=" ")
  stop(paste("get.neighborhood.active is returning wrong neighborhood in tests",
             bad.tests))
}

expect_equal(get.neighborhood.active(network.initialize(0),v=1),integer(0))

cat("ok\n")

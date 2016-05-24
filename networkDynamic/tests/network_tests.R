########################################################################
# This file contains the testing suite for the "network" methods, i.e.:
#
#   network.extract              network.edgecount.active
#   network.dynamic.check        network.naedgecount.active
#   network.dyadcount.active     network.size.active   
#########################################################################

require(networkDynamic)
require(testthat) # helper package for unit tests


#------------------------ NETWORK.EXTRACT ----------------------------
# Notes:
#  - I've tested combinations of retain.all.vertices with bipartite
#    and non-bipartite networks, directed and undirected networks,
#    networks w/ and w/o loops or multi-edges.  I have not tested
#    what would happen if users added onto the basic network structure,
#    or what happens if edges/nodes are missing.
#--------------------------------------------------------------------

cat("testing network.extract ... ")

# a network with no active anything and no edges
anet0 <- network.initialize(7, directed=FALSE)
bnet0=network.extract(anet0, at=10)
b0 = (bnet0%n%"n" == 7)
# a network with no active anything
anet1 <- network.initialize(7, directed=FALSE)
add.edges(anet1, c(1,2,3,4,4,4),
                 c(2,4,4,6,7,5))
anet.copy <- anet1
deactivate.edges(anet1)
deactivate.vertices(anet1)
bnet1=network.extract(anet1, at=10)
bnet1.t=anet1%t%10
b1 = all.equal(bnet1, as.networkDynamic(network.initialize(0,directed=FALSE)))
b2 = identical(bnet1, bnet1.t)
# a network with no active edges, but active vertices, retain.all=F
anet2 <- anet.copy
activate.vertices(anet2, 10, 20, v=c(1,3,4))
deactivate.vertices(anet2, v=c(2,4,5,7))
deactivate.edges(anet2)
bnet2=network.extract(anet2, at=10)
bnet2.t=anet2%t%10
b2 = (bnet2%n%"n" == 3)
b3 = (bnet2%n%"mnext" == 1)
b4 = identical(bnet2, bnet2.t)
# a network with no active edges, but active vertices, retain.all=T
anet3 <- anet2
bnet3=network.extract(anet3, at=10, retain.all.vertices=TRUE)
b5 = (bnet3%n%"n" == 7)
b6 = (bnet3%n%"mnext" == 1)
# a network with no active vertices, but active edges
anet4 <- anet.copy
activate.vertices(anet4, 10, Inf)
activate.edges(anet4, -Inf, 20)
bnet4=network.extract(anet4, -10, 0)
b7 = all.equal(bnet4, as.networkDynamic(network.initialize(0,directed=FALSE)))
# networks with active edges and vertices
#  a non-bipartite simple network, retain.all=F
anet5 <- anet.copy
activate.vertices(anet5, 10, 40, v=c(1,3,4,5,7))
deactivate.vertices(anet5, v=c(2,6))
activate.edges(anet5, at=30, e=c(1,3,4,6))
deactivate.edges(anet5, e=c(2,5))
bnet5=network.extract(anet5, at=30)
bnet5.t=anet5%t%30
el5 = as.matrix(bnet5, matrix.type="edgelist")
b8 = (bnet5%n%"n" == 5)
b9 = (bnet5%n%"mnext" == 3)
b10 = all(el5==matrix(c(2,3,3,4),2,2))
b11 = identical(bnet5, bnet5.t)
#  a non-bipartite simple network, retain.all=T
anet6 <- anet5
bnet6=network.extract(anet6, at=30, retain.all=TRUE)
el6 = as.matrix(bnet6, matrix.type="edgelist")
b12 = (bnet6%n%"n" == 7)
b13 = (bnet6%n%"mnext" == 3)
b14 = all(el6==matrix(c(3,4,4,5),2,2))
# bipartite simple network, retain.all = F
anet7 <- network.initialize(10, bipartite=6)
add.edges(anet7, c(1,2,3,4,4,6), c(7,8,7,8,9,10))
activate.vertices(anet7, -Inf, Inf, v=c(1,3,6,7,9,10))
activate.edges(anet7, 10, Inf, e=c(1,3,5))
bnet7=network.extract(anet7, 0, 20, active.default=FALSE)
el7 = as.matrix(bnet7, matrix.type="edgelist")
b15 = (bnet7%n%"n" == 6)
b16 = (bnet7%n%"mnext" == 3)
b17 = (bnet7%n%"bipartite"== 3)
b18 = all(el7==matrix(c(1,2,4,4),2,2))
# bipartite simple network, retain.all = T
anet8 <- anet7
bnet8=network.extract(anet8, 0, 20, active.default=FALSE, retain.all.vertices=T)
el8 = as.matrix(bnet8, matrix.type="edgelist")
b19 = (bnet8%n%"n" == 10)
b20 = (bnet8%n%"mnext" == 3)
b21 = (bnet8%n%"bipartite"== 6)
b22 = all(el8==matrix(c(1,3,7,7),2,2))
# a non-bipartite network, with loops, retain.all=F
anet9 <- anet.copy
add.edges(anet9, c(2,5,6), c(2,5,6))
activate.vertices(anet9, v=c(1,2,4,6))
deactivate.vertices(anet9, v=c(3,5,7))
activate.edges(anet9, e=c(1,2,6,7,8))
deactivate.edges(anet9, e=c(3,4,5,9))
bnet9=network.extract(anet9, at=-Inf)
bnet9.t=anet9%t%-Inf
el9 = as.matrix(bnet9, matrix.type="edgelist")
b23 = (bnet9%n%"n" == 4)
b24 = (bnet9%n%"mnext" == 4)
b25 = all(el9==matrix(c(1,2,2,2,3,2),3,2))
b26 = identical(bnet9, bnet9.t)
# a non-bipartite network, with loops, retain.all=F
anet10 <- anet9
bnet10=network.extract(anet10, at=Inf, retain.all.vertices=T)
el10 = as.matrix(bnet10, matrix.type="edgelist")
b27 = (bnet10%n%"n" == 7)
b28 = (bnet10%n%"mnext" == 4)
b29 = all(el10==matrix(c(1,2,2,2,4,2),3,2))
# a non-bipartite network, with multi-edges, retain.all=F
anet11 <- network.initialize(8)
add.edges(anet11, c(1,1,1,1,1,4,4,5,5,6,6,8),
                  c(3,3,3,6,2,2,5,4,8,7,7,5))
activate.vertices(anet11, 10,20, v=c(1,3,4,5,6,7))
deactivate.vertices(anet11, v=c(2,8))
activate.edges(anet11, 10, 15, e=c(1:3,7,9,11))
deactivate.edges(anet11, e=c(4,5,6,8,10,12))
bnet11=network.extract(anet11, 11,14)
el11 = as.matrix(bnet11, matrix.type="edgelist")
b30 = (bnet11%n%"n" == 6)
b31 = (bnet11%n%"mnext" == 6)
b32 = all(el11==matrix(c(1,1,1,3,5,2,2,2,4,6),5,2))
# a non-bipartite network, with multi-edges, retain.all=T
anet12 <- anet11
bnet12=network.extract(anet12, 14,19, retain.all.vertices=T)
el12 = as.matrix(bnet12, matrix.type="edgelist")
b33 = (bnet12%n%"n" == 8)
b34 = (bnet12%n%"mnext" == 6)
b35 = all(el12==matrix(c(1,1,1,4,6,3,3,3,5,7),5,2))
# pavel's test with nullified edges
data(flo)
anet13 <- network(flo)
anet13[15,] = 0   # remove all out edges from node 15
el13 <- as.matrix(anet13, matrix.type="edgelist")
deactivate.edges(anet13)
activate.edges(anet13, e=seq(2,40,2))
bnet13 <- anet13 %t% 5
el13.b1 <- el13[c(2,4,6,8,9,seq(12,24,2),seq(28,36,2)),]
el13.b2 <- as.matrix(bnet13, matrix.type="edgelist")
b36 = (bnet13%n%"n" == 16)
b37 = all(el13.b1==el13.b2)


b.tests = paste("b", seq(0,37), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  stop(paste("network.extract is incorrectly extracting networks in tests",
             bad.tests))
}

nd<-network.initialize(7)
add.edges.active(nd,onset=0,terminus=5,tail=1,head=2)
add.edges.active(nd,onset=4,terminus=6,tail=2,head=3)
add.edges.active(nd,onset=0,terminus=5,tail=3,head=4)
add.edges.active(nd,onset=4,terminus=8,tail=4,head=5)
activate.edge.attribute(nd,"letters","a",onset=0,terminus=5)
activate.vertex.attribute(nd,"numbers","one",onset=0,terminus=5)
activate.network.attribute(nd,"squashes","pumpkin",onset=0,terminus=5)
nd1<-network.extract(nd,onset=1,terminus=3,trim.spells=TRUE)
if(!all(get.change.times(nd1) == c(1,3))){
  stop("network.extract did not trim edge activity spells to query range as expected")
}
if(!all(get.edge.value.active(nd1,"letters",onset=-Inf,terminus=Inf,return.tea=TRUE)[[1]][[2]]==c(1,3))){
  stop("network.extract did not trim edge attribute activity spells to query range as expected")
}
if(!all(get.vertex.attribute.active(nd1,"numbers",onset=-Inf,terminus=Inf,return.tea=TRUE)[[1]][[2]]==c(1,3))){
  stop("network.extract did not trim vertex attribute activity spells to query range as expected")
}
if(!all(get.network.attribute.active(nd1,"squashes",onset=-Inf,terminus=Inf,return.tea=TRUE)[[2]]==c(1,3))){
  stop("network.extract did not trim network attribute activity spells to query range as expected")
}

# test trim with "at" version
nd3<-network.extract(nd,at=3,trim.spells=TRUE)
if (!all(nd3$mel[[1]]$atl$letters.active[[2]] == c(3,3), nd3$mel[[1]]$atl$letters.active[[1]]=='a')){
  stop("network.extract did not trim edge attribute spells to 'at' query range as expected")
}

if (!all(nd3$val[[1]]$numbers.active[[2]] == c(3,3),nd3$val[[1]]$numbers.active[[1]]=='one')){
  stop("network.extract did not trim vertex attribute spells to 'at' query range as expected")
}

if (!all(nd3$gal$squashes.active[[2]] == c(3,3),nd3$gal$squashes.active[[1]]=='pumpkin')){
  stop("network.extract did not trim network attribute spells to 'at' query range as expected")
}

if (!all(nd3$mel[[1]]$atl$active==c(3,3))){
  stop("network.extract did not trim edge activity spells to 'at' query range as expecrted")
}

if (!all(nd3$val[[1]]$active==c(3,3))){
  stop("network.extract did not trim vertex activity spells to 'at' query range as expecrted")
}

# check that vertices included by retain.all=TRUE are marked inactive
net<-network.initialize(3)
activate.vertices(net,onset=0,terminus=2,v=2:3)
net1<-network.extract(net,at=1,active.default=FALSE,retain.all.vertices=TRUE,trim.spells=TRUE)
expect_true(is.null(get.vertex.activity(net1)[[1]]),info="check that vertices included by retain.all=TRUE are marked inactive by trim.spells")
expect_equal(unlist(get.vertex.activity(net1)),c(1,1,1,1),info="check that vertices included by retain.all=TRUE are marked inactive by trim.spells")

# check if edges to inactive vertices are retained when retain.all
net<-network.initialize(3)
activate.vertices(net,onset=0,terminus=2,v=2:3)
net[1,2]<-1
expect_equal(network.edgecount(network.extract(net,at=1,retain.all.vertices=TRUE,active.default=FALSE)),0,info="check if edges to inactive vertices are retained when retain.all.vertices=TRUE")

test_that("net.obs.period updated appropriately",{
  net<-network.initialize(3)
  activate.vertices(net,onset=0,terminus=31)
  net%n%'net.obs.period'<-list(observations=list(c(-3,-1),c(0,24),c(25,26),c(27,31),c(32,32)),mode="discrete", time.increment=1,time.unit="day")
  netex<-network.extract(net,onset=20,terminus=30)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(unlist(unlist(net.obs$observations)),c(20,24,25,26,27,30))
  
  netex<-network.extract(net,onset=25,terminus=26)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(unlist(net.obs$observations),c(25,26))
  
  netex<-network.extract(net,onset=25.5,terminus=25.5)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(unlist(net.obs$observations),c(25.5,25.5))
  
  netex<-network.extract(net,onset=26.5,terminus=26.5)  # query inside gap in range
  net.obs<-netex%n%'net.obs.period'
  #expect_equal(unlist(net.obs$observations),c(26.5,26.5))
  expect_equal(unlist(net.obs$observations),c(Inf,Inf))
  
  netex<-network.extract(net,onset=25,terminus=26.5)  # query overlap gap in range
  net.obs<-netex%n%'net.obs.period'
  expect_equal(unlist(net.obs$observations),c(25,26))
  
  netex<-network.extract(net,onset=32,terminus=32)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(unlist(net.obs$observations),c(32,32))

  
  netex<-network.extract(net,onset=32,length=1)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(unlist(net.obs$observations),c(32,32))
  
  netex<-network.extract(net,onset=33,terminus=34)
  net.obs<-netex%n%'net.obs.period'
  #expect_equal(unlist(net.obs$observations),c(33,34)) # queries outside of range construct new range
  expect_equal(unlist(net.obs$observations),c(Inf,Inf))
  
  netex<-network.extract(net,onset=27,length=1)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(unlist(net.obs$observations),c(27,28))
  
  netex<-network.extract(net,terminus=28,length=1)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(unlist(net.obs$observations),c(27,28))
  
  netex<-network.extract(net,at=26)  # at query at open end of interval
  net.obs<-netex%n%'net.obs.period'
  #expect_equal(unlist(net.obs$observations),c(26,26))
  expect_equal(unlist(net.obs$observations),c(Inf,Inf))
  
  netex<-network.extract(net,onset=-4,terminus=0)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(unlist(net.obs$observations),c(-3,-1))
  
  netex<-network.extract(net,onset=-4,terminus=1)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(unlist(net.obs$observations),c(-3,-1,0,1))
  
  # make sure the time range is not expanded
  netex<-network.extract(net,onset=-Inf,terminus=Inf)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(net.obs$observations,list(c(-3,-1),c(0,24),c(25,26),c(27,31),c(32,32)))
  
  netex<-network.extract(net,onset=-10,terminus=40)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(net.obs$observations,list(c(-3,-1),c(0,24),c(25,26),c(27,31),c(32,32)))
  
  # test with at spell at begining of spell #575
  test<-network.initialize(1)
  test%n%'net.obs.period'<-list(observations=list(c(0,1),c(2,3)),mode="discrete", time.increment=1,time.unit="day")
  netex<-network.extract(test,onset=0,terminus=0)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(net.obs$observations,list(c(0,0)))
  
  netex<-network.extract(test,onset=2,terminus=2)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(net.obs$observations,list(c(2,2)))
  
  
  # test for 'null' observation spell
  test<-network.initialize(1)
  test%n%'net.obs.period'<-list(observations=list(c(Inf,Inf)),mode="discrete", time.increment=1,time.unit="day")
  netex<-network.extract(test,onset=0,terminus=0)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(net.obs$observations,list(c(Inf,Inf)))
  
  test<-network.initialize(1)
  test%n%'net.obs.period'<-list(observations=list(c(1,1)),mode="discrete", time.increment=1,time.unit="day")
  netex<-network.extract(test,onset=1,terminus=1)
  net.obs<-netex%n%'net.obs.period'
  expect_equal(net.obs$observations,list(c(1,1)))
               
})

# check extraction of network size 0
nd<-activate.vertices(network.initialize(5),onset=1,terminus=Inf)
net0<-network.extract(nd,at=0)
expect_true(is.network(net0),info="check network object returned for net of zero vertices")
expect_equal(network.size(net0),0,info="check network object returned for net of zero vertices has 0")

# check ordering of elements after extraction
nd<-network.initialize(5)
network.vertex.names(nd)<-LETTERS[5:1]
n1<-network.extract(nd,at=0)

# does order match if no verts removed?
expect_equal(network.vertex.names(n1),LETTERS[5:1])

deactivate.vertices(nd,onset=-1,terminus=2,v=2:4)
n2<-network.extract(nd,at=0)
expect_equal(network.vertex.names(n2),c("E","A"))

# check that network attributes respected by extraction
net<-network.initialize(3,directed=FALSE,bipartite=0,hyper=TRUE,loops=TRUE,multiple=TRUE)
activate.vertices(net,onset=0,terminus=3)
netout<-network.extract(net,onset=1,terminus=2)
expect_equal(is.bipartite(netout),TRUE)
expect_equal(netout%n%'bipartite',0)
expect_equal(is.multiplex(netout),TRUE)
expect_equal(is.hyper(netout),TRUE)
expect_equal(is.directed(netout),FALSE)



#-------------------- NETWORK.DYNAMIC.CHECK ------------------------
# Notes:
#  - 
#--------------------------------------------------------------------

cat("testing network.dynamic.check ... ")

# a network with no activity specified at all
cnet0 <- network.initialize(5)
add.edges(cnet0, c(1,3,5), c(4,4,2))
check0 = network.dynamic.check(cnet0, verbose=F)
c0 = (sum(!check0$vertex.checks)==0 &&
      sum(!check0$edge.checks)==0)

# a network with complete inactivity specified
cnet1 <- cnet0
deactivate.edges(cnet1)
deactivate.vertices(cnet1)
check1 = network.dynamic.check(cnet1, verbose=F)
c1 = (sum(!check1$vertex.checks)==0 &&
      sum(!check1$edge.checks)==0)

# a network that is A-okay
cnet2 <- cnet1
activate.edges(cnet2, 10,20, e=1:2)
activate.vertices(cnet2, 0,30)
check2 = network.dynamic.check(cnet2, verbose=F)
c2 = (sum(!check2$vertex.checks)==0 &&
      sum(!check2$edge.checks)==0 &&
      sum(!check2$dyad.checks)==0)

# a network with illegal node activity matrices
cnet3 <- cnet2
cnet3$val[[2]]$active <- matrix(c(10,0,20,5),2,2)
cnet3$val[[3]]$active <- matrix(c(10,15,20,25),2,2)
cnet3$val[[4]]$active <- matrix(c(10,40,20,30),2,2)
check3 = network.dynamic.check(cnet3, verbose=F)
c3 = (all(check3$vertex.checks==c(T,F,F,F,T)) &&
      sum(!check3$edge.checks)==0)

# a network with illegal edge activity matrices
cnet4 <- cnet2
cnet4$mel[[1]]$atl$active <- matrix(c(10,0,20,5),2,2)
cnet4$mel[[2]]$atl$active <- matrix(c(10,15,20,25),2,2)
cnet4$mel[[3]]$atl$active <- matrix(c(10,40,20,30),2,2)
check4 = network.dynamic.check(cnet4, verbose=F)
c4 = (sum(!check4$vertex.checks)==0 &&
      sum(check4$edge.checks)==0)

# a network with active edges that have inactive end points
cnet5 <- cnet1
activate.edges(cnet5,10,20)
activate.vertices(cnet5,0,3, v=c(1,3,4))
activate.vertices(cnet5,0,30, v=c(2,5))
check5 = network.dynamic.check(cnet5, verbose=F)
c5 = all(check5$dyad.checks==c(F,F,T))


# a network with multiple problems
cnet6 <- cnet5
cnet6$val[[2]]$active <- matrix(c(10,0,20,5),2,2)
cnet6$mel[[3]]$atl$active <- matrix(c(10,15,20,25),2,2)
check6 = network.dynamic.check(cnet6, verbose=F)
c6 = (check6$vertex.checks[2]==F &&
      check6$edge.checks[3]==F)


# a network with a point-activated edge
cnet7 <- network.initialize(2)
cnet7[1,2]<-1
activate.edges(cnet7,at=1)
check7<-network.dynamic.check(cnet7)
c7 <- (check7$vertex.checks==c(T,T) &&
       check7$edge.checks==T &&
       check7$dyad.checks==T)

c.tests = paste("c", seq(1,7), sep="")
c.results= sapply(c.tests, function(x){eval(parse(text=x))})
if(any(!c.results)){
  bad.tests = paste("c", which(!c.results), sep="", collapse=" ")
  stop(paste("network.dynamic.check is returning incorrect results in tests",
             bad.tests))
}

# a network with no edges
nd <-network.initialize(3)
if(!length(network.dynamic.check(nd,verbose=FALSE)$edge.checks)==0){
  stop("edge check returned wrong value for empty network")
}

# an vertex attribute named active that is not a list
nd <-network.initialize(3)
set.vertex.attribute(nd,"letters.active","a")
if(any(network.dynamic.check(nd,verbose=FALSE)$vertex.tea.checks)){
  stop("network.dynamic.check failed to flag bad vertex TEA attributes")
}

# multiple spells
nd <-network.initialize(3)
activate.vertex.attribute(nd,"letters","a",onset=1,terminus=2)
activate.vertex.attribute(nd,"letters","b",onset=2,terminus=3)
if(!all(network.dynamic.check(nd,verbose=FALSE)$vertex.tea.checks)){
  stop("network.dynamic.check failed on perfectly good vertex TEA attributes")
}

# missing attributes
nd <-network.initialize(3)
activate.vertex.attribute(nd,"letters","a",at=1)
if(!all(network.dynamic.check(nd,verbose=FALSE)$vertex.tea.checks)){
  stop("network.dynamic.check failed on perfectly good vertex TEA attributes")
}

# mangled spell matrix
nd <-network.initialize(3)
activate.vertex.attribute(nd,"letters","a",onset=1,terminus=2)
activate.vertex.attribute(nd,"letters","b",onset=2,terminus=3)
nd$val[[2]]$letters.active[[2]][1,1]<-5
if(network.dynamic.check(nd,verbose=FALSE)$vertex.tea.checks[2]){
  stop("network.dynamic.check failed to flag bad vertex TEA attributes")
}

# test for edges
nd <-network.initialize(3)
add.edges(nd,1,2)
add.edges(nd,2,3)
add.edges(nd,3,1)
set.edge.attribute(nd,"letters.active","a")
if(any(network.dynamic.check(nd,verbose=FALSE)$edge.tea.checks)){
  stop("network.dynamic.check failed to flag bad edge TEA attributes")
}

# test multiple spells
nd <-network.initialize(3)
add.edges(nd,1,2)
add.edges(nd,2,3)
add.edges(nd,3,1)
activate.edge.attribute(nd,"letters","a",onset=1,terminus=2)
activate.edge.attribute(nd,"letters","b",onset=2,terminus=3)
activate.edge.attribute(nd,"letters","b",onset=3,terminus=4)
if(!all(network.dynamic.check(nd,verbose=FALSE)$edge.tea.checks)){
  stop("network.dynamic.check incorrectly flaged good edge TEA attributes")
}

# mangled spell matrix
nd$mel[[2]]$atl$letters.active[[2]][1,1]<-100
if(network.dynamic.check(nd,verbose=FALSE)$edge.tea.checks[2]){
  stop("network.dynamic.check failed to flag bad edge TEA attributes")
}

# check network attrs
nd <-network.initialize(3)
set.network.attribute(nd,"letters.active","a")
if(network.dynamic.check(nd,verbose=FALSE)$network.tea.checks){
  stop("network.dynamic.check failed to flag bad network TEA attributes")
}

# check multiple spls network attrs
nd <-network.initialize(3)
activate.network.attribute(nd,'letters',"a",onset=1,terminus=2)
activate.network.attribute(nd,'letters',"b",onset=2,terminus=3)
activate.network.attribute(nd,'letters',"c",onset=3,terminus=4)
if(!network.dynamic.check(nd,verbose=FALSE)$network.tea.checks){
  stop("network.dynamic.check incorrectly flaged good network TEA attributes")
}

# check malformed spls
nd$gal$letters.active[[2]][1,1]<-Inf
if(network.dynamic.check(nd,verbose=FALSE)$network.tea.checks){
  stop("network.dynamic.check failed to flag bad network TEA attributes")
}

# check net.obs.period
net<-network.initialize(3)
add.edges.active(net,tail=1,head=2,onset=0,terminus=3)
checks<-network.dynamic.check(net,verbose=FALSE)
expect_true(is.null(checks$net.obs.period.check))

# check for bad net obs period
nop <- list(observations=list(c("a",2)), mode="discrete", time.increment=1,time.unit="step")
set.network.attribute(net,'net.obs.period',nop)
checks<-network.dynamic.check(net,verbose=FALSE)
expect_false(checks$net.obs.period.check)

# check for bad net obs period
nop <- list(observations=list(c(1,2),c(1,2,3)), mode="discrete", time.increment=1,time.unit="step")
set.network.attribute(net,'net.obs.period',nop)
checks<-network.dynamic.check(net,verbose=FALSE)
expect_false(checks$net.obs.period.check)

# check for good net obs period
nop <- list(observations=list(c(0,4)), mode="discrete", time.increment=1,time.unit="step")
set.network.attribute(net,'net.obs.period',nop)
checks<-network.dynamic.check(net,verbose=FALSE)
expect_true(checks$net.obs.period.check)

# check for net obs period conflicting with network time range
nop <- list(observations=list(c(1,2)), mode="discrete", time.increment=1,time.unit="step")
set.network.attribute(net,'net.obs.period',nop)
checks<-network.dynamic.check(net,verbose=TRUE)
expect_false(checks$net.obs.period.check)

# check network size 0
expect_equal(network.dynamic.check(network.initialize(0))$vertex.check,logical(0))

cat("ok\n")




#---------------- NETWORK.EDGECOUNT.ACTIVE --------------------------
# Notes:
#  - this function also relies primarily on 'is.active', so the testing
#    done here is minimal
#----------------------------------------------------------------------

cat("testing network.edgecount.active ... ")
data(flo)
net0 <- network.initialize(3)
net1 <- network(flo)
net2 <- net1
for (i in 1:5)
  net1$mel[[i]]$atl$na = TRUE
for (i in 1:8)
  net2$mel[[i]]$atl$na = TRUE
deactivate.edges(net1)
deactivate.edges(net2)
activate.edges(net1, 10,20, e=3:7)
activate.edges(net2, 10,20, e=3:7)
b1=(network.edgecount.active(net1, 5,15)==2)            # a positive count, na.omit=T
b2=(network.edgecount.active(net1, 5,15, na.omit=F)==5) # a positive count, na.omit=F
b3=(network.edgecount.active(net1, 40, Inf)==0)         # a 0 count b/c of non-activity
b4=(network.edgecount.active(net2, 5,15)==0)            # a 0 count b/c of missing edges
b5=(network.edgecount.active(net0, 5,15)==0)            # a 0 count b/c of no edges

b.tests = paste("b", seq(1,5), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  error(paste("network.edgecount.active is incorrectly counting edges in tests",
             bad.tests))
}

# a network with edge activity spells, some of which are missing
net7<-network(flo)
deactivate.edges(net7)
activate.edges(net7, 10,20, e=3:7)
delete.edges(net7,eid=1)
expect_equal(network.edgecount.active(net7,at=10),5)

expect_equal(network.edgecount(network.initialize(0)),0)

cat("ok\n")


#---------------- NETWORK.NAEDGECOUNT.ACTIVE --------------------------
# Notes:
#  - this function also relies primarily on 'is.active', so the testing
#    done here is minimal
#----------------------------------------------------------------------
cat("testing network.naedgecount.active ... ")
data(flo)
# a network with missing edges
net3 <- network(flo)
net4 <- net3
for (i in 1:5)
  net3$mel[[i]]$atl$na = TRUE
deactivate.edges(net3)
activate.edges(net3, 10,20, e=3:7)
# a network without missing edges
deactivate.edges(net4)
activate.edges(net4, 10,20, e=3:7)
# a network without missing edges, but with nullified edges
net5 <- net4
net5$mel[[1]] <- NULL
net5$mel[[2]] <- NULL
# a network with missing edges and with nullified edges
net6 <- net3
net6$mel[[1]] <- NULL
net6$mel[[2]] <- NULL

# a network with edge activity spells, missing edge,  and some deleted edges
net7<-net3
delete.edges(net7,eid=1)
expect_equal(network.naedgecount.active(net7,at=10),3)

b1=(network.naedgecount.active(net3, 5,15)==3)     # a positive count
b2=(network.naedgecount.active(net6, 5,15)==2)     # a positive count
b3=(network.naedgecount.active(net3, 40, Inf)==0)  # a 0 count b/c of non-activity
b4=(network.naedgecount.active(net4, 40, Inf)==0)  # a 0 count b/c of non-activity
b5=(network.naedgecount.active(net5, 40, Inf)==0)  # a 0 count b/c of non-activity
b6=(network.naedgecount.active(net6, 40, Inf)==0)  # a 0 count b/c of non-activity
b7=(network.naedgecount.active(net4, 5,15)==0)     # a 0 count b/c of no missing edges
b8=(network.naedgecount.active(net5, 5,15)==0)     # a 0 count b/c of no missing edges
b9=(network.naedgecount.active(net0, 5,15)==0)     # a 0 count b/c of no edges

b.tests = paste("b", seq(1,9), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  error(paste("network.naedgecount.active is incorrectly counting edges in tests",
             bad.tests))
}

expect_equal(network.naedgecount.active(network.initialize(0)),0)

cat("ok\n")

#-------------------- NETWORK.SIZE.ACTIVE --------------------------
# Notes:
#  - this function contains a single line and hinges on 'is.active'
#    returning the right thing, so in the interest of time, I am 
#    running 2 simple tests.
#--------------------------------------------------------------------
cat("testing network.size.active ... ")
data(flo)
net6 <- network(flo)
deactivate.vertices(net6)
activate.vertices(net6,-Inf,20, v=1:10)
b1=(network.size.active(net6, 5,15)==10)    # a positive count
b2=(network.size.active(net6, 40, Inf)==0)  # a 0 count b/c of non-activity

b.tests = paste("b", seq(1,2), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  error(paste("network.size.active is incorrectly counting nodes in tests",
             bad.tests))
}

expect_equal(network.size.active(network.initialize(0),at=0),0)

cat("ok\n")


#---------------- NETWORK.DYADCOUNT.ACTIVE --------------------------
# Notes:
#  - this function also relies primarily on 'is.active', so again,
#    we'll test using a white box approach.
#----------------------------------------------------------------------

cat("testing network.dyadcount.active ... ")

# a directed network, no missing edges
net7 <- network.initialize(6)
activate.vertices(net7, 10, 20, v=1:4)
activate.vertices(net7, 0, 10, v=4:5)
activate.vertices(net7, 30, 40, v=6)
b1 = (network.dyadcount.active(net7, 50, Inf)==0)  # a 0 count, b/c of nonactivity
b2 = (network.dyadcount.active(net7, at=35)==0)    # a 0 count, b/c of only 1 active node
b3 = (network.dyadcount.active(net7,  0, 10)==2)   # a 2 count
b4 = (network.dyadcount.active(net7, 10, 15)==12)  # a >2 count
# an directed bipartitie network, no missing edges
net8 <- network.initialize(10, bipartite=7)
activate.vertices(net8, 10, 40, v=1:6)
activate.vertices(net8, 0, 10, v=7:8)
activate.vertices(net8, 10, 30, v=9)
activate.vertices(net8, 50, 60, v=10)
b5 = (network.dyadcount.active(net8, 80, 90)==0)   # a 0 count, b/c of nonactivity
b6 = (network.dyadcount.active(net8, 50, 55)==0)   # a 0 count, b/c of only 1 active node
b7 = (network.dyadcount.active(net8, 0, 10)==2)    # a 2 count
b8 = (network.dyadcount.active(net8, 10, 15)==12)  # a >2 count
# an undirected bipartitie network, no missing edges
net9 <- network.initialize(10, bipartite=7, directed=FALSE)
activate.vertices(net9, 10, 40, v=1:6)
activate.vertices(net9, 0, 10, v=7:8)
activate.vertices(net9, 10, 30, v=9)
activate.vertices(net9, 10, 40, v=10)
b9 = (network.dyadcount.active(net9, 40, 50)==0)    # a 0 count
b10 = (network.dyadcount.active(net9, 0, 10)==1)    # a 1 count
b11 = (network.dyadcount.active(net9, 10, 15)==12)  # a >1 count
# an undirected non-bipartitie network, no missing edges
net10 <- network(flo, directed=FALSE)
activate.vertices(net10, 20, Inf, v=1:8)
activate.vertices(net10, 10, 10, v=9:14)
b12 = (network.dyadcount.active(net10, at=0, active.default=F)==0)   # a 0 count
b13 = (network.dyadcount.active(net10, at=11)==1)                    # a 1 count
b14 = (network.dyadcount.active(net10, at=30)==45)                   # a >1 count
# a directed non-bipartite network with missing edges
add.edges(net7, c(4,5,1), c(5,4,4))
for (i in 1:3)
  net7$mel[[i]]$atl$na=TRUE
b15 = (network.dyadcount.active(net7, at=15)==11)            # a >1 count, na.omit=T
b16 = (network.dyadcount.active(net7, at=15, na.omit=F)==12) # a >1 count, na.omit=F
b17 = (network.dyadcount.active(net7, at=5)==0)              # a 0 count b/c of missing edges
# a directed bipartite network with missing edges
add.edges(net8, c(3,8,7), c(9,7,8))
for (i in 1:3)
  net8$mel[[i]]$atl$na=TRUE
b18 = (network.dyadcount.active(net8, 25,30)==11)            # a >1 count, na.omit=T
b19 = (network.dyadcount.active(net8, 25,30, na.omit=F)==12) # a >1 count, na.omit=F
b20 = (network.dyadcount.active(net8, 2,3)==0)               # a 0 count b/c of missing edges
# an undirected bipartitie network with missing edges
add.edges(net9, c(4,5,7), c(10,10,8))
for (i in 1:3)
  net9$mel[[i]]$atl$na=TRUE
b21 = (network.dyadcount.active(net9, 31,32)==4)            # a >1 count, na.omit=T
b22 = (network.dyadcount.active(net9, 31,32, na.omit=F)==6) # a >1 count, na.omit=F
b23 = (network.dyadcount.active(net9, 2,3)==0)              # a 0 count b/c of missing edges
# an undirected non-bipartitie network with missing edges
net11 <- network.initialize(16, directed=FALSE)
activate.vertices(net11, 20, Inf, v=1:8)
activate.vertices(net11, 10, 10, v=9:14)
add.edges(net11, c(15, 9), c(16,11))
activate.edges(net11)
for (i in 1:2)
  net11$mel[[i]]$atl$na=TRUE
b24 = (network.dyadcount.active(net11, at=10, active.default=F)==14)  # a >1 count, na.omit=T
b25 = (network.dyadcount.active(net11, at=10, na.omit=F)==28)         # a >1 count, na.omit=F
b26 = (network.dyadcount.active(net11, -5,0)==0)                      # a 0 count b/c of missing edges

b.tests = paste("b", seq(1,26), sep="")
b.results= sapply(b.tests, function(x){eval(parse(text=x))})
if(any(!b.results)){
  bad.tests = paste("b", which(!b.results), sep="", collapse=" ")
  stop(paste("network.dyadcount.active is incorrectly counting dyads in tests",
             bad.tests))
}

expect_equal(network.dyadcount.active(network.initialize(0)),0)

# check for case when first mode of bipartite = 0
test<-network.initialize(3,bipartite=0)
network.dyadcount.active(test,at=1)

cat("ok\n")

# ------ network.collapse tests----
cat("testing network.collapse ...")
test<-network.initialize(5)
add.edges.active(test, tail=c(1,2,3), head=c(2,3,4),onset=0,terminus=1)
activate.edges(test,onset=3,terminus=5)
activate.edges(test,onset=-2,terminus=-1)

# test edge.count
net <-network.collapse(test)

# should not be nD
if (is.networkDynamic(net)){
  stop("network.collapse did not correctly remove networkDynamic class attributes from result")
}

# activity attributes removed?
if ('active'%in%list.vertex.attributes(net)){
  stop("network.collapse did not correctly remove vertex activity data from result")
}

if ('active'%in%list.edge.attributes(net)){
  stop("network.collapse did not correctly remove edge activity data from result")
}

# count and druration not added
if('activity.count'%in%list.vertex.attributes(net) | 'activity.duration'%in%list.vertex.attributes(net)) {
  stop("network.collpase included activity summary attributes when rm.time.info=TRUE")
}

net <-network.collapse(test,rm.time.info=FALSE)

# count and druration not added
if(!'activity.count'%in%list.vertex.attributes(net) | !'activity.duration'%in%list.vertex.attributes(net)) {
  stop("network.collpase did not include activity summary attributes when rm.time.info=FALSE")
}

if (!all(get.edge.value(net,'activity.count')==c(3,3,3))){
  stop("network.collapse did not calculate edge activity.count correctly")
}
    
if (!all(get.vertex.attribute(net,'activity.count')==c(1,1,1,1,1))){
  stop("network.collapse did not calculate vertex activity.count correctly when no vertex spells defined")
}

if (!all(get.edge.value(net,'activity.duration')==c(4,4,4))){
  stop("network.collapse did not calculate edge activity.duration correctly")
}

if (!all(get.vertex.attribute(net,'activity.duration')==c(Inf,Inf,Inf,Inf,Inf))){
  stop("network.collapse did not calculate vertex activity.duration correctly when to vertex spells defined")
}


# test no edges example
test<-network.initialize(5)
activate.vertices(test,onset=0,terminus=5)
net <-network.collapse(test,rm.time.info=FALSE)

if (!all(get.vertex.attribute(net,'activity.duration')==c(5,5,5,5,5))){
  stop("network.collapse did not calculate vertex activity.duration correctly")
}

# test version where only some edges have activity
test<-network.initialize(5)
activate.vertices(test)
add.edges(test, tail=c(1,2,3), head=c(2,3,4))
activate.edges(test,onset=1,terminus=2,e=2)
net <-network.collapse(test,onset=1,terminus=4,rm.time.info=FALSE)

if(!all(get.edge.value(net,'activity.count')==c(1,1,1),get.edge.value(net,'activity.duration')==c(3,1,3))){
  stop("network.collapse did not calculate edge count and duration correctly when only some edges have activity spells")
}

# test version where only some vertices have activity
test<-network.initialize(5)
activate.vertices(test,onset=0,terminus=5,v=2:3)
net <-network.collapse(test,rm.time.info=FALSE)
if (!all(get.vertex.attribute(net,'activity.duration')==c(Inf,5,5,Inf,Inf))){
  stop("network.collapse did not calculate vertex duration correctly when only some vertices have activity")
}



# test time parameters
test<-network.initialize(5)
add.edges.active(test, tail=c(1,2,3), head=c(2,3,4),onset=0,terminus=1)
activate.edges(test,onset=3,terminus=5)
activate.edges(test,onset=-2,terminus=-1)
activate.edge.attribute(test,'weight',5,onset=3,terminus=4)
activate.edge.attribute(test,'weight',3,onset=4,terminus=5)
activate.vertex.attribute(test,'stuff',3, onset=3,terminus=4)
activate.vertex.attribute(test,'moreliststuff',list(list(x=1,z=2)), onset=3,terminus=4)
activate.network.attribute(test,'morestuff',3, onset=3,terminus=4)
activate.network.attribute(test,'liststuff',list(a=1,b=2), onset=3,terminus=4)

# are attributes flattened correctly?
net3 <-network.collapse(test,onset=3,terminus=4)
expect_equal(get.edge.value(net3,'weight'),c(5,5,5))
expect_equal(get.network.attribute(net3,'morestuff'),3)
expect_equal(get.network.attribute(net3,'liststuff'),list(a=1,b=2))
expect_equal(get.vertex.attribute(net3,'stuff'),c(3,3,3,3,3))
expect_equal(get.vertex.attribute(net3,'moreliststuff',unlist=FALSE)[[1]],list(x=1,z=2))



net4 <-network.collapse(test,onset=4,terminus=5)
if(!all(all(get.edge.value(net3,'weight')==c(5,5,5),get.edge.value(net4,'weight')==c(3,3,3)))){
  stop('network.collapse did not flatten attributes as expected')
}

# if it includes multiple attribute values, should produce warning
expect_that(net <-network.collapse(test,onset=3,terminus=5), 
            gives_warning("Multiple attribute values matched query spell"))

# test 'at' version
net <-network.collapse(test,at=4)
if (!all(get.edge.value(net,'weight')==c(3,3,3))){
  stop("network.collapse did not flatten attributes as expected for 'at' param")
}


# test operator version
net <-test%k%4
if (!all(get.edge.value(net,'weight')==c(3,3,3))){
  stop("network.collapse did not flatten attributes as expected for '%k%' operator")
}

# test for case with null spell Inf,Inf
net <- network.initialize(3)
deactivate.vertices(net,onset=-Inf,terminus=Inf,v=1)
net1<-net%k%1

# test retain.all.vertices
net <- network.initialize(5)
set.vertex.attribute(net,'letters',LETTERS[1:5])
activate.vertices(net,onset=c(0,10,0,0,10),terminus=20)

expect_equal(network.size(network.collapse(net,at=5)),3)
expect_equal(network.size(network.collapse(net,at=5,retain.all.vertices=TRUE)),5)
expect_equal(get.vertex.attribute(network.collapse(net,at=5),'letters'),c('A','C','D'))

# test active default
net<-network.initialize(3)
activate.vertices(net,v=2:3)
expect_equal(network.size(network.collapse(net,at=1,active.default=FALSE)),2,info="check active.default param")


# check extraction of network size 0
nd<-activate.vertices(network.initialize(5),onset=1,terminus=Inf)
net0<-network.collapse(nd,at=0)
expect_true(is.network(net0),info="check network object returned for net of zero vertices")
expect_equal(network.size(net0),0,info="check network object returned for net of zero vertices has 0")

# ----- check that the earliest and latest rules work----
testD<-network.initialize(4)
testD[1,2:4]<-1
activate.vertex.attribute(testD,'color','red',onset=0,terminus=1)
activate.vertex.attribute(testD,'color','green',onset=1,terminus=2)
activate.vertex.attribute(testD,'color','blue',onset=2,terminus=3)
activate.edge.attribute(testD,'letter',"a",onset=0,terminus=1)
activate.edge.attribute(testD,'letter',"b",onset=1,terminus=2)
activate.edge.attribute(testD,'letter',"c",onset=2,terminus=3)

testEarly<-network.collapse(testD,rule='earliest')
expect_equal(testEarly%v%'color',rep('red',4))
expect_equal(testEarly%e%'letter',rep('a',3))

testLate<-network.collapse(testD,rule='latest')
expect_equal(testLate%v%'color',rep('blue',4))
expect_equal(testLate%e%'letter',rep('c',3))


test_that("net.obs.period is appropriately removed or retained",{
  nD<-activate.vertices(network.initialize(3))
  obs<-list(observations=list(c(0,100)),mode="discrete", time.increment=1,time.unit="step")
  nD%n%'net.obs.period'<-obs
  n<-network.collapse(nD)
  expect_equal(n%n%'net.obs.period',NULL)
  n<-network.collapse(nD,rm.time.info=FALSE)
  expect_equal(n%n%'net.obs.period',obs)
})

# test for form of dynamic network attribute, bug #1184
net <- network.initialize(5)
x <- matrix(0,2,2) # a matrix
activate.network.attribute(net, "dyadic", x, onset=1, length=1)

# query it and get back a marix
activeQ<-get.network.attribute.active(net, "dyadic", at=1) # a matrix
net1 <- (net %k% 1) # Collapse at t=1.
expect_equal(net1 %n% "dyadic",activeQ) # a list containing a matrix

net <- network.initialize(5)
x <- matrix(0,2,2) # a matrix
y <- matrix(1,2,2) # a matrix
activate.network.attribute(net, "dyadic", x, onset=1, length=1)
#activate.network.attribute(net, "dyadic2", y, onset=1, length=1)

# query it and get back a marix
activeQ<-get.network.attribute.active(net, "dyadic", at=1) # a matrix
net1 <- (net %k% 1) # Collapse at t=1.
expect_equal(net1 %n% "dyadic",activeQ) # a list containing a matrix

#----- get.networks -------

#test for get.networks

test <- network.initialize(5)
add.edges.active(test, tail=c(1,2,3), head=c(2,3,4),onset=0,terminus=1)
activate.edges(test,onset=3,terminus=5)
activate.edges(test,onset=-2,terminus=-1)
activate.edge.attribute(test,'weight',5,onset=3,terminus=4)
activate.edge.attribute(test,'weight',3,onset=4,terminus=5,e=1:2)

# uset start and end args
netlist <- get.networks(test,start=0, end=5)
expect_equal(length(netlist),5,info='get.networks tets')
expect_equal(network.edgecount(netlist[[1]]),3,info='get.networks tets')
expect_equal(network.edgecount(netlist[[2]]),0,info='get.networks tets')
expect_equal(network.edgecount(netlist[[4]]),3,info='get.networks tets')
expect_equal(network.edgecount(netlist[[5]]),3,info='get.networks tets')

# use start end and increment
expect_warning(netlist <- get.networks(test,start=0, end=5,time.increment=5),'Multiple attribute values matched')
expect_equal(length(netlist),1)
# warning-free extraction using 'latest' rule
netlist <- get.networks(test,start=0, end=5,time.increment=5,rule='latest')
expect_equal(netlist[[1]]%e%'weight',c(3,3,5))


# test error cases
expect_error(get.networks(test,start=1),"Unable to infer appropriate onsets and term")
expect_error(get.networks(test,end=1),"Unable to infer appropriate onsets and term")
expect_error(get.networks(test,start=0,end=1,onsets=1),"onsets & termini cannot be specified with start & end arguments")
expect_error(get.networks(test,onsets=1,termini=1:5),"onsets and termini must have the same number of elements")
expect_error(get.networks(test,onsets=1:3,termini=3:1),"Onset times must precede terminus times")
expect_error(get.networks(test,start=-Inf,end=1),"start and end values must be finite")

#use onset and terminus
netlist<-get.networks(test,onsets=c(0,1,2),termini=c(1,2,3))
expect_equal(length(netlist),3)

# use net obs period to infer params
test%n%'net.obs.period'<-list(observations=list(c(-1,5)),mode='discrete',time.increment=1,time.unit='step')
netlist<-get.networks(test)
expect_equal(length(netlist),6)

# params should override net.obs.period
netlist<-get.networks(test,onsets=c(0,1,2),termini=c(1,2,3))
expect_equal(length(netlist),3)

# check network size changes and argument passing
test<-network.initialize(5)
activate.vertices(test,onset=1:5,terminus=c(5,5,5,5,5))
netlist<-get.networks(test,start=1,end=5)
expect_equal(sapply(netlist,network.size),1:4)

netlist<-get.networks(test,start=1,end=5,retain.all.vertices=TRUE)
expect_equal(sapply(netlist,network.size),c(5,5,5,5))


cat("ok\n")
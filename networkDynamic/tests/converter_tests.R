#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2013 the statnet development team
######################################################################

#test networkDynamic conversion functionality

require(networkDynamic)
require(testthat)


# ------ get.vertex.activity test ------

net<-network.initialize(5)
activate.vertices(net,onset=c(1,2,3,-Inf,100),terminus=c(10,20,Inf,30,Inf))
activate.vertices(net,at=1000)
expect_equal(unlist(get.vertex.activity(net)),c(1,1000,10,1000,2,1000,20,1000,3,Inf,-Inf,1000,30,1000,100,Inf),info='check if spell list has correct elements in get.vertex.activity')

expect_equivalent(get.vertex.activity(net,as.spellList=TRUE),data.frame(onset=c(1,1000,2,1000,3,-Inf,1000,100),terminus=c(10,1000,20,1000,Inf,30,1000,Inf),vertex.id=c(1,1,2,2,3,4,4,5),onset.censored=c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE),terminus.censored=c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),duration=c(9,0,18,0,Inf,Inf,0,Inf)),info="check if get.vertex.activity(spellList) returns correct data")

# check sort order
expect_equal(get.vertex.activity(net,as.spellList=TRUE)$vertex.id,c(1,1,2,2,3,4,4,5),info="check sort order of get.vertex.activity vertex.id")

expect_equal(get.vertex.activity(net,as.spellList=TRUE)$onset,c(1,1000,2,1000,3,-Inf,1000,100),info="check sort order of get.vertex.activity onset")
expect_equal(get.vertex.activity(net,as.spellList=TRUE)$terminus,c(10,1000,20,1000,Inf,30,1000,Inf),info="check sort order of get.vertex.activity terminus")

expect_equal(get.vertex.activity(network.initialize(0)),list())




# check that net.obs.period triggers censoring
nop <- list(observations=list(c(-5,1001)), mode="discrete", time.increment=1,time.unit="step")
set.network.attribute(net,'net.obs.period',nop)
expect_equivalent(get.vertex.activity(net,as.spellList=TRUE),data.frame(onset=c(1,1000,2,1000,3,-5,1000,100),terminus=c(10,1000,20,1000,1001,30,1000,1001),vertex.id=c(1,1,2,2,3,4,4,5),onset.censored=c(FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE),terminus.censored=c(FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE),duration=c(9,0,18,0,998,35,0,901)),info="check that net.obs.period triggers censoring in get.vertex.activity(spellList)")


# check for bug with no activity #192
nd<-network.initialize(2)
nd<-add.edges.active(nd,at=1,tail=1,head=2)

expect_equivalent(get.vertex.activity(nd,as.spellList=TRUE),data.frame(onset=c(-Inf,-Inf),terminus=c(Inf,Inf),vertex.id=c(1,2),onset.censored=c(TRUE,TRUE),terminus.censored=c(TRUE,TRUE),duration=c(Inf,Inf)),info="checking with no vertex activity")

# test if activity is missing for an element
activate.vertices(nd,onset=1,terminus=2,v=2)
expect_equivalent(get.vertex.activity(nd,as.spellList=TRUE),data.frame(onset=c(-Inf,1),terminus=c(Inf,2),vertex.id=c(1,2),onset.censored=c(TRUE,FALSE),terminus.censored=c(TRUE,FALSE),duration=c(Inf,1)),info="test when only some verticies have activity")

# test when default activity set false
expect_equivalent(
  get.vertex.activity(nd,as.spellList=TRUE,active.default=FALSE),
  data.frame(onset=1,terminus=2,vertex.id=2,onset.censored=FALSE,terminus.censored=FALSE,duration=1),info="test when only some verticies have activity and active.default=FALSE")

# test if v specified for spellmatrix
nd <- network.initialize(5)
activate.vertices(nd,onset=0:4,terminus=1:5)
expect_equivalent(
 get.vertex.activity(nd,as.spellList=TRUE,v=2:4)[,1:3],
 data.frame(onset=1:3,terminus=2:4,vertex.id=2:4),info="test when when v argument specified with spellmatrix")

expect_equivalent(
  unlist(get.vertex.activity(nd,,v=2:4)),c(1,2,2,3,3,4)
  ,info="test when when v argument specified")


# test if called with ordinary network
expect_equal(unlist(get.vertex.activity(network.initialize(3))),c(-Inf,Inf,-Inf,Inf,-Inf,Inf),info="calling with ordinary network")

# test when called with something else
expect_error(get.vertex.activity(list()),'requires an argument of class network')

# test with null spell
net<- network.initialize(2)
deactivate.vertices(net,onset=-Inf,terminus=Inf)
expect_equal(sapply(get.vertex.activity(net),is.null),c(TRUE,TRUE),info="testing for null spell")


#test for active.default behavior
expect_equal(get.vertex.activity(network.initialize(1))[[1]],matrix(c(-Inf,Inf),ncol=2),info='check get.vertex.activity active default behavior')

expect_equal(get.vertex.activity(network.initialize(1),active.default=FALSE)[[1]],NA,info='check get.vertex.activity active.default=FALSE behavior')

expect_equal(nrow(get.vertex.activity(network.initialize(0),as.spellList=TRUE)),0)

# test with deleted spell (was failing in v 0.4)
net<-network.initialize(3)
delete.vertex.activity(net,v=2)
expect_equal(get.vertex.activity(net,as.spellList=TRUE)$onset,c(-Inf,-Inf,-Inf))

# test with deactiviate spell (was failing in v 0.4)
net<-network.initialize(3)
deactivate.vertices(net,v=2)
expect_equal(get.vertex.activity(net,as.spellList=TRUE)$onset,c(-Inf,-Inf))
expect_equal(nrow(get.vertex.activity(net,as.spellList=TRUE,active.default=FALSE)),0)



# --------------- networkDynamic() conversion test -----

# ---- networkDynamic() vertex spells -----
vert.spls <-matrix(
  c(1,2,1,
  2,3,2,
  3,4,3,
  4,5,4),ncol=3,byrow=TRUE)

nd <-networkDynamic(vertex.spells=vert.spls)
if (network.size(nd)!=4){
  stop("networkDynamic did not create network of appropriate size from vertex.spells")
}
if(!all(get.vertex.activity(nd,as.spellList=TRUE)[1:3]==vert.spls)){
  stop("networkDynamic did not create network with appropriate vertex.spells from vertex.spells")
}
   
# vertex.spells  impute missing vertex?   
vert.spls <-matrix(
     c(1,2,1,
       2,3,2,
       3,4,3,
       4,5,9),ncol=3,byrow=TRUE)
nd <-networkDynamic(vertex.spells=vert.spls)
if (network.size(nd)!=9){
  stop("networkDynamic did not create network of appropriate imputed size from vertex.spells")
}
   
# vertex.spells  bad spell input
vert.spls <-matrix(
     c(1,2,1,
       5,3,2,  #uh oh, this spell is backwards
       3,4,3,
       4,5,4),ncol=3,byrow=TRUE)
expect_error(
  networkDynamic(vertex.spells=vert.spls) #throws error, but doesn't say what line
,"Onset times must precede terminus times in activate.vertices")
   
# vertex.spells empty network
vert.spls <-matrix(numeric(0),ncol=3,nrow=0)
nd<-networkDynamic(vertex.spells=vert.spls)
expect_equal(network.size(nd),0,info="checking creating zero vertex network with networkDynamic converter")
   
#  ---- networkDynamic() edge spells ---------------
edge.spls <-matrix(
  c(1,2,1,2,
    2,3,2,3,
    3,4,3,4,
    4,5,4,5),ncol=4,byrow=TRUE)
nd <-networkDynamic(edge.spells=edge.spls)
if(network.edgecount(nd)!=4){
  stop("networkDynamic() did not create appropriate number of edges from edge.spells")
}
if(network.size(nd)!=5){
  stop("networkDynamic() did not create appropriate number of vertices from edge.spells")
}

# edge spells - 0 edges
# this generates 'mysterious' warning
# https://statnet.csde.washington.edu/trac/ticket/189
net<-network.initialize(5)  #will crash if no vertices
edge.spls <-matrix(0,ncol=4,nrow=0,byrow=TRUE)
nd <-networkDynamic(base.net=net,edge.spells=edge.spls)


# edge spells - Loops
edge.spls <-matrix(
  c(1,2,1,1,
    1,2,2,2),ncol=4,byrow=TRUE)
nd <-networkDynamic(edge.spells=edge.spls)


# NA values  
edge.spls <-matrix(
  c(1,2,1,NA,
    1,2,"a",2),ncol=4,byrow=TRUE)
expect_error(networkDynamic(edge.spells=edge.spls),"must be numeric")

# edge spells - infer network size
edge.spls <-matrix(
  c(1,2,1,1,
    1,2,2,9),ncol=4,byrow=TRUE)
nd <-networkDynamic(edge.spells=edge.spls)
if(network.size(nd)!=9){
  stop("networkDynamic() did not infer network size from edge ids as expected")
}

# network properties from base net preserved
net <-network.initialize(5,directed=FALSE,loops=TRUE)
edge.spls <-matrix(
  c(1,2,1,2,
    2,3,2,3,
    3,4,3,4,
    4,5,4,5),ncol=4,byrow=TRUE)
nd <-networkDynamic(base.net=net,edge.spells=edge.spls)
if (is.directed(nd) & !has.loops(nd)){
  stop("networkDynamic did not preserve basic network settings from base.net argument")
}

# constsruct with both edges and vertices
edge.spls <-matrix(
  c(1,2,1,2,
    2,3,2,3,
    3,4,3,4,
    4,5,4,5),ncol=4,byrow=TRUE)
vert.spls <-matrix(
  c(1,2,1,
    2,3,2,
    3,4,3,
    4,5,4),ncol=3,byrow=TRUE)
nd<-networkDynamic(vertex.spells=vert.spls, edge.spells=edge.spls)

# inconsistant edges and vertices


#  ---- networkDynamic() vertex toggles - no base net ---------------
vrt.tog <-matrix(
  c(1,1,
    2,2,
    3,2),ncol=2,byrow=TRUE)
nd <-networkDynamic(vertex.toggles=vrt.tog)
if(!all(get.vertex.activity(nd,as.spellList=TRUE)[,1:3]==matrix(c(-Inf,1,1, -Inf,2,2, 3,Inf,2),ncol=3,byrow=TRUE))){
  stop("networkDynamic() did not produduce expected output for vertex.toggles ")
}

# check the list version with the inf values
expect_equal(unlist(get.vertex.activity(nd)),c(-Inf,1,-Inf,3,2,Inf))

# vertex toggles - with base net
net <-network.initialize(3)
vrt.tog <-matrix(
  c(1,1,
    2,2,
    3,2),ncol=2,byrow=TRUE)
nd <-networkDynamic(base.net=net,vertex.toggles=vrt.tog)
if(!all(get.vertex.activity(nd,as.spellList=TRUE)[,1:3]==matrix(c(-Inf,1,1,
                                                                  -Inf,2,2,
                                                                  3,Inf,2,
                                                                  -Inf,Inf,3),ncol=3,byrow=TRUE))){
  stop("networkDynamic() did not produce expected output for vertex.toggles and base.net")
}

# check the list version with the inf values
expect_equal(unlist(get.vertex.activity(nd)),c(-Inf,1,-Inf,3,2,Inf,-Inf,Inf))


#  ---- networkDynamic() edge toggles ------------
edge.tog <- matrix(
  c(1,1,2,
    2,2,3,
    3,2,3),ncol=3,byrow=TRUE)
nd<-networkDynamic(edge.toggles=edge.tog)
if(!all(get.edge.activity(nd,as.spellList=TRUE)[,1:4]==matrix(c(1,Inf,1,2, 2,3,2,3),ncol=4,byrow=TRUE))){
  stop("networkDynamic() did not produce expected output for edge.toggles")
}
expect_equal(unlist(get.edge.activity(nd)),c(1,Inf,2,3))

net<-network.initialize(4)
net[3,4]<-1
net[1,4]<-1
edge.tog <- matrix(
  c(1,1,2,
    2,2,3,
    3,2,3,
    0,1,4),ncol=3,byrow=TRUE)
nd<-networkDynamic(base.net=net,edge.toggles=edge.tog)
if(!all(get.edge.activity(nd,as.spellList=TRUE)[,1:4]==matrix(c(-Inf,Inf,3,4, 
                                                                -Inf,0,1,4, 
                                                                1,Inf,1,2, 
                                                                2,3,2,3),ncol=4,byrow=TRUE))){
  stop("networkDynamic() did not produce expected output for edge.toggles with base.net")
}
expect_equal(unlist(get.edge.activity(nd)),c(-Inf,Inf,-Inf,0,1,Inf,2,3))

#  ---- networkDynamic() vertex changes ----------
vrt.cng <-matrix(
  c(1,1,1,
    2,2,1,
    3,2,0),ncol=3,byrow=TRUE)
nd <-networkDynamic(vertex.changes=vrt.cng)
if (!all(get.vertex.activity(nd,as.spellList=TRUE)[,1:3]==matrix(c(1,Inf,1,2, 2,3,2,3),ncol=4,byrow=TRUE) )){
  stop("networkDynamic() did not produce expected output for vertex.changes argument")
}

expect_equal(unlist(get.vertex.activity(nd)),c(1,Inf,2,3))

#  ---- networkDynamic() edge changes -----------
edge.cng <-matrix(
  c(1,1,2,1,
    2,2,3,1,
    3,2,3,0),ncol=4,byrow=TRUE)
nd <-networkDynamic(edge.changes=edge.cng)
if (!all(get.edge.activity(nd,as.spellList=TRUE)[,1:4]==matrix(c(1,Inf,1,2, 2,3,2,3),ncol=4,byrow=TRUE))){
  stop("networkDynamic() did not produce expected output for edge.changes argument")
}

expect_equal(unlist(get.edge.activity(nd)),c(1,Inf,2,3))
# check net.obs.period defaults
expect_equal((nd%n%'net.obs.period')$'observations'[[1]],c(1,Inf))

# edge changes - activate edge allready active ignored
edge.cng <-matrix(
  c(2,2,3,1,
    3,2,3,1),ncol=4,byrow=TRUE)
nd <-networkDynamic(edge.changes=edge.cng)
expect_equivalent(as.numeric(get.edge.activity(nd,as.spellList=TRUE)[1:4]),c(2,Inf,2,3))

# check net.obs.period defaults
expect_equal((nd%n%'net.obs.period')$'observations'[[1]],c(2,Inf))




#  ---- networkDynamic() list of networks ---------- 
#try converting the newcomb panel data (working 9/3)
data(newcomb)
newDyn <- networkDynamic(network.list=newcomb[1:3])
#does it pass a consistency check?
check <- network.dynamic.check(newDyn) 
if (!all(sapply(check, all)))
  stop("newcomb network.list conversion did not pass network.dynamic.check")
#is the matrix equal to the input matrix
for (k in 1:3) {
  if (!all(as.sociomatrix(newcomb[[k]]) == as.sociomatrix(network.extract(newDyn,onset=k-1,terminus=k)))){
    stop("FAIL: networkDynamic conversion: 1st input matrix does not match crosssection from time 0 to 1 for newcomb example")
  }
}

# try converting a list that includes different size networks. (working 9/10)
# note that beach[[25]] is NA (missing)
data(windsurferPanels)
beach<-beach[1:7]
# should return error
dynBeach=NULL
expect_error(dynBeach<-networkDynamic(network.list = beach),
             "vertex.pid must be specified")
                   
if (is.network(dynBeach)) stop("did not ask for vertex.pid when network.list have different size networks")

dynBeach<-networkDynamic(network.list=beach, vertex.pid="vertex.names")

#data level node indicies are stored as the vertex names
expect_equal(get.network.attribute(dynBeach,'vertex.pid'),'vertex.names',info='check vertex.pid name set correctly')

#check if the neighborhood is the same for both.
for (i in 1:length(beach)) {
  if (!identical(beach[[i]], NA)) {
    beach[[i]]<-set.network.attribute(beach[[i]],'vertex.pid','vertex.names')
    for (j in network.vertex.names(beach[[i]])) {
      ng1 <-get.neighborhood(beach[[i]], v=get.vertex.id(beach[[i]], j))
      ng2 <- get.neighborhood.active(dynBeach, onset=i-1, terminus=i, v=get.vertex.id(dynBeach, j))
      # the following line is much much slower
      #ng2 <-get.neighborhood(network.extract(dynBeach,onset=i-1,terminus=i,retain.all.vertices=T),v=get.vertex.id(dynBeach, j))
      # need to check the vertex names, not the ids which are changed when converting to a networkDynamic object
      names1 <- sort(network.vertex.names(beach[[i]])[ng1])
      names2 <- sort(network.vertex.names(dynBeach)[ng2])
      # print these if necessary
      # print(paste('============ ', i, '=========='))
      # print(names1)
      # print(names2)
      if (!identical(names1, names2)) {
        print(paste("FAIL: networkDynamic(): neigborhoods do not match for variable sized network list example (windsurfers)",
                   " at time", i, 'vertex', j))
        print(names1)
        print(names2)
      }
    }
    
  }
}

# try a better (but truncated to make test fast) representation that preserves edge times and gaps

data(windsurferPanels)
beach<-beach[c(24,26)]
dynBeach<-networkDynamic(network.list=beach, vertex.pid="vertex.names",onsets=c(24,26),termini=c(25,27))


# make sure day 25 really is missing
if (any(is.active(dynBeach, at=25,v=1:network.size(dynBeach)))){
  stop("onsets and termini did not correctly omit day 25 in windsurfer example")
}

# check net.obs.period

expect_equal(do.call(rbind,(dynBeach%n%'net.obs.period')$observations),cbind(c(24,26),c(25,27)),info='was net.obs.period created by default by networkDynamic() with onsets and termini did not have expected range')




#  ---- networkDynamic()  network list TEAs -----------

#try a reduced newcomb version that has edge weights
newRankDyn <-networkDynamic(network.list=newcomb.rank[1:2],create.TEAs=TRUE)
# check that it matches original
expect_equal(as.matrix(network.collapse(newRankDyn,at=0),attrname='rank'),as.matrix(newcomb.rank[[1]],attrname='rank'))
expect_equal(as.matrix(network.collapse(newRankDyn,at=1),attrname='rank'),as.matrix(newcomb.rank[[2]],attrname='rank'))

# test vertex TEA from list of odd-sized networks
netlist<-list(network.initialize(3),network.initialize(1),network.initialize(2),network.initialize(2))
netlist[[1]]<-set.vertex.attribute(netlist[[1]],'id',c('a','b','c'))
netlist[[2]]<-set.vertex.attribute(netlist[[2]],'id','b')
netlist[[3]]<-set.vertex.attribute(netlist[[3]],'id',c('c','d'))
netlist[[4]]<-set.vertex.attribute(netlist[[4]],'id',c('c','d'))
netlist[[1]]<-set.vertex.attribute(netlist[[1]],'val',c('a1','b1','c1'))
netlist[[2]]<-set.vertex.attribute(netlist[[2]],'val','b2')
netlist[[3]]<-set.vertex.attribute(netlist[[3]],'val',c('c3','d3'))
dyn<-networkDynamic(network.list=netlist,vertex.pid='id',create.TEAs=TRUE)
expect_equal(get.vertex.attribute.active(dyn,'val',at=0),c('a1','b1','c1',NA))
expect_equal(get.vertex.attribute.active(dyn,'val',at=1),c(NA,'b2',NA,NA))
expect_equal(get.vertex.attribute.active(dyn,'val',at=2),c(NA,NA,'c3','d3'))
expect_true(all(is.na(get.vertex.attribute.active(dyn,'val',at=3))))

# test edge TEA from list of odd-sized networks
netlist[[1]]<-add.edges(netlist[[1]],tail=1:2,head=2:3)
netlist[[3]]<-add.edges(netlist[[3]],tail=1,head=2)
netlist[[4]]<-add.edges(netlist[[4]],tail=1,head=2)
netlist[[1]]<-set.edge.attribute(netlist[[1]],'eid','ab1',e=1)
netlist[[1]]<-set.edge.attribute(netlist[[1]],'eid','bc1',e=2)
netlist[[3]]<-set.edge.attribute(netlist[[3]],'eid','cd3',e=1)
netlist[[4]]<-set.edge.attribute(netlist[[4]],'eid','cd4',e=1)
dyn<-networkDynamic(network.list=netlist,vertex.pid='id',create.TEAs=TRUE)
expect_equal(get.edge.attribute.active(dyn,'eid',at=0),c("ab1","bc1" ,NA   ))
expect_equal(get.edge.attribute.active(dyn,'eid',at=1),c(NA,NA ,NA))
expect_equal(get.edge.attribute.active(dyn,'eid',at=2),c(NA,NA,"cd3"))
expect_equal(get.edge.attribute.active(dyn,'eid',at=3),c(NA,NA,"cd4"))

# test network TEA from list of odd-sized networks
netlist[[1]]<-set.network.attribute(netlist[[1]],'netname','first')
netlist[[2]]<-set.network.attribute(netlist[[2]],'netname','second')
netlist[[3]]<-set.network.attribute(netlist[[3]],'netname','third')
netlist[[4]]<-set.network.attribute(netlist[[4]],'netname','forth')
dyn<-networkDynamic(network.list=netlist,vertex.pid='id',create.TEAs=TRUE)
expect_equal(unlist(get.network.attribute.active(dyn,'netname',onset=-Inf,terminus=Inf,return.tea=TRUE)[[1]]),c("first",  "second", "third",  "forth"))


# test TEAs with same size networks, no vertex.pid
netlist<-list(network.initialize(4),network.initialize(4),network.initialize(4),network.initialize(4))
template<-netlist[[1]]
template<-set.vertex.attribute(template,'id',c('a','b','c','d'))
netlist[[1]]<-set.vertex.attribute(netlist[[1]],'val',c('a1','b1','c1'),v=1:3)
netlist[[2]]<-set.vertex.attribute(netlist[[2]],'val','b2',v=2)
netlist[[3]]<-set.vertex.attribute(netlist[[3]],'val',c('c3','d3'),v=3:4)
dyn<-networkDynamic(network.list=netlist,base.net=template,create.TEAs=TRUE)
expect_equal(get.vertex.attribute(dyn,'id'),c('a','b','c','d'))
expect_equal(get.vertex.attribute.active(dyn,'val',at=0),c('a1','b1','c1',NA))
expect_equal(get.vertex.attribute.active(dyn,'val',at=1),c(NA,'b2',NA,NA))
expect_equal(get.vertex.attribute.active(dyn,'val',at=2),c(NA,NA,'c3','d3'))
expect_true(all(is.na(get.vertex.attribute.active(dyn,'val',at=3))))

# test edge TEA from list of same-sized networks
netlist[[1]]<-add.edges(netlist[[1]],tail=1:2,head=2:3)
netlist[[3]]<-add.edges(netlist[[3]],tail=3,head=4)
netlist[[4]]<-add.edges(netlist[[4]],tail=3,head=4)
netlist[[1]]<-set.edge.attribute(netlist[[1]],'eid','ab1',e=1)
netlist[[1]]<-set.edge.attribute(netlist[[1]],'eid','bc1',e=2)
netlist[[3]]<-set.edge.attribute(netlist[[3]],'eid','cd3',e=1)
netlist[[4]]<-set.edge.attribute(netlist[[4]],'eid','cd4',e=1)
dyn<-networkDynamic(network.list=netlist,create.TEAs=TRUE)
expect_equal(get.edge.attribute.active(dyn,'eid',at=0),c("ab1","bc1" ,NA   ))
expect_equal(get.edge.attribute.active(dyn,'eid',at=1),c(NA,NA ,NA))
expect_equal(get.edge.attribute.active(dyn,'eid',at=2),c(NA,NA,"cd3"))
expect_equal(get.edge.attribute.active(dyn,'eid',at=3),c(NA,NA,"cd4"))

# test network TEA from list of same-sized networks
netlist[[1]]<-set.network.attribute(netlist[[1]],'netname','first')
netlist[[2]]<-set.network.attribute(netlist[[2]],'netname','second')
netlist[[3]]<-set.network.attribute(netlist[[3]],'netname','third')
netlist[[4]]<-set.network.attribute(netlist[[4]],'netname','forth')
dyn<-networkDynamic(network.list=netlist,create.TEAs=TRUE)
expect_equal(unlist(get.network.attribute.active(dyn,'netname',onset=-Inf,terminus=Inf,return.tea=TRUE)[[1]]),c("first",  "second", "third",  "forth"))


# -------- networkDynamic() edge spell tea tests ---------

# test dyad eid lookup
test<-network.initialize(6,loops=TRUE)
add.edges(test,tail=1:3,head=2:4)
add.edges(test,tail=5,head=5)

# dyads with no edge should return NA
expect_equal(networkDynamic:::get.dyads.eids(test,1,1),NA) 
# get eids back
expect_equal(networkDynamic:::get.dyads.eids(test,1:3,2:4),1:3)
# self loops work
expect_equal(networkDynamic:::get.dyads.eids(test,5,5),4)  
# error if lengths of head and tails differ
expect_error(networkDynamic:::get.dyads.eids(test,1,1:3),regexp = 'length of the tails and heads parameters must be the same') 

# test for multiplex throws warning
add.edges(test,tail=1,head=2)
expect_warning(expect_equal(networkDynamic:::get.dyads.eids(test,1,2),1), regexp = 'only smallest eid returned')

# ok, now test edge eids


# create an edge spell matrix where the last two columns will be TEA vals
testnumbers<-matrix(c(1,2,1,2, 1, 0.5,
                      3,4,1,2, 2, 0.1,
                      4,5,1,2, 0, 0.1,
                      5,7,2,3, 3, -1,
                      5,7,3,4, 4, 1.5),ncol=6,byrow=TRUE)

testnet<-networkDynamic(edge.spells=testnumbers,create.TEAs = TRUE,edge.TEA.names = c('value','weight'))


expect_true('value.active'%in%list.edge.attributes(testnet))
expect_true('weight.active'%in%list.edge.attributes(testnet))

expect_equal(get.edge.attribute.active(testnet,'value',at=1),c(1,NA,NA))
expect_equal(get.edge.attribute.active(testnet,'value',at=4),c(0,NA,NA))
expect_equal(get.edge.attribute.active(testnet,'value',at=5),c(NA,3,4))
expect_equal(get.edge.attribute.active(testnet,'weight',at=5),c(NA,-1,1.5))

# test for mismatch between number of cols and number of names
expect_error(testnet<-networkDynamic(edge.spells=testnumbers,create.TEAs = TRUE,edge.TEA.names = c('value','weight','foo')),regexp = 'edge.TEA.names must match the number of remaining columns in edge')

# now try with a non-numeric value
testletters<-data.frame(onset=c(1,2,5,5),
                    terminus=c(2,4,7,7),
                        head=c(1,1,2,3),
                        tail=c(2,2,3,4),
                      value=c('A','B','C','D'),stringsAsFactors=FALSE)
testnet<-networkDynamic(edge.spells=testletters,create.TEAs = TRUE,edge.TEA.names = c('value'))

# careful, these tests fail if character vector in data frame is converted to a factor
expect_equal(get.edge.attribute.active(testnet,'value',at=1),c('A',NA,NA))
expect_equal(get.edge.attribute.active(testnet,'value',at=2),c('B',NA,NA))
expect_equal(get.edge.attribute.active(testnet,'value',at=5),c(NA,'C','D'))

# test guessing col names from data.frame
testnet<-networkDynamic(edge.spells=testletters,create.TEAs = TRUE,edge.TEA.names = c('value'))
expect_true('value.active'%in%list.edge.attributes(testnet))

testnet<-networkDynamic(edge.spells=testletters,create.TEAs = FALSE)
expect_false('value.active'%in%list.edge.attributes(testnet))

# run test on a smallish realistic dataset
vertexData <-read.table(system.file('extdata/cls33_10_16_96_vertices.tsv', 
                                    package='networkDynamic'),header=TRUE,stringsAsFactors=FALSE)
edgeData <-read.table(system.file('extdata/cls33_10_16_96_edges.tsv', 
                                  package='networkDynamic'),header=TRUE,stringsAsFactors=FALSE)
classDyn <- networkDynamic(vertex.spells=vertexData[,c(3,4,1)],
                             edge.spells=edgeData[,c(3,4,1,2,5,6)],
                             create.TEAs=TRUE,edge.TEA.names=c('weight','type'))
expect_equal(c('weight.active','type.active')%in%list.edge.attributes(classDyn),c(TRUE,TRUE))
# check that it actually stored some values
expect_equal(get.edge.attribute.active(classDyn,'weight',onset=0,terminus=5,rule='earliest'),c(1, 1, 1, 1, 1, 1, 1, 1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,  0.2, 0.2, 0.2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, NA, NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA, NA, NA, NA, NA, NA))

expect_equal(get.edge.attribute.active(classDyn,'type',onset=0,terminus=5,rule='earliest'),c("social", "social", "sanction", "sanction", "sanction", "sanction",  "task", "task", "task", "task", "task", "task", "task", "task",  "task", "task", "task", "task", "task", "task", "task", "task",  "task", "task", "task", "task", "task", "task", "task", "task",  "task", "task", "task", "task", "task", "task", "task", "task",  "task", "task", "task", "task", "task", "task", "social", "social",  "social", "social", "social", "social", "social", "social", "social",  "social", "social", "social", NA, NA, NA, NA, NA, NA, NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,  NA))


# check a perverse case with unsorted values
testnumbers<-matrix(c(3,4,1,2, 2, 0.1,
                      1,2,1,2, 1, 0.5,
                      4,5,1,2, 0, 0.1),ncol=6,byrow=TRUE)

testnet<-networkDynamic(edge.spells=testnumbers,create.TEAs = TRUE,edge.TEA.names = c('value','weight'))
expect_equal(testnet$mel[[1]]$atl$value.active[[2]], structure(c(1, 3, 4, 2, 4, 5), .Dim = c(3L, 2L)))

# check a perverse case with intersecting values
testnumbers<-matrix(c(3,9,1,2, 2, 0.1,
                      1,4,1,2, 1, 0.5,
                      4,6,1,2, 0, 0.1),ncol=6,byrow=TRUE)

expect_warning(networkDynamic(edge.spells=testnumbers,create.TEAs = TRUE,edge.TEA.names = c('value','weight')),regexp = 'invalid spell matrix for edge TEA')


#-------- networkDynamic() vertex spell tea tests ---------
  
# create an vertex spell matrix where the last two columns will be TEA vals
testnumbers<-matrix(c(1,2,1, 1, 0.5,
                      3,4,1, 2, 0.1,
                      4,5,1, 0, 0.1,
                      5,7,2, 3, -1,
                      5,7,3, 4, 1.5),ncol=5,byrow=TRUE)

testnet<-networkDynamic(vertex.spells=testnumbers,create.TEAs = TRUE,vertex.TEA.names = c('value','weight'))

expect_true('value.active'%in%list.vertex.attributes(testnet))
expect_true('weight.active'%in%list.vertex.attributes(testnet))

expect_equal(get.vertex.attribute.active(testnet,'value',at=1),c(1,NA,NA))
expect_equal(get.vertex.attribute.active(testnet,'value',at=4),c(0,NA,NA))
expect_equal(get.vertex.attribute.active(testnet,'value',at=5),c(NA,3,4))
expect_equal(get.vertex.attribute.active(testnet,'weight',at=5),c(NA,-1,1.5))

# test for mismatch between number of cols and number of names
expect_error(testnet<-networkDynamic(vertex.spells=testnumbers,create.TEAs = TRUE,vertex.TEA.names = c('value','weight','foo')),regexp = 'vertex.TEA.names must match the number of remaining columns in vertex')

# now try with a non-numeric value
testletters<-data.frame(onset=c(1,2,5,5),
                        terminus=c(2,4,7,7),
                        vid=c(1,1,2,3),
                        value=c('A','B','C','D'),stringsAsFactors=FALSE)
testnet<-networkDynamic(vertex.spells=testletters,create.TEAs = TRUE,vertex.TEA.names = c('value'))

# careful, these tests fail if character vector in data frame is converted to a factor
expect_equal(get.vertex.attribute.active(testnet,'value',at=1),c('A',NA,NA))
expect_equal(get.vertex.attribute.active(testnet,'value',at=2),c('B',NA,NA))
expect_equal(get.vertex.attribute.active(testnet,'value',at=5),c(NA,'C','D'))

# test guessing col names from data.frame
testnet<-networkDynamic(vertex.spells=testletters,create.TEAs = TRUE,vertex.TEA.names = c('value'))
expect_true('value.active'%in%list.vertex.attributes(testnet))

# make sure it doesn't create if instructed not to
testnet<-networkDynamic(vertex.spells=testletters,create.TEAs = FALSE)
expect_false('value.active'%in%list.vertex.attributes(testnet))


# check a perverse case with unsorted values
testnumbers<-matrix(c(3,4,1, 2, 0.1,
                      1,2,1, 1, 0.5,
                      4,5,1, 0, 0.1),ncol=5,byrow=TRUE)

testnet<-networkDynamic(vertex.spells=testnumbers,create.TEAs = TRUE,vertex.TEA.names = c('value','weight'))
expect_equal(testnet$val[[1]]$value.active[[2]], structure(c(1, 3, 4, 2, 4, 5), .Dim = c(3L, 2L)))

# check a perverse case with intersecting values
testnumbers<-matrix(c(3,9,1, 2, 0.1,
                      1,4,1, 1, 0.5,
                      4,6,1, 0, 0.1),ncol=5,byrow=TRUE)

expect_warning(networkDynamic(vertex.spells=testnumbers,create.TEAs = TRUE,vertex.TEA.names = c('value','weight')),regexp = 'invalid spell matrix for vertex TEA')

# ----------- as.networkDynamic.data.frame tests -------

# check correct spells printed for edges
test <- network.initialize(3)
test[1,2]<-1
test[2,3]<-1
activate.edges(test,onset=1,terminus=2)
ref <- matrix(c(1,2,1,2,0,0,1,1, 1,2,2,3,0,0,1,2),ncol=8,byrow=TRUE)
if (!all(as.data.frame(test) == ref)){
  stop("unexpected output for edge spells from as.networkDynamic.data.frame") 
}

test <- network.initialize(3)
test[1,2]<-1
test[2,3]<-1
activate.edges(test,at=3)
ref <- matrix(c(3,3,1,2,0,0,0,1, 3,3,2,3,0,0,0,2),ncol=8,byrow=TRUE)
if (!all(as.data.frame(test) == ref)){
  stop("unexpected output for edge spells of zero length from as.networkDynamic.data.frame") 
}

# test for missing activity attribute (only set on one edge)
test <- network.initialize(3)
test[1,2]<-1
test[2,3]<-1
activate.edges(test,at=3,e=1)
as.data.frame(test)
tryCatch(
  as.data.frame(test), error = function(e){ warning(paste("error in as.networkDynamic.data.frame  for edge with missing activity attribute",e))} )


#check for duration for funny length spells
test <- network.initialize(3)
test[1,2]<-1
test[2,3]<-1
activate.edges(test,at=1,e=1)
activate.edges(test,onset=2.7,terminus=5, e=2)
ref <- matrix(c(1.0,1,1,2,0,0,0,1, 2.7,5,2,3,0,0,2.3,2),ncol=8,byrow=TRUE)
if (!all(as.data.frame(test) == ref)){
  stop("unexpected output for edge spells from as.networkDynamic.data.frame") 
}

# check multiple spells per edge
test <- network.initialize(3)
test[1,2]<-1
activate.edges(test,onset=1,terminus=2)
activate.edges(test,onset=3,terminus=4)
ref <- matrix(c(1,2,1,2,0,0,1,1, 3,4,1,2,0,0,1,1),ncol=8,byrow=TRUE)
if (!all(as.data.frame(test) == ref)){
  stop("unexpected output for multiple spells per edge for as.networkDynamic.data.frame") 
}

# check censoring arguments when passed in
test <- network.initialize(3)
test[1,2]<-1
activate.edges(test,onset=-Inf,terminus=10)
ref <- matrix(c(5,10,1,2,1,0,5,1),ncol=8,byrow=TRUE)
if (!all(as.data.frame(test,start=5) == ref)){
  stop("unexpected output for 'start' left censoring argument from as.networkDynamic.data.frame") 
}

test <- network.initialize(3)
test[1,2]<-1
activate.edges(test,onset=0,terminus=Inf)
ref <- matrix(c(0,5,1,2,0,1,5,1),ncol=8,byrow=TRUE)
if (!all(as.data.frame(test,end=5) == ref)){
  stop("unexpected output for 'end' left censoring argument from as.networkDynamic.data.frame") 
}

test <- network.initialize(3)
test[1,2]<-1
activate.edges(test,onset=-Inf,terminus=Inf)
ref <- matrix(c(-Inf,Inf,1,2,1,1,Inf, 1),ncol=8,byrow=TRUE)
if (!all(as.data.frame(test) == ref)){
  warning("unexpected output for as.networkDynamic.data.frame: Inf and -Inf times not treated as censored") 
}

# check censoring arguments when set on input object using attr
test <- network.initialize(3)
test[1,2]<-1
activate.edges(test,onset=-Inf,terminus=Inf)
attr(test,"start")<-5
expect_warning(as.data.frame(test),'has been deprecated',info='specifying start and end using attrs deprecated')
#skye: this seems to work, but I want to remove the feature in favor of net.obs.period


# test using net.obs.period
test <- network.initialize(3)
test[1,2]<-1
activate.edges(test,onset=-Inf,terminus=Inf)
set.network.attribute(test,'net.obs.period',list(observations=list(c(5,10)),mode='discrete',time.increment=1,time.unit='step'))
expect_equivalent(as.numeric(as.data.frame(test)[,1:2]),c(5,10),info='net.obs.list provides censoring info')


# =============== TESTING networkDynamic ===========================

# working 9/3/2012. Just compare the data frame output with edgetimes
# a really crude edgelist example
edgetimes <- as.data.frame(matrix( c(1,2,1,2, 1,2,3,4,  2,3,1,3 ),ncol=4,byrow=TRUE))
edgetimetest<-networkDynamic(edge.spells = edgetimes)
# do the edges and spells match when spit back out?
if (!all(as.data.frame(edgetimetest)[,1:4]==edgetimes)){
  stop("output spell matrix does not match input for networkDynamic()")
}

#does the internal representation match?
if( !all(as.vector(edgetimetest$mel[[1]]$atl$active) == c(1,2))){
  stop("networkDynamic() gave unexpected internal representation spells")
}


# combining multiple spells (should combine the spells for edge between v1 and v2)
edgetimes <- as.data.frame(matrix( c(1,2,1,2, 2,3,1,2,  1,2,3,4 ),ncol=4,byrow=TRUE))
edgetimetest<-networkDynamic(edge.spells = edgetimes)
if (!all(as.data.frame(edgetimetest)[,1:4]==matrix(c(1,3,1,2, 1,2,3,4),ncol=4,byrow=TRUE))){
  stop("output spell matrix did not merge input spells as expected for networkDynamic()")
}


# with censoring
# Skye:  why does input of Inf mean that it is right censored?
edgetimes <- as.data.frame(matrix( c(1,Inf,1,2, 2,3,2,3),ncol=4,byrow=TRUE))
edgetimetest<-networkDynamic(edge.spells = edgetimes)
if (!all(as.data.frame(edgetimetest)[,1:4]==edgetimes)){
  stop("output spell matrix did not merge input spells as expected for networkDynamic()")
}

# with missing node (should fill in the missing node)
edgetimes <- as.data.frame(matrix( c(1,2,1,2, 2,4,1,2,  1,2,1,4 ),ncol=4,byrow=TRUE))
edgetimetest<-networkDynamic(edge.spells = edgetimes)
if (network.size(edgetimetest)!=4){
  stop("networkDynamic() did not create network with implied size of 4")
}


# create with vertex and edge dynamics specified
nodetimes <-as.data.frame(matrix( c(1,1,2,1, 2,2,3,2,  3,3,4,3 ),ncol=4,byrow=TRUE))
edgetimes <- as.data.frame(matrix( c(1,2,1,2, 2,4,1,2,  1,2,1,4 ),ncol=4,byrow=TRUE))
nd<-networkDynamic(vertex.spells=nodetimes,edge.spells=edgetimes)

# check for net.obs.period
expect_equal((nd%n%'net.obs.period')$observations[[1]],c(1,4),info='net.obs.period created by default by networkDynamic(vertex.spells,edge.spells) did not have expected range')

# ---- networkDynamic changes conversion ----
#[time,tail,head,direction]
echange<-matrix( c(1,1,2,1, 
                   2,1,2,0,
                   2,2,3,1,
                   3,1,3,0),ncol=4,byrow=TRUE)

nd<-networkDynamic(edge.changes=echange)
# the expected output
spls<-data.frame(onset=c(1,2,-Inf),terminus=c(2,Inf,3),tail=c(1,2,1),head=c(2,3,3),onset.censored=c(FALSE,FALSE,TRUE),terminus.censored=c(FALSE,TRUE,FALSE),duration=c(1,Inf,Inf),edge.id=c(1,2,3))

# check that spells constructed correctly from edge.changes
expect_equivalent(as.data.frame(nd),spls,info='comparing networkDynamic(edge.changes) to resulting spell list')

# check net.obs.period construction
expect_equal((nd%n%'net.obs.period')$observations[[1]],c(-Inf,Inf),info='net.obs.period created by default by networkDynamic(edge.changes) did not have expected range')

expect_equal((nd%n%'net.obs.period')$mode,'discrete',info='net.obs.period created by default by networkDynamic(edge.changes) did not have expected mode')

# [time,vertex.id,direction,red_herring]
vchange<-matrix( c(1,1,1, 
                   2,1,0,
                   2,2,1,
                   3,3,0),ncol=3,byrow=TRUE)

nd<-networkDynamic(vertex.changes=vchange)
expect_equivalent(get.vertex.activity(nd,as.spellList=TRUE),data.frame(onset=c(1,2,-Inf),terminus=c(2,Inf,3),vertex.id=c(1,2,3),onset.censored=c(FALSE,FALSE,TRUE),terminus.censored=c(FALSE,TRUE,FALSE),duration=c(1,Inf,Inf)),info='networkDynamic(vertex.change) did not give expected spell list result')


# [time,vertex.id,direction,red_herring]
# checks for internal bug in column name assignments
vchange<-matrix( c(1,1,1,5, 
                   2,1,0,6,
                   2,2,1,7,
                   3,3,0,8),ncol=4,byrow=TRUE)
nd<-networkDynamic(vertex.changes=vchange)
expect_equivalent(get.vertex.activity(nd,as.spellList=TRUE),data.frame(onset=c(1,2,-Inf),terminus=c(2,Inf,3),vertex.id=c(1,2,3),onset.censored=c(FALSE,FALSE,TRUE),terminus.censored=c(FALSE,TRUE,FALSE),duration=c(1,Inf,Inf)),info='networkDynamic(vertex.change) did not give expected spell list result')


# ----- networkDynamic network.list conversion ----

# does it create a dynamic network
d1 <- network.initialize(3)
d1[1,2]<-1
d2 <- network.initialize(3)
d2[1,2]<-1
d2[2,3]<-1
d3 <- network.initialize(3)
d3[3,1]<-1

# default timing
ddyn <- networkDynamic(network.list = list(d1,d2,d3))

if(!is.networkDynamic(ddyn)){
  stop("as.networkDynamic.list didn't create a dynamic network from ")
}

# check that correct spells with unit lengths were created
ref <- matrix(c(0,2,1,2,0,0,2, 1,2,2,3,0,0,1, 2,3,3,1,0,0,1),ncol=7,byrow=TRUE)
if (!all(as.data.frame(ddyn)[,1:7]==ref)){
  stop("correct unit length spells were not created for list input networks in networkDynamic()")
}

# check that default net.obs.period created
expect_equal((ddyn%n%'net.obs.period')$observations[[1]],c(0,3),info='was net.obs.period created by default by networkDynamic() did not have expected range')

# does it preserve network attributes of passed in network
d1 <- network.initialize(2,directed=F,bipartite=T,multiple=T,loops=T)
d2 <- network.initialize(2,directed=F,bipartite=T,multiple=T,loops=T)
dlist <- list(d1,d2)
ddyn <- networkDynamic(network.list = dlist)

if (is.directed(ddyn) != FALSE){
  stop("'directed' argument if initial network in list not respected in dynamic version")
}
if (is.bipartite(ddyn) != TRUE){
  stop("'bipartite' argument if initial network in list not respected in dynamic version")
}
if (is.multiplex(ddyn) != TRUE){
  stop("'multiple' argument if initial network in list not respected in dynamic version")
}
if (has.loops(ddyn) != TRUE){
  stop("'loops' argument if initial network in list not respected in dynamic version")
}


# does it warn if network attributes of passed in networks do not match
# working 9/10
d1 <- network.initialize(2,directed=F)
d2 <- network.initialize(2,directed=T)
dlist <- list(d1,d2)
expect_warning(networkDynamic(network.list = dlist), "have different network properties",info="different network attributes in network.list did not result in a warning")

# are vertex attributes included?
d1 <- network.initialize(2)
d2 <- network.initialize(2)
dbase<-network.initialize(2)
set.vertex.attribute(dbase,"test","one")
set.network.attribute(dbase,'another',"two")
dlist <- list(d1,d2)
dnet<-networkDynamic(base.net=dbase,network.list = dlist)

expect_true('test'%in%list.vertex.attributes(dnet),info="vertex attributes from base.net network not copied by networkDynamic()")

expect_true('another'%in%list.network.attributes(dnet),info="user specified network attributes from base.net network not copied by networkDynamic()")

dlist[[1]]<-set.vertex.attribute(dlist[[1]],"test","one")
dlist[[1]]<-set.network.attribute(dlist[[1]],'another',"two")
dnet<-networkDynamic(network.list = dlist)
expect_true('test'%in%list.vertex.attributes(dnet),info="vertex attributes from first item of network list network not copied by networkDynamic()")

expect_true('another'%in%list.network.attributes(dnet),info="user specified network attributes from base.net network not copied by networkDynamic()")

dnet<-networkDynamic(network.list = dlist,create.TEAs=TRUE)
expect_true('test.active'%in%list.vertex.attributes(dnet),info="vertex attributes from first item of network list network not copied as TEA by networkDynamic()")

expect_true('another.active'%in%list.network.attributes(dnet),info="user specified network attributes from base.net network not copied as TEA by networkDynamic()")


# are edge attributes included

# specify net.obs.period
d1 <- network.initialize(3)
d1[1,2]<-1
d2 <- network.initialize(3)
d2[1,2]<-1
d2[2,3]<-1
d3 <- network.initialize(3)
d3[3,1]<-1
nop = list(observations=list(c(3,5)), mode="discrete", time.increment=1,time.unit="step")
ddyn <- networkDynamic(network.list = list(d1,d2,d3), net.obs.period=nop)

expect_false(is.null(ddyn%n%'net.obs.period'),info="net.obs.period argument in input to networkDynamic reproduced in output")



# check for error when edge spell list has only one row issue #190
edge.spls <-matrix( c(1,2,1,2),ncol=4,byrow=TRUE)
nd <-networkDynamic(edge.spells=edge.spls)
expect_equivalent(as.matrix(get.edge.activity(nd,as.spellList=TRUE)[,1:4]),edge.spls,info="check when edge spell list has only one row")


# check for char data in input tables edge spells
net <-network.initialize(4)
edgetimes <- as.data.frame( list(c("1",1,2),c(2,2,3),c(1,3,1),c(2,4,3)))
expect_error(edgetimetest<-networkDynamic(base.net=net,edge.spells = edgetimes),"the onset time column of the edge.spells argument to networkDynamic must be numeric",info="testing non-numeric input to networkDynamic edgespells onset")

edgetimes <- as.data.frame( list(c(1,1,2),c(2,"2",3),c(1,3,1),c(2,4,3)))
expect_error(edgetimetest<-networkDynamic(base.net=net,edge.spells = edgetimes),"the terminus time column of the edge.spells argument to networkDynamic must be numeric",info="testing non-numeric input to networkDynamic edgespells terminus")

edgetimes <- as.data.frame( list(c(1,1,2),c(2,2,3),c(1,"3",1),c(2,4,3)))
expect_error(edgetimetest<-networkDynamic(base.net=net,edge.spells = edgetimes),"must be a numeric",info="testing non-numeric input to networkDynamic edgespells terminus")

edgetimes <- as.data.frame( list(c(1,1,2),c(2,2,3),c(1,3,1),c(2,4,"3")))
expect_error(edgetimetest<-networkDynamic(base.net=net,edge.spells = edgetimes),"must be a numeric",info="testing non-numeric input to networkDynamic edgespells terminus")

# check for char data in input tables edge toggles
edgetimes <- as.data.frame( list(c(1,1,2),c(2,"2",3),c(1,0,1)))
expect_error(edgetimetest<-networkDynamic(base.net=net,edge.toggles = edgetimes),"must be a numeric",info="testing non-numeric input to networkDynamic edge.toggles tail")

edgetimes <- as.data.frame( list(c(1,"1",2),c(2,2,3),c(1,0,1)))
expect_error(edgetimetest<-networkDynamic(base.net=net,edge.toggles = edgetimes),"must be numeric",info="testing non-numeric input to networkDynamic edge.toggles time")

edgetimes <- as.data.frame( list(c(1,1,2),c(2,2,3),c(1,"0",1)))
expect_error(edgetimetest<-networkDynamic(base.net=net,edge.toggles = edgetimes),"must be a numeric",info="testing non-numeric input to networkDynamic edge.toggles head")


# check for char data in input tables edge changes
edgetimes <- as.data.frame( list(c(1,"1",2),c(2,2,3),c(2,4,3),c(1,0,1)))
expect_error(edgetimetest<-networkDynamic(base.net=net,edge.changes = edgetimes),"must be numeric",info="testing non-numeric input to networkDynamic edge.changes time")

edgetimes <- as.data.frame( list(c(1,1,2),c(2,2,"3"),c(2,4,3),c(1,0,1)))
expect_error(edgetimetest<-networkDynamic(base.net=net,edge.changes = edgetimes),"must be a numeric",info="testing non-numeric input to networkDynamic edge.changes tail")

edgetimes <- as.data.frame( list(c(1,1,2),c(2,2,3),c(2,4,"3"),c(1,0,1)))
expect_error(edgetimetest<-networkDynamic(base.net=net,edge.changes = edgetimes),"must be a numeric",info="testing non-numeric input to networkDynamic edge.changes head")

edgetimes <- as.data.frame(list(c(1,1,2),c(2,2,3),c(2,4,3),c(1,0,"1")))
expect_error(edgetimetest<-networkDynamic(base.net=net,edge.changes = edgetimes),"must be numeric",info="testing non-numeric input to networkDynamic edge.changes direction")


# check for char data in input tables vertex toggles

verttimes <- as.data.frame(list(c(1,"1",2),c(1,2,3)))
expect_error(networkDynamic(base.net=net,vertex.toggles = verttimes),"be numeric",info="testing non-numeric input to networkDynamic vertex.toggles time")

verttimes <- as.data.frame(list(c(1,1,2),c(1,2,"3")))
expect_error(networkDynamic(base.net=net,vertex.toggles = verttimes),"be numeric",info="testing non-numeric input to networkDynamic vertex.toggles vertex.id")

# check for char data in input tables vertex spells
verttimes <- as.data.frame(list(c(1,"1",2),c(1,1,2),c(1,2,3)))
expect_error(networkDynamic(base.net=net,vertex.spells = verttimes),"be numeric",info="testing non-numeric input to networkDynamic vertex.spells onset")

verttimes <- as.data.frame(list(c(1,1,2),c(1,"1",2),c(1,2,3)))
expect_error(networkDynamic(base.net=net,vertex.spells = verttimes),"be numeric",info="testing non-numeric input to networkDynamic vertex.spells terminus")

verttimes <- as.data.frame(list(c(1,1,2),c(1,1,2),c(1,2,"3")))
expect_error(networkDynamic(base.net=net,vertex.spells = verttimes),"be numeric",info="testing non-numeric input to networkDynamic vertex.spells vertex.id")

# check for char data in input tables vertex changes

verttimes <- as.data.frame(list(c(1,"1",2),c(1,1,2),c(1,0,1)))
expect_error(networkDynamic(base.net=net,vertex.changes = verttimes),"be numeric",info="testing non-numeric input to networkDynamic vertex.changes time")

verttimes <- as.data.frame(list(c(1,1,2),c(1,"1",2),c(1,0,1)))
expect_error(networkDynamic(base.net=net,vertex.changes = verttimes),"be numeric",info="testing non-numeric input to networkDynamic vertex.changes vertex.id")

verttimes <- as.data.frame(list(c(1,1,2),c(1,1,2),c(1,0,"1")))
expect_error(networkDynamic(base.net=net,vertex.changes = verttimes),"be numeric",info="testing non-numeric input to networkDynamic vertex.changes direction")

# network list with network size 0

expect_equal(network.size(networkDynamic(network.list=list(network.initialize(0)))),0)

# check for edge.spells plus base.net with pre-existing edges bug #556
test<-network.initialize(3)
test[1,2]<-1
networkDynamic(base.net=test,edge.spells=matrix(c(0,1,1,2),ncol=4))


# =============== TESTING as.data.frame.networkDynamic ====
message('testing as.data.frame.networkDynamic\n')
#lets try passing in a noncensored and duration, see if we get the same things back
edgetimes <- as.data.frame(matrix( c(1,2,1,2,0,0, 3,5,1,2,0,0,1,2,3,4,0,0, 3,4,2,4,0,0 ),ncol=6,byrow=TRUE))
colnames(edgetimes)<-c("onset","terminus","tail","head","onset.censored","terminus.censored")
testnet <- network.initialize(4)
add.edges.active(testnet,onset=1,terminus=2,tail=1,head=2)
add.edges.active(testnet,onset=1,terminus=2,tail=3,head=4)
activate.edges(testnet,onset=3,terminus=5,e=get.edgeIDs(testnet, v=1,alter=2))
add.edges.active(testnet,onset=3,terminus=4,tail=2,head=4)
# these should match
if (!all(as.data.frame(testnet)[,1:6] == edgetimes)){
  stop("FAIL: output data.frame from as.data.frame.networkDynamic did not match input")
}

#check column names
if(!all(names(as.data.frame(testnet))==c("onset","terminus","tail","head","onset.censored","terminus.censored","duration","edge.id"))){
  stop("Unexpected column names returned by as.data.frame.networkDynamic")
  
}

# censoring should set the appropriate start or end to Inf
# skye: Is this the behavior we want for as.data.frame?
edgetimes <- as.data.frame(matrix( c(0,2,1,2,1,0,2, 3,5,1,2,0,0,2,  1,2,3,4,0,0,1, 3,6,2,4,0,1,3 ),ncol=7,byrow=TRUE))
colnames(edgetimes)<-c("onset","terminus","tail","head","onset.censored","terminus.censored","duration")
testnet <- network.initialize(4)
add.edges.active(testnet,onset=-Inf,terminus=2,tail=1,head=2)
add.edges.active(testnet,onset=1,terminus=2,tail=3,head=4)
activate.edges(testnet,onset=3,terminus=5,e=get.edgeIDs(testnet, v=1,alter=2))
add.edges.active(testnet,onset=3,terminus=Inf,tail=2,head=4)

# these should match
if(!all(as.data.frame(testnet,start=0,end=6)[-8]==edgetimes)){
  stop("as.data.frame.networkDynamic gave unexpected censored spell matrix output")
}

# check censoring for non-Inf edges when start and end are set narrower
tel<-matrix(c(40,  72,  10, 4,
         214, 247, 1,  11,  
         224, 256, 7,10),ncol=4,byrow=TRUE)
test<-networkDynamic(edge.spells=tel) 
result<-as.data.frame(test,start=50,end=60)
expect_equal(nrow(result),1)
expect_equal(as.numeric(result[,1:4]),c(50,60,10,4))
expect_equal(as.logical(result[,5:6]),c(TRUE,TRUE))



# properly handle edges with no spell activity
test <- network.initialize(3)
test[1,2]<-1
test[2,3]<-1
activate.edges(test,at=3,e=1)
temp = as.data.frame(test,active.default=FALSE)
if (!all(temp == c(3,3,1,2,F,F,0,1))) stop('did not handle edges without spell activity')

expect_equivalent(as.data.frame(test,active.default=TRUE),
data.frame(onset=c(3,-Inf),terminus=c(3,Inf),
          tail=c(1,2),head=c(2,3),onset.censored=c(FALSE,TRUE),
           terminus.censored=c(FALSE,TRUE),duration=c(0,Inf),edge.id=c(1,2)),
                  info="test active.default=TRUE for edge with no spell")


# properly handle nD object with no edges at all
# active default means it should still return two edges
net <-network.initialize(3)
activate.vertices(net, onset=1, terminus=Inf)
temp<-as.data.frame(net)

if (nrow(temp) != 0) {
  stop("as.data.frame.networkDynamic() did not handle an object without any edge activity")
}

# check for crash with network size 0
expect_equal(nrow(get.edge.activity(network.initialize(0),as.spellList=TRUE)),0)

# check active.default
net<-network.initialize(3)
add.edges(net,tail=1:3,head=c(2,3,1))
expect_equal(get.edge.activity(net,as.spellList=TRUE)$onset,c(-Inf,-Inf,-Inf))
expect_equal(nrow(get.edge.activity(net,as.spellList=TRUE,active.default=FALSE)),0)

# check for 'null' spell
net<-network.initialize(3)
add.edges(net,tail=1:3,head=c(2,3,1))
deactivate.edges(net,e=2)
spls<-as.data.frame(net,as.spellList=TRUE)
expect_equal(spls$onset,c(-Inf,-Inf),info='check for null spell and active.default')
expect_equal(spls$terminus,c(Inf,Inf),info='check for null spell and active.default')
expect_equal(spls$edge.id,c(1,3),info='check for null spell and active.default')

expect_equal(nrow(as.data.frame(net,as.spellList=TRUE,active.default=FALSE)),0,info='check for null spell and active.default=FALSE')

# check for multiplex case

nd <-network.initialize(2,multiple=TRUE)
add.edges.active(nd,onset=12,terminus=12.1,tail=1,head=2)
add.edges.active(nd,onset=12,terminus=12.5,tail=1,head=2)
as.data.frame(nd)
expect_equal(nrow(as.data.frame(nd)),2,info="test as.data.frame with multiplex network")

# check for correct sort ordering  by onset,terminus, edge.id
m<-matrix(c(10,11,1,2,  12,13,1,2, 5,6,2,3, 12,13,2,3, 12,12.5,1,3  ),ncol=4,byrow=TRUE)
nd<-networkDynamic(network.initialize(3,multiple=TRUE),edge.spells=m)
# add in an edge to test multiplex case
add.edges.active(nd,onset=12,terminus=12.1,tail=1,head=3)
ndmat<-as.data.frame(nd)
expect_equal(ndmat$onset,c(5,12,10,12,12,12),info='check as.data.frame.networkDynamic onset sorting')
expect_equal(ndmat$terminus,c(6,13,11,13,12.5,12.1),info='check as.data.frame.networkDynamic terminus sorting')
expect_equal(ndmat$edge.id,c(1,1,2,2,3,4),info='check as.data.frame.networkDynamic edge.id sorting')

# check that e argument removes appropriate edges' spells
ndmat<-as.data.frame(nd,e=c(1,4))
expect_equal(ndmat$edge.id,c(1,1,4),info='check that e argument removes edge spells from as.data.frame.networkDynamic')

# check for behavior with deleted edges
m<-matrix(c(10,11,1,2,  12,13,1,2, 5,6,2,3, 12,13,2,3, 12,12.5,1,3  ),ncol=4,byrow=TRUE)
nd<-networkDynamic(network.initialize(3,multiple=TRUE),edge.spells=m)
# add in an edge to test multiplex case
add.edges.active(nd,onset=12,terminus=12.1,tail=1,head=3)
delete.edges(nd,eid=2)
ndmat<-as.data.frame(nd)
expect_equal(ndmat$edge.id,c(1,1,3,4),info='check that deleted edge doesnt cause problem for as.data.frame.networkDynamic')

# check for appropriate exclusion of terminating edges
test<-network.initialize(3)
add.edges.active(test,tail=c(1,2,3),head=c(2,3,1),onset=c(0,1,2),terminus=c(1,2,3))
# first edge should be excluded because it lies entirly outside query range
# second edge could be excluded because it terminates at onset of query range
expect_equal(as.data.frame(test,start=2,end=3)$edge.id,3)
# also check for appropriate inclusion at other boundry
expect_equal(as.data.frame(test,start=1,end=2)$edge.id,2)

# and for at spell query
test<-network.initialize(3)
add.edges.active(test,tail=c(1,2,3,3),head=c(2,3,1,1),onset=c(0,1,2,1),terminus=c(1,2,3,1))
expect_equal(as.data.frame(test,start=1,end=1)$edge.id,c(2,4))


# ----- get.edge.activity ----
# some of this is not tested becaue it actually calls as.data.frame internally

# check subsetting with e argument
net<-network.initialize(5)
add.edges.active(net,tail=1:4,head=2:5,onset=1:4,terminus=2:5)
expect_equal(unlist(get.edge.activity(net)),c(1,2,2,3,3,4,4,5))

expect_equivalent( get.edge.activity(net,as.spellList=TRUE,e=2:3),
                    data.frame(onset=2:3,terminus=3:4,tail=2:3,head=3:4,
                    onset.censored=c(FALSE,FALSE),terminus.censored=c(FALSE,FALSE),
                               duration=c(1,1),edge.id=2:3) ,
                   info="test edge activity subsetting with e")

expect_equal(unlist(get.edge.activity(net,e=2:3)),c(2,3,3,4))


expect_equal(get.edge.activity(network.initialize(0)),list())

# test with ordinary net and default activity
net<-network.initialize(3)
add.edges(net,c(1,2,3),c(2,3,1))

expect_equal(sapply(get.edge.activity(net,active.default=FALSE),is.null),
             c(TRUE,TRUE,TRUE),info='check get.edge.activity with ordinary net and active.default=FALSE')

expect_equal(unlist(get.edge.activity(net,active.default=TRUE)),
             c(-Inf,Inf,-Inf,Inf,-Inf,Inf),info='check get.edge.activity with ordinary net and active.default=TRUE')

# can it distinguish from deleted edges?

delete.edges(net,e=2)
spls<-get.edge.activity(net,active.default=TRUE)
expect_true(is.null(spls[[2]]),info='check get.edge.activity with ordinary net and active.default=TRUE distinguish deleted edge')
expect_equal(spls[[3]],matrix(c(-Inf,Inf),ncol=2),info='check get.edge.activity with ordinary net and active.default=TRUE distinguish deleted edge')

# test `null` spell
net<-network.initialize(3)
add.edges(net,c(1,2,3),c(2,3,1))
deactivate.edges(net,e=2)
expect_equal(sapply(get.edge.activity(net,active.default=FALSE),is.null),
             c(TRUE,TRUE,TRUE),info='check get.edge.activity with ordinary net and null spell active.default=FALSE')

spls<-get.edge.activity(net,active.default=TRUE)
expect_true(is.null(spls[[2]]),info='check get.edge.activity with ordinary net and active.default=TRUE distinguish deleted edge')
expect_equal(spls[[3]],matrix(c(-Inf,Inf),ncol=2),info='check get.edge.activity with ordinary net and active.default=TRUE and null spell')

# check that it doesn't censor by default with net.obs.period
net<-network.initialize(2)
add.edges.active(net,tail=1,head=2,onset=-Inf,terminus=Inf)
net%n%'net.obs.period'<-list(observations=list(c(0,100)),mode="discrete", time.increment=1,time.unit="step")
expect_equal(as.numeric(get.edge.activity(net,as.spellList=TRUE)[1:2]),c(-Inf,Inf))




# ------- networkDynamic from toggles (function used by TERGM) ------------

test <- network.initialize(3)
test[1,3]<-1
tog <- matrix(c(1,1,2, 1,2,3, 2,1,2, 4,1,3, 4,1,2), ncol=3, byrow=TRUE)
net<-networkDynamic(base.net=test,edge.toggles=tog)
spells <-as.data.frame(net)[1:6]
# first spell should be onset censored because edge was in original net
# all edges in original net are considered active before toggles
if (!all(spells[1,]==c(-Inf,4,1,3,1,0))){
  stop("networkDynamic() did not record initial toggle correctly")
}
# 2nd spell toggles twice
if (!all(spells[2,]==c(1,2,1,2,0,0))){
  stop("networkDynamic() did not record double toggle correctly")
}
if (!all(spells[4,]==c(1,Inf,2,3,0,1)) | !all(spells[3,]==c(4,Inf,1,2,0,1))){
  stop("networkDynamic() did not record toggles correctly")
}

# check for net.obs.period
expect_equal((net%n%'net.obs.period')$observations[[1]],c(-Inf,Inf),info='net.obs.period created by default by networkDynamic(edge.toggles) did not have expected range')

expect_equal((net%n%'net.obs.period')$mode,'discrete',info='net.obs.period created by default by networkDynamic(edge.toggles) did not have expected mode')




# ==================== TESTING duration.matrix
# this function is used internally by as.networkDynamic.network() and should not be called by user
net <-network.initialize(3)
net[1,2]<-1;
net[2,3]<-1;
net[1,3]<-1;
# toggle list: time, tail, head
tog<-matrix(c(1,1,2, 1,2,3, 2,1,2, 4,1,3, 4,1,2), ncol=3, byrow=TRUE)
# we expect this matrix
ref <- matrix(c(0,1,1,2,1,0,1, 0,1,2,3,1,0,1, 0,4,1,3,1,0,4, 2,4,1,2,0,0,2),ncol=7,byrow=TRUE)
if (!all(networkDynamic:::duration.matrix(net, changes=tog, start=0, end=5)==ref)){
  stop("duration.matrix returned an unexpected spell list for its input toggles")
}
# testing start and end
ref1 <- matrix(c(1,1,1,2,1,0,0, 1,1,2,3,1,0,0, 1,4,1,3,1,0,3, 2,4,1,2,0,0,2),ncol=7,byrow=TRUE)
if (!all(networkDynamic:::duration.matrix(net, changes=tog, start=1, end=5)==ref1)){
  stop("duration.matrix returned an unexpected spell list for its input toggles")
}
# testing start and end
ref2 <- matrix(c(0,1,1,2,1,0,1, 0,1,2,3,1,0,1, 0,4,1,3,1,0,4, 2,4,1,2,0,0,2),ncol=7,byrow=TRUE)
if (!all(networkDynamic:::duration.matrix(net, changes=tog, start=0, end=8)==ref2)){
  stop("duration.matrix returned an unexpected spell list for its input toggles")
}

networkDynamic:::duration.matrix(network.initialize(0),changes=tog,start=0,end=1)



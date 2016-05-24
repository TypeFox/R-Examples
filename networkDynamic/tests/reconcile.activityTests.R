#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2013 the statnet development team
######################################################################
# tests for the reconcile activity functions
require(networkDynamic)
require(testthat)

# ---- reconcile.vertex.activity tests ------------
# test when called with wrong object
expect_error(reconcile.vertex.activity("net"), 'only be applied to networkDynamic objects')

# test with isolate and inactive edge
nd<-network.initialize(6)
add.edges.active(nd,tail=1:3,head=2:4,onset=1,terminus=3)
add.edges.active(nd,tail=4,head=1,onset=5,terminus=7)
add.edge(nd,tail=1,head=6)
deactivate.edges(nd, onset=-Inf, terminus=Inf, e=get.edgeIDs(nd, v=1, alter=6))
add.edge(nd,tail=1, head=5)  # default activated
nd2<-reconcile.vertex.activity(nd,mode='match')
spls<-get.vertex.activity(nd2,as.spellList=TRUE)
expect_equal(spls$onset,c(-Inf, 1,1,1,5,-Inf),info='test reconcile.vertex.activity w isolate and inactive')
expect_equal(spls$terminus,c(Inf,3,3,3,7,Inf))
expect_equal(spls$vertex.id,c(1,2,3,4,4,5))

# test modifiy in place
nd<-network.initialize(6)
add.edges.active(nd,tail=1:3,head=2:4,onset=1,terminus=3)
add.edges.active(nd,tail=4,head=1,onset=5,terminus=7)
add.edge(nd,tail=1,head=6)
deactivate.edges(nd, onset=-Inf, terminus=Inf, e=get.edgeIDs(nd, v=1, alter=6))
add.edge(nd,tail=1, head=5)  # default activated
reconcile.vertex.activity(nd,mode='match')
spls<-get.vertex.activity(nd,as.spellList=TRUE)
expect_equal(spls$onset,c(-Inf, 1,1,1,5,-Inf),info='reconcile.vertex.activity modify in place')
expect_equal(spls$terminus,c(Inf,3,3,3,7,Inf),info='reconcile.vertex.activity modify in place')
expect_equal(spls$vertex.id,c(1,2,3,4,4,5),info='reconcile.vertex.activity modify in place')

# test with edge active default FALSE
nd<-network.initialize(6)
add.edges.active(nd,tail=1:3,head=2:4,onset=1,terminus=3)
add.edges.active(nd,tail=4,head=1,onset=5,terminus=7)
add.edge(nd,tail=1,head=6)
deactivate.edges(nd, onset=-Inf, terminus=Inf, e=get.edgeIDs(nd, v=1, alter=6))
add.edge(nd,tail=1, head=5)  # default activated
nd2<-reconcile.vertex.activity(nd, mode='match',edge.active.default=FALSE)
spls<-get.vertex.activity(nd2,as.spellList=TRUE)
expect_equal(spls$onset,c(1, 5, 1, 1, 1, 5))
expect_equal(spls$terminus,c(3, 7, 3, 3, 3, 7))
expect_equal(spls$vertex.id,c(1, 1, 2, 3, 4, 4))

# test with bad mode
expect_error(reconcile.vertex.activity(nd,mode='foobar'),"'arg' should be one of")

# ---- reconcile.vertex.activity tests (expand.to.edges mode) ------------

nd<-network.initialize(6)
deactivate.vertices(nd, onset=-Inf, terminus=Inf)
activate.vertices(nd, v=2, onset=2, terminus=5)
add.edges.active(nd,tail=1:3,head=2:4,onset=1,terminus=3)
add.edges.active(nd,tail=4,head=1,onset=5,terminus=7)
add.edge(nd,tail=1,head=6)
deactivate.edges(nd, onset=-Inf, terminus=Inf, e=get.edgeIDs(nd, v=1, alter=6))
add.edge(nd,tail=1, head=5)  # default activated
nd2<-reconcile.vertex.activity(nd, mode='expand.to.edges')
spls<-get.vertex.activity(nd2,as.spellList=TRUE)
expect_equal(spls$onset,c(-Inf, 1,1,1,5,-Inf))
expect_equal(spls$terminus,c(Inf,5,3,3,7,Inf))
expect_equal(spls$vertex.id,c(1,2,3,4,4,5))

nd<-network.initialize(6)
deactivate.vertices(nd, onset=-Inf, terminus=Inf)
activate.vertices(nd, v=2, onset=5, terminus=7)
add.edges.active(nd,tail=1:3,head=2:4,onset=1,terminus=3)
add.edges.active(nd,tail=4,head=1,onset=5,terminus=7)
add.edge(nd,tail=1,head=6)
deactivate.edges(nd, onset=-Inf, terminus=Inf, e=get.edgeIDs(nd, v=1, alter=6))
add.edge(nd,tail=1, head=5)  # default activated
nd2<-reconcile.vertex.activity(nd, mode='expand.to.edges')
spls<-get.vertex.activity(nd2,as.spellList=TRUE)
expect_equal(spls$onset,c(-Inf, 1,5,1,1,5,-Inf))
expect_equal(spls$terminus,c(Inf,3,7,3,3,7,Inf))
expect_equal(spls$vertex.id,c(1,2,2,3,4,4,5))

# test active default

nd<-network.initialize(3,directed=FALSE)
activate.vertices(nd,at=0)
nd[,]<-1
reconcile.vertex.activity(nd,mode='expand',edge.active.default=TRUE)
spls<-get.vertex.activity(nd,as.spellList=TRUE)
expect_equal(spls$onset,c(-Inf,-Inf,-Inf))
expect_equal(spls$terminus,c(Inf,Inf,Inf))

nd<-network.initialize(3,directed=FALSE)
activate.vertices(nd,at=0)
nd[,]<-1
reconcile.vertex.activity(nd,mode='expand',edge.active.default=FALSE)
spls<-get.vertex.activity(nd,as.spellList=TRUE)
expect_equal(spls$onset,c(0,0,0))
expect_equal(spls$terminus,c(0,0,0))


# ---- reconcile.vertex.activity tests (encompass.edges mode) ------------

nd<-network.initialize(3)
activate.vertices(nd,v=1,onset=0,terminus=2)
activate.vertices(nd,v=1,onset=3,terminus=4)
add.edges.active(nd,tail=1,head=2,onset=-1,terminus=5)
reconcile.vertex.activity(nd, mode='encompass.edges')


nd<-network.initialize(6)
deactivate.vertices(nd, onset=-Inf, terminus=Inf)
activate.vertices(nd, v=2, onset=2, terminus=5)
activate.vertices(nd, v=2, onset=8, terminus=10)
add.edges.active(nd,tail=1:3,head=2:4,onset=1,terminus=3)
add.edges.active(nd,tail=4,head=1,onset=5,terminus=7)
add.edge(nd,tail=1,head=6)
deactivate.edges(nd, onset=-Inf, terminus=Inf, e=get.edgeIDs(nd, v=1, alter=6))
add.edge(nd,tail=1, head=5)  # default activated
nd2<-reconcile.vertex.activity(nd, mode='encompass.edges')
spls<-get.vertex.activity(nd2,as.spellList=TRUE)
# combines the spells for vertex 4
expect_equal(spls$onset,c(-Inf, 1,1,1,-Inf))
expect_equal(spls$terminus,c(Inf,3,3,7,Inf))
expect_equal(spls$vertex.id,c(1,2,3,4,5))

# testing 0-duration spells
nd<-network.initialize(6)
deactivate.vertices(nd, onset=-Inf, terminus=Inf)
activate.vertices(nd, v=2, onset=2, terminus=5)
add.edges.active(nd,tail=1:2,head=2:3,onset=1,terminus=3)

add.edges.active(nd,tail=4,head=1,onset=9,terminus=9)
add.edge(nd,tail=1,head=6)
deactivate.edges(nd, onset=-Inf, terminus=Inf, e=get.edgeIDs(nd, v=1, alter=6))
add.edge(nd,tail=1, head=5)  # default activated
activate.edges(nd, e=2, at=3)
#get.edge.activity(nd, as.spell=T)
#get.vertex.activity(nd, as.spell=T)
nd2<-reconcile.vertex.activity(nd, mode='encompass.edges')
spls<-get.vertex.activity(nd2,as.spellList=TRUE)
# combines the spells for vertex 4
expect_equal(spls$onset,c(-Inf, 1,1,9,-Inf))
expect_equal(spls$terminus,c(Inf,4,4,9,Inf))
expect_equal(spls$vertex.id,c(1,2,3,4,5))


# ---- reconcile.edge.activity tests (reduce.to.vertices mode) ------------
nd<- network.initialize(6)
activate.vertices(nd, onset=1, terminus=2)
add.edges.active(nd,tail=1:3,head=2:4,onset=1,terminus=3)
add.edges.active(nd,tail=4,head=1,onset=5,terminus=7)
add.edge(nd,tail=1,head=6)
deactivate.edges(nd, onset=-Inf, terminus=Inf, e=get.edgeIDs(nd, v=1, alter=6))
add.edge(nd,tail=1, head=5)  # default activated
#get.edge.activity(nd, as.spell=T)
#get.vertex.activity(nd, as.spell=T)
nd2<-reconcile.edge.activity(nd, mode='reduce.to.vertices')
spls<-get.edge.activity(nd2,as.spellList=TRUE)
expect_equal(spls$onset,c(1,1,1,1))
expect_equal(spls$terminus,c(2,2,2,2))
expect_equal(spls$edge.id,c(1,2,3,6))

# censored spells
nd<- network.initialize(6)
activate.vertices(nd, onset=2, terminus=Inf)
add.edges.active(nd,tail=1:3,head=2:4,onset=1,terminus=3)
add.edges.active(nd,tail=4,head=1,onset=5,terminus=7)
add.edge(nd,tail=1,head=6)
deactivate.edges(nd, onset=-Inf, terminus=Inf, e=get.edgeIDs(nd, v=1, alter=6))
add.edge(nd,tail=1, head=5)  # default activated
#get.edge.activity(nd, as.spell=T)
#get.vertex.activity(nd, as.spell=T)
nd2<-reconcile.edge.activity(nd, mode='reduce.to.vertices')
spls<-get.edge.activity(nd2,as.spellList=TRUE)
expect_equal(spls$onset,c(2,2,2,5,2))
expect_equal(spls$terminus,c(3,3,3,7,Inf))
expect_equal(spls$edge.id,c(1,2,3,4,6))

# 0-duration spells
nd<- network.initialize(6)
activate.vertices(nd, onset=1, terminus=2)
activate.vertices(nd, at=6)
add.edges.active(nd,tail=1:3,head=2:4,onset=1,terminus=3)
add.edges.active(nd,tail=4,head=1,onset=5,terminus=7)
add.edge(nd,tail=1,head=6)
deactivate.edges(nd, onset=-Inf, terminus=Inf, e=get.edgeIDs(nd, v=1, alter=6))
add.edge(nd,tail=1, head=5)  # default activated
#get.edge.activity(nd, as.spell=T)
#get.vertex.activity(nd, as.spell=T)
nd2<-reconcile.edge.activity(nd, mode='reduce.to.vertices')
spls<-get.edge.activity(nd2,as.spellList=TRUE)
expect_equal(spls$onset,c(1,1,1,6,1,6))
expect_equal(spls$terminus,c(2,2,2,6,2,6))
expect_equal(spls$edge.id,c(1,2,3,4,6,6))


# ---- reconcile.edge.activity tests (match.to.vertices mode) ------------
nd<- network.initialize(6)
activate.vertices(nd, onset=1, terminus=2)
add.edges.active(nd,tail=1:3,head=2:4,onset=1,terminus=3)
add.edges.active(nd,tail=4,head=1,onset=5,terminus=7)
add.edge(nd,tail=1,head=6)
deactivate.edges(nd, onset=-Inf, terminus=Inf, e=get.edgeIDs(nd, v=1, alter=6))
add.edge(nd,tail=1, head=5)  # default activated
get.edge.activity(nd, as.spell=T)
get.vertex.activity(nd, as.spell=T)
nd2<-reconcile.edge.activity(nd, mode='match.to.vertices')
spls<-get.edge.activity(nd2,as.spellList=TRUE)
expect_equal(spls$onset,c(1,1,1,1,1,1))
expect_equal(spls$terminus,c(2,2,2,2,2,2))
expect_equal(spls$edge.id,c(1,2,3,4,5,6))

# 0-duration spells
nd<- network.initialize(6)
activate.vertices(nd, onset=1, terminus=2, v=c(1,2,3))
deactivate.vertices(nd, v=4, onset=-Inf, terminus=Inf)
activate.vertices(nd, at=6, v=c(1,6))
add.edges.active(nd,tail=1:3,head=2:4,onset=1,terminus=3)
add.edges.active(nd,tail=4,head=1,onset=5,terminus=7)
add.edge(nd,tail=1,head=6)
add.edge(nd,tail=1, head=5)  # default activated
#get.edge.activity(nd, as.spell=T)
#get.vertex.activity(nd, as.spell=T)
nd2<-reconcile.edge.activity(nd, mode='match.to.vertices')
spls<-get.edge.activity(nd2,as.spellList=TRUE)
expect_equal(spls$onset,c(1,1,6,1,6))
expect_equal(spls$terminus,c(2,2,6,2,6))
expect_equal(spls$edge.id,c(1,2,5,6,6))

# test error from bad mode
expect_error(nd2<-reconcile.edge.activity(nd, mode='destroy.everything!'),"'arg' should be one of")

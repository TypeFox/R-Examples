# tests for the path distance functions

require(tsna)
require(testthat)
require(networkDynamicData)
linegraph<-network.initialize(10)
add.edges(linegraph,tail=1:9,head=2:10)

data(concurrencyComparisonNets)



# ----- tests for tPath ----



test_that('tPath basic tests',{
  line<-network.initialize(4)
  add.edges.active(line,tail=1:3,head=2:4,onset=0:2,terminus=1:3)
  # check return format
  expect_equal(names(tPath(line,v=1)),c('tdist','previous','gsteps','start','end','direction','type'))
  expect_is(tPath(line,v=1),class = 'tPath')
  
  # check args
  expect_error(tPath(line,v=1,type='foo'))
  expect_error(tPath(line,v=1,direction='foo'))
  
  # check unimplemented
  expect_error(tPath(line,v=1,type='latest.depart'),regexp='method is not yet implemented')
  
  # check basic line with default earliest arriving fwd
  expect_equal(tPath(line,v=1)$tdist,c(0, 0, 1, 2))
  expect_equal(tPath(line,v=2)$tdist,c(Inf,0,1,2))
  
  # test starting and ending flags
  expect_equal(tPath(line,v=1,start=0.5)$tdist, c(0,0,0.5,1.5))
  expect_equal(tPath(line,v=1,start=2)$tdist, c(0,Inf,Inf,Inf))
  expect_equal(tPath(line,v=1,end=2)$tdist, c(0,0,1,Inf))
  
  line<-network.initialize(4)
  add.edges.active(line,tail=1:3,head=2:4,onset=c(2,1,3),terminus=c(3,2,4))
  expect_equal(tPath(line,v=1)$tdist,c(0,1,Inf,Inf))
  
  # test active default
  test<-as.networkDynamic(network.initialize(4))
  add.edges(test,1:3,2:4)
  expect_equal(tPath(test,v=1,start=0)$tdist,c(0,0,0,0))
  expect_equal(tPath(test,v=1,active.default=FALSE,start=0)$tdist,c(0,Inf,Inf,Inf))
  
  test<-network.initialize(4)
  add.edges(test,1:3,3:4)
  activate.edges(test,e=1,at=2)
  
  # test start message
  test<-as.networkDynamic(network.initialize(4))
  expect_message(tPath(test,v=1),regexp="'start' time parameter for paths was not specified")
  
  # test wrong object
  expect_error(tPath(network.initialize(3)),regexp='first argument must be a networkDynamic object')
  
  # test no v specified
  expect_error(tPath(as.networkDynamic(network.initialize(2)),regexp='argument with valid vertex ids was not given'))
  
  # test on network size 0
  expect_equal(tPath(as.networkDynamic(network.initialize(0)),start=0,v=numeric(0))$tdist,numeric(0))
})


# ----- tests for forward earliest paths ----

test_that("path in large base network matches",{
  fwdDFS<-tPath(base,v=24)
  expect_equal(sum(fwdDFS$tdist<Inf),772) # should find 772 vertices, because that is what we found with BFS search
  
  # check if we can find the same set as the 'infected'
  infset<-which(get.vertex.attribute.active(base,'status',at=102)>0)
  pathset<-which(tPath(base,v=24,graph.step.time=1)$tdist<Inf)
})

data(moodyContactSim)


# tests with moody's example network
test_that("test of moody's example network",{
  
  paths<-tPath(moodyContactSim,v=10)
  
  expect_equal(paths$tdist,c(543, 454, 594,   0, 672, 661, 184, 679, 634,   0, 709, 581, 413, 625, 669, 535))
  expect_equal(paths$previous,c(16,13,13,10,13,16,10,13,1,0,8,1,4,4,2,2))
  expect_equal(paths$gsteps,c(5, 3, 3, 1, 3, 5, 1, 3, 6, 0, 4, 6, 2, 2, 4, 4))
  
  # render a pdf for visual inspection of correctness
  # tree<-create_tree(paths)
  # pdf(file="MoodyTestNetTree.pdf",width=10,height=5)
  # par(mfcol=c(1,2))
  # plot(moodyContactSim,displaylabels=TRUE,
  #      edge.label=lapply(get.edge.activity(moodyContactSim),
  #                         function(spl){
  #                           paste("(",spl[,1],"-",spl[,2],")",sep=
  #                                   ''
  #                           )
  #                         }),
  #      edge.label.col='blue',
  #      edge.label.cex=0.6,
  #      main="moody example net")
  # 
  # plot(tree,
  #      coord=layout.normalize(network.layout.animate.Graphviz(tree,layout.par=list(gv.engine='dot')),keep.aspect.ratio=FALSE),
  #      displaylabels=TRUE,
  #      jitter=FALSE,
  #      label.pos=2,
  #      main='earliest paths from v10',
  #      edge.label=lapply(get.edge.activity(tree),
  #                        function(spl){
  #                          paste("(",spl[,1],")",sep=
  #                                  ''
  #                          )
  #                        }),
  #      edge.label.col='blue',
  #      edge.label.cex=0.6)
  # par(mfcol=c(1,1))
  # dev.off()
})



test_that("test on network with two components",{
  test<-network.initialize(10)
  activate.vertices(test)
  test[1:5,5:1]<-1
  test[6:10,10:6]<-1
  expect_equal(which(tPath(test,v=1)$tdist!=Inf),1:5)
  expect_equal(which(tPath(test,v=6)$tdist!=Inf),6:10)
})



# test path distance
test_that("graph step time param works",{
  test<-network.initialize(4)
  add.edges.active(test,tail=1:3,head=2:4,onset=0:2,terminus=1:3)
  # count each geodesic step as instantaneous
  expect_equal(tPath(test,v=1,graph.step.time=0)$tdist,c(0, 0, 1, 2))
  # count each geodesic step as something less than 1
  expect_equal(tPath(test,v=1,graph.step.time=0.5)$tdist,c(0, 0.5, 1.5, 2.5))
  # count each geodesic step as 1
  expect_equal(tPath(test,v=1,graph.step.time=1)$tdist,c(0, 1, 2, 3))
  # count each geodesic step as 2
  expect_equal(tPath(test,v=1,graph.step.time=2)$tdist,c(0, Inf, Inf, Inf))
  
  # test with always active edges
  test<-network.initialize(4)
  add.edges.active(test,tail=1:3,head=2:4,onset=0,terminus=10)
  # count each geodesic step as 1
  expect_equal(tPath(test,v=1,graph.step.time=1)$tdist,c(0, 1, 2, 3))
  # count each geodesic step as 2
  expect_equal(tPath(test,v=1,graph.step.time=2)$tdist,c(0, 2, 4, 6))
  # count each geodesic step as 0
  expect_equal(tPath(test,v=1,graph.step.time=0)$tdist,c(0, 0, 0, 0))
  
  test<-as.networkDynamic(network.initialize(4))
  add.edges(test,tail=1:3,head=2:4)
  # count each geodesic step as 1
  expect_equal(tPath(test,v=1,graph.step.time=1)$tdist,c(0, 1, 2, 3))
  # count each geodesic step as 2
  expect_equal(tPath(test,v=1,graph.step.time=2)$tdist,c(0, 2, 4, 6))
  
  # test with an edge with multiple activity spells, some later
  test<-network.initialize(4)
  add.edges.active(test,tail=1:3,head=2:4,onset=0:2,terminus=1:3)
  activate.edges(test,e=1,onset=5,terminus=10)
  expect_equal(tPath(test,v=1,graph.step.time=2)$tdist,c(0, 7, Inf, Inf))
  
  # test with combination of start and end values
  test<-network.initialize(10)
  add.edges.active(test,tail=1:9,head=2:10,onset=0:9,terminus=10)
  tPath(test,v=1,start=5,graph.step.time=2)$tdist
  
})

# ----- tests for shortest geo fwd path ----
# construct a network in which distinguishes paths (
#
# tel<-matrix(c(2,3,1,2,
#               3,4,2,5,
#               4,5,5,3,
#               0,1,1,3,
#               5,6,1,4,
#               7,8,4,3),byrow=TRUE,ncol=4)
# # shortest path from 1 to 3 is the direct route, arriving at time
# test<-networkDynamic(edge.spells=tel)
# plot(test,displaylabels=TRUE,edge.label=paste(get.edge.attribute(test,'active',unlist=FALSE)),label.col='blue',edge.label.cex=0.7)
# 
# #path starting at time 0
# 
# path0<-tPath(test,v=1,type='fewest.steps')
# expect_equal(path0$tdist,c(0, 2, 0, 5, 3))
# expect_equal(path0$previous,c(0, 1, 1, 1, 2))
# expect_equal(path0$gsteps,c(0, 1, 1, 1, 2))
# 
# # path starting at v1 time 1
# # (direct path to v3 eliminated because it is too early)
# #   shortest geodesic path to v3 is  1-3
# #   earliest temporal path to v3 is 1-2-5-3
# #   shortest geodesic temporal path to v3 is 1-4-3
# path1<-tPath(test,v=1,start=1,type='fewest.steps')
# expect_equal(path1$tdist,c(0, 1, 6, 4, 2))
# expect_equal(path1$previous,c(0, 1, 4, 1, 2))
# expect_equal(path1$gsteps,c(0, 1, 2, 1, 2))
# 
# 
# # test for a later-leaving path arriving earlier


# ----- tests for paths.bkwd.latest -----



# reverse-ordered edge spells
test<-network.initialize(10)
add.edges(test,tail=1:9,head=2:10)
activate.edges(test,onset=10:0,terminus=11:1)
results<-tPath(test,v=5,direction='bkwd',type='latest.depart')
expect_equal(results$tdist,c(Inf, Inf, Inf,   3,   0, Inf, Inf, Inf, Inf, Inf))
expect_equal(results$previous,c(0, 0, 0, 5, 0, 0, 0, 0, 0, 0))
expect_equal(results$gsteps,c(Inf, Inf, Inf, 1, 0, Inf, Inf, Inf, Inf, Inf))

# forward-ordred edge spells
test<-network.initialize(10)
add.edges(test,tail=1:9,head=2:10)
activate.edges(test,onset=0:10,terminus=1:11)
results<-tPath(test,v=10,direction='bkwd',type='latest.depart')
expect_equal(results$tdist,c(8,7,6,5,4,3,2,1,0,0))
expect_equal(results$previous,c(2,3,4,5,6,7,8,9,10,0))
expect_equal(results$gsteps,c(9, 8, 7, 6, 5, 4, 3, 2, 1, 0))

# moody sim
results<-tPath(moodyContactSim,v=10,direction='bkwd',type='latest.depart')
expect_equal(results$tdist,c(Inf, Inf, Inf, 723, Inf, Inf, 539, Inf, Inf,   0, Inf, Inf, Inf, Inf, Inf, Inf))
expect_equal(results$previous,c(0,  0,  0, 10,  0,  0, 10,  0,  0,  0,  0,  0,  0,  0,  0,  0))

results<-tPath(moodyContactSim,v=16,direction='bkwd',type='latest.depart')
expect_equal(results$tdist,c(180, 196, Inf,  13, Inf,  62, Inf, Inf, Inf, 723, 548, Inf, 271, 103, Inf,   0))
expect_equal(results$previous,c(16, 16,  0, 16,  0, 16,  0,  0,  0,  4,  1,  0,  2,  4,  0,  0))

test_that("graph step time param works for bakward path",{
  test<-network.initialize(4)
  add.edges.active(test,tail=1:3,head=2:4,onset=0:2,terminus=1:3)
  # count each geodesic step as instantaneous
  expect_equal(tPath(test,v=4,graph.step.time=0,direction='bkwd',type='latest.depart')$tdist,c(2, 1, 0, 0))
  # count each geodesic step as something less than 1
  expect_equal(tPath(test,v=4,graph.step.time=0.5,direction='bkwd',type='latest.depart')$tdist,c(2.5, 1.5, 0.5, 0.0))
  # count each geodesic step as 1
  expect_equal(tPath(test,v=4,graph.step.time=1,direction='bkwd',type='latest.depart')$tdist,c(3, 2, 1, 0))
  # count each geodesic step as 2
  expect_equal(tPath(test,v=4,graph.step.time=2,direction='bkwd',type='latest.depart')$tdist,c( Inf, Inf, Inf,0))
  
  # test with always active edges
  test<-network.initialize(4)
  add.edges.active(test,tail=1:3,head=2:4,onset=0,terminus=10)
  # count each geodesic step as 1
  expect_equal(tPath(test,v=4,graph.step.time=1,direction='bkwd',type='latest.depart')$tdist,c(3, 2, 1, 0))
  # count each geodesic step as 2
  expect_equal(tPath(test,v=4,graph.step.time=2,direction='bkwd',type='latest.depart')$tdist,c(6, 4, 2, 0))
  # count each geodesic step as 0
  expect_equal(tPath(test,v=4,graph.step.time=0,direction='bkwd',type='latest.depart')$tdist,c(0, 0, 0, 0))
  
  test<-as.networkDynamic(network.initialize(4))
  add.edges(test,tail=1:3,head=2:4)
  # count each geodesic step as 1
  expect_equal(tPath(test,v=4,graph.step.time=1,direction='bkwd',type='latest.depart')$tdist,c(3, 2, 1, 0))
  # count each geodesic step as 2
  expect_equal(tPath(test,v=4,graph.step.time=2,direction='bkwd',type='latest.depart')$tdist,c(6, 4, 2, 0))
  
  # test with an edge with multiple activity spells, some later
  test<-network.initialize(4)
  add.edges.active(test,tail=1:3,head=2:4,onset=0:2,terminus=1:3)
  activate.edges(test,e=1,onset=5,terminus=10)
  expect_equal(tPath(test,v=4,graph.step.time=2,direction='bkwd',type='latest.depart')$tdist,c(Inf, Inf, 9, 0))
  
})


# --------- tests for tsna:::paths.fwd.latest ---------
# two paths, does it 
test<-network.initialize(2)
add.edges.active(test,tail=1,head=2,onset=0,terminus=1)
activate.edges(test,onset=2,terminus=3)
tPath(test,v=2,start=0,end=3)
tsna:::paths.fwd.latest(test,v=2,start=0,end=3)


# create a network in which the latest-starting path and
# the latest ending path are not the same
test<-network.initialize(5,direct=FALSE)
add.edges(test,tail=c(1,1,2,4),head=c(3,2,4,3))
activate.edges(test,at=c(1,2,3,4))
plot(test,displaylabels=TRUE,edge.label=get.edge.activity(test))
# latest starting path v1 to v4 should be at time 2 (via edge 2)
# latest ending path v1 to v4 should be at time 4 (via edge 3)
tsna:::paths.fwd.latest(test,v=1,start=0,end=4)


# create a network in which the latest-starting path and
# the latest ending path are not the same
test<-network.initialize(5,direct=FALSE)
add.edges(test,tail=c(1,1,2,4),head=c(3,2,4,3))
activate.edges(test,at=c(1,2,3,4))
plot(test,displaylabels=TRUE,edge.label=get.edge.activity(test))
tPath(test,v=1,start=0)
tsna:::paths.fwd.latest(test,v=1,start=0)

# create a network in which an early-leaving path arrives latest
# the latest path from v1 to v3 should arrive at t4 via v4
test<-network.initialize(4,directed=FALSE)
add.edges(test,tail=c(1,1,2,4),head=c(3,2,4,3))
activate.edges(test,at=c(2,1,3,4))
plot(test,displaylabels=TRUE,edge.label=get.edge.activity(test))



# test network illustrating problems with implementation of graph.step.time
test<-network.initialize(4)
add.edges.active(test,1:3,2:4,at=1)
tPath(test,v=1,start=0)$tdist  # all vertices rechable at t=1
tPath(test,v=1,start=0,graph.step.time = 1)$tdist  # no vertices are reachable because no path is open long enough for transmission to occur



# the network below illustrates the various possible paths
# tests ability to distinguish paths
# howver, it is not a great complex test case since there are no
# indirect paths 
paths5<-network.initialize(7)
network.vertex.names(paths5)<-LETTERS[1:7]
add.edges.active(paths5,tail=c(1,2),head=c(2,7),onset=c(1,4),terminus=c(2,5))
add.edges.active(paths5,tail=c(1,3),head=c(3,7),onset=c(0,6),terminus=c(2,7))
add.edges.active(paths5,tail=c(1,4),head=c(4,7),onset=c(4,5),terminus=c(5,6))
add.edges.active(paths5,tail=c(1,5),head=c(5,7),onset=c(6,9),terminus=c(7,10))
add.edges.active(paths5,tail=c(1,6),head=c(6,7),onset=c(4,10),terminus=c(5,11))
plot(paths5, mode='circle',displaylabels=TRUE,edge.label=get.edge.activity(paths5),edge.label.col='blue',edge.label.cex=0.6)
as.data.frame(paths5)
# FORWARDS
# earliest leaving ACG @ 6
# earliest arriving ABG @ 4
res2<-tPath(paths5,v=1)
expect_equal(res2$tdist[7],4)
# latest leaving AEG @ 10


# quickest ADG @ 5
# latest ariving  AFG @ 11

#BACKWARDS


# ----- tests for paths.fwd.approx ---
# data(moodyContactSim)
# # hard to test because random algorithm
# 
# split<-network.initialize(2)
# add.edges.active(split,tail=1,head=2,onset=0,terminus=1)
# arrivePercent<-paths.fwd.approx(split,v=1,tries=1000)
# # would expect 1/2,1/2
# 
# split<-network.initialize(2,directed=FALSE)
# add.edges.active(split,tail=1,head=2,onset=0,terminus=1)
# arrivePercent<-paths.fwd.approx(split,v=1,tries=1000)
# # would expect 1/2,1/2
# 
# split<-network.initialize(3)
# add.edges.active(split,tail=c(1,1),head=2:3,onset=0,terminus=1)
# arrivePercent<-paths.fwd.approx(split,v=1,tries=1000)
# # would expect 1/3,1/3,1/3
# 
# split<-network.initialize(5)
# add.edges.active(split,tail=c(1,1),head=2:3,onset=c(0,0.5),terminus=c(1,1))
# add.edges.active(split,tail=3,head=4,onset=0.5,terminus=1)
# # 
# plot(split,displaylabels=TRUE,edge.label=get.edge.activity(split))
# # 1/4 paths should remain on v1
# # 1/2 paths remain v2
# # 1/8 paths remain v3
# # 1/8 paths remain v4
# # 0 paths reach v5
# paths.fwd.approx(split,v=1,tries=100)

# for these values, set a rng seed
#set.seed(123)
#fwdProbs<-paths.fwd.approx(moodyContactSim,v=1)
#expect_equal(fwdProbs,c(0.858750, 0.000000, 0.000000, 0.000625, 0.000000, 0.001250, 0.000000, 0.001875, 0.036250, 0.000000, 0.023125, 0.038750, 0.000000, 0.000000, 0.000000, 0.039375))
#fwdProbs<-paths.fwd.approx(moodyContactSim,v=10,tries=100)
#expect_equal(fwdProbs,c(0, 0, 0, 0.01, 0, 0, 0.01, 0, 0, 0.98))


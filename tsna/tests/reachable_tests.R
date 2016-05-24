# test the forward and backward componetnt functions
require(tsna)
require(testthat)
require(networkDynamicData)
linegraph<-network.initialize(10)
add.edges(linegraph,tail=1:9,head=2:10)

# ---- forward.reachable tests -----
# test non timed version error
test<-linegraph
expect_error(forward.reachable(test,v=1),'must be a networkDynamic')

activate.edges(test,onset=0,terminus=3)

expect_equal(forward.reachable(test,v=1,per.step.depth=Inf),1:10)

# reverse-ordered edge spells
test<-linegraph
activate.edges(test,onset=10:0,terminus=11:1)
expect_equal(forward.reachable(test,v=5,per.step.depth=Inf),5:6)
expect_equal(forward.reachable(test,v=10,per.step.depth=Inf),10)

# forward-ordred edge spells
test<-linegraph
activate.edges(test,onset=0:10,terminus=1:11)
expect_equal(forward.reachable(test,v=5,per.step.depth=Inf),5:10)
expect_equal(forward.reachable(test,v=10,per.step.depth=Inf),10)



# test with two seeds
expect_equal(forward.reachable(test,v=c(1,5)),c(1,5,2,3,4,6,7,8,9,10))

# test on undirected case
test<-linegraph
set.network.attribute(test,'directed',FALSE)
activate.edges(test,onset=0,terminus=3)
expect_equal(forward.reachable(test,v=5,per.step.depth=Inf),c(5,4,6,3,7,2,8,1,9,10))

# test on network with bounded time
test<-linegraph
activate.edges(test,onset=0:10,terminus=1:11)
expect_equal(forward.reachable(test,v=5,start=4,end=6),c(5,6,7))
expect_equal(forward.reachable(test,v=1,end=6,per.step.depth=1),c(1,2,3,4,5,6,7))


# test with finite depthtest<-linegraph
test<-linegraph
activate.edges(test,onset=0,terminus=10)
expect_equal(forward.reachable(test,v=1,per.step.depth=2,end=3),1:7)
expect_equal(forward.reachable(test,v=1,per.step.depth=1.5,end=3),1:5)

test_that("test on network with two components",{
  test<-network.initialize(10)
  activate.vertices(test)
  test[1:5,5:1]<-1
  test[6:10,10:6]<-1
  expect_equal(forward.reachable(test,v=1),1:5)
  expect_equal(forward.reachable(test,v=6),6:10)
})


test_that("network with at spell durations",{
  test<-linegraph
  activate.edges(test,onset=0:10,terminus=0:10)
  expect_equal(forward.reachable(test,v=5,per.step.depth=Inf),5:10)
  expect_equal(forward.reachable(test,v=10,per.step.depth=Inf),10)
  
})



# test on network with net.obs.period set

# test on network with way too many time steps
data(hospital_contact)
#forward.reachable(hospital,v=1,start=120, end=347640,interval=300,per.step.depth=1)


# test that it matches infction in example network

data(concurrencyComparisonNets)
which(get.vertex.attribute.active(base,'status',at=1)==1)
# size of forward component in base from v 24
sum(get.vertex.attribute.active(base,'status',at=102)==1)

epiFound<-which(get.vertex.attribute.active(base,'status',at=102)==1)
fwdFound<-forward.reachable(base,v=24,per.step.depth=1,end=100)
expect_equal(length(setdiff(fwdFound,epiFound)),0)

# --- profiling -----#
# fiveRuns1<-function(){
#   forward.reachable1(base,v=1,per.step.depth=1)
#   forward.reachable1(base,v=2,per.step.depth=1)
#   forward.reachable1(base,v=3,per.step.depth=1)
#   forward.reachable1(base,v=4,per.step.depth=1)
#   forward.reachable1(base,v=5,per.step.depth=1)
# }
# fiveRuns2<-function(){
#   forward.reachable2(base,v=1,per.step.depth=1)
#   forward.reachable2(base,v=2,per.step.depth=1)
#   forward.reachable2(base,v=3,per.step.depth=1)
#   forward.reachable2(base,v=4,per.step.depth=1)
#   forward.reachable2(base,v=5,per.step.depth=1)
# }
# 
# # memory profiling
# Rprof(filename='fwdReachable.before')
# fiveRuns2()
# Rprof(NULL)
# summaryRprof(filename='fwdReachable.before')

# # time profiling
# library(microbenchmark)
# timing<-microbenchmark(fiveRuns1(),fiveRuns2(),times=1)

# ------ tests for the reachable set sizes ----

data(moodyContactSim)
expect_equal(sort(tReach(moodyContactSim)),c(2,  2,  3,  3,  3,  3,  5,  5,  6,  8, 10, 10, 14, 15, 16, 16))

expect_equal(sort(tReach(moodyContactSim,direction='bkwd')),c(3, 4, 5, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9))




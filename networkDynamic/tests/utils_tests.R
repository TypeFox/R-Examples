# tests for utilities functions
require(networkDynamic)
require(testthat)
# ------- test for adjust.activity -----

test<-network.initialize(3)
activate.vertices(test,onset=0,terminus=3,v=1:2)
add.edges.active(test,tail=1:2,head=2:3,onset=0,terminus=3)
add.edge(test,tail=3,head=1)
activate.vertex.attribute(test,'fruit','apple',v=1:2,onset=0,terminus=3)
activate.edge.attribute(test,'veggie','carrot',e=1:2,onset=0,terminus=3)
activate.network.attribute(test,'meat','pork',onset=0,terminus=3)
test%n%'net.obs.period'<-list(observations=list(c(0,1),c(1,2),c(2,3)),mode="discrete", time.increment=1,time.unit="step")

# test offset and return value
test2<-adjust.activity(test,offset=5)

expect_equal(test2$val[[1]]$active,matrix(c(5,8),ncol=2))
expect_equal(test2$val[[3]]$active,NULL)
expect_equal(test2$mel[[1]]$atl$active,matrix(c(5,8),ncol=2))
expect_equal(test2$mel[[3]]$atl$active,NULL)
expect_equal(test2$val[[1]]$'fruit.active'[[2]],matrix(c(5,8),ncol=2))
expect_equal(test2$val[[3]]$'fruit.active'[[2]],NULL)
expect_equal(test2$mel[[1]]$atl$'veggie.active'[[2]],matrix(c(5,8),ncol=2))
expect_equal(test2$mel[[3]]$atl$'veggie.active'[[2]],NULL)
expect_equal(test2$gal$'meat.active'[[2]],matrix(c(5,8),ncol=2))
expect_equal(unlist((test2%n%'net.obs.period')$observations),c(5,6,6,7,7,8))
expect_equal((test2%n%'net.obs.period')$time.increment,1)

# test factor and modify-in-place
adjust.activity(test,factor=.5)

expect_equal(test$val[[1]]$active,matrix(c(2.5,4),ncol=2))
expect_equal(test$val[[3]]$active,NULL)
expect_equal(test$mel[[1]]$atl$active,matrix(c(2.5,4),ncol=2))
expect_equal(test$mel[[3]]$atl$active,NULL)
expect_equal(test$val[[1]]$'fruit.active'[[2]],matrix(c(2.5,4),ncol=2))
expect_equal(test$val[[3]]$'fruit.active'[[2]],NULL)
expect_equal(test$mel[[1]]$atl$'veggie.active'[[2]],matrix(c(2.5,4),ncol=2))
expect_equal(test$mel[[3]]$atl$'veggie.active'[[2]],NULL)
expect_equal(test$mel[[1]]$atl$'veggie.active'[[2]],matrix(c(2.5,4),ncol=2))
expect_equal(unlist((test%n%'net.obs.period')$observations),c(2.5,3,3,3.5,3.5,4))
expect_equal((test%n%'net.obs.period')$time.increment,0.5)

# ---- test for add.vertices.active -----
net<-network.initialize(3)
# test for adding zero vertices
add.vertices.active(net,nv=0)
expect_equal(network.size(net),3)
expect_true(is.networkDynamic(net))

add.vertices.active(net,nv=2, onset=1,terminus=2)
expect_equal(network.size(net),5)
expect_true(is.networkDynamic(net))
expect_equal(unlist(get.vertex.activity(net,as.spellList=TRUE)[4:5,1:2]),c(1,1,2,2),check.names=FALSE)


# ---- tests for get.dyads.active ----
test_that("get.dyads.active works",{
  
  expect_error( get.dyads.active(network.initialize(3,hyper=TRUE),at=1),regexp="does not currently support hypergraphic",info="error on hyper")
  expect_equal( nrow(get.dyads.active(network.initialize(0),at=1)),0,info="network size zero case")
  expect_equal( nrow(get.dyads.active(network.initialize(3),at=1)),0,info="zero edges case")
  
  test<-network.initialize(5)
  add.edges.active(test,tail=1,head=2,onset=0,terminus=1)
  add.edges.active(test,tail=2,head=3,onset=1,terminus=2)
  add.edges.active(test,tail=3,head=4,onset=1,terminus=3)
  activate.edges(test,e=1,onset=2,terminus=3)
  as.data.frame(test)
  expect_equal(get.dyads.active(test,at=0),rbind(1:2))
  expect_equal(get.dyads.active(test,at=1),cbind(2:3,3:4))
  expect_equal(get.dyads.active(test,onset=0,terminus=4),rbind(1:2,2:3,3:4))
  
  # test with no (default) dynamics
  test2<-network.initialize(3)
  test2[1,2]<-1
  expect_equal(get.dyads.active(test2,onset=0,terminus=4),rbind(1:2))
  # test active default
  expect_equal(get.dyads.active(test2,onset=0,terminus=4,active.default=FALSE),cbind(list(),list()), info='test active default arg')
  
  # deleted edges
  test2<-network.initialize(3)
  test2[1,2]<-1
  test2[2,3]<-1
  test2[3,1]<-1
  delete.edges(test2,eid=2)
  expect_equal(get.dyads.active(test2,at=1),rbind(1:2,c(3,1)),info='deleted edge case')
  
})


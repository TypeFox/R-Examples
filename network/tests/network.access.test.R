library(network)

binet = network.initialize(10, bipartite = 6)
set.vertex.attribute(binet, 'myval', paste('b1', 1:6), v=1:6)
set.vertex.attribute(binet, 'myval', paste('b2', 1:4), v=7:10)

check <- vector()
check[1] <- all(get.vertex.attribute(binet, 'myval') == c("b1 1", "b1 2", "b1 3", "b1 4", "b1 5", "b1 6", "b2 1", "b2 2", "b2 3" ,"b2 4"))

# check for distinction between bipartite=FALSE and bipartite=0
testA<-network.initialize(3,bipartite=0)
if(!is.bipartite(testA)){
  stop('failed test of is.bipartite for bipartite=0')
}

testB<-network.initialize(3,bipartite=FALSE)
if(is.bipartite(testB)){
  stop('failed test of is.bipartite for bipartite=FALSE')
}

testC<-network.initialize(3,bipartite=TRUE)
if(!is.bipartite(testC)){
  stop('failed test of is.bipartite for bipartite=TRUE')
}

if(!is.bipartite(binet)){
  stop('failed test of is.bipartite for bipartite=6')
}

# add vertices to bipartite graphs
g = binet; add.vertices(g, 5, last.mode=F)
check[2] <- network.size(g) == 15
check[3] <- get.network.attribute(g, 'bipartite') == 11
check[4] <- identical(get.vertex.attribute(g, 'myval'),
  c("b1 1", "b1 2", "b1 3", "b1 4", "b1 5", "b1 6", NA,NA,NA,NA,NA,"b2 1","b2 2","b2 3","b2 4"))

test<-network.initialize(3,bipartite=0)
test%v%'letters'<-LETTERS[1:3]
add.vertices(test,nv=1,last.mode=FALSE)
if(!identical(test%v%'letters',c(NA,"A","B","C"))){
  stop("Error adding vertices to first mode of network with biparite=0")
}

test<-network.initialize(3,bipartite=0)
test%v%'letters'<-LETTERS[1:3]
add.vertices(test,nv=1,last.mode=TRUE)
if(!identical(test%v%'letters',c("A","B","C",NA))){
  stop("Error adding vertices to last mode of network with biparite=0")
}


g = binet
add.vertices(g, 5, last.mode=T)
check[5] <- network.size(g) == 15
check[6] <- get.network.attribute(g, 'bipartite') == 6
check[7] <- identical(get.vertex.attribute(g, 'myval'),
                  c("b1 1", "b1 2", "b1 3", "b1 4", "b1 5", "b1 6","b2 1","b2 2","b2 3","b2 4", NA,NA,NA,NA,NA))

# replacement operators should always replace
y <- network.initialize(4,dir=FALSE) # This network can have at most 1 edge.
y[1,2] <- NA # Assign NA to (1,2)
y[1,2] <- NA
check[8] <- network.edgecount(y) == 0
check[9] <- network.edgecount(y, na.omit=F) == 1

y[,] <- 1
check[10] <- network.edgecount(y) == 6
y[,] <- NA
check[11] <- network.edgecount(y) == 0
check[12] <- network.edgecount(y, na.omit=F) == 6
y[,] <- 0
check[13] <- network.edgecount(y, na.omit=F) == 0


# ------ test valid.eids function
net<-network.initialize(4)
net[,]<-1
delete.edges(net,eid=4:6)
if(!all(valid.eids(net)==c(1,2,3,7,8,9,10,11,12))){
  stop('valid.eids did not return correct ids for non-null elements of network')
}

#If everything worked, check is TRUE
if(!all(check)){                                               #Should be TRUE
  stop(paste("network package test failed on test(s):",which(!check)))
}
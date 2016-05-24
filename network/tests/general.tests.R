#The following battery of tests is intended to verify the functionality of
#the network library

library(network)

# ----- check assigning multiple attribute values in a single call ------
test<-network.initialize(3)
set.vertex.attribute(test,c('a','b'),c(1,2))
if(!all(test%v%'a'==c(1,1,1) & test%v%'b'==c(2,2,2))){
  stop('setting multiple attribute values with set.vertex.attribute failed')
}

test<-network.initialize(3)
set.vertex.attribute(test,list('a','b'),c(1,2))
if(!all(test%v%'a'==c(1,1,1) & test%v%'b'==c(2,2,2))){
  stop('setting multiple attribute values with set.vertex.attribute failed')
}

test<-network.initialize(3)
set.vertex.attribute(test,c('a','b'),list(c(1,2,3),c(4,5,6)))
if(!all(test%v%'a'==c(1,2,3) & test%v%'b'==c(4,5,6))){
  stop('setting multiple attribute values with set.vertex.attribute failed')
}

test<-network.initialize(3)
set.vertex.attribute(test,c('a','b'),list(list(1,2,3),list(4,5,6)))
if(!all(test%v%'a'==c(1,2,3) & test%v%'b'==c(4,5,6))){
  stop('setting multiple attribute values with set.vertex.attribute failed')
}

test<-network.initialize(3)
obj<-list(one='a complex object',two=c('with muliple','parts'))
set.vertex.attribute(test,c('a','b'),list(list(as.list(obj)),list(as.list(obj))))
if(!all(all.equal(get.vertex.attribute(test,'a',unlist=FALSE)[[1]],obj) & all.equal(get.vertex.attribute(test,'b',unlist=FALSE)[[1]],obj))){
  stop('setting multiple attribute values with list values in set.vertex.attribute failed')
}

# check assignment to list of networks
net <- network.initialize(2)
netlist <- list(net)
set.network.attribute(netlist[[1]],"test","a value")
if (!"test" %in% list.network.attributes(netlist[[1]]))
  stop('assignment to list of networks failed')

# test multiple assignment for network

test<-network.initialize(3)
set.network.attribute(test,c("a","b"),1:2)
if (!all(test%n%'a'==1,test%n%'b'==2)){
  stop('mulltiple attribute assignment failed for set.network.attribute')
}

test<-network.initialize(3)
set.network.attribute(test,list("a","b"),as.list(1:2))
if (!all(test%n%'a'==1,test%n%'b'==2)){
  stop('mulltiple attribute assignment failed for set.network.attribute')
}



#  test multiple assignment for edges 

test<-network.initialize(3)
add.edges(test,tail=1:3,head=c(2,3,1))
net<-test
set.edge.attribute(net,c("a","b"),1:2)
if (!all(net%n%'a'==1,net%n%'b'==2)){
  stop('mulltiple attribute assignment failed for set.edge.attribute')
}

net<-test
set.edge.attribute(net,c('a','b'),list(c(1,2,3),c(4,5,6)))
if(!all(net%e%'a'==c(1,2,3) & net%e%'b'==c(4,5,6))){
  stop('setting multiple attribute values with set.edge.attribute failed')
}

net<-test
set.edge.attribute(net,c('a','b'),list(list(1,2,3),list(4,5,6)))
if(!all(net%e%'a'==c(1,2,3) & net%e%'b'==c(4,5,6))){
  stop('setting multiple attribute values with set.edge.attribute failed')
}

net<-test
obj<-list(one='a complex object',two=c('with muliple','parts'))
set.edge.attribute(net,c('a','b'),list(list(as.list(obj)),list(as.list(obj))))
if(!all(all.equal(get.edge.attribute(net,'a',unlist=FALSE)[[1]],obj) & all.equal(get.edge.attribute(net,'b',unlist=FALSE)[[1]],obj))){
  stop('setting multiple attribute values with list values in set.edge.attribute failed')
}



# ---- checks for  get.edge.attribute overloading and omit args ----
net<-network.initialize(3)
add.edges(net,c(1,2,3),c(2,3,1))
set.edge.attribute(net,'test',"a")
if(!all(get.edge.attribute(net,'test')==c("a","a","a"))){stop("overloading of get.edge.attribute to get.edge.value not working correctly ")}

# check list output of get.edge.attribute with deleted.edges.omit
delete.edges(net,2)
set.edge.attribute(net,'foo','bar',1)
if(!identical(list('bar',NULL,NULL),get.edge.attribute(net,'foo',unlist=FALSE,  deleted.edges.omit = FALSE))){
  stop("deleted.edges.omit argument causing bad return values in get.edge.attribute ")
}
if(!identical(list('bar',NULL),get.edge.attribute(net,'foo',unlist=FALSE,  deleted.edges.omit = TRUE))){
  stop("deleted.edges.omit argument causing bad return values in get.edge.attribute ")
}

# check unlisted output of get.edge.attribute with na.omit and deleted.edges.omit
if(!identical(c('bar'),get.edge.attribute(net,'foo',unlist=TRUE,deleted.edges.omit=TRUE))){
  stop("omission argument causing bad return values in get.edge.attribute")
}
if(!identical(c('bar'),get.edge.attribute(net,'foo',unlist=TRUE,deleted.edges.omit=TRUE))){
  stop("omission  arguments causing bad return values in get.edge.attribute")
}

# check for null.na recoding of non-set attributes
if(!identical(c('bar'),get.edge.attribute(net,'foo',unlist=TRUE,null.na=FALSE))){
  stop("null.na  arguments causing bad return values in get.edge.attribute")
}
if(!identical(c('bar',NA),get.edge.attribute(net,'foo',unlist=TRUE,null.na=TRUE))){
  stop("null.na  arguments causing bad return values in get.edge.attribute")
}
if(!identical(list('bar',NULL,NULL),get.edge.attribute(net,'foo',unlist=FALSE,null.na=FALSE))){
  stop("null.na  arguments causing bad return values in get.edge.attribute")
}
if(!identical(list('bar',NULL,NA),get.edge.attribute(net,'foo',unlist=FALSE,null.na=TRUE))){
  stop("null.na  arguments causing bad return values in get.edge.attribute")
}



#mark an edge as missing to test na.omit
set.edge.attribute(net,'na',TRUE,e=1)

# check that values corresponding to missing edges are ommited
if(!identical(list('bar',NULL,NULL),get.edge.attribute(net,'foo',unlist=FALSE,na.omit=FALSE))){
  stop("na.omit argument causing bad return values in get.edge.attribute")
}
if(!identical(list(NULL,NULL),get.edge.attribute(net,'foo',unlist=FALSE,na.omit=TRUE))){
  stop("na.omit argument causing bad return values in get.edge.attribute")
}

if(!identical(c('bar'),get.edge.attribute(net,'foo',unlist=TRUE,na.omit=FALSE))){
  stop("na.omit argument causing bad return values in get.edge.attribute")
}
if(!identical(NULL,get.edge.attribute(net,'foo',unlist=TRUE,na.omit=TRUE))){
  stop("na.omit argument causing bad return values in get.edge.attribute")
}
# check for behavior when querying the 'na' attribute
if(!identical(c(TRUE,FALSE),get.edge.attribute(net,'na',na.omit=FALSE))){
  stop("get.edge.attribute did not return correct values for 'na' attribute with na.omit=FALSE")
}
if(!identical(c(FALSE),get.edge.attribute(net,'na',na.omit=TRUE))){
  stop("get.edge.attribute did not return correct values for 'na' attribute with na.omit=TRUE")
}

# check behavior on a network with no edges
if(!identical(list(),get.edge.attribute(network.initialize(3),'foo',unlist=FALSE))){
  stop("get.edge.attribute did not return correct values network with no edges")
}

if(!identical(NULL,get.edge.attribute(network.initialize(3),'foo',unlist=TRUE))){
  stop("get.edge.attribute did not return correct values network with no edges")
}

if(!identical(NULL,get.edge.attribute(net,'bar'))){
  stop("get.edge.attribute did not return correct values for attribute that does not exist")
}


# check for behavior of attribute values explicitly set to null
net<-network.initialize(3)
net[1,2]<-1
net[1,3]<-1
set.edge.attribute(net,'nullval',list(NULL))

# expect NULL,NULL
if(!identical(list(NULL,NULL),get.edge.attribute(net,'nullval',unlist=FALSE,null.na=FALSE))){
  stop("get.edge.attribute not returning NULL values stored as edge attribute correctly")
}

# expect that this should return NULL values, which will dissappear on unlisting
# do NOT want to see NA,NA
if(!identical(NULL,get.edge.attribute(net,'nullval',null.na=FALSE))){
  stop("get.edge.attribute not returning NULL values stored as edge attribute correctly")
}
if(!identical(NULL,get.edge.attribute(net,'nullval',null.na=TRUE))){
  stop("get.edge.attribute not returning NULL values stored as edge attribute correctly")
}


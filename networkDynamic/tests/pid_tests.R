#  Part of the statnet package, http://statnetproject.org
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) in
#    http://statnetproject.org/attribution
#
#  Copyright 2013 the statnet development team
######################################################################

require(networkDynamic)
require(testthat)

#Create a network with three edges
m<-matrix(0,3,3)
m[1,2]<-1; m[2,3]<-1; m[3,1]<-1
g<-network(m)

#Create a matrix of values corresponding to edges
mm<-m
mm[1,2]<-7; mm[2,3]<-4; mm[3,1]<-2

#Assign some attributes
set.edge.attribute(g,"myeval",3:5)
set.edge.value(g,"myeval2",mm)
set.network.attribute(g,"mygval","boo")
set.vertex.attribute(g,"myvval",letters[1:3])
network.vertex.names(g) <- LETTERS[1:10]

set.network.attribute(g, 'vertex.pid', 'vertex.names')


# check if some pid missing

#------ get.vertex.id checks ------

expect_equal(get.vertex.id(g, 'A'),1)
expect_equal(get.vertex.id(g, 'B'),2)
# returns NA if not found
expect_true(is.na(get.vertex.id(g, 'D')))

# check multiple works
expect_equal(get.vertex.id(g, c('B','C','D')),c(2,3,NA))

expect_error(get.vertex.id(network.initialize(5),"does not have a 'vertex.pid' attribute"))

expect_error(get.vertex.id(network.initialize(0),"does not have a 'vertex.pid' attribute"))


#------- get.vertex.pid checks ---

expect_equal(get.vertex.pid(g, 2),'B')
# returns NA if not found
expect_true(is.na(get.vertex.pid(g, 5)))

# multiple works
expect_equal(get.vertex.pid(g,c(2,3,4)),c("B","C",NA))

# check when no vertex pid speced
expect_error(get.vertex.pid(network.initialize(5),"does not have a 'vertex.pid' attribute"))

expect_error(get.vertex.pid(network.initialize(0),"does not have a 'vertex.pid' attribute"))


# ---------- add.vertices  checks -----

# test calling original function, direct assignment
net <- as.networkDynamic(network.initialize(1))
net <-add.vertices(net,nv=3)
expect_equal(network.size(net),4,info='add.vertices direct assignment')

# test calling original function, modify inplace
net <- as.networkDynamic(network.initialize(1))
add.vertices(net,nv=3)
expect_equal(network.size(net),4,info='add.vertices modify in place')


net <- as.networkDynamic(network.initialize(1))
set.network.attribute(net,'vertex.pid','data_id')
set.vertex.attribute(net,'data_id','one')

# adding wrong number of ids gives error
expect_error(add.vertices(net,4,vertex.pid=c('two','three','four')), info='does not match number of new vertices')

# adding duplicate ids gives error

expect_error(add.vertices(net,4,vertex.pid=c('two','three','three','three')), info='vertex.pid values must be unique')

# error did not modify network
expect_equal(network.size(net),1)

# adding correctly 
add.vertices(net,4,vertex.pid=c('two','three','four','five'))
expect_equal(network.size(net),5, info='add.vertices check verts added')
expect_equal(net%v%'data_id',c("one","two","three","four","five" ),info='add.vertices added vertex.pids')

# adding with no specified pid
add.vertices(net,3)
expect_equal(anyDuplicated(get.vertex.attribute(net,'data_id')),0)

# adding with pid disabled
net<-as.networkDynamic(network.initialize(3))
set.network.attribute(net,'vertex.pid',NULL)
expect_equal(network.size(add.vertices(net,3)),6)

# adding to net of size 0

expect_equal(network.size(add.vertices(network.initialize(0),1)),1)

# ------------ add.edges checks ----

# no pid defined, modify in place
nd<-as.networkDynamic(network.initialize(3))
add.edges(nd,tail=1:3,head=c(2,3,1))
expect_equal(network.edgecount(nd),3)

# direct assignement
nd<-as.networkDynamic(network.initialize(3))
nd2<-add.edges(nd,tail=1:3,head=c(2,3,1))
expect_equal(network.edgecount(nd2),3)

# pid defined
nd<-as.networkDynamic(network.initialize(3))
set.network.attribute(nd,'edge.pid','myFavoriteId')
add.edges(nd,tail=1:3,head=c(2,3,1))
expect_true(nd%n%'edge.pid'=='myFavoriteId')
expect_equal(length(get.edge.attribute(nd,'myFavoriteId')),3,info='check add.edges created pids for edges')
expect_true('myFavoriteId'%in%list.edge.attributes(nd),info='check add.edges created edge.pid with correct name')

# adding to net with edges, and passing 
nd<-as.networkDynamic(network.initialize(3))
add.edges(nd,tail=1:3,head=c(2,3,1))
set.network.attribute(nd,'edge.pid','edge.pid')
set.edge.attribute(nd,'edge.pid',c("A","B","C"))
add.edges(nd,tail=3,head=1,edge.pid="D")
add.edges(nd,tail=3,head=2)

expect_equal(length(get.edge.attribute(nd,'edge.pid')),5)
expect_equal(get.edge.attribute(nd,'edge.pid')[1:4],LETTERS[1:4])
             
# ------------ add.edge checks ----
             
# no pid defined, modify in place
nd<-as.networkDynamic(network.initialize(3))
add.edge(nd,tail=1,head=2)
expect_equal(network.edgecount(nd),1)
             
# direct assignement
nd<-as.networkDynamic(network.initialize(3))
nd2<-add.edge(nd,tail=1,head=2)
expect_equal(network.edgecount(nd2),1)
             
# pid defined
nd<-as.networkDynamic(network.initialize(3))
set.network.attribute(nd,'edge.pid','myFavoriteId')
add.edge(nd,tail=1,head=2)
expect_true(nd%n%'edge.pid'=='myFavoriteId')
expect_equal(length(get.edge.attribute(nd,'myFavoriteId')),1,info='check add.edge created pids for edges')
expect_true('myFavoriteId'%in%list.edge.attributes(nd),info='check add.edge created edge.pid with correct name')
             
# adding to net with edges, and passing 
nd<-as.networkDynamic(network.initialize(3))
add.edges(nd,tail=1:3,head=c(2,3,1))
set.network.attribute(nd,'edge.pid','edge.pid')
set.edge.attribute(nd,'edge.pid',c("A","B","C"))
add.edge(nd,tail=3,head=1,edge.pid="D")
add.edge(nd,tail=3,head=2)
             
expect_equal(length(get.edge.attribute(nd,'edge.pid')),5)
expect_equal(get.edge.attribute(nd,'edge.pid')[1:4],LETTERS[1:4])             
             
             
# check error for non-unique
nd<-as.networkDynamic(network.initialize(3))
set.network.attribute(nd,'edge.pid','edge.pid')
expect_error(add.edge(nd,tail=1,head=2,edge.pid=c("A","A","A")), 'Only one edge.pid can be specified')   
             
# check for errror from existign non-unique
nd<-as.networkDynamic(network.initialize(3))
set.network.attribute(nd,'edge.pid','edge.pid')  
add.edges(nd,tail=1:3,head=c(2,3,1))
set.edge.attribute(nd,'edge.pid',"A")             
expect_error(add.edge(nd,tail=3,head=1,edge.pid="B"),"edge.pid attribute must be specified and unique for each edge")
             
             

# ---- intitialize.pids ----

test<-as.networkDynamic(network.initialize(30))
add.edges(test,1:29,2:30)
initialize.pids(test)
expect_equal(anyDuplicated(get.vertex.attribute(test,'vertex.pid')),0)
expect_equal(anyDuplicated(get.edge.attribute(test,'edge.pid')),0)

initialize.pids(network.initialize(0))

# ----- get.edge.id ----------------
net<-as.networkDynamic(network.initialize(5))
add.edges(net,1:4,2:5)
set.edge.attribute(net,'data_id',LETTERS[1:4])
set.network.attribute(net,'edge.pid','data_id')
expect_equal(get.edge.id(net,c("B","D")),c(2,4))

# error if not defined
expect_error(get.edge.id(network.initialize(4)),"does not have an 'edge.pid' attribute")

# NA if not existing

expect_true(is.na(get.edge.id(net,"L")))

# ----- get.edge.pid -----
expect_equal(get.edge.pid(net,c(1,4)),c("A","D"))

# NA if out of range
expect_true(identical(get.edge.pid(net,c(4,5)),c("D",NA)))

# error if not defined
expect_error(get.edge.pid(network.initialize(4)),"does not have an 'edge.pid' attribute")

expect_error(get.edge.pid(network.initialize(0)),"does not have an 'edge.pid' attribute")
             
             


#----- edge.pid.check  checks ---------
nd <-as.networkDynamic(network.initialize(5))
add.edges(nd,1:4,2:5)
set.edge.attribute(nd,"myId",LETTERS[1:4])
set.network.attribute(nd,'edge.pid','myId')

expect_true(edge.pid.check(nd))

# missing
delete.edge.attribute(nd,"myId")
expect_error(edge.pid.check(nd),'Missing edge.pids')

# partially missing
set.edge.attribute(nd,"myId",LETTERS[1:3],e=1:3)
expect_error(edge.pid.check(nd),'must be specified and unique')


# not unique
set.edge.attribute(nd,"myId","a")
expect_error(edge.pid.check(nd),'must be specified and unique')

# not defined
expect_warning(edge.pid.check(network.initialize(2)),"does not have an 'edge.pid' attribute")

expect_warning(edge.pid.check(network.initialize(0)),"does not have an 'edge.pid' attribute")



# ----- vertex.pid.check checks ------
nd <-as.networkDynamic(network.initialize(5))
set.vertex.attribute(nd,"myId",LETTERS[1:5])
set.network.attribute(nd,'vertex.pid','myId')
expect_true(vertex.pid.check(nd),info='checking correctly formatted edge.pid')

nd <-as.networkDynamic(network.initialize(5))
set.vertex.attribute(nd,"myId",LETTERS[1:4],v=1:4)
set.network.attribute(nd,'vertex.pid','myId')
expect_error(vertex.pid.check(nd),info='error for mis-formatted vertex.pid')

expect_warning(vertex.pid.check(network.initialize(3)),"does not have a 'vertex.pid' attribute")

expect_warning(vertex.pid.check(network.initialize(0)),"does not have a 'vertex.pid' attribute")

# ------- extraction check ----
nd <-as.networkDynamic(network.initialize(5))
set.vertex.attribute(nd,"myId",LETTERS[1:5])
set.network.attribute(nd,'vertex.pid','myId')
activate.vertices(nd,onset=c(1,2,3,4,5),terminus=c(3,4,5,6,7))
n3<-network.extract(nd,at=3)
expect_equal(get.vertex.attribute(n3,'myId'),c("B","C"))
expect_equal(get.vertex.id(n3,c("A","B","C")),c(NA,1,2))


# find vertex corresponding to extracted vertex
haystack<-network.initialize(30)
activate.vertices(haystack,v=10:20)
# hide a needle somewhere in the haystack
set.vertex.attribute(haystack,'needle',TRUE,v=10)
initialize.pids(haystack)
# some hay is removed over time ...
newstack<-network.extract(haystack,at=100,active.default=FALSE)
# we find the needle!
needleId <-which(get.vertex.attribute(newstack,'needle'))

# which vertex is the corresponding one in original stack?
oldId<-get.vertex.id(haystack,get.vertex.pid(newstack,needleId))

expect_true(get.vertex.attribute(haystack,'needle')[oldId],info="find vertex corresponding to extracted vertex")





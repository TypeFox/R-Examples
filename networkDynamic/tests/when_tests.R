# tests for the when functions
require(testthat)
require(networkDynamic)

#------ when.vertex.attrs.match tests ----- 

test1<-network.initialize(5)
test1<-activate.vertex.attribute(test1,'myValue',1:5,onset=1,terminus=2)
test1<-activate.vertex.attribute(test1,'myValue',10:15,onset=2,terminus=3)
test1<-activate.vertex.attribute(test1,'myValue',-5:0,onset=5,terminus=12)

# basic test
expect_equal(when.vertex.attrs.match(test1,attrname='myValue',value=3),c(Inf,Inf,1,Inf,Inf))
expect_equal(when.vertex.attrs.match(test1,attrname='myValue',value=-5),c(5,Inf,Inf,Inf,Inf))
expect_equal(when.vertex.attrs.match(test1,attrname='myValue',value=200),c(Inf,Inf,Inf,Inf,Inf))

# changing no.match value
expect_equal(when.vertex.attrs.match(test1,attrname='myValue',value=200,no.match=NA),c(NA,NA,NA,NA,NA))

# changing operation
expect_equal(when.vertex.attrs.match(test1,attrname='myValue',value=3,match.op='>'),c(2,2,2,1,1))

# complex operation, complex value
expect_equal(when.vertex.attrs.match(test1,attrname='myValue',value=c(1,2,3),match.op='%in%'),c(1,1,1,Inf,Inf))


# bad type of operation
expect_error(when.vertex.attrs.match(test1,attrname='myValue',value=3,match.op='min'),regexp='provide Logical results for every attribute value')

# not a networkDynamic
expect_error(when.vertex.attrs.match(5,attrname='myValue',value=3),regexp='argument to be a networkDynamic object')

# non-existing attrname
expect_equal(when.vertex.attrs.match(test1,attrname='foo',value=3),c(NA,NA,NA,NA,NA))

# missing attrname
expect_error(when.vertex.attrs.match(test1,value=3),regexp="providing an 'attrname' argument")

# missing value
expect_error(when.vertex.attrs.match(test1,attrname='myValue'),regexp="requires providing a 'value' argument")

# use v argument to query specific vertices
expect_equal(when.vertex.attrs.match(test1,attrname='myValue',value=3,match.op='>',v=3:4),c(2,1))

# bad rule
expect_error(when.vertex.attrs.match(test1,attrname='myValue',value=3,rule='guess'), regexp="no matching methods implemented for rule 'guess'")

# latest rule
expect_equal(when.vertex.attrs.match(test1,attrname='myValue',value=3, match.op='>',rule='latest'),c(3,3,3,3,3))


# what if attributes are missing from some elements

test2<-network.initialize(5)
test2<-activate.vertex.attribute(test2,'myValue',c(1,2,5),onset=1,terminus=2,v=c(1,2,5))
test2<-activate.vertex.attribute(test2,'myValue',c(10,11,15),onset=2,terminus=3,v=c(1,2,5))
test2<-activate.vertex.attribute(test2,'myValue',c(-5,-4,-1),onset=5,terminus=12,v=c(1,2,5))
expect_equal(when.vertex.attrs.match(test2,attrname='myValue',value=3,match.op='<'),c(1,1,NA,NA,5))


#------ when.edge.attrs.match tests ----- 

test1<-network.initialize(5)
add.edges(test1,tail=1:5,head=c(2,3,4,5,1))
test1<-activate.edge.attribute(test1,'myValue',1:5,onset=1,terminus=2)
test1<-activate.edge.attribute(test1,'myValue',10:15,onset=2,terminus=3)
test1<-activate.edge.attribute(test1,'myValue',-5:0,onset=5,terminus=12)

# basic test
expect_equal(when.edge.attrs.match(test1,attrname='myValue',value=3),c(Inf,Inf,1,Inf,Inf))
expect_equal(when.edge.attrs.match(test1,attrname='myValue',value=-5),c(5,Inf,Inf,Inf,Inf))
expect_equal(when.edge.attrs.match(test1,attrname='myValue',value=200),c(Inf,Inf,Inf,Inf,Inf))

# changing no.match value
expect_equal(when.edge.attrs.match(test1,attrname='myValue',value=200,no.match=NA),c(NA,NA,NA,NA,NA))

# changing operation
expect_equal(when.edge.attrs.match(test1,attrname='myValue',value=3,match.op='>'),c(2,2,2,1,1))

# complex operation, complex value
expect_equal(when.edge.attrs.match(test1,attrname='myValue',value=c(1,2,3),match.op='%in%'),c(1,1,1,Inf,Inf))


# bad type of operation
expect_error(when.edge.attrs.match(test1,attrname='myValue',value=3,match.op='min'),regexp='provide Logical results for every attribute value')

# not a networkDynamic
expect_error(when.edge.attrs.match(5,attrname='myValue',value=3),regexp='argument to be a networkDynamic object')

# non-existing attrname
expect_equal(when.edge.attrs.match(test1,attrname='foo',value=3),c(NA,NA,NA,NA,NA))

# missing attrname
expect_error(when.edge.attrs.match(test1,value=3),regexp="providing an 'attrname' argument")

# missing value
expect_error(when.edge.attrs.match(test1,attrname='myValue'),regexp="requires providing a 'value' argument")

# use e argument to query specific vertices
expect_equal(when.edge.attrs.match(test1,attrname='myValue',value=3,match.op='>',e=3:4),c(2,1))

# bad rule
expect_error(when.edge.attrs.match(test1,attrname='myValue',value=3,rule='guess'), regexp="no matching methods implemented for rule 'guess'")

# latest rule
expect_equal(when.edge.attrs.match(test1,attrname='myValue',value=3, match.op='>',rule='latest'),c(3,3,3,3,3))

# what if some edges are deleted
test1<-delete.edges(test1,e=c(2,3))
expect_equal(when.edge.attrs.match(test1,attrname='myValue',value=3),c(Inf,NA,NA,Inf,Inf))

# what if attributes are missing from some elements
test2<-network.initialize(5)
add.edges(test2,tail=1:5,head=c(2,3,4,5,1))
test2<-activate.edge.attribute(test2,'myValue',c(1,2,5),onset=1,terminus=2,e=c(1,2,5))
test2<-activate.edge.attribute(test2,'myValue',c(10,11,15),onset=2,terminus=3,e=c(1,2,5))
test2<-activate.edge.attribute(test2,'myValue',c(-5,-4,-1),onset=5,terminus=12,e=c(1,2,5))
expect_equal(when.edge.attrs.match(test2,attrname='myValue',value=3,match.op='<'),c(1,1,NA,NA,5))

# ---- when.next.edge.change ----
# test<-network.initialize(5)
# add.edges.active(test,tail=1:5,head=c(2,3,4,5,1),onset=1:5,terminus=2:6)
# activate.edges(test,e=2,onset=6,terminus=100)
# add.edge(test,tail=4,head=5)
# expect_equal(when.next.edge.change(test,at=0),1)
# expect_equal(when.next.edge.change(test,at=2),3)
# expect_equal(when.next.edge.change(test,at=2,v=1),Inf)
# # test a query in the gap
# expect_equal(when.next.edge.change(test,at=3,v=2),6)  




#The following battery of tests is intended to verify the functionality of
#the network library

library(network)
#These functions are intended to mimic functionality from the sna package.
#Said package is not required to use network, but was used in creating this
#battery of tests.
rgraph<-function(n){
  m<-matrix(rbinom(n*n,1,0.5),n,n)
  diag(m)<-0
  m
}
degree<-function(d,cmode = "freeman")
{
  n <- dim(d)[1]
  diag(d) <- NA
  switch(cmode, indegree = apply(d, 2, sum, na.rm = TRUE),
    outdegree = apply(d, 1, sum, na.rm = TRUE), freeman = apply(d,
    2, sum, na.rm = TRUE) + apply(d, 1, sum, na.rm = TRUE))
}
#gctorture(TRUE)     #Uncomment to perform a more intensive (SLOW) test

# ---- Check assignment, deletion, and adjacency for dyadic graphs ----
check<-vector()
temp<-network(matrix(0,5,5))
temp[1,2]<-1                 #Add edge
check[1]<-temp[1,2]==1       #Check adjacency
check[2]<-get.network.attribute(temp,"mnext")==2  #Check count
temp[1,2]<-1                 #Should have no effect
check[3]<-get.network.attribute(temp,"mnext")==2  #Check count
temp[1,1]<-1                 #Should have no effect
check[4]<-temp[1,1]==0       #Shouldn't be present
check[5]<-get.network.attribute(temp,"mnext")==2  #Check count
temp[,2]<-1                  #Should add 3 edges
check[6]<-get.network.attribute(temp,"mnext")==5  #Check count
check[7]<-all(temp[,2]==c(1,0,1,1,1))  #Verify row
temp[3:4,3:4]<-1             #Should add 2 edges
check[8]<-get.network.attribute(temp,"mnext")==7  #Check count
temp[,]<-0                   #Delete edges
check[9]<-all(temp[,]==matrix(0,5,5))  #Verify that edges were removed
temp[1:2,3:5]<-1             #Add new edges
check[10]<-sum(temp[,])==6   #Check edge sum
temp<-add.vertices(temp,3)   #Add vertices
check[11]<-network.size(temp)==8
check[12]<-sum(temp[,])==6   #Edges should still be there
check[13]<-all(temp[,5]==c(1,1,0,0,0,0,0,0))
temp[8,]<-1                  #Add edges to new vertex
check[14]<-all(temp[8,]==c(1,1,1,1,1,1,1,0))  #Verify
temp<-delete.vertices(temp,c(7,8))  #Remove vertices
check[15]<-network.size(temp)==6  #Verify removal
check[16]<-sum(temp[,])==6   #Check edge sum
check[17]<-!any(c(7,8)%in%c(sapply(temp$mel,"[[","inl"),sapply(temp$mel,"[[","outl")))  #Make sure they're really gone!
temp<-network(matrix(0,5,5),directed=FALSE,loops=TRUE)  #Create undir graph
check[18]<-is.directed(temp)==FALSE    #Some simple gal tests
check[19]<-has.loops(temp)==TRUE
temp[1,]<-1
check[20]<-all(temp[,1]==temp[1,])   #Verify edges
temp<-permute.vertexIDs(temp,5:1)       #Permute 
check[21]<-all(temp[1,]==c(0,0,0,0,1))  #Verify permutation
check[22]<-all(temp[,5]==rep(1,5))
check[23]<-all(get.neighborhood(temp,1)%in%c(5,1)) #Check neighborhoods
check[24]<-all(sort(get.neighborhood(temp,5))==1:5)
check[25]<-length(get.edges(temp,5))==5            #Check get.edges
check[26]<-length(get.edges(temp,5,2))==1
g<-rgraph(10)
temp<-network(g)
check[27]<-all(g==temp[,])                         #Yet more edge checkage
check[28]<-all(g[3:1,-(4:3)]==temp[3:1,-(4:3)])
temp[,,,names.eval="newval"]<-matrix(1:100,10,10)  #Edge value assignment
check[29]<-all(as.sociomatrix(temp,"newval")==matrix(1:100,10,10)*g)
check[30]<-all(apply(as.matrix.network.incidence(temp),1,sum)==(degree(g,cmode="indegree")-degree(g,cmode="outdegree")))  #Check incidence matrix
check[31]<-all(dim(as.matrix.network.incidence(temp))==c(10,sum(g)))
check[32]<-all(apply(as.matrix.network.incidence(temp,"newval"),1,sum)==(degree(matrix(1:100,10,10)*g,cmode="indegree")-degree(matrix(1:100,10,10)*g,cmode="outdegree")))
check[33]<-all(as.matrix.network.edgelist(temp,"newval")==cbind(row(g)[g>0],col(g)[g>0],matrix(1:100,10,10)[g>0]))
temp[1:3,1:5,names.eval="newval"]<-matrix(1:15,3,5)
check[34]<-all(as.sociomatrix(temp,"newval")[1:3,1:5]==g[1:3,1:5]*matrix(1:15,3,5))
temp[,,"na"]<-TRUE                         #Verify NA filtering
check[35]<-sum(temp[,,na.omit=TRUE])==0
check[36]<-sum(is.na(temp[,,na.omit=FALSE]))==sum(g)

#---- Check assignment, deletion, and adjacency for hypergraphs ----
temp<-network.initialize(10,directed=F,hyper=T,loops=T)
check[37]<-sum(temp[,])==0
temp<-add.edge(temp,1:4,1:4,"value",list(5))
temp<-add.edge(temp,3:5,3:5,"value",list(6))
temp<-add.edge(temp,4:7,4:7,"value",list(7))
temp<-add.edge(temp,6:10,6:10,"value",list(8))
check[38]<-all(as.matrix.network.incidence(temp)==cbind(c(1,1,1,1,0,0,0,0,0,0),c(0,0,1,1,1,0,0,0,0,0),c(0,0,0,1,1,1,1,0,0,0),c(0,0,0,0,0,1,1,1,1,1)))
check[39]<-all(as.matrix.network.incidence(temp,"value")==cbind(5*c(1,1,1,1,0,0,0,0,0,0),6*c(0,0,1,1,1,0,0,0,0,0),7*c(0,0,0,1,1,1,1,0,0,0),8*c(0,0,0,0,0,1,1,1,1,1)))
check[40]<-all(temp[,]==((as.matrix.network.incidence(temp)%*%t(as.matrix.network.incidence(temp)))>0))

#---- Check coercion and construction methods ----
g<-rgraph(10)
temp<-network(g)
check[41]<-all(temp[,]==g)
temp<-as.network(g*matrix(1:100,10,10),names.eval="value",ignore.eval=FALSE)
check[42]<-all(as.sociomatrix(temp,"value")==g*matrix(1:100,10,10))
temp<-as.network.matrix(as.matrix.network.edgelist(temp,"value"),matrix.type="edgelist",names.eval="value",ignore.eval=FALSE)
check[43]<-all(as.sociomatrix(temp,"value")==g*matrix(1:100,10,10))
temp<-as.network.matrix(as.matrix.network.incidence(temp,"value"),matrix.type="incidence",names.eval="value",ignore.eval=FALSE)
check[44]<-all(as.sociomatrix(temp,"value")==g*matrix(1:100,10,10))

# check functioning of na.rm argument #922
plain<-as.network.matrix(matrix(c(0,1,NA,NA),ncol=2),na.rm=TRUE)
if (network.naedgecount(plain) != 0){
  stop('problem with na values in adjacency matrix coericon')
}
plain<-as.network.matrix(matrix(c(0,1,NA,NA),ncol=2),na.rm=FALSE)
if (network.naedgecount(plain) != 1){
  stop('problem with na values in adjacnecy matrix coericon')
}



# check creating of network using dataframe with named cols
edata <-data.frame(
  tails=c(1,2,3),
  heads=c(2,3,1),
  love=c('yes','no','maybe'),
  hate=c(3,0,2)
  )

temp<-as.network(edata,matrix.type="edgelist",ignore.eval=FALSE)
if(!all(list.edge.attributes(temp)==c('hate','love','na'))){
  stop("problem with network edgelist coercion from data frame")
}
if(!all(temp%e%'hate'==c(3,0,2))){
  stop("problem with network edgelist coercion from data frame")
}
   
# ditto, but with passed in names for attributes
temp<-as.network(edata,matrix.type="edgelist",ignore.eval=FALSE,names.eval=c('hello','goodbye'))
if(!all(list.edge.attributes(temp)==c('goodbye','hello','na'))){
     stop("problem with network edgelist coercion from data frame")
}
if(!all(temp%e%'goodbye'==c(3,0,2))){
    stop("problem with network edgelist coercion from data frame")
}   
# test for as.matrix.network edgelist bug #935
x <- network.initialize(10)
add.edge(x,1,2)
add.edge(x,2,3)
set.edge.attribute(x,'foo','bar',e=2) # i.e. the edge from 2 to 3
if (!identical(as.matrix.network.edgelist(x,'foo'),structure(c("1", "2", "2", "3", NA, "bar"), .Dim = 2:3, n = 10, vnames = 1:10))){
  stop("problem with as.matrix.network.edgelist with attribute and deleted edge")
}

   

#---- Check attribute assignment/access ----
g<-rgraph(10)
temp<-network(g)
temp<-set.vertex.attribute(temp,"value",1:10)
check[45]<-all(get.vertex.attribute(temp,"value")==1:10)
temp<-delete.vertex.attribute(temp,"value")
check[46]<-all(is.na(get.vertex.attribute(temp,"value")))
temp<-set.vertex.attribute(temp,"value",1:5,c(2,4,6,8,10))
check[47]<-all(get.vertex.attribute(temp,"value")[c(2,4,6,8,10)]==1:5)
temp<-set.network.attribute(temp,"value","pork!")
check[48]<-get.network.attribute(temp,"value")=="pork!"
temp<-delete.network.attribute(temp,"value")
check[49]<-is.null(get.network.attribute(temp,"value"))
temp<-set.edge.attribute(temp,"value",5)
check[50]<-all(get.edge.attribute(temp$mel,"value")==5)
temp<-delete.edge.attribute(temp,"value")
check[51]<-all(is.null(get.edge.attribute(temp$mel,"value")))
temp<-set.edge.value(temp,"value",g*matrix(1:100,10,10))
check[52]<-all(get.edge.value(temp,"value")==(g*matrix(1:100,10,10))[g>0])
check[53]<-all(as.sociomatrix(temp,"value")==(g*matrix(1:100,10,10)))


#---- Check additional operators ----
g<-rgraph(10)
temp<-network(g,names.eval="value",ignore.eval=FALSE)
temp2<-network(g*2,names.eval="value",ignore.eval=FALSE)
check[54]<-all(g==as.sociomatrix(temp+temp2))
check[55]<-all(g*3==as.sociomatrix(sum(temp,temp2,attrname="value"),"value"))
check[56]<-all(g==as.sociomatrix(temp*temp2))
check[57]<-all(g*2==as.sociomatrix(prod(temp,temp2,attrname="value"),"value"))
check[58]<-all(0==as.sociomatrix(temp-temp2))
check[59]<-all(-g==as.sociomatrix(sum(temp,-as.sociomatrix(temp2,"value"),attrname="value"),"value"))
check[60]<-all(((g%*%g)>0)==as.sociomatrix("%c%.network"(temp,temp2)))
check[61]<-all(((g%*%g)>0)==as.sociomatrix(temp%c%temp2))
check[62]<-all(((!temp)[,]==!g)[diag(10)<1])
check[63]<-all((temp|temp2)[,]==g)
check[64]<-all((temp&temp2)[,]==g)
temp%v%"value"<-1:10
check[65]<-all(temp%v%"value"==1:10)
temp%n%"value"<-"pork!"
check[66]<-temp%n%"value"=="pork!"

# ---- Check to ensure that in-place modification is not producing side effects ----
g<-network.initialize(5); checkg<-g; add.vertices(g,3)
check[67]<-(network.size(checkg)==5)&&(network.size(g)==8)
g<-network.initialize(5); checkg<-g; delete.vertices(g,2)
check[68]<-(network.size(checkg)==5)&&(network.size(g)==4)
g<-network.initialize(5); checkg<-g; add.edge(g,2,3)
check[69]<-(sum(checkg[,])==0)&&(sum(g[,])==1)
g<-network.initialize(5); checkg<-g; add.edges(g,c(2,2,2),c(1,3,4))
check[70]<-(sum(checkg[,])==0)&&(sum(g[,])==3)
g<-network.initialize(5); checkg<-g; g%v%"boo"<-1:5
check[71]<-all(is.na(checkg%v%"boo"))&&all(g%v%"boo"==1:5)
g<-network.initialize(5); checkg<-g; g%n%"boo"<-1:5
check[72]<-is.null(checkg%n%"boo")&&all(g%n%"boo"==1:5)
g<-network.initialize(5); g[1,]<-1; checkg<-g; g%e%"boo"<-col(matrix(0,5,5))
check[73]<-is.null(checkg%e%"boo")&&all(g%e%"boo"==2:5)
g<-network.initialize(5); checkg<-g; permute.vertexIDs(g,5:1)
check[74]<-all(checkg%v%"vertex.names"==1:5)&&all(g%v%"vertex.names"==5:1)
g<-network.initialize(5); temp<-(function(){add.vertices(g,3); network.size(g)})()
check[75]<-(network.size(g)==5)&&(temp==8)
g<-network.initialize(5); (function(){g<-network.initialize(4); add.vertices(g,3)})()
check[76]<-(network.size(g)==5)

# check for operators with undirected edge error ticket #279
# nw1 is assigned tail<head
nw1<-network.initialize(3,directed=FALSE)
nw1[1,2]<-1

# nw2 is assigned tail>head
nw2<-network.initialize(3,directed=FALSE)
nw2[2,1]<-1

# Which, the binary network operators don't take into account:
check[77]<-network.edgecount(nw1-nw2)==0 # Should have 0, has 1.
check[78]<-network.edgecount(nw1|nw2)==1 # Should have 1, has 2 (1->2 and 2->1).
check[79]<-network.edgecount(nw1&nw2)==1 # Should have 1, has 0 (since it treats 1->2 and 2->1 differently).
check[80]<-network.edgecount(!nw1)==2 # Should have choose(3,2)-1=2, has 3.
check[81]<-network.edgecount(!nw2)==2 # Should have choose(3,2)-1=2, has 2.

#If everything worked, check is TRUE
if(!all(check)){                                               #Should be TRUE
  stop(paste("network package test failed on test(s):",which(!check)))
}

# ----- checks for network edgecount ------ 
require(testthat)
test<-network.initialize(4)
# directed
expect_equal(network.dyadcount(test),12)
# undirected
test%n%'directed'<-FALSE
expect_equal(network.dyadcount(test),6)

# loops allowed
test%n%'loops'<-TRUE
#undirected
expect_equal(network.dyadcount(test),10)
# directed
test%n%'directed'<-TRUE
expect_equal(network.dyadcount(test),16)

# directed bipartite
test%n%'loops'<-FALSE
test%n%'bipartite'<-1
expect_equal(network.dyadcount(test),6)

# undirected bipartite
test%n%'directed'<-FALSE
expect_equal(network.dyadcount(test),3)

# NA values
test[1,2]<-NA
expect_equal(network.dyadcount(test,na.omit = TRUE),2)


# ----- checks for dyads eids -----

data(emon)
el<-as.matrix.network.edgelist(emon[[1]])
expect_equal(get.dyads.eids(emon[[1]],el[,1],el[,2]),as.list(1:83))
expect_equal(get.dyads.eids(emon[[1]],el[5:10,1],el[5:10,2]),as.list(5:10))
expect_error(get.dyads.eids(emon[[1]],1,2:3),regexp = 'heads and tails vectors must be the same length')
expect_error(get.dyads.eids(network.initialize(0),1,2),regexp = 'invalid vertex id in heads or tails vector')

mult<-network.initialize(5,multi=TRUE)
add.edges(mult,1,2)
add.edges(mult,1,2)
expect_warning(expect_true(is.na(get.dyads.eids(mult,1,2)[[1]])),regexp = 'multiple edge ids for dyad')

expect_equal(get.dyads.eids(network.initialize(0),numeric(0),numeric(0)), list())
expect_equal(get.dyads.eids(network.initialize(5),tails=1:2,heads=3:4),list(numeric(0),numeric(0)))

# check oposite matching for undirected nets
undir<-network.initialize(3,directed=FALSE)
undir[1,2]<-1
expect_equal(get.dyads.eids(undir,2,1),list(1))
expect_equal(get.dyads.eids(undir,1,2),list(1))


undir%n%'directed'<-TRUE
expect_equal(get.dyads.eids(undir,2,1),list(integer(0)))
expect_equal(get.dyads.eids(undir,1,2),list(1))

expect_equal(get.dyads.eids(undir,2,1,neighborhood='in'),list(1))
expect_equal(get.dyads.eids(undir,1,2,neighborhood='in'),list(integer(0)))



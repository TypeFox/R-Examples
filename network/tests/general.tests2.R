# additional tests of misc network functionality split off from general.tests.R to avoid speed warnings
library(network)

# ----- check memory saftey with a big assignment ---
net<-network.initialize(100000)
net<-add.edges(net,1:99999,2:100000)
set.edge.attribute(net,'LETTERS',LETTERS)

# --- tests for get.induced.subgraph additions --
data(emon)
# extract the network of responders in MtStHelens network with interaction Frequency of 4
subG4<-get.inducedSubgraph(emon$MtStHelens,eid=which(emon$MtStHelens%e%'Frequency'==4))
if(network.size(subG4)!=24){
  stop('wrong size eid induced subgraph')
}

if (any(subG4%e%'Frequency'!=4)){
  stop('bad edges in eid induced subgraph')
}

# checks for error conditions
# can't specify eid with v or alter
# get.inducedSubgraph(v=1:2,emon$MtStHelens,eid=which(emon$MtStHelens%e%'Frequency'==4))
# get.inducedSubgraph(alter=1:2,emon$MtStHelens,eid=which(emon$MtStHelens%e%'Frequency'==4))
# get.inducedSubgraph(emon$MtStHelens,eid=200:300)


# ---- tests for specific bugs/edgecases -----

# ticket #180 (used to throw error if no edges exist)
set.edge.attribute(network.initialize(3),"test","a")

# check for network of zero size --used to give error ticket #255
set.vertex.attribute(network.initialize(0),'foo','bar')


# check for is.na.network problems #619
x2<-network.initialize(3)
x2[1,2]<-NA
if(is.na.network(x2)[1,2]!=1){
  stop('problem iwth is.na.netowrk')
}

# check for na problems in which.matrix.type #926
mat <- matrix(rbinom(200, 1, 0.2), nrow = 20)
naIndices <- sample(1:200, 20)
mat[naIndices] <- NA
nw <- network(mat)

# ---- check for undirected loops getID cases #327 #609 -----
net<-network.initialize(2,loops=TRUE,directed=FALSE)
net[1,1]<-1
net[1,2]<-1
net[2,2]<-1
if(get.edgeIDs(net,v=1,alter=1)!=1){
  stop("problem with get.edgeIDs on undirected network with loops")
}
if(get.edgeIDs(net,v=2,alter=2)!=3){
  stop("problem with get.edgeIDs on undirected network with loops")
}

net<-network.initialize(2,loops=TRUE,directed=FALSE)
net[1,2]<-1
if(length(get.edgeIDs(net,v=2,alter=2))>0){
  stop("problem with get.edgeIDs on undirected network with loops")
}

# check for problem with as.network.edgelist with zero edges #1138
result <- as.matrix.network.edgelist(network.initialize(5),as.sna.edgelist = TRUE)
if (nrow(result) != 0){
  stop('as.matrix.network.edgelist did not return correct value for net with zero edges')
}
as.matrix.network.edgelist(network.initialize(0),as.sna.edgelist = TRUE)
result2<-as.matrix.network.adjacency(network.initialize(5))
if(nrow(result2) != 5 & ncol(result2) != 5){
  stop('as.matrix.network.adjacency did not return matrix with correct dimensions')
}
result3<-as.matrix.network.adjacency(network.initialize(0))
if(nrow(result3) != 0 & ncol(result3) != 0){
  stop('as.matrix.network.adjacency did not return matrix with correct dimensions')
}
result4<-as.matrix.network.incidence(network.initialize(5))
if(nrow(result4) != 5 & ncol(result4) != 0){
  stop('as.matrix.network.incidence did not return matrix with correct dimensions')
}
result5<-as.matrix.network.incidence(network.initialize(0))
if(nrow(result5) != 0 & ncol(result5) != 0){
  stop('as.matrix.network.incidence did not return matrix with correct dimensions')
}

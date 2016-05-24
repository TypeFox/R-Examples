require(network)

# --------- test list.vertex.attributes ---

net<-network.initialize(3)

list.vertex.attributes(net)

if(!all(list.vertex.attributes(net)==c('na','vertex.names'))){
  stop('list.vertex.attribute did not report default attributes corrrectly')
}

set.vertex.attribute(net,'letters',c("a","b","c"))

if(!all(list.vertex.attributes(net)==c('letters','na','vertex.names'))){
  stop('list.vertex.attribute did not report added attributes corrrectly')
}


# ----- test list.edge.attributes ----

net<-network.initialize(3)
if(length(list.edge.attributes(net))!=0){
  stop("list.edge.attributes did not return empty list for network with no edges")
}

add.edges(net,1,2)
add.edges(net,2,3)
if(list.edge.attributes(net)!='na'){
  stop("list.edge.attributes did not return 'na' for network with only edges")
}

set.edge.attribute(net,'letter',c("a","b"))
if(!all(list.edge.attributes(net)==c('letter','na'))){
  stop("list.edge.attributes did not return attribute names for network with edges")
}
   
delete.edges(net,eid=1)   
if(!all(list.edge.attributes(net)==c('letter','na'))){
   stop("list.edge.attributes did not return attribute names for network deleted edge")
}
   
# ---- test list.network.attributes ----
net<-network.initialize(3)   
if(!all(list.network.attributes(net)==c("bipartite", "directed",  "hyper","loops","mnext",     "multiple","n" ))){
  stop("list.network.attributes returned unexpected values for default attributes of a network")
} 
   
set.network.attribute(net,'letter',"a")   
   if(!all(list.network.attributes(net)==c("bipartite", "directed",  "hyper","letter","loops","mnext",     "multiple","n" ))){
     stop("list.network.attributes returned unexpected values for network with attribute added")
   } 

# ----- tests for printing function for edges cases ------
net<-network.initialize(100)
net%n%'a_matrix'<-matrix(1:100,nrow=10,ncol=10)
net%n%'a_null'<-NULL
net%n%'a_list'<-list(part1=list(c("A","B")),part2=list("c"))
net%n%'a_desc_vec'<-numeric(rep(100,1))
net%n%'a_net'<-network.initialize(5)
print.network(net)

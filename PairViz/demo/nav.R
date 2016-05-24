# This file constructs some of the graphs used in "Graphs as navigational infrastructure for high dimensional data spaces."

require(PairViz)
require(Rgraphviz)
#----------------
# Versions of Fig 1,2,3
n <- colnames(iris)[1:4] <- c("SL","SW","PL","PW")
g <- mk_complete_graph(n)
plot(g,"circo")

space2 <- kspace_graph(n,2)
plot(space2,"circo")


lg <- mk_line_graph(g)
plot(lg,"circo")



lgc <- complement(lg)
plot(lgc)

space21 <- kspace_graph(n,2,1) # same as lg

space20 <- kspace_graph(n,2,0) # same as lgc

space32 <- kspace_graph(n,3,2)
plot(space32,"circo")



#----------------
# Versions of Fig 7,8

n <- LETTERS[1:5]
g <- mk_complete_graph(n)
lg <- mk_line_graph(g,sep="")
plot(lg,"circo")
petersen <- complement(lg)
plot(petersen,"circo")

space31 <- kspace_graph(n,2,1) # same as lg
space30 <-  kspace_graph(n,2,0) # same as petersen

plot(space32,"circo")

space31 <- kspace_graph(n,3,1)

plot(space31,"circo")



#----------------
# Versions of Fig 9,10
require(gclus)
data(ozone)
n <- colnames(ozone)

g <- new("graphNEL", nodes=n)
g <- addEdge(n[1], n[2:5],g)
plot(g,"circo")
 
lg <- mk_line_graph(g)
plot(lg,"circo")

    
d <- abs(cor(ozone, method="spearman"))
g <- mk_complete_graph(d)
q <- quantile(as.dist(d),.9)
g1 <- dn_graph(g,q, `>=`)
plot(g1)
lg1 <- mk_line_graph(g1)
plot(lg1,"circo")
plot(complement(lg1),"circo")


#----------------
# Versions of Fig 14

n <- LETTERS[1:5]
space32 <-  kspace_graph(n,3,2) 

plot(space32,"circo")


space31 <- kspace_graph(n,3,1) 

plot(space31,"circo")

space32c <- complement(space32) # same as space31


#----------------
# Versions of Fig 15-19

xn <-paste("X",1:3,sep="")
yn <- paste("Y",1:2,sep="")
g <- bipartite_graph(xn,yn)
plot(g)
lg <- mk_line_graph(g,sep="")
plot(lg)
plot(complement(lg))

space3d <- mk_line_graph(lg)
plot(space3d,"circo")

space3d <- iterated_line_graph(g,sep=".")
plot(space3d,"circo") # this version does compression, not needed here


u <-  new("graphNEL", nodes=paste("U",1:3,sep=""))
u <- addEdge("U2",c("U1","U3"),u)
plot(u)
v <-  new("graphNEL", nodes=paste("V",1:2,sep=""))
v <- addEdge("V2","V1",v)
plot(v)

g1 <- graph_product(u,v,"cartesian")
plot(g1)


g2 <- graph_product(u,v,"tensor")
plot(g2)

g3<- graph_product(u,v,"strong")
plot(g3)


 
g4<- graph_compose(u,v)
plot(g4)

g5<- graph_compose(v,u)
plot(g5)
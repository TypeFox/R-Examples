## ---- echo=FALSE, results="hide"-----------------------------------------
	library(rnetcarto)
	require("igraph")

## ---- echo=TRUE----------------------------------------------------------
 # Generate a simple random network
 a = matrix(as.integer(runif(100)<.3), ncol=10) 
 a[lower.tri(a)] = 0
 rownames(a) = c('a','b','b','c','d','e','f','g','h','i')
 colnames(a) = rownames(a)
 # Find an optimal partition for modularity using netcarto.
 #  The output consists in a table containing node properties,
 #  and the modularity value of the partition.
 netcarto(a)

## ---- echo=TRUE----------------------------------------------------------
    input = matrix(0,3,3)
    input[1,2] = 1
    input[2,3] = 1
    input[3,1] = 1
	input[2,1] = 1
    input[3,2] = 1
    input[1,3] = 1
    rownames(input) = c("A","B","C")
    colnames(input) = rownames(input)
	print(input)

## ---- echo=TRUE----------------------------------------------------------
    # import from rnetcarto matrix format to igraph:
    G = igraph::graph.adjacency(input,weighted=TRUE,mode="undirected")
    # Export to a matrix compatible with netcarto:
	input = igraph::get.adjacency(G,sparse=FALSE)

## ---- echo=FALSE---------------------------------------------------------
	plot(G, layout = igraph::layout.circle, ,
       vertex.size = 60,
       vertex.color="red",
       vertex.frame.color= "white",
       vertex.label.color = "white",
       vertex.label.family = "sans",
       edge.width=1,
       edge.color="black")

## ---- echo=TRUE----------------------------------------------------------
    input = matrix(0,7,7)
    input[1,2] = 10
    input[2,3] = 10
    input[3,1] = 10
    input[4,5] = 10
    input[5,6] = 10
    input[6,4] = 10
    rownames(input) = c("A","B","C","D","E","F","G")
    colnames(input) = rownames(input)

## ---- echo=FALSE---------------------------------------------------------
    input = matrix(0,6,6)
    input[1,2] = 10
    input[2,3] = 10
    input[3,1] = 10
    input[4,5] = 10
    input[5,6] = 10
    input[6,4] = 10
	input = input+t(input)-diag(input)
    rownames(input) = c("A","B","C","D","E","F")
    colnames(input) = rownames(input)
	print(input)

## ---- echo=FALSE---------------------------------------------------------
    G = igraph::graph.adjacency(input,weighted=TRUE,mode="undirected")
	plot(G, layout = layout.circle, ,
       vertex.size = 60,
       vertex.color="red",
       vertex.frame.color= "white",
       vertex.label.color = "white",
       vertex.label.family = "sans",
       edge.width=1,
       edge.color="black")

## ---- echo=TRUE----------------------------------------------------------
    input = matrix(0,6,2)
    input[1,1] = 1
    input[2,1] = 1
    input[3,1] = 1
    input[4,2] = 1
    input[5,2] = 1
    input[6,2] = 1
    rownames(input) = c("A","B","C","D","E","F")
    colnames(input) = c("Team 1", "Team 2")
	print(input)

## ---- echo=TRUE----------------------------------------------------------
    nd1 = c("A","B","C","D","E","F","C")
    nd2 = c("B","C","A","E","F","D","D")
	web = list(nd1,nd2,weights)
    print(list(nd1,nd2))

## ---- echo=TRUE----------------------------------------------------------
    nd1 = c("A","B","C","D","E","F","C","A")
    nd2 = c("B","C","A","E","F","D","D","D")
	weights = c(10,10,10,10,10,10,10,10,1)
	web = list(nd1,nd2,weights)
    print(web)

## ---- echo=TRUE----------------------------------------------------------
    nd1 = c("A","B","C","D","E","F","C","A")
    nd2 = c("Team1","Team2","Team1","Team1","Team2","Team1","Team1","Team2")
	bipartite = list(nd1,nd2)
    print(bipartite)

## ---- echo=TRUE----------------------------------------------------------
    netcarto(igraph::get.adjacency(G,sparse=FALSE))

## ---- echo=TRUE----------------------------------------------------------
   netcarto(bipartite, bipartite=TRUE)


library(CRF)

nNodes <- 4
nStates <- 2

adj <- matrix(0, nrow=nNodes, ncol=nNodes)
for (i in 1:(nNodes-1))
{
	adj[i,i+1] <- 1
	adj[i+1,i] <- 1
}
adj[nNodes,1] <- 1
adj[1,nNodes] <- 1

crf <- make.crf(adj, nStates)

crf$node.pot[1,] <- c(1, 3)
crf$node.pot[2,] <- c(9, 1)
crf$node.pot[3,] <- c(1, 3)
crf$node.pot[4,] <- c(9, 1)

for (i in 1:crf$n.edges)
{
   crf$edge.pot[[i]][1,] <- c(2, 1)
   crf$edge.pot[[i]][2,] <- c(1, 2)
}

answer <-
structure(list(decode = c(1L, 1L, 1L, 1L), node.bel = structure(c(0.485633270321361, 
0.859168241965974, 0.485633270321361, 0.859168241965974, 0.514366729678639, 
0.140831758034026, 0.514366729678639, 0.140831758034026), .Dim = c(4L, 
2L)), edge.bel = list(structure(c(0.455954631379962, 0.403213610586011, 
0.0296786389413989, 0.111153119092628), .Dim = c(2L, 2L)), structure(c(0.455954631379962, 
0.403213610586011, 0.0296786389413989, 0.111153119092628), .Dim = c(2L, 
2L)), structure(c(0.455954631379962, 0.0296786389413989, 0.403213610586011, 
0.111153119092628), .Dim = c(2L, 2L)), structure(c(0.455954631379962, 
0.403213610586011, 0.0296786389413989, 0.111153119092628), .Dim = c(2L, 
2L))), logZ = 8.57357352485234), .Names = c("decode", "node.bel", 
"edge.bel", "logZ"))


Loop <- list()
Loop$crf <- crf
Loop$answer <- answer
save(Loop, file="Loop.RData")

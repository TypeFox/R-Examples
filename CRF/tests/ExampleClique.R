library(CRF)

nNodes <- 4
nStates <- 2

adj <- matrix(1, nrow=nNodes, ncol=nNodes)
for (i in 1:(nNodes-1))
{
	adj[i,i] <- 0
}

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
structure(list(decode = c(1L, 1L, 1L, 1L), node.bel = structure(c(0.518774157923799, 
0.892048591938156, 0.518774157923799, 0.892048591938156, 0.481225842076201, 
0.107951408061844, 0.481225842076201, 0.107951408061844), .Dim = c(4L, 
2L)), edge.bel = list(structure(c(0.504417448923247, 0.387631143014909, 
0.0143567090005522, 0.0935946990612921), .Dim = c(2L, 2L)), structure(c(0.368028713418001, 
0.150745444505798, 0.150745444505798, 0.330480397570403), .Dim = c(2L, 
2L)), structure(c(0.504417448923247, 0.387631143014909, 0.0143567090005522, 
0.0935946990612921), .Dim = c(2L, 2L)), structure(c(0.504417448923247, 
0.0143567090005522, 0.387631143014909, 0.0935946990612921), .Dim = c(2L, 
2L)), structure(c(0.827443401435671, 0.0646051905024848, 0.0646051905024848, 
0.0433462175593595), .Dim = c(2L, 2L)), structure(c(0.504417448923247, 
0.387631143014909, 0.0143567090005522, 0.0935946990612921), .Dim = c(2L, 
2L))), logZ = 9.58107599956325), .Names = c("decode", "node.bel", 
"edge.bel", "logZ"))


Clique <- list()
Clique$crf <- crf
Clique$answer <- answer
save(Clique, file="Clique.RData")

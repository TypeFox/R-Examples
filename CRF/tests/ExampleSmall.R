library(CRF)

nNodes <- 4
nStates <- 2

adj <- matrix(0, nrow=nNodes, ncol=nNodes)
for (i in 1:(nNodes-1))
{
	adj[i,i+1] <- 1
	adj[i+1,i] <- 1
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
structure(list(decode = c(2L, 1L, 1L, 1L), node.bel = structure(c(0.359630606860158, 
0.843007915567282, 0.486279683377309, 0.881002638522427, 0.640369393139842, 
0.156992084432718, 0.513720316622691, 0.118997361477573), .Dim = c(4L, 
2L)), edge.bel = list(structure(c(0.337203166226913, 0.505804749340369, 
0.0224274406332454, 0.134564643799472), .Dim = c(2L, 2L)), structure(c(0.451187335092348, 
0.0350923482849604, 0.391820580474934, 0.121899736147757), .Dim = c(2L, 
2L)), structure(c(0.460686015831135, 0.420316622691293, 0.0255936675461741, 
0.0934036939313984), .Dim = c(2L, 2L))), logZ = 8.24012129807647), .Names = c("decode", 
"node.bel", "edge.bel", "logZ"))


Small <- list()
Small$crf <- crf
Small$answer <- answer
save(Small, file="Small.RData")

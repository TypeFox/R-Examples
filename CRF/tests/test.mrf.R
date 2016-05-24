library(CRF)

data(Rain)
rain <- Rain$rain
months <- Rain$months

p.rain <- sum(rain == 2) / length(rain)
nll <- log(p.rain) * sum(rain == 2) + log(1-p.rain) * sum(rain == 1)

n.nodes <- dim(rain)[2]
adj <- matrix(0, n.nodes, n.nodes)
for (i in 1:(n.nodes-1))
{
	adj[i, i+1] <- 1
}

crf <- make.crf(adj, 2)
crf <- make.features(crf)
crf <- make.par(crf, 2)

crf$node.par[,1,] <- 1
for (i in 1:crf$n.edges)
{
	crf$edge.par[[i]][1,1,] <- 2
	crf$edge.par[[i]][2,2,] <- 2
}

crf$par.stat <- mrf.stat(crf, rain)

nll <- mrf.nll(crf$par, crf, rain)
crf$nll
crf$gradient

crf <- train.mrf(crf, rain)
crf$par
crf$node.pot[1,]
crf$edge.pot[[1]]

##################################################

crf <- make.par(crf, 4)

for (i in 1:crf$n.edges)
{
	crf$edge.par[[i]][1,1,] <- 2
	crf$edge.par[[i]][1,2,] <- 3
	crf$edge.par[[i]][2,1,] <- 4
	crf$edge.par[[i]][2,2,] <- 0
}

crf <- train.mrf(crf, rain)
crf$par
crf$node.pot[1,]
crf$edge.pot[[1]]

##################################################

crf <- make.par(crf, 5)

crf$node.par[1,1,] <- 5
crf$node.par[crf$n.nodes,1,] <- 5

crf <- train.mrf(crf, rain)
crf$par
crf$node.pot[1,]
crf$edge.pot[[1]]

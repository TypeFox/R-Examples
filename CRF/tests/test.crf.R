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
crf <- make.par(crf, 4)

crf$node.par[,1,] <- 1
for (i in 1:crf$n.edges)
{
	crf$edge.par[[i]][1,1,] <- 2
	crf$edge.par[[i]][1,2,] <- 3
	crf$edge.par[[i]][2,1,] <- 4
}

node.fea <- lapply(1:dim(rain)[1], function(i) matrix(1, crf$n.nf, crf$n.nodes))
edge.fea <- lapply(1:dim(rain)[1], function(i) matrix(1, crf$n.ef, crf$n.edges))

nll <- crf.nll(crf$par, crf, rain, node.fea, edge.fea)
crf$nll
crf$gradient

crf <- train.crf(crf, rain, node.fea, edge.fea)
crf$par
crf$node.pot[1,]
crf$edge.pot[[1]]

##################################################

crf <- make.features(crf, 0, 0)
crf <- make.par(crf, 4)

node.ext <- list()
edge.ext <- list()

node.ext[[1]] <- list()
node.ext[[1]][[1]] <- matrix(0, crf$n.nodes, crf$max.state)
node.ext[[1]][[1]][,1] <- 1
node.ext[[1]][2:crf$n.par] <- NaN

edge.ext[[1]] <- list()
edge.ext[[1]][[1]] <- NaN
edge.ext[[1]][[2]] <- list()
edge.ext[[1]][[2]][[1]] <- matrix(0, 2, 2)
edge.ext[[1]][[2]][[1]][1,1] <- 1
edge.ext[[1]][[2]][2:crf$n.edges] <- edge.ext[[1]][[2]][1]
edge.ext[[1]][[3]] <- list()
edge.ext[[1]][[3]][[1]] <- matrix(0, 2, 2)
edge.ext[[1]][[3]][[1]][1,2] <- 1
edge.ext[[1]][[3]][2:crf$n.edges] <- edge.ext[[1]][[3]][1]
edge.ext[[1]][[4]] <- list()
edge.ext[[1]][[4]][[1]] <- matrix(0, 2, 2)
edge.ext[[1]][[4]][[1]][2,1] <- 1
edge.ext[[1]][[4]][2:crf$n.edges] <- edge.ext[[1]][[4]][1]
node.ext[2:dim(rain)[1]] <- node.ext[1]
edge.ext[2:dim(rain)[1]] <- edge.ext[1]

nll <- crf.nll(crf$par, crf, rain, NULL, NULL, node.ext, edge.ext)
crf$nll
crf$gradient

crf <- train.crf(crf, rain, NULL, NULL, node.ext, edge.ext)
crf$par
crf$node.pot[1,]
crf$edge.pot[[1]]

##################################################

crf <- make.features(crf, 13)
crf <- make.par(crf, 16)

for (i in 1:crf$n.nf)
{
	crf$node.par[,1,i] <- i
}
for (i in 1:crf$n.edges)
{
	crf$edge.par[[i]][1,1,] <- crf$n.nf + 1
	crf$edge.par[[i]][1,2,] <- crf$n.nf + 2
	crf$edge.par[[i]][2,1,] <- crf$n.nf + 3
}

node.fea <- lapply(1:dim(rain)[1], function(i) matrix(0, crf$n.nf, crf$n.nodes))
for (i in 1:dim(rain)[1])
{
	node.fea[[i]][1,] <- 1
	node.fea[[i]][months[i]+1,] <- 1
}
edge.fea <- lapply(1:dim(rain)[1], function(i) matrix(1, crf$n.ef, crf$n.edges))

nll <- crf.nll(crf$par, crf, rain, node.fea, edge.fea)
crf$nll
crf$gradient

crf <- train.crf(crf, rain, node.fea, edge.fea)
crf$par
crf$node.pot[1,]
crf$edge.pot[[1]]

##################################################

crf <- make.par(crf, 5)

crf$node.par[1,1,] <- 5
crf$node.par[crf$n.nodes,1,] <- 5

node.fea <- lapply(1:dim(rain)[1], function(i) matrix(1, crf$n.nf, crf$n.nodes))
edge.fea <- lapply(1:dim(rain)[1], function(i) matrix(1, crf$n.ef, crf$n.edges))

crf <- train.crf(crf, rain, node.fea, edge.fea)
crf$par
crf$node.pot[1,]
crf$edge.pot[[1]]

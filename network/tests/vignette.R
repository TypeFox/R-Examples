require("network")
set.seed(1702)

results = NULL

data("flo")
data("emon")

net <- network.initialize(5)
net

nmat <- matrix(rbinom(25, 1, 0.5), nr = 5, nc = 5)
net <- network(nmat, loops = TRUE)
net

summary(net)
results[1] = all(nmat == net[,])

net <- as.network(nmat, loops = TRUE)
results[2] = all(nmat == net[,])

nflo <- network(flo, directed = FALSE)
nflo

results[3] = all(nflo[9,] == c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1))
results[4] = nflo[9,1] == 1
results[5] = nflo[9,4] == 0
results[6] = is.adjacent(nflo, 9, 1) == TRUE
results[7] = is.adjacent(nflo, 9, 4) == FALSE

results[8] = network.size(nflo) ==  16
results[9] = network.edgecount(nflo) == 20
results[10] = network.density(nflo) == 1/6
results[11] = has.loops(nflo) == FALSE
results[12] = is.bipartite(nflo) == FALSE
results[13] = is.directed(nflo) == FALSE
results[14] = is.hyper(nflo) == FALSE
results[15] = is.multiplex(nflo) == FALSE

as.sociomatrix(nflo)

results[16] = all(nflo[,] == as.sociomatrix(nflo))
results[17] = all(as.matrix(nflo) == as.sociomatrix(nflo))
as.matrix(nflo,matrix.type = "edgelist")

net <- network.initialize(5, loops = TRUE)
net[nmat>0] <- 1
results[18] = all(nmat == net[,])

net[,] <- 0
net[,] <- nmat
results[19] = all(nmat == net[,])

net[,] <- 0
for(i in 1:5)
  for(j in 1:5)
    if(nmat[i,j])
      net[i,j] <- 1
results[20] = all(nmat == net[,])

net[,] <- 0
add.edges(net, row(nmat)[nmat>0], col(nmat)[nmat>0])
results[21] = all(nmat == net[,])

net[,] <- as.numeric(nmat[,])
results[22] = all(nmat == net[,])

net <- network.initialize(5)
add.edge(net, 2, 3)
net[,]
results[23] = net[2,3] == 1

add.edges(net, c(3, 5), c(4, 4))
net[,]
results[24] = (net[3,4] == 1 && net[5,4] == 1)

net[,2] <- 1
net[,]
results[25] = net[2,2] == 0

delete.vertices(net, 4)
results[26] = all(net[,] == matrix(c(0,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0), byrow=T, nrow=4))

add.vertices(net, 2)
net[,]

get.edges(net, 1)
get.edges(net, 2, neighborhood = "in")
get.edges(net, 1, alter = 2)

results[27] = get.edgeIDs(net, 1) == 4
results[28] = all(get.edgeIDs(net, 2, neighborhood = "in") == c(7, 5, 4))
results[29] = get.edgeIDs(net, 1, alter = 2) == 4

results[30] = get.neighborhood(net, 1) == 2
results[31] = all(get.neighborhood(net, 2, type = "in") == c(4, 3, 1))

net[2,3] <- 0
results[32] = net[2,3] == 0

delete.edges(net, get.edgeIDs(net, 2, neighborhood = "in"))
results[33] = all(net[,] == matrix(0, 6,6))

net <- network.initialize(5)
set.network.attribute(net, "boo", 1:10)
net %n% "hoo" <- letters[1:7]

results[34] = 'boo' %in% list.network.attributes(net)
results[35] = 'hoo' %in% list.network.attributes(net)

results[36] = all(get.network.attribute(net, "boo") == c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
results[37] = all(net %n% "hoo" == c("a", "b", "c", "d", "e", "f", "g"))

delete.network.attribute(net, "boo")
results[38] = 'boo' %in% list.network.attributes(net) == FALSE

set.vertex.attribute(net, "boo", 1:5)
net %v% "hoo" <- letters[1:5]
results[39] = 'boo' %in% list.vertex.attributes(net)
results[40] = 'hoo' %in% list.vertex.attributes(net)

results[41] = all(get.vertex.attribute(net, "boo") == 1:5)
results[42] = all(net %v% "hoo" == letters[1:5])
delete.vertex.attribute(net, "boo")
results[43] = 'boo' %in% list.vertex.attributes(net) == FALSE

net <- network(nmat)
set.edge.attribute(net, "boo", sum(nmat):1)
set.edge.value(net, "hoo", matrix(1:25, 5, 5))
net %e% "woo" <- matrix(rnorm(25), 5, 5)
net[,, names.eval = "zoo"] <- nmat * 6
results[44] = 'boo' %in% list.edge.attributes(net)
results[45] = 'hoo' %in% list.edge.attributes(net)

results[46] = all(get.edge.attribute(get.edges(net, 1), "boo") == c(3,7))
results[47] = all(get.edge.value(net, "hoo") == c(2, 3, 11, 14, 17, 18, 21))
net %e% "woo"
as.sociomatrix(net, "zoo")
delete.edge.attribute(net, "boo")
results[48] = 'boo' %in% list.edge.attributes(net) == FALSE

MtSHloc <- emon$MtStHelens %v% "Location"
MtSHimat <- cbind(MtSHloc %in% c("L", "B"), MtSHloc %in% c("NL", "B"))
MtSHbyloc <- network(MtSHimat, matrix = "incidence", hyper = TRUE,
                        directed = FALSE, loops = TRUE)
MtSHbyloc %v% "vertex.names" <- emon$MtStHelens %v% "vertex.names"
MtSHbyloc

plot(nflo, displaylabels = TRUE, boxed.labels = FALSE)
plot(nflo, displaylabels = TRUE, mode = "circle")
plot(emon$MtSi)

if (!all(results)) {
  stop(paste('The following tests in vignette.R failed:', which(results==FALSE)))
}


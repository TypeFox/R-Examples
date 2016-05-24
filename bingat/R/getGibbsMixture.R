getGibbsMixture <-
function(data, type, desiredGroups, maxIter=50, digits=3){ 
if(missing(data) || missing(desiredGroups) || missing(type))
stop("data, type, and/or desiredGroups is missing.")

if(maxIter <= 0)
stop("maxIter must be an integer greater than 0.")

nodes <- getNumNodes(data, type)
edges <- getNumEdges(nodes, type)
groups <- sample(1:desiredGroups, ncol(data), replace=TRUE)

gstars <- list()
weights <- NULL
taus <- NULL

gcheck <- 0
tcheck <- 0
wcheck <- 1
group <- NULL
converge <- NA

#Find the starting gstars/tau/weights for the given groups
for(j in 1:desiredGroups){
weights[j] <- sum(groups==j) / length(groups)
gstars[[j]] <- estGStar(data[, groups==j, drop=FALSE]) 
taus[j] <- estTau(data[, groups==j, drop=FALSE], type, gstars[[j]])
}

for(iter in 1:maxIter){
#Calculate distances to the gstars
distances <- matrix(NA, ncol(data), desiredGroups)
for(j in 1:desiredGroups)
distances[,j] <- apply(data, 2, function(x){calcDistance(x, gstars[[j]], type)})

t <- matrix(NA, ncol(data), desiredGroups)
pij <- matrix(NA, ncol(data), desiredGroups)
for(i in 1:ncol(data)){
for(j in 1:desiredGroups)
t[i, j] <- -edges*log(1 + exp(-taus[j])) - taus[j]*distances[i, j] + log(weights[j])
for(j in 1:desiredGroups)
pij[i, j] <- 1/sum(exp(t[i,] - t[i, j]))
}

#Check for covergence and quit if true
if(gcheck + tcheck + wcheck == 0){
for(i in 1:ncol(data))
group <- c(group, which(pij[i,] == max(pij[i,])))

converge <- TRUE
break
}

#Recompute gstars/taus/weights for new groups
gstars2 <- list()
weights2 <- NULL
taus2 <- NULL
for(j in 1:desiredGroups){
weights2[j] <- sum(pij[,j]) / ncol(data)
gnum <- apply(t(data)*pij[,j], 2, sum)
gdem <- sum(pij[,j])
gstars2[[j]] <- ifelse(gnum/gdem > .5, 1, 0)  
distj <- apply(data, 2, function(x){calcDistance(x, gstars2[[j]], type)})

tnum <- sum(pij[,j] * distj)
tdem <- sum(pij[,j] * edges) - tnum
tfrac <- tnum/tdem
tfrac <- ifelse(tfrac==0, .Machine$double.xmin, tfrac)

taus2[j] <- -log(tfrac) 
}

#Recalculate the checks for convergence
gcheck <- 0
tcheck <- 0
wcheck <- 0
for(j in 1:desiredGroups){
gcheck <- gcheck + sum(gstars[[j]] != gstars2[[j]])
tcheck <- tcheck + as.numeric(round(taus[j], digits) != round(taus2[j], digits))
wcheck <- wcheck + as.numeric(round(weights[j], digits) != round(weights2[j], digits))
}

if(iter == maxIter){
warning(sprintf("EM algorithm did not converge on %s groups with %s iterations", as.character(desiredGroups), as.character(maxIter)))
converge <- FALSE
}

gstars <- gstars2
taus <- taus2
weights <- weights2
}

return(list(weights=weights2, gstars=gstars2, taus=taus2, converge=converge, 
iterations=iter, numgroups=desiredGroups, type=type, pij=pij, group=group))
}

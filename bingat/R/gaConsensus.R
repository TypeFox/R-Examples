gaConsensus <-
function(data, groups, iters=10, nRuns=1, popSize=200, method="manhattan", parallel=FALSE, cores=3){
if(missing(data) || missing(groups))
stop("data and/or groups is missing.")

if(length(unique(groups)) != 2)
stop("There must be exactly two groups.")

### Get our groups in the right order
x1 <- sum(groups==unique(groups)[1])
x2 <- sum(groups==unique(groups)[2])
if(x1 <= 0 || x2 <= 0)
stop("Each group must have at least 1 subject.")

grp1 <- data[,groups==unique(groups)[1]]
grp2 <- data[,groups==unique(groups)[2]]
groupSize1 <- ncol(grp1)
groupSize2 <- ncol(grp2)
data <- cbind(grp1, grp2)

### Set up our function for the GA
b <- dist(c(rep(0, groupSize1), rep(1, groupSize2)), method)
myScore <- function(indices) {
if(sum(indices) == 0)
return(0)
edges <- which(indices==1)

a <- Reduce("+", dists[edges])
mycor <- cor(a, b)
return(-mycor)
}

### Precompute all our distances
dists <- vector("list", nrow(data))
for(i in 1:nrow(data))
dists[[i]] <- dist(t(data[i,]), method)

### Run the GA
if(parallel){
cl <- parallel::makeCluster(cores) 
doParallel::registerDoParallel(cl)

res <- foreach::foreach(i=1:nRuns, .combine=rbind, .inorder=FALSE, .multicombine=TRUE, .packages=c("genalg")) %dopar%{
tempGa <- genalg::rbga.bin(size=nrow(data), iters=iters, popSize=popSize, evalFunc=myScore, verbose=TRUE)
return(tempGa$population[which(tempGa$evaluation==min(tempGa$evaluation)),])
}
parallel::stopCluster(cl) 
}else{
res <- NULL
for(i in 1:nRuns){
tempGa <- genalg::rbga.bin(size=nrow(data), iters=iters, popSize=popSize, evalFunc=myScore, verbose=TRUE)
res <- rbind(res, tempGa$population[which(tempGa$evaluation==min(tempGa$evaluation)),])
}
}

### Clean up our results
res <- t(res)
res <- res[,!duplicated(t(res)), drop=FALSE]
corrs <- apply(res, 2, myScore)
res <- res[,order(corrs), drop=FALSE]

colnames(res) <- paste("Solution", 1:ncol(res))
corrs <- corrs[order(corrs)] * -1

return(list(solutions=res, corrs=corrs))
}

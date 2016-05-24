kullbackLeiber <-
function(data, plot=TRUE, parallel=FALSE, cores=3){
if(missing(data))
stop("data missing.")

numData <- length(data)
nameData <- names(data)

data <- lapply(data, function(x) x+1)  # Add 1 so we don't ever get an all 0 comparison

if(!parallel){
cl <- parallel::makeCluster(min(cores, numData)) 
doParallel::registerDoParallel(cl)

results <- foreach::foreach(i=1:numData, .combine=list, .multicombine=TRUE, .inorder=TRUE, .packages=c("dirmult")) %dopar%{
mle.param <- dirmult::dirmult(data[[i]], trace=FALSE)
return(mle.param)
}
parallel::stopCluster(cl)
}else{
results <- vector("list", numData)
for(i in 1:numData)
results[[i]] <- dirmult::dirmult(data[[i]], trace=FALSE)
}

alpha <- lapply(results, function(x) x$gamma)
names(alpha) <- nameData
LL.list <- mapply(function(x, y) LLDM(x, y), x=data, y=alpha)

KLmat <- matrix(0, numData, numData)
for(i in 1:numData){
for(j in 1:numData){
ll <- LLDM(data[[i]], alpha[[j]])
KLmat[i, j] <- LL.list[i]- ll
}
} 
colnames(KLmat) <- nameData
rownames(KLmat) <- nameData

if(plot){
KLdist <- dist(KLmat, method="euclidean")

gplots::heatmap.2(as.matrix(KLdist), dendrogram="both", Rowv=TRUE, Colv=TRUE, trace="none", symm=FALSE,
main="Kullback-Leibler Divergences", margins=c(12,9), density.info="none")
}

return(KLmat)
}

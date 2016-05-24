pairedCompareTwoDataSets <-
function(data1, data2, numPerms=1000, parallel=FALSE, cores=3){
if(missing(data1) || missing(data2))
stop("data is missing.")

if(numPerms <= 0)
stop("The number of permutations must be an integer greater than 0.")

groupSize <- ncol(data1)
if(groupSize != ncol(data2))
stop("Groups must be the same size.")

data <- cbind(data1, data2)
gstar1 <- getMLEandLoglike(data1)$mleTree
gstar2 <- getMLEandLoglike(data2)$mleTree
gstarDistance <- sqrt(sum((gstar1-gstar2)^2))

if(parallel){
cl <- parallel::makeCluster(cores) 
doParallel::registerDoParallel(cl)

permDistances <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
samps <- sample(0:1, groupSize, replace=TRUE)
samps <- samps*groupSize + 1:groupSize

gstar1 <- getMLEandLoglike(data[,samps])$mleTree
gstar2 <- getMLEandLoglike(data[,-samps])$mleTree
return(sqrt(sum((gstar1-gstar2)^2)))
}
parallel::stopCluster(cl) 
}else{
permDistances <- NULL
for(i in 1:numPerms){ 
samps <- sample(0:1, groupSize, replace=TRUE)
samps <- samps*groupSize + 1:groupSize

gstar1 <- getMLEandLoglike(data[,samps])$mleTree
gstar2 <- getMLEandLoglike(data[,-samps])$mleTree

permDistances <- c(permDistances, sqrt(sum((gstar1-gstar2)^2)))
}
}

Pvalue <- sum(permDistances>=gstarDistance)/numPerms 
return(Pvalue)
}

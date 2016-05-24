pairedPvalue <-
function(data, type, groups, numPerms=10, parallel=FALSE, cores=3){
if(missing(data) || missing(type) || missing(groups))
stop("data, type, and/or groups is missing.")

if(numPerms <= 0)
stop("The number of permutations must be an integer greater than 0.")

grp1 <- data[,groups==unique(groups)[1]]
grp2 <- data[,groups==unique(groups)[2]]
groupSize <- ncol(grp1)
if(groupSize != ncol(grp2))
stop("Groups must be the same size.")

data <- cbind(grp1, grp2)
gstarDistance <- calcDistance(estGStar(grp1), estGStar(grp2), type) 

if(parallel){
cl <- parallel::makeCluster(cores) 
doParallel::registerDoParallel(cl)

permDistances <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
samps <- sample(0:1, groupSize, replace=TRUE)
samps <- samps*groupSize + 1:groupSize

gstar1 <- estGStar(data[,samps])
gstar2 <- estGStar(data[,-samps])
return(calcDistance(gstar1, gstar2))
}
parallel::stopCluster(cl) 
}else{
permDistances <- NULL
for(i in 1:numPerms){ 
samps <- sample(0:1, groupSize, replace=TRUE)
samps <- samps*groupSize + 1:groupSize

gstar1 <- estGStar(data[,samps])
gstar2 <- estGStar(data[,-samps])

permDistances <- c(permDistances, calcDistance(gstar1, gstar2))
}
}

Pvalue <- sum(permDistances>=gstarDistance)/numPerms 
return(Pvalue)
}

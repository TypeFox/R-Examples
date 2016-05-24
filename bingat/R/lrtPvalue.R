lrtPvalue <-
function(data, type, groups, numPerms=10, parallel=FALSE, cores=3){
if(missing(data) || missing(type) || missing(groups))
stop("data, type, and/or groups is missing.")

if(numPerms <= 0)
stop("The number of permutations must be an integer greater than 0.")

if(length(unique(groups)) != 2)
stop("There must be exactly two groups.")

x1 <- sum(groups==unique(groups)[1])
x2 <- sum(groups==unique(groups)[2])
if(x1 <= 0 || x2 <= 0)
stop("Each group must have at least 1 subject.")

grp1 <- data[,groups==unique(groups)[1]]
grp2 <- data[,groups==unique(groups)[2]]
groupSize <- ncol(grp1)

gstar <- estGStar(data)
gstar1 <- estGStar(grp1)
gstar2 <- estGStar(grp2)

tau <- estTau(data, type, gstar)
tau1 <- estTau(grp1, type, gstar1)
tau2 <- estTau(grp2, type, gstar2)

ll <- estLogLik(data, type, gstar, tau)
ll1 <- estLogLik(grp1, type, gstar1, tau1)
ll2 <- estLogLik(grp2, type, gstar2, tau2)

e10raw <- -2 * (ll-ll1-ll2)

if(parallel){
cl <- parallel::makeCluster(cores) 
doParallel::registerDoParallel(cl)

lambda  <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
g <- sample(1:ncol(data), groupSize)
grp1 <- data[,g]
grp2 <- data[,-g]

gstar1 <- estGStar(grp1)
gstar2 <- estGStar(grp2)
tau1 <- estTau(grp1, type, gstar1)
tau2 <- estTau(grp2, type, gstar2)
ll1 <- estLogLik(grp1, type, gstar1, tau1)
ll2 <- estLogLik(grp2, type, gstar2, tau2)

e10 <- -2 * (ll-ll1-ll2)
return(e10)
}
parallel::stopCluster(cl) 
}else{
lambda <- NULL
for(i in 1:numPerms){
g <- sample(1:ncol(data), groupSize)
grp1 <- data[,g]
grp2 <- data[,-g]

gstar1 <- estGStar(grp1)
gstar2 <- estGStar(grp2)
tau1 <- estTau(grp1, type, gstar1)
tau2 <- estTau(grp2, type, gstar2)
ll1 <- estLogLik(grp1, type, gstar1, tau1)
ll2 <- estLogLik(grp2, type, gstar2, tau2)

e10 <- -2 * (ll-ll1-ll2)
lambda <- c(lambda, e10)
}
}
pValue <- sum(lambda > e10raw)/numPerms

return(pValue)
}

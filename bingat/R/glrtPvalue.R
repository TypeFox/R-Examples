glrtPvalue <-
function(data, type, groups, numPerms=10, parallel=FALSE, cores=3){
if(missing(data) || missing(type) || missing(groups))
stop("data, type, and/or groups is missing.")

if(numPerms <= 0)
stop("The number of permutations must be an integer greater than 0.")

if(length(unique(groups)) != 2)
stop("There must be exactly two groups.")

if(max(groups) != 1 || min(groups) != 0)
stop("'groups' must use 0 and 1 to denote groups.")

x1 <- sum(groups==unique(groups)[1])
x2 <- sum(groups==unique(groups)[2])

if(x1 <= 0 || x2 <= 0)
stop("Each group must have at least 1 subject.")

reg <- glmReg(data, type, groups)
glrt <- list(glrtReg(data, type, groups))
names(glrt) <- "GLRT"

if(parallel){
cl <- parallel::makeCluster(cores) 
doParallel::registerDoParallel(cl)

GLRTpermut  <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .packages=c("bingat")) %dopar%{
p1 <- sample(x1+x2, x1, replace=FALSE)
grps <- rep(0, x1+x2)
grps[p1] <- 1
glrt.mcnp <- glrtReg(data, type, grps)
return(glrt.mcnp)
}
parallel::stopCluster(cl) 
}else{
GLRTpermut <- apply(as.matrix(1:numPerms), 1, function(x, x1, x2, data, type){
p1 <- sample(x1+x2, x1, replace=FALSE)
grps <- rep(0, x1+x2)
grps[p1] <- 1
glrt.mcnp <- glrtReg(data, type, grps)
return(glrt.mcnp)
}, x1=x1, x2=x2, data=data, type=type)
}

pvalue <- list(sum(GLRTpermut >= glrt)/numPerms)
names(pvalue) <- "pvalue"

#Label our output list
l <- list(reg, glrt, pvalue)
keys <- unique(unlist(lapply(l, names)))
results<- setNames(do.call(mapply, c(FUN=c, lapply(l, `[`, keys))), keys)

return(results)
}

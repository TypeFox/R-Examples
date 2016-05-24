compareTwoDataSets <-
function(data1, data2, numBootStraps=1000, enableMC=FALSE, cores=3){
if(missing(data1) || missing(data2))
stop("Two valid data sets are required.")
if(numBootStraps <= 0)
numBootStraps <- 1

dataComb <- merge(data1, data2, by=0, all=TRUE)
rownames(dataComb) <- dataComb[,1]
dataComb <- dataComb[,-1]
dataComb[is.na(dataComb)] <- 0

data1 <- dataComb[,1:ncol(data1), drop=FALSE]
data2 <- dataComb[,(ncol(data1)+1):ncol(dataComb), drop=FALSE]

fit1 <- getMLEandLoglike(data1)
fit2 <- getMLEandLoglike(data2)
fitComb <- getMLEandLoglike(dataComb)
LRTobs <- -2*(fitComb$Loglik-fit1$Loglik-fit2$Loglik)
if(LRTobs == 0) #exactly the same
return(1)

bootSample <- 1:numBootStraps
x1 <- ncol(data1)
x2 <- ncol(data2)

if(enableMC){
cl <- parallel::makeCluster(cores) 
doParallel::registerDoParallel(cl)

LRTbootstrap <- foreach::foreach(i=1:numBootStraps, .combine=cbind, .inorder=FALSE, .multicombine=TRUE, .export=c("getMLEandLoglike")) %dopar%{
p1 <- sample(x1+x2, replace=TRUE)
data1b <- dataComb[,p1[c(1:x1)], drop=FALSE]
data2b <- dataComb[,p1[c((x1+1):(x1+x2))], drop=FALSE]
datacombb <- cbind(data1b, data2b)

fit1 <- getMLEandLoglike(data1b)
fit2 <- getMLEandLoglike(data2b)
fitComb <- getMLEandLoglike(datacombb)
LRT <- -2*(fitComb$Loglik-fit1$Loglik-fit2$Loglik)
}
parallel::stopCluster(cl) 
}else{ 
LRTbootstrap <- apply(as.matrix(bootSample), 1, function(x, x1, x2, datacomb){
p1 <- sample(x1+x2, replace=TRUE)
data1b <- datacomb[,p1[c(1:x1)], drop=FALSE]
data2b <- datacomb[,p1[c((x1+1):(x1+x2))], drop=FALSE]
datacombb <- cbind(data1b, data2b)

fit1 <- getMLEandLoglike(data1b)
fit2 <- getMLEandLoglike(data2b)
fitcomb <- getMLEandLoglike(datacombb)
LRT <- -2*(fitcomb$Loglik-fit1$Loglik-fit2$Loglik)
}, x1=x1, x2=x2, datacomb=dataComb)
}

pValue <- sum(LRTbootstrap > LRTobs)/numBootStraps
return(pValue)
}

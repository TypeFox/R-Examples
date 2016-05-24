glmReg <-
function(data, type, groups){
if(dim(table(groups)) != 2){
gstar <- estGStar(data)
b0 <- gstar
b0b1 <- NULL
b1 <- NULL
hammingError <- NULL
tau <- estTau(data, type, gstar)
loglik <- estLogLik(data, type, gstar, tau)
}else{
#Estimate the gstars for each group
b0 <- estGStar(data[,groups==unique(groups)[1]])  
b0b1 <- estGStar(data[,groups==unique(groups)[2]])  
b1 <- xor(b0, b0b1)*1

index <- as.matrix(1:length(groups))
hammingError <- apply(index, 1, function(x, b0, b1, grps, datap, typ){
calcDistance(datap[,x], xor(b0, b1*grps[x]), typ)
}, b0=b0, b1=b1, grps=groups, datap=data, typ=type)

#Get the tau and loglik values
gstar <- matrix(0, nrow(data), ncol(data))
for(i in 1:length(groups))
gstar[,i] <- as.matrix(xor(b0, b1*groups[i])) 

tau <- estTau(data, type, gstar)
loglik <- estLogLik(data, type, gstar, tau)
}

results <- list(b0, b1, b0b1, hammingError, loglik, tau)
names(results) <- c("b0.covs0", "b1.Differences", "b0b1.covs1", "hammingError", "loglik", "tau")

return(results)
}

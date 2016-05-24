Xmcupo.sevsample <-
function(group.data){
if(missing(group.data))
stop("group.data is missing.")

K <- ncol(group.data[[1]])
numGrps <- length(group.data)

groupParameterEstimated <- vector("list", numGrps)
for(i in 1:numGrps){
data <- group.data[[i]]
nreads <- rowSums(data)
pi <- colSums(data)/sum(colSums(data))
theta <- weirMoM(data, pi)
groupParameterEstimated[[i]] <- c(pi, theta, nreads)
}

Xmcupo <- Xmcupo.statistics(groupParameterEstimated, K)
pvalue <- 1-pchisq(q=Xmcupo, df=(numGrps-1)*(K-1), ncp=0, lower.tail=TRUE)

ret <- list(Xmcupo, pvalue)
names(ret) <- c("Xmcupo statistics", "p value")

return(ret)
}

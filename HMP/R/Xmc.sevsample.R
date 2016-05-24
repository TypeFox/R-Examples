Xmc.sevsample <-
function(group.data, pi0){
if(missing(group.data) || missing(pi0))
stop("group.data and/or pi0 missing.")

numGroups <- length(group.data)
group.parameter.estimated <- vector("list", numGroups)

if(ncol(group.data[[1]]) != length(pi0))
stop("group.data and pi0 need to be the same length.")

for(x in 1:numGroups){
data <- group.data[[x]]
nreads.data <- as.matrix(rowSums(data))

pi.MoM <- colSums(data)/sum(colSums(data))
theta.MoM <- weirMoM(data, pi.MoM)
group.parameter.estimated[[x]] <- c(pi.MoM, theta.MoM, t(nreads.data))
}

rank.Bj <- sum(pi0>0)-1 
Xmc <- Xmc.statistics(group.parameter.estimated, pi0)
p.value <- 1-pchisq(q=Xmc, df=numGroups*rank.Bj, ncp=0, lower.tail=TRUE)

sevRAD.mean.test <- list(Xmc, p.value)
names(sevRAD.mean.test) <- c("Xmc statistics", "p value")

return(sevRAD.mean.test)
}

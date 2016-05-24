#calculate mantel tests
wassermann <- function(gen.mat, cost.mats, eucl.mat=NULL, plot=TRUE, nperm=999)
{
mats <- cost.mats
if (!is.null(eucl.mat)) mats$Euclidean <- eucl.mat

n.mats <- length(mats)

if (n.mats<2) stop("There are not enough resistance matrices for this approach")
mantel.tab <- data.frame(model=NA, r=NA, p=NA)
cc<- 1


for (i in 1:(n.mats-1))
{

for (j in 2:n.mats)

{

if (i<j)
{
mant1 <- mantel.partial(gen.mat, mats[[i]], mats[[j]],permutations=nperm)
mant2 <- mantel.partial(gen.mat, mats[[j]], mats[[i]],permutations=nperm)


mantel.tab[cc,] <- c(paste("Gen ~",names(mats)[i]," | ", names(mats)[j],sep=""), round(mant1$statistic,4), round(mant1$signif,4))
mantel.tab[cc+1,] <- c(paste("Gen ~",names(mats)[j]," | ", names(mats)[i],sep=""), round(mant2$statistic,4), round(mant2$signif,4))
cc<-cc+2

if (plot==TRUE)
{
x<- densityplot(permustats(mant1),main=paste("Gen ~",names(mats)[i]," | ", names(mats)[j],sep=""), cex.main=0.6)
print(x)

x<- densityplot(permustats(mant2), main=paste("Gen ~",names(mats)[j]," | ", names(mats)[i],sep=""), cex.main=0.6)
print(x)
}

}
}
}

mantel.tab <- mantel.tab[order(mantel.tab$r, decreasing = TRUE),]

return(list(mantel.tab=mantel.tab))
}

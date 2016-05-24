ckappa <-
function(r)
{

r <- na.omit(r)
r1 <- r[,1]
r2 <- r[,2]
n1 <- as.character(r1)
n2 <- as.character(r2)
lev <- levels(as.factor(c(n1,n2)))
p <- length(lev)

tab <- matrix(nrow=p,ncol=p)
dimnames(tab) <- list(levels(as.factor(c(n1,n2))),levels(as.factor(c(n1,n2))))
dim1 <- dimnames(tab)[[1]]
dim2 <- dimnames(tab)[[2]]

tabi <- table(n1,n2)
dimi1 <- dimnames(tabi)[[1]]
dimi2 <- dimnames(tabi)[[2]]

for(i in 1:p)for(j in 1:p)
{
	if((sum(dim1[i]==dimi1)==1)&(sum(dim2[j]==dimi2)==1))
	tab[i,j] <- tabi[dim1[i],dim2[j]]
	else
	tab[i,j] <- 0
}

tsum <- sum(tab)
ttab <- tab/tsum
tm1 <- apply(ttab, 1, sum)
tm2 <- apply(ttab, 2, sum)
agreeP <- sum(diag(ttab))
chanceP <- sum(tm1 * tm2)
kappa2 <- (agreeP - chanceP)/(1 - chanceP)
result <- list("table"=tab,"kappa"=kappa2)
result
}

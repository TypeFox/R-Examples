wkappa <-
function(r,weights="squared") 
{
    r <- na.omit(r)
    r1 <- r[, 1]
    r2 <- r[, 2]
    n1 <- as.character(r1)
    n2 <- as.character(r2)
    lev <- levels(as.factor(c(n1, n2)))
    p <- length(lev)
    if (weights != "squared") weights <- "absolute"
    tab <- matrix(nrow = p, ncol = p)
    weight <- matrix(nrow = p, ncol = p)
    rvect <- as.matrix(r)
    dim(rvect) <- c(dim(r)[1]*dim(r)[2],1)
    if (is.numeric(rvect))
    {
	tabi <- table(r[,1],r[,2])
	tabname <- as.character(sort(as.numeric(levels(as.factor(c(dimnames(tabi)[[1]],dimnames(tabi)[[2]]))))))
        dimnames(tab) <- list(tabname,tabname)
	dim1 <- tabname
	dim2 <- tabname
        dimi1 <- as.character(dimnames(tabi)[[1]])
        dimi2 <- as.character(dimnames(tabi)[[2]])
    }else
    {
    dimnames(tab) <- list(levels(as.factor(c(n1, n2))), levels(as.factor(c(n1, n2))))
    dim1 <- dimnames(tab)[[1]]
    dim2 <- dimnames(tab)[[2]]
    tabi <- table(n1, n2)
    dimi1 <- dimnames(tabi)[[1]]
    dimi2 <- dimnames(tabi)[[2]]
    }
    for (i in 1:p) for (j in 1:p)
    {
        if ((sum(dim1[i] == dimi1) == 1) & (sum(dim2[j] == dimi2) == 1)) tab[i, j] <- tabi[dim1[i], dim2[j]]
        else tab[i, j] <- 0
	if (weights == "squared") weight[i,j] <- 1 - (i - j)^2/(p - 1)^2
        else weight[i,j] <- 1 - abs(i - j)/abs(p - 1)
    }
    tsum <- sum(tab)
    ttab <- tab/tsum
    agreeP <- sum(ttab*weight)
    tm1 <- apply(ttab, 1, sum)
    tm2 <- apply(ttab, 2, sum)
    ttabchance <- tm1%*%t(tm2)
    chanceP <- sum(ttabchance*weight)
    kappa2 <- (agreeP - chanceP)/(1 - chanceP)
    result <- list(table = tab, weights=weights, kappa = kappa2)
    result
}

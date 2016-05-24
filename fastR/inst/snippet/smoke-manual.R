row.sum <- apply(smokeTab,1,sum)
col.sum <- apply(smokeTab,2,sum)
grandTotal <- sum(smokeTab)
e <- outer(row.sum,col.sum)/grandTotal; e
o <- smokeTab
stat <- sum ( (e-o)^2/e); stat
pval <- 1 - pchisq(stat,df=2); pval

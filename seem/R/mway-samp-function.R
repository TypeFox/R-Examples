mway.samp <- function(pnom,plab,perc.v,runs){ 

# latin hypercubes from a uniform distribution
# how many parameters
np <- length(pnom)

# intervals
n <- runs
# param.range has min, max for each parameter
min <- pnom-perc.v*pnom/100; max <- pnom+perc.v*pnom/100

# nominal are means
# divide interval in n equal parts
 pequal <- (max - min)/n
 interval <- matrix(ncol=np, nrow=n)
# sample without replacement
 for(i in 1:np) interval[,i] <- sample(n, replace=F)
# sample from each interval 
 pval <- interval
 for (i in 1:n) {
    for(k in 1:np){
        j <- interval[i,k]
        minj <- pequal[k]*(j-1) + min[k]
        maxj <- pequal[k]*j + min[k]
        pval[i,k] <- runif(1, minj, maxj)
    }
 }
 pval=rbind(pnom, pval,deparse.level=0)

 return(list(plab=plab, pval=round(pval,4), pnom=pnom, fact=F))
}

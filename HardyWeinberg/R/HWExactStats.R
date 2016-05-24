HWExactStats <- function(X,x.linked=FALSE,...) {
nmarkers <- nrow(X)
pvalvec <- numeric(nmarkers)
for (i in 1:nmarkers) {
   pvalvec[i] <- HWExact(X[i,],x.linked=x.linked,...)$pval
}
return(pvalvec)
}


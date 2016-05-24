clx <-
function(fm,cluster){
dfcw=1
M <- length(unique(cluster))
N <- length(cluster)
dfc <- (M/(M-1))*((N-1)/(N-fm$rank))
u <- apply(estfun(fm),2,
function(x) tapply(x, cluster, sum))
vcovCL <- dfc*sandwich(fm, meat.=crossprod(u)/N)*dfcw
vcovCL
coeftest(fm, vcovCL) 
}

`flexmonoreg` <-
function(y, x, lambda=0, alpha.location=1, alpha.adjacency=0.5){
f1 <- flexisoreg(y, x, lambda, alpha.location, alpha.adjacency)
f2 <- flexisoreg(y, 0-x, lambda, alpha.location, alpha.adjacency)
if(f2$statistic>f1$statistic){
m <- max(f2$groups)
return( list(groups=m+1-f2$groups, estimates=f2$estimates, statistic=0-f2$statistic) )
}else{
return( f1 )
}
}

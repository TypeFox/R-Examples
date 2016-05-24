`flexmonoreg.stat` <-
function(y, x, lambda=0, alpha.location=1, alpha.adjacency=0.5){
f1 <- flexisoreg(y, x, lambda, alpha.location, alpha.adjacency)$statistic
f2 <- flexisoreg(y, 0-x, lambda, alpha.location, alpha.adjacency)$statistic
return(ifelse(f2>f1, 0-f2, f1))
}

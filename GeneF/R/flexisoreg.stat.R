`flexisoreg.stat` <-
function(y, x, lambda=0, alpha.location=1, alpha.adjacency=0.5){
return(flexisoreg(y, x, lambda, alpha.location, alpha.adjacency)$statistic)
}

ZT.statistics.Hnull.Ha <-
function(Nrs, fit, type){
if(tolower(type) == "hnull"){
data <- Multinomial(Nrs, probs=fit$pi)
}else if(tolower(type) == "ha"){
data <- Dirichlet.multinomial(Nrs, shape=fit$gamma)
}else{
stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
}

ZT <- c(Z.statistics(data), T.statistics(data))

return(ZT)
}

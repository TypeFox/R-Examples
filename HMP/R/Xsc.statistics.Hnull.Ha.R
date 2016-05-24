Xsc.statistics.Hnull.Ha <-
function(Nrs, fit, type, pi0){
data <- Dirichlet.multinomial(Nrs, shape=fit$gamma)
fit.null <- DM.MoM(data)

if(tolower(type) == "hnull"){
Xsc <- Xsc.statistics(pi1=fit.null$pi, theta=fit.null$theta, Nrs, pi0=fit$pi)
}else if(tolower(type) == "ha"){
Xsc <- Xsc.statistics(pi1=fit.null$pi, theta=fit.null$theta, Nrs, pi0=pi0)
}else{
stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
}

return(Xsc)
}

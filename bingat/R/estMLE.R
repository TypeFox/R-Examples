estMLE <-
function(data, type){
gstar <- estGStar(data)
tau <- estTau(data, type, gstar)

return(list(gstar=gstar, tau=tau))
}

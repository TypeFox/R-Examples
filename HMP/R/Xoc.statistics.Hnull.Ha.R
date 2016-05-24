Xoc.statistics.Hnull.Ha <-
function(Nrs, group.alphap, n.groups, type){
group.data.null <- vector("list", n.groups)

if(tolower(type) == "hnull"){
for(x in 1:n.groups)
group.data.null[[x]] <- Dirichlet.multinomial(Nrs[[x]], shape=group.alphap)
}else if(tolower(type) == "ha"){
for(x in 1:n.groups)
group.data.null[[x]] <- Dirichlet.multinomial(Nrs[[x]], shape=group.alphap[x,]) 
}else{
stop(sprintf("Type '%s' not found. Type must be 'ha' for power or 'hnull' for size.\n", as.character(type)))
}

Xoc <- Xoc.statistics(group.data.null)

return(Xoc)
}

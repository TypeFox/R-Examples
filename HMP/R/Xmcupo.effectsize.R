Xmcupo.effectsize <-
function(group.data){
if(missing(group.data))
stop("data.groups missing.")

par.groups <- lapply(group.data, function(x){
p <- DM.MoM(x)
p$reads <- rowSums(x)
return(p)
})

N.group <- length(par.groups)
Kc <- length(par.groups[[1]]$pi)
N.total <- do.call(sum, lapply(par.groups, function(x){sum(x$reads)}))

group.parameter.estimated <- lapply(par.groups, function(x){c(x$pi, x$theta, x$reads)})
Xmcupo <- Xmcupo.statistics(group.parameter.estimated, K=Kc)

if(Kc >= N.group){
pi.groups <- diag(rep(1, N.group))
}else{
stop("The number of taxa must be greater than the number of groups.")
}

gp.max <- lapply(as.list(1:N.group), function(x,gp.1, pi.groups){
gpp <- gp.1[[x]]
ret <- c(pi.groups[x,], gpp$theta, gpp$reads)
return(ret)
}, gp.1=par.groups, pi.groups=pi.groups)

Xmcupo.gp <- Xmcupo.statistics(gp.max, K=N.group)
CramerV<- sqrt(Xmcupo[1]/(N.total*min(Kc-1, N.group-1)))
Mod.CramerV <- sqrt(Xmcupo[1]/(Xmcupo.gp[1]*min(Kc-1, N.group-1)))

result <- c(Xmcupo[1],CramerV, Mod.CramerV)
names(result) <- c("Chi-Squared", "Cramer Phi", "Modified-Cramer Phi")
return(result)
}

Xdc.statistics <-
function(group.data, epsilon=10^(-4)){
fit <- lapply(group.data, function(x, epsilon){
dirmult::dirmult(x, initscalar=DM.MoM(x)$theta, epsilon=epsilon, trace=FALSE)
}, epsilon=epsilon)

logliks <- unlist(lapply(fit, function(x){x$loglik}))
groupData <- NULL

for(i in 1:length(group.data))
groupData <- rbind(groupData, group.data[[i]])
fit.group <- dirmult::dirmult(groupData, initscalar=DM.MoM(groupData)$theta, epsilon=epsilon, trace=FALSE)

Xdc <- -2*(fit.group$loglik-sum(logliks))

return(Xdc)
}

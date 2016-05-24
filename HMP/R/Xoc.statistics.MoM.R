Xoc.statistics.MoM <-
function(group.data){
fit <- lapply(group.data, DM.MoM)
logliks <- unlist(lapply(fit, function(x){x$loglik}))
groupData <- NULL
numGrps <- length(group.data)

for(i in 1:numGrps)
groupData <- rbind(groupData, group.data[[i]])
fit.group <- DM.MoM(groupData)

pigroups <- lapply(fit, function(x){x$pi})
thetagroup <- fit.group$theta
indexp <- as.matrix(1:numGrps)
equal.theta.loglik <- sum(apply(indexp, 1, 
function(x){loglikDM(group.data[[x]], pigroups[[x]]*(1-thetagroup)/thetagroup)}))

Xoc <- -2*(equal.theta.loglik-sum(logliks))

return(Xoc)
}

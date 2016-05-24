Xdc.statistics.MoM <-
function(group.data){
fit <- lapply(group.data, DM.MoM)

logliks <- unlist(lapply(fit, function(x){x$loglik}))
groupData <- NULL

for(i in 1:length(group.data))
groupData <- rbind(groupData, group.data[[i]])
fit.group <- DM.MoM(groupData)

Xdc <- -2*(fit.group$loglik-sum(logliks))

return(Xdc)
}

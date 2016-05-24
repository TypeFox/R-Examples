Xoc.sevsample <-
function(group.data, epsilon=10^(-4)){
if(missing(group.data))
stop("group.data missing.")

n.groups <- length(group.data)
Xoc <- Xoc.statistics(group.data, epsilon)

p.value <- 1-pchisq(q=Xoc, df=(n.groups-1), ncp=0, lower.tail=TRUE)

sev.overd.test <- list(Xoc, p.value)
names(sev.overd.test) <- c("Xoc statistics", "p value")

return(sev.overd.test)
}

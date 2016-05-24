`ordicluster2` <-
function(ordiplot,cluster,mingroups=1,maxgroups=nrow(ordiplot$sites),...){
    mingroups <- max(1,mingroups)
    maxgroups <- min(maxgroups,nrow(ordiplot$sites))
    n <- nrow(ordiplot$sites)
    for (i in mingroups:maxgroups) {
        groups <- cutree(cluster,k=i)
        ordihull(ordiplot,groups,...)
    }
}


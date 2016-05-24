change.order.clusters <- function(clustering, new.order){
    n.cluster=length(new.order)
    stopifnot(1:n.cluster==new.order[order(new.order)])
    stopifnot(max(clustering$cluster)==n.cluster)

    #move all cluster numbers above actual range
    clustering$cluster <- clustering$cluster+n.cluster
    centers.new <- clustering$centers
    #store old memberships
    membership <- clustering$membership

    #assign new order
    for(i in 1:n.cluster){
        clustering$cluster[clustering$cluster==i+n.cluster] <- new.order[i] 
        centers.new[new.order[i],] <- clustering$centers[i,]
        clustering$membership[,new.order[i]] <- membership[,i]
    }

    clustering$centers <- centers.new 
    return(clustering)
}

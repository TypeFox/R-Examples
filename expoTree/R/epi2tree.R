epi2tree <- function(epi) {
  make.internal <- function(epi,node,parent,cur.id) {
    node.i <- which(epi$id==node)
    dt <- epi$dtimes[node.i]
    # get offspring
    offspring <- epi$parent==node
    if (sum(offspring) > 0) {
      off.id <- epi$id[offspring]
      # create internal node for each branching event
      inf.t <- epi$itimes[offspring]
      dt <- epi$dtimes[node.i]
      intervals <- c(epi$itimes[node.i],inf.t) - c(inf.t,dt)
      ages <- c(inf.t,dt)
      new.id <- 1:length(inf.t) + cur.id
      cur.id <- max(new.id)
      parents <- c(parent,new.id)
      edge.list <- cbind(parents,c(new.id,node),intervals,ages)
      for (j in 1:length(off.id)) {
        out <- make.internal(epi=epi,node=off.id[j],parent=new.id[j],
                             cur.id=cur.id)
        edge.list <- rbind(edge.list,out[[1]])
        cur.id <- out[[2]]
      }
      return(list(edge.list,cur.id))
    } else {
      edge.list <- c(parent,node,epi$itimes[node.i]-dt,dt)
      return(list(edge.list,cur.id))
    }
  }
  edge.list <- make.internal(epi,0,-1,max(epi$id))
  edges <- matrix(edge.list[[1]][-1,1:2]+1,ncol=2)
  root.edge <- edge.list[[1]][1,3]
  edge.lengths <- edge.list[[1]][-1,3]
  tip.label <- as.character(epi$id)
  coords <- epi.coords(epi)[[2]]
  coords[,1] <- coords[,1]+1
  tree <- list(edge=edges,edge.length=edge.lengths,tip.label=tip.label,
               Nnode=max(edge.list[[1]])-max(epi$id),root.edge=root.edge,
               coords=coords)
  class(tree) <- "phylo"
  tree
}

